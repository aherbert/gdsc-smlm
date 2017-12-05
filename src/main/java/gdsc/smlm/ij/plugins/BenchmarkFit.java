package gdsc.smlm.ij.plugins;

import java.awt.Rectangle;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;

import gdsc.core.ij.Utils;
import gdsc.core.utils.Statistics;
import gdsc.core.utils.StoredDataStatistics;
import gdsc.core.utils.TextUtils;
import gdsc.smlm.data.config.CalibrationReader;
import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.PSFProtos.PSFType;
import gdsc.smlm.data.config.PSFProtosHelper;
import gdsc.smlm.engine.FitConfiguration;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolver;
import gdsc.smlm.fitting.Gaussian2DFitter;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.ij.IJImageSource;
import gdsc.smlm.ij.plugins.CreateData.BenchmarkParameters;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.IJImageConverter;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.ExtendedGenericDialog;
import ij.plugin.PlugIn;
import ij.plugin.WindowOrganiser;
import ij.text.TextWindow;

/**
 * Fits the benchmark image created by CreateData plugin.
 */
public class BenchmarkFit implements PlugIn
{
	private static final String TITLE = "Benchmark Fit";

	private static int regionSize = 3;
	private static double lastS = 0;
	private static boolean offsetFitting = true;
	private static double startOffset = 0.5;
	private static boolean comFitting = false;
	private static boolean backgroundFitting = true;
	private static boolean signalFitting = true;
	private static boolean showHistograms = false;
	private static boolean saveRawData = false;
	private static String rawDataDirectory = "";
	private static int histogramBins = 100;

	private static TextWindow summaryTable = null, analysisTable = null;

	private static final String[] NAMES = new String[] { "dB (photons)", "dSignal (photons)", "dAngle (deg)", "dX (nm)",
			"dY (nm)", "dSx (nm)", "dSy (nm)", "Time (ms)", "dActualSignal (photons)", "dSax (nm)", "dSay (nm)" };
	private static final int TIME = 7;
	private static final int ACTUAL_SIGNAL = 8;
	private static final int ADJUSTED_X_SD = 9;
	private static final int ADJUSTED_Y_SD = 10;
	private static boolean[] displayHistograms = new boolean[NAMES.length];
	static
	{
		for (int i = 0; i < displayHistograms.length; i++)
			displayHistograms[i] = true;
	}

	private FitConfiguration fitConfig;
	private ImagePlus imp;
	private CreateData.BenchmarkParameters benchmarkParameters;
	private double[] answer = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
	private Rectangle region = null;
	private AtomicInteger comValid = new AtomicInteger();

	// Used to store all the results for cross-method comparison
	private double[][] results;
	private long[] resultsTime;

	public class BenchmarkResult
	{
		/**
		 * The parameters used to create the data
		 */
		final BenchmarkParameters benchmarkParameters;
		/**
		 * The actual parameters (with XY position adjusted for the region size)
		 */
		final double[] answer;
		/**
		 * The string description of the parameters used to create and then fit the data
		 */
		final String parameters;
		/**
		 * Results conversion factors
		 */
		final double[] convert;
		/**
		 * The results of fitting the data. Results are only stored if fitting was successful.
		 */
		final double[][] results;
		/**
		 * The time for fitting the data.
		 */
		final long[] resultsTime;

		public BenchmarkResult(BenchmarkParameters benchmarkParameters, double[] answer, String parameters,
				double[] convert, double[][] results, long[] resultsTime)
		{
			this.benchmarkParameters = benchmarkParameters;
			this.answer = answer;
			this.parameters = parameters;
			this.convert = convert;
			this.results = results;
			this.resultsTime = resultsTime;
		}
	}

	/**
	 * Store all the results from fitting on the same benchmark dataset
	 */
	public static LinkedList<BenchmarkResult> benchmarkResults = new LinkedList<BenchmarkFit.BenchmarkResult>();

	/**
	 * Used to allow multi-threading of the fitting method
	 */
	private class Worker implements Runnable
	{
		volatile boolean finished = false;
		final BlockingQueue<Integer> jobs;
		final Statistics[] stats = new Statistics[NAMES.length];
		final ImageStack stack;
		final Rectangle region;
		final double[][] xy;
		final FitConfiguration fitConfig;
		final double sa;

		double[] data = null;
		private double[] lb, ub = null;
		private double[] lc, uc = null;

		public Worker(BlockingQueue<Integer> jobs, ImageStack stack, Rectangle region, FitConfiguration fitConfig)
		{
			this.jobs = jobs;
			this.stack = stack;
			this.region = region;
			this.fitConfig = fitConfig.clone();
			this.xy = getStartPoints();

			for (int i = 0; i < stats.length; i++)
				stats[i] = (showHistograms || saveRawData) ? new StoredDataStatistics() : new Statistics();
			sa = getSa();

			createBounds();
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		public void run()
		{
			try
			{
				while (true)
				{
					Integer job = jobs.take();
					if (job == null || job.intValue() < 0)
						break;
					if (!finished)
						// Only run if not finished to allow queue to be emptied
						run(job.intValue());
				}
			}
			catch (InterruptedException e)
			{
				System.out.println(e.toString());
				throw new RuntimeException(e);
			}
			finally
			{
				finished = true;
			}
		}

		private void run(int frame)
		{
			if (Utils.isInterrupted())
			{
				finished = true;
				return;
			}

			showProgress();

			// Extract the data
			data = IJImageConverter.getDoubleData(stack.getPixels(frame + 1), stack.getWidth(), stack.getHeight(),
					region, data);

			final int size = region.height;
			final int totalFrames = benchmarkParameters.frames;

			// Subtract the bias
			final double bias = benchmarkParameters.bias;
			for (int i = 0; i < data.length; i++)
				data[i] -= bias;
			// Remove the gain
			double gain = 1;
			if (!fitConfig.isFitCameraCounts())
			{
				gain = benchmarkParameters.gain;
				for (int i = 0; i < data.length; i++)
					data[i] /= gain;
			}

			// Get the background and signal estimate
			final double b = (backgroundFitting) ? getBackground(data, size, size)
					: answer[Gaussian2DFunction.BACKGROUND] * gain;
			final double signal = (signalFitting) ? getSignal(data, b)
					//: benchmarkParameters.p[frame];
					: answer[Gaussian2DFunction.SIGNAL] * gain;

			// Find centre-of-mass estimate
			if (comFitting)
			{
				getCentreOfMass(data, size, size, xy[xy.length - 1]);

				double dx = xy[xy.length - 1][0] - answer[Gaussian2DFunction.X_POSITION];
				double dy = xy[xy.length - 1][1] - answer[Gaussian2DFunction.Y_POSITION];
				if (dx * dx + dy * dy < startOffset * startOffset)
					comValid.incrementAndGet();
			}

			double[] initialParams = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
			initialParams[Gaussian2DFunction.BACKGROUND] = b;
			initialParams[Gaussian2DFunction.SIGNAL] = signal;
			initialParams[Gaussian2DFunction.X_SD] = fitConfig.getInitialXSD();
			initialParams[Gaussian2DFunction.Y_SD] = fitConfig.getInitialYSD();

			double[][] bounds = null;
			double[][] result = new double[xy.length][];
			long[] time = new long[xy.length];
			int c = 0;
			int resultPosition = frame;
			for (double[] centre : xy)
			{
				long start = System.nanoTime();

				// Do fitting
				final double[] params = initialParams.clone();
				params[Gaussian2DFunction.X_POSITION] = centre[0];
				params[Gaussian2DFunction.Y_POSITION] = centre[1];
				fitConfig.initialise(1, size, size, params);
				FunctionSolver solver = fitConfig.getFunctionSolver();
				if (solver.isBounded())
					bounds = setBounds(solver, initialParams, bounds);
				else if (solver.isConstrained())
					setConstraints(solver);

				final FitStatus status = solver.fit(data, null, params, null);
				if (isValid(status, params, size))
				{
					// TODO - Check this is OK for the MLE camera model
					// That estimates counts. We require an estimate of photons.
					if (fitConfig.isFitCameraCounts())
					{
						// Update all the parameters to be in photons
						params[Gaussian2DFunction.BACKGROUND] /= gain;
						params[Gaussian2DFunction.SIGNAL] /= gain;
					}
					result[c] = params;
					time[c] = System.nanoTime() - start;
					// Store all the results for later analysis
					results[resultPosition] = params;
					resultsTime[resultPosition] = time[c];
					c++;
				}
				else
				{
					//System.out.println(status);
				}
				resultPosition += totalFrames;
			}

			addResults(stats, answer, benchmarkParameters.p[frame], sa, time, result, c);
		}

		/**
		 * Set background using the average value of the edge in the data
		 * 
		 * @param data
		 * @param maxx
		 * @param maxy
		 * @return The background
		 */
		private double getBackground(double[] data, int maxx, int maxy)
		{
			return Gaussian2DFitter.getBackground(data, maxx, maxy, 1);
		}

		private double[][] setBounds(FunctionSolver solver, double[] params, double[][] bounds)
		{
			if (bounds == null)
			{
				double[] lower = null;
				double[] upper = null;
				// Check the bounds 
				if (params[Gaussian2DFunction.BACKGROUND] < lb[Gaussian2DFunction.BACKGROUND])
				{
					lower = lb.clone();
					lower[Gaussian2DFunction.BACKGROUND] = params[Gaussian2DFunction.BACKGROUND] -
							Math.abs(lb[Gaussian2DFunction.BACKGROUND] - params[Gaussian2DFunction.BACKGROUND]);
					if (fitConfig.requireStrictlyPositiveFunction() && lower[Gaussian2DFunction.BACKGROUND] < 0)
						lower[Gaussian2DFunction.BACKGROUND] = 0;
				}
				if (params[Gaussian2DFunction.BACKGROUND] > ub[Gaussian2DFunction.BACKGROUND])
				{
					upper = ub.clone();
					upper[Gaussian2DFunction.BACKGROUND] = params[Gaussian2DFunction.BACKGROUND] +
							Math.abs(ub[Gaussian2DFunction.BACKGROUND] - params[Gaussian2DFunction.BACKGROUND]);
				}
				if (params[Gaussian2DFunction.SIGNAL] < lb[Gaussian2DFunction.SIGNAL])
				{
					if (lower == null)
						lower = lb.clone();
					lower[Gaussian2DFunction.SIGNAL] = params[Gaussian2DFunction.SIGNAL] -
							Math.abs(lb[Gaussian2DFunction.SIGNAL] - params[Gaussian2DFunction.SIGNAL]);
					if (fitConfig.requireStrictlyPositiveFunction() && lower[Gaussian2DFunction.SIGNAL] < 0)
						lower[Gaussian2DFunction.SIGNAL] = 0;
				}
				if (params[Gaussian2DFunction.SIGNAL] > ub[Gaussian2DFunction.SIGNAL])
				{
					if (upper == null)
						upper = ub.clone();
					upper[Gaussian2DFunction.SIGNAL] = params[Gaussian2DFunction.SIGNAL] +
							Math.abs(ub[Gaussian2DFunction.SIGNAL] - params[Gaussian2DFunction.SIGNAL]);
				}
				if (lower == null)
					lower = lb;
				if (upper == null)
					upper = ub;
				bounds = new double[][] { lower, upper };
			}
			solver.setBounds(bounds[0], bounds[1]);
			return bounds;
		}

		private void createBounds()
		{
			if (ub == null)
			{
				ub = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
				lb = new double[ub.length];

				// Background could be zero so always have an upper limit
				ub[Gaussian2DFunction.BACKGROUND] = Math.max(0, 2 * benchmarkParameters.getBackground());
				double signal = benchmarkParameters.getSignal();
				lb[Gaussian2DFunction.SIGNAL] = signal * 0.5;
				ub[Gaussian2DFunction.SIGNAL] = signal * 2;
				ub[Gaussian2DFunction.X_POSITION] = 2 * regionSize + 1;
				ub[Gaussian2DFunction.Y_POSITION] = 2 * regionSize + 1;
				lb[Gaussian2DFunction.ANGLE] = -Math.PI;
				ub[Gaussian2DFunction.ANGLE] = Math.PI;
				lb[Gaussian2DFunction.Z_POSITION] = Double.NEGATIVE_INFINITY;
				ub[Gaussian2DFunction.Z_POSITION] = Double.POSITIVE_INFINITY;
				double wf = 1.5;
				double s = benchmarkParameters.s / benchmarkParameters.a;
				lb[Gaussian2DFunction.X_SD] = s / wf;
				ub[Gaussian2DFunction.X_SD] = s * wf;
				lb[Gaussian2DFunction.Y_SD] = s / wf;
				ub[Gaussian2DFunction.Y_SD] = s * wf;
				if (fitConfig.isFitCameraCounts())
				{
					final double gain = benchmarkParameters.gain;
					// Update all the parameters affected by gain
					lb[Gaussian2DFunction.BACKGROUND] *= gain;
					lb[Gaussian2DFunction.SIGNAL] *= gain;
					ub[Gaussian2DFunction.BACKGROUND] *= gain;
					ub[Gaussian2DFunction.SIGNAL] *= gain;
				}
			}
		}

		private void setConstraints(FunctionSolver solver)
		{
			if (uc == null)
			{
				lc = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
				uc = new double[lc.length];
				Arrays.fill(lc, Float.NEGATIVE_INFINITY);
				Arrays.fill(uc, Float.POSITIVE_INFINITY);
				lc[Gaussian2DFunction.BACKGROUND] = 0;
				lc[Gaussian2DFunction.SIGNAL] = 0;
			}
			solver.setConstraints(lc, uc);
		}

		private boolean isValid(FitStatus status, double[] params, int size)
		{
			if (status != FitStatus.OK)
			{
				return false;
			}

			// Reject fits that are outside the bounds of the data
			if (params[Gaussian2DFunction.SIGNAL] < 0 || params[Gaussian2DFunction.X_POSITION] < 0 ||
					params[Gaussian2DFunction.Y_POSITION] < 0 || params[Gaussian2DFunction.X_POSITION] > size ||
					params[Gaussian2DFunction.Y_POSITION] > size)
			{
				return false;
			}

			// Q. Should we do width bounds checking?
			if (fitConfig.isXSDFitting())
			{
				if (params[Gaussian2DFunction.X_SD] < lb[Gaussian2DFunction.X_SD] ||
						params[Gaussian2DFunction.X_SD] > ub[Gaussian2DFunction.X_SD])
				{
					return false;
				}
			}
			if (fitConfig.isYSDFitting())
			{
				if (params[Gaussian2DFunction.Y_SD] < lb[Gaussian2DFunction.Y_SD] ||
						params[Gaussian2DFunction.Y_SD] > ub[Gaussian2DFunction.Y_SD])
				{
					return false;
				}
			}

			return true;
		}

	}

	/**
	 * Add the results to the statistics
	 * 
	 * @param stats
	 * @param answer
	 * @param photons
	 * @param sa
	 * @param time
	 * @param result
	 * @param c
	 *            Count of the number of results
	 */
	private static void addResults(Statistics[] stats, double[] answer, double photons, double sa, long[] time,
			double[][] result, int c)
	{
		// Store the results from each run
		for (int i = 0; i < c; i++)
		{
			addResult(stats, answer, photons, sa, result[i], time[i]);
		}
	}

	/**
	 * Add the given results to the statistics
	 * 
	 * @param stats
	 * @param answer
	 * @param photons
	 * @param sa
	 * @param result
	 * @param time
	 */
	private static void addResult(Statistics[] stats, double[] answer, double photons, double sa, double[] result,
			long time)
	{
		for (int j = 0; j < result.length; j++)
		{
			stats[j].add(result[j] - answer[j]);
		}
		stats[7].add(time);
		stats[8].add(result[Gaussian2DFunction.SIGNAL] - photons);
		stats[9].add(result[Gaussian2DFunction.X_SD] - sa);
		stats[10].add(result[Gaussian2DFunction.Y_SD] - sa);
	}

	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if ("analysis".equals(arg))
		{
			if (benchmarkResults.isEmpty())
			{
				IJ.error(TITLE,
						"No benchmark results in memory.\n \n" + TextUtils.wrap(
								"Run the Fit Benchmark Data plugin and results will be stored for comparison analysis.",
								60));
				return;
			}
			runAnalysis();
		}
		else
		{
			if (CreateData.benchmarkParameters == null)
			{
				IJ.error(TITLE,
						"No benchmark parameters in memory.\n \n" + TextUtils.wrap(
								"Run the " + CreateData.TITLE +
										" plugin in benchmark mode with a fixed number of photons per localisation.",
								60));
				return;
			}
			benchmarkParameters = CreateData.benchmarkParameters;
			imp = CreateData.getImage();
			if (imp == null || imp.getStackSize() != benchmarkParameters.frames)
			{
				IJ.error(TITLE, "No benchmark image to match the parameters in memory");
				return;
			}

			if (!showDialog())
				return;

			run();
		}
	}

	private boolean showDialog()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		final double sa = getSa();
		gd.addMessage(
				String.format("Fits the benchmark image created by CreateData plugin.\nPSF width = %s, adjusted = %s",
						Utils.rounded(benchmarkParameters.s / benchmarkParameters.a), Utils.rounded(sa)));

		FitEngineConfiguration config = SettingsManager.readFitEngineConfiguration(0);
		fitConfig = config.getFitConfiguration();
		fitConfig.setNmPerPixel(benchmarkParameters.a);

		// For each new benchmark width, reset the PSF width to the square pixel adjustment
		if (lastS != benchmarkParameters.s)
		{
			lastS = benchmarkParameters.s;
			fitConfig.setInitialPeakStdDev(sa);
		}

		gd.addSlider("Region_size", 2, 20, regionSize);
		PeakFit.addPSFOptions(gd, fitConfig);
		gd.addChoice("Fit_solver", SettingsManager.getFitSolverNames(), fitConfig.getFitSolver().ordinal());
		gd.addCheckbox("Offset_fit", offsetFitting);
		gd.addNumericField("Start_offset", startOffset, 3);
		gd.addCheckbox("Include_CoM_fit", comFitting);
		gd.addCheckbox("Background_fitting", backgroundFitting);
		gd.addMessage("Signal fitting can be disabled for " + PSFProtosHelper.getName(PSFType.ONE_AXIS_GAUSSIAN_2D) +
				" function");
		gd.addCheckbox("Signal_fitting", signalFitting);
		gd.addCheckbox("Show_histograms", showHistograms);
		gd.addCheckbox("Save_raw_data", saveRawData);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		regionSize = (int) Math.abs(gd.getNextNumber());
		fitConfig.setPSFType(PeakFit.getPSFTypeValues()[gd.getNextChoiceIndex()]);
		fitConfig.setFitSolver(gd.getNextChoiceIndex());
		offsetFitting = gd.getNextBoolean();
		startOffset = Math.abs(gd.getNextNumber());
		comFitting = gd.getNextBoolean();
		backgroundFitting = gd.getNextBoolean();
		signalFitting = gd.getNextBoolean();
		showHistograms = gd.getNextBoolean();
		saveRawData = gd.getNextBoolean();

		gd.collectOptions();

		if (!comFitting && !offsetFitting)
		{
			IJ.error(TITLE, "No initial fitting positions");
			return false;
		}

		if (regionSize < 1)
			regionSize = 1;

		if (gd.invalidNumber())
			return false;

		// Initialise the correct calibration
		CalibrationWriter calibration = CalibrationWriter.create(fitConfig.getCalibration());
		calibration.setNmPerPixel(benchmarkParameters.a);
		calibration.setCountPerPhoton(benchmarkParameters.gain);
		calibration.setQuantumEfficiency(benchmarkParameters.qe);
		calibration.setBias(benchmarkParameters.bias);
		calibration.setCameraType(benchmarkParameters.cameraType);
		calibration.setReadNoise(benchmarkParameters.readNoise);
		calibration.setExposureTime(1000);
		fitConfig.setCalibration(calibration.getCalibration());

		if (!PeakFit.configureFitSolver(config, IJImageSource.getBounds(imp), 0))
			return false;

		if (showHistograms)
		{
			gd = new ExtendedGenericDialog(TITLE);
			gd.addMessage("Select the histograms to display");
			gd.addNumericField("Histogram_bins", histogramBins, 0);

			double[] convert = getConversionFactors();

			for (int i = 0; i < displayHistograms.length; i++)
				if (convert[i] != 0)
					gd.addCheckbox(NAMES[i].replace(' ', '_'), displayHistograms[i]);
			gd.showDialog();
			if (gd.wasCanceled())
				return false;
			histogramBins = (int) Math.abs(gd.getNextNumber());
			for (int i = 0; i < displayHistograms.length; i++)
				if (convert[i] != 0)
					displayHistograms[i] = gd.getNextBoolean();
		}

		return true;
	}

	private double getSa()
	{
		final double sa = PSFCalculator.squarePixelAdjustment(benchmarkParameters.s, benchmarkParameters.a) /
				benchmarkParameters.a;
		return sa;
	}

	/** The total progress. */
	int progress, stepProgress, totalProgress;

	/**
	 * Show progress.
	 */
	private synchronized void showProgress()
	{
		if (progress % stepProgress == 0)
		{
			if (Utils.showStatus("Frame: " + progress + " / " + totalProgress))
				IJ.showProgress(progress, totalProgress);
		}
		progress++;
	}

	private void run()
	{
		// Initialise the answer.
		answer[Gaussian2DFunction.BACKGROUND] = benchmarkParameters.getBackground();
		answer[Gaussian2DFunction.SIGNAL] = benchmarkParameters.getSignal();
		answer[Gaussian2DFunction.X_POSITION] = benchmarkParameters.x;
		answer[Gaussian2DFunction.Y_POSITION] = benchmarkParameters.y;
		answer[Gaussian2DFunction.X_SD] = benchmarkParameters.s / benchmarkParameters.a;
		answer[Gaussian2DFunction.Y_SD] = benchmarkParameters.s / benchmarkParameters.a;

		// Set up the fit region. Always round down since 0.5 is the centre of the pixel.
		int x = (int) benchmarkParameters.x;
		int y = (int) benchmarkParameters.y;
		region = new Rectangle(x - regionSize, y - regionSize, 2 * regionSize + 1, 2 * regionSize + 1);
		if (!new Rectangle(0, 0, imp.getWidth(), imp.getHeight()).contains(region))
		{
			// Check if it is incorrect by only 1 pixel
			if (region.width <= imp.getWidth() + 1 && region.height <= imp.getHeight() + 1)
			{
				Utils.log("Adjusting region %s to fit within image bounds (%dx%d)", region.toString(), imp.getWidth(),
						imp.getHeight());
				region = new Rectangle(0, 0, imp.getWidth(), imp.getHeight());
			}
			else
			{
				IJ.error(TITLE, "Fit region does not fit within the image");
				return;
			}
		}

		// Adjust the centre & account for 0.5 pixel offset during fitting
		x -= region.x;
		y -= region.y;
		answer[Gaussian2DFunction.X_POSITION] -= (region.x + 0.5);
		answer[Gaussian2DFunction.Y_POSITION] -= (region.y + 0.5);

		// Configure for fitting
		fitConfig.setBackgroundFitting(backgroundFitting);
		fitConfig.setNotSignalFitting(!signalFitting);
		fitConfig.setComputeDeviations(false);

		final ImageStack stack = imp.getImageStack();

		// Create a pool of workers
		int nThreads = Prefs.getThreads();
		BlockingQueue<Integer> jobs = new ArrayBlockingQueue<Integer>(nThreads * 2);
		List<Worker> workers = new LinkedList<Worker>();
		List<Thread> threads = new LinkedList<Thread>();
		for (int i = 0; i < nThreads; i++)
		{
			Worker worker = new Worker(jobs, stack, region, fitConfig);
			Thread t = new Thread(worker);
			workers.add(worker);
			threads.add(t);
			t.start();
		}

		final int totalFrames = benchmarkParameters.frames;

		// Store all the fitting results
		results = new double[totalFrames * getNumberOfStartPoints()][];
		resultsTime = new long[results.length];

		// Fit the frames
		totalProgress = totalFrames;
		stepProgress = Utils.getProgressInterval(totalProgress);
		progress = 0;
		for (int i = 0; i < totalFrames; i++)
		{
			// Only fit if there were simulated photons
			if (benchmarkParameters.p[i] > 0)
			{
				put(jobs, i);
			}
		}
		// Finish all the worker threads by passing in a null job
		for (int i = 0; i < threads.size(); i++)
		{
			put(jobs, -1);
		}

		// Wait for all to finish
		for (int i = 0; i < threads.size(); i++)
		{
			try
			{
				threads.get(i).join();
			}
			catch (InterruptedException e)
			{
				e.printStackTrace();
			}
		}
		threads.clear();

		if (comFitting)
			Utils.log(TITLE + ": CoM within start offset = %d / %d (%s%%)", comValid.intValue(), totalFrames,
					Utils.rounded((100.0 * comValid.intValue()) / totalFrames));

		IJ.showProgress(1);
		IJ.showStatus("Collecting results ...");

		// Collect the results
		Statistics[] stats = new Statistics[NAMES.length];
		for (int i = 0; i < workers.size(); i++)
		{
			Statistics[] next = workers.get(i).stats;
			for (int j = 0; j < next.length; j++)
			{
				if (stats[j] == null)
					stats[j] = next[j];
				else
					stats[j].add(next[j]);
			}
		}
		workers.clear();

		// Show a table of the results
		summariseResults(stats);

		// Optionally show histograms
		if (showHistograms)
		{
			IJ.showStatus("Calculating histograms ...");

			int[] idList = new int[NAMES.length];
			int count = 0;
			double[] convert = getConversionFactors();

			boolean requireRetile = false;
			for (int i = 0; i < NAMES.length; i++)
			{
				if (displayHistograms[i] && convert[i] != 0)
				{
					// We will have to convert the values...
					double[] tmp = ((StoredDataStatistics) stats[i]).getValues();
					for (int j = 0; j < tmp.length; j++)
						tmp[j] *= convert[i];
					StoredDataStatistics tmpStats = new StoredDataStatistics(tmp);
					idList[count++] = Utils.showHistogram(TITLE, tmpStats, NAMES[i], 0, 0, histogramBins,
							String.format("%s +/- %s", Utils.rounded(tmpStats.getMean()),
									Utils.rounded(tmpStats.getStandardDeviation())));
					requireRetile = requireRetile || Utils.isNewWindow();
				}
			}

			if (count > 0 && requireRetile)
			{
				idList = Arrays.copyOf(idList, count);
				new WindowOrganiser().tileWindows(idList);
			}
		}

		if (saveRawData)
		{
			String dir = Utils.getDirectory("Data_directory", rawDataDirectory);
			if (dir != null)
				saveData(stats, dir);
		}

		IJ.showStatus("");
	}

	private void saveData(Statistics[] stats, String dir)
	{
		rawDataDirectory = dir;
		for (int i = 0; i < NAMES.length; i++)
		{
			saveStatistics((StoredDataStatistics) stats[i], NAMES[i]);
		}
	}

	private void saveStatistics(StoredDataStatistics stats, String title)
	{
		String filename = rawDataDirectory + title.replace(" ", "_") + ".txt";

		BufferedWriter out = null;
		try
		{
			FileOutputStream fos = new FileOutputStream(filename);
			out = new BufferedWriter(new OutputStreamWriter(fos, "UTF-8"));
			//out.write(title);
			//out.newLine();
			double[] data = stats.getValues();
			Arrays.sort(data);
			for (double d : data)
			{
				//out.write(Utils.rounded(d, 4)); // rounded
				out.write(Double.toString(d));
				out.newLine();
			}
		}
		catch (Exception e)
		{
		}
		finally
		{
			if (out != null)
			{
				try
				{
					out.close();
				}
				catch (IOException e)
				{
				}
			}
		}
	}

	private void put(BlockingQueue<Integer> jobs, int i)
	{
		try
		{
			jobs.put(i);
		}
		catch (InterruptedException e)
		{
			throw new RuntimeException("Unexpected interruption", e);
		}
	}

	/**
	 * @return The starting points for the fitting
	 */
	private double[][] getStartPoints()
	{
		double[][] xy = new double[getNumberOfStartPoints()][];
		int ii = 0;

		if (offsetFitting)
		{
			if (startOffset == 0)
			{
				xy[ii++] = new double[] { answer[Gaussian2DFunction.X_POSITION],
						answer[Gaussian2DFunction.Y_POSITION] };
			}
			else
			{
				// Fit using region surrounding the point. Use -1,-1 : -1:1 : 1,-1 : 1,1 directions at 
				// startOffset pixels total distance 
				final double distance = Math.sqrt(startOffset * startOffset * 0.5);

				for (int x = -1; x <= 1; x += 2)
					for (int y = -1; y <= 1; y += 2)
					{
						xy[ii++] = new double[] { answer[Gaussian2DFunction.X_POSITION] + x * distance,
								answer[Gaussian2DFunction.Y_POSITION] + y * distance };
					}
			}
		}
		// Add space for centre-of-mass at the end of the array
		if (comFitting)
			xy[ii++] = new double[2];
		return xy;
	}

	private int getNumberOfStartPoints()
	{
		int n = (offsetFitting) ? 1 : 0;
		if (startOffset > 0)
			n *= 4;
		return (comFitting) ? n + 1 : n;
	}

	/**
	 * Sum the intensity above background to estimate the signal
	 * 
	 * @param data
	 * @param b
	 *            background
	 * @return The signal
	 */
	public static double getSignal(double[] data, double b)
	{
		double s = 0;
		for (double d : data)
			s += d;
		// Subtract the background per pixel and ensure at least 1 photon in the signal
		return Math.max(1, s - b * data.length);
	}

	/**
	 * Get the centre of mass of the data
	 * 
	 * @param data
	 * @param maxx
	 * @param maxy
	 * @param com
	 *            The centre-of-mass
	 */
	public static void getCentreOfMass(double[] data, int maxx, int maxy, double[] com)
	{
		com[0] = com[1] = 0;
		double sum = 0;
		for (int y = 0, index = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++, index++)
			{
				final double value = data[index];
				sum += value;
				com[0] += x * value;
				com[1] += y * value;
			}
		}

		for (int i = 2; i-- > 0;)
		{
			com[i] /= sum;
		}
	}

	private void summariseResults(Statistics[] stats)
	{
		createTable();

		StringBuilder sb = new StringBuilder();

		// Create the benchmark settings and the fitting settings
		sb.append(benchmarkParameters.getMolecules()).append('\t');
		sb.append(Utils.rounded(benchmarkParameters.getSignal())).append('\t');
		sb.append(Utils.rounded(benchmarkParameters.s)).append('\t');
		sb.append(Utils.rounded(benchmarkParameters.a)).append('\t');
		sb.append(Utils.rounded(getSa() * benchmarkParameters.a)).append('\t');
		// Report XY in nm from the pixel centre
		sb.append(Utils.rounded(distanceFromCentre(benchmarkParameters.x))).append('\t');
		sb.append(Utils.rounded(distanceFromCentre(benchmarkParameters.y))).append('\t');
		sb.append(Utils.rounded(benchmarkParameters.a * benchmarkParameters.z)).append('\t');
		sb.append(Utils.rounded(benchmarkParameters.gain)).append('\t');
		sb.append(Utils.rounded(benchmarkParameters.readNoise)).append('\t');
		sb.append(Utils.rounded(benchmarkParameters.getBackground())).append('\t');
		sb.append(Utils.rounded(benchmarkParameters.noise)).append('\t');

		sb.append(Utils.rounded(benchmarkParameters.getSignal() / benchmarkParameters.noise)).append('\t');
		sb.append(Utils.rounded(benchmarkParameters.precisionN)).append('\t');
		sb.append(Utils.rounded(benchmarkParameters.precisionX)).append('\t');
		sb.append(Utils.rounded(benchmarkParameters.precisionXML)).append('\t');
		sb.append(region.width).append("x");
		sb.append(region.height).append('\t');
		sb.append(Utils.rounded(fitConfig.getInitialPeakStdDev() * benchmarkParameters.a)).append('\t');
		sb.append(PSFProtosHelper.getName(fitConfig.getPSF().getPsfType()));
		if (fitConfig.isFixedPSF())
		{
			// Only fixed fitting can ignore the signal
			if (!signalFitting)
				sb.append("NS");
		}
		if (!backgroundFitting)
			sb.append("NB");
		sb.append(":").append(PeakFit.getSolverName(fitConfig));
		if (fitConfig.isModelCameraMLE())
		{
			sb.append(":Camera\t");

			// Add details of the noise model for the MLE
			CalibrationReader r = new CalibrationReader(fitConfig.getCalibration());
			sb.append("EM=").append(r.isEMCCD());
			sb.append(":G=").append(r.getCountPerPhoton());
			sb.append(":N=").append(r.getReadNoise());
		}
		else
			sb.append('\t');

		// Convert to units of the image (ADUs and pixels)		
		double[] convert = getConversionFactors();

		// Store the results for fitting on this benchmark dataset
		BenchmarkResult benchmarkResult = new BenchmarkResult(benchmarkParameters, answer, sb.toString(), convert,
				this.results, this.resultsTime);
		if (!benchmarkResults.isEmpty())
		{
			// Clear the results if the benchmark has changed
			if (benchmarkResults.getFirst().benchmarkParameters.id != benchmarkParameters.id)
				benchmarkResults.clear();
		}
		benchmarkResults.add(benchmarkResult);

		// Now output the actual results ...		
		sb.append('\t');
		final double recall = (stats[0].getN() / (double) getNumberOfStartPoints()) /
				benchmarkParameters.getMolecules();
		sb.append(Utils.rounded(recall));

		for (int i = 0; i < stats.length; i++)
		{
			if (convert[i] != 0)
				sb.append('\t').append(Utils.rounded(stats[i].getMean() * convert[i], 6)).append('\t')
						.append(Utils.rounded(stats[i].getStandardDeviation() * convert[i]));
			else
				sb.append("\t0\t0");
		}
		summaryTable.append(sb.toString());
	}

	/**
	 * Get the factors to convert the fitted units into calibrated photons and nm units. Set
	 * the conversion to zero if the function does not fit the specified statistic.
	 * 
	 * @return The conversion factors
	 */
	private double[] getConversionFactors()
	{
		final double[] convert = new double[NAMES.length];
		convert[Gaussian2DFunction.BACKGROUND] = (fitConfig.isBackgroundFitting()) ? 1 : 0;
		convert[Gaussian2DFunction.SIGNAL] = (fitConfig.isNotSignalFitting() && fitConfig.isFixedPSF()) ? 0 : 1;
		convert[Gaussian2DFunction.ANGLE] = (fitConfig.isAngleFitting()) ? 180.0 / Math.PI : 0;
		convert[Gaussian2DFunction.X_POSITION] = benchmarkParameters.a;
		convert[Gaussian2DFunction.Y_POSITION] = benchmarkParameters.a;
		convert[Gaussian2DFunction.X_SD] = (fitConfig.isXSDFitting()) ? benchmarkParameters.a : 0;
		convert[Gaussian2DFunction.Y_SD] = (fitConfig.isYSDFitting()) ? benchmarkParameters.a : 0;
		convert[TIME] = 1e-6;
		convert[ACTUAL_SIGNAL] = convert[Gaussian2DFunction.SIGNAL];
		convert[ADJUSTED_X_SD] = convert[Gaussian2DFunction.X_SD];
		convert[ADJUSTED_Y_SD] = convert[Gaussian2DFunction.Y_SD];
		return convert;
	}

	private double distanceFromCentre(double x)
	{
		x -= 0.5;
		final int i = (int) Math.round(x);
		x = x - i;
		return x * benchmarkParameters.a;
	}

	private void createTable()
	{
		if (summaryTable == null || !summaryTable.isVisible())
		{
			summaryTable = new TextWindow(TITLE, createHeader(false), "", 1000, 300);
			summaryTable.setVisible(true);
		}
	}

	private void createAnalysisTable()
	{
		if (analysisTable == null || !analysisTable.isVisible())
		{
			analysisTable = new TextWindow(TITLE + " Combined Analysis", createHeader(true), "", 1000, 300);
			analysisTable.setVisible(true);
		}
	}

	private String createHeader(boolean extraRecall)
	{
		StringBuilder sb = new StringBuilder(createParameterHeader() + "\tRecall");
		if (extraRecall)
			sb.append("\tOrigRecall");
		for (int i = 0; i < NAMES.length; i++)
		{
			sb.append('\t').append(NAMES[i]).append("\t+/-");
		}
		return sb.toString();
	}

	private String createParameterHeader()
	{
		return "Molecules\tN\ts (nm)\ta (nm)\tsa (nm)\tX (nm)\tY (nm)\tZ (nm)\tGain\tReadNoise (ADUs)\tB (photons)\tNoise (photons)\tSNR\tLimit N\tLimit X\tLimit X ML\tRegion\tWidth\tMethod\tOptions";
	}

	private void runAnalysis()
	{
		benchmarkParameters = benchmarkResults.getFirst().benchmarkParameters;
		final double sa = getSa();

		// The fitting could have used centre-of-mass or not making the number of points different.
		// Find the shortest array (this will be the one where the centre-of-mass was not used)
		int length = Integer.MAX_VALUE;
		for (BenchmarkResult benchmarkResult : benchmarkResults)
			if (length > benchmarkResult.results.length)
				length = benchmarkResult.results.length;

		// Build a list of all the frames which have results
		int[] valid = new int[length];
		int j = 0;
		int[] count = new int[benchmarkResults.size()];
		for (BenchmarkResult benchmarkResult : benchmarkResults)
		{
			int c = 0;
			for (int i = 0; i < valid.length; i++)
				if (benchmarkResult.results[i] != null)
				{
					c++;
					valid[i]++;
				}
			count[j++] = c;
		}

		final int target = benchmarkResults.size();

		// Check that we have data
		if (!validData(valid, target))
		{
			IJ.error(TITLE, "No frames have fitting results from all methods");
			return;
		}

		// Get the number of start points valid for all the results 
		final int totalFrames = benchmarkParameters.frames;
		final double numberOfStartPoints = (double) (length / totalFrames);

		createAnalysisTable();

		// Create the results using only frames where all the fitting methods were successful
		j = 0;
		for (BenchmarkResult benchmarkResult : benchmarkResults)
		{
			final double[] answer = benchmarkResult.answer;

			Statistics[] stats = new Statistics[NAMES.length];
			for (int i = 0; i < stats.length; i++)
				stats[i] = new Statistics();

			for (int i = 0; i < valid.length; i++)
			{
				if (valid[i] < target)
					continue;

				addResult(stats, answer, benchmarkParameters.p[i % totalFrames], sa, benchmarkResult.results[i],
						benchmarkResult.resultsTime[i]);
			}

			StringBuilder sb = new StringBuilder(benchmarkResult.parameters);

			// Now output the actual results ...		
			sb.append('\t');
			final double recall = (stats[0].getN() / numberOfStartPoints) / benchmarkParameters.getMolecules();
			sb.append(Utils.rounded(recall));
			// Add the original recall
			sb.append('\t');
			final double recall2 = (count[j++] / numberOfStartPoints) / benchmarkParameters.getMolecules();
			sb.append(Utils.rounded(recall2));

			// Convert to units of the image (ADUs and pixels)		
			final double[] convert = benchmarkResult.convert;

			for (int i = 0; i < stats.length; i++)
			{
				if (convert[i] != 0)
					sb.append('\t').append(Utils.rounded(stats[i].getMean() * convert[i], 6)).append('\t')
							.append(Utils.rounded(stats[i].getStandardDeviation() * convert[i]));
				else
					sb.append("\t0\t0");
			}
			analysisTable.append(sb.toString());
		}
	}

	private boolean validData(int[] valid, int target)
	{
		for (int i = 0; i < valid.length; i++)
			if (valid[i] == target)
				return true;
		return false;
	}
}
