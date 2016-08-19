package gdsc.smlm.ij.plugins;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import gdsc.core.ij.Utils;
import gdsc.core.match.ClassificationResult;
import gdsc.core.match.Coordinate;
import gdsc.core.match.FractionClassificationResult;
import gdsc.core.threshold.AutoThreshold;
import gdsc.core.threshold.FloatHistogram;
import gdsc.core.threshold.Histogram;
import gdsc.core.utils.ImageExtractor;
import gdsc.core.utils.NoiseEstimator.Method;
import gdsc.core.utils.Statistics;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.engine.FitWorker;
import gdsc.smlm.filters.Spot;
import gdsc.smlm.fitting.FitConfiguration;
import gdsc.smlm.fitting.Gaussian2DFitter;
import gdsc.smlm.ij.plugins.BenchmarkSpotFilter.FilterResult;
import gdsc.smlm.ij.plugins.BenchmarkSpotFilter.ScoredSpot;
import gdsc.smlm.ij.plugins.ResultsMatchCalculator.PeakResultPoint;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.ImageConverter;
import gdsc.smlm.results.MemoryPeakResults;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.plugin.PlugIn;
import ij.text.TextWindow;

/**
 * Attempt to classify the spot candidates into those that do match a result and those that do not. This smart ranking
 * can be used to ensure all the good candidates are processed per frame before the fail limit takes effect.
 */
public class BenchmarkSmartSpotRanking implements PlugIn
{
	private static final String TITLE = "Smart Spot Ranking";

	static FitConfiguration fitConfig;
	private static FitEngineConfiguration config;
	static
	{
		fitConfig = new FitConfiguration();
		config = new FitEngineConfiguration(fitConfig);
		config.setNoiseMethod(Method.QUICK_RESIDUALS_LEAST_MEAN_OF_SQUARES);
	}

	private static AutoThreshold.Method[] thresholdMethods;
	private static double[] snrLevels;
	private static boolean[] thresholdMethodOptions;
	private static String[] thresholdMethodNames;
	static
	{
		snrLevels = new double[100];
		int i = 0;
		for (int snr = 20; snr <= 70; snr += 5)
			snrLevels[i++] = snr;
		snrLevels = Arrays.copyOf(snrLevels, i);

		thresholdMethods = AutoThreshold.Method.values();
		thresholdMethodOptions = new boolean[thresholdMethods.length + snrLevels.length];
		thresholdMethodNames = new String[thresholdMethodOptions.length];
		thresholdMethodOptions[AutoThreshold.Method.NONE.ordinal()] = true;
		for (i = 0; i < thresholdMethods.length; i++)
		{
			thresholdMethodNames[i] = thresholdMethods[i].name;
			thresholdMethodOptions[i] = true;
		}
		for (int j = 0; i < thresholdMethodNames.length; i++, j++)
		{
			thresholdMethodNames[i] = "SNR" + snrLevels[j];
			thresholdMethodOptions[i] = true;
		}

		// Turn some off
		// These often fail to converge
		thresholdMethodOptions[AutoThreshold.Method.INTERMODES.ordinal()] = false;
		thresholdMethodOptions[AutoThreshold.Method.MINIMUM.ordinal()] = false;
		// These are slow
		thresholdMethodOptions[AutoThreshold.Method.SHANBHAG.ordinal()] = false;
		thresholdMethodOptions[AutoThreshold.Method.RENYI_ENTROPY.ordinal()] = false;
	}
	private AutoThreshold.Method[] methods = null;
	private double[] levels = null;
	private String[] methodNames = null;

	private static double fractionPositives = 100;
	private static double fractionNegativesAfterAllPositives = 50;
	private static int negativesAfterAllPositives = 10;
	private static boolean selectMethods = true;
	private static int compactBins = 1024;
	private static String[] SORT = new String[] { "(None)", "tp", "fp", "tn", "fn", "Precision", "Recall", "F0.5", "F1",
			"F2", "Jaccard", "MCC" };
	private static int sortIndex = SORT.length - 3; // F2 to favour recall
	private static boolean useFractionScores = true;
	private static boolean showOverlay = false;

	private boolean extraOptions = false;

	private static TextWindow summaryTable = null;

	private MemoryPeakResults results;
	private CreateData.SimulationParameters simulationParameters;

	private static HashMap<Integer, ArrayList<Coordinate>> actualCoordinates = null;
	private static HashMap<Integer, FilterCandidates> filterCandidates;
	private static double fP, fN;
	private static int nP, nN;

	static int lastId = -1, lastFilterId = -1;
	private static double lastFractionPositives = -1;
	private static double lastFractionNegativesAfterAllPositives = -1;
	private static int lastNegativesAfterAllPositives = -1;

	private ImagePlus imp;

	// Allow other plugins to access the results
	static int rankResultsId = 0;
	static HashMap<Integer, RankResults> rankResults;
	static double distanceInPixels;
	static double lowerDistanceInPixels;
	static double candidateTN, candidateFN;

	private class FilterCandidates
	{
		// Integer counts of positives (matches) and negatives
		final int p, n;
		// Double sums of the fractions match score and antiscore 
		final double np, nn;
		final ScoredSpot[] spots;

		public FilterCandidates(int p, int n, double np, double nn, ScoredSpot[] spots)
		{
			this.p = p;
			this.n = n;
			this.np = np;
			this.nn = nn;
			this.spots = spots;
		}
	}

	private static final byte TP = (byte) 1;
	private static final byte FP = (byte) 2;
	private static final byte TN = (byte) 3;
	private static final byte FN = (byte) 4;

	private class RankResult
	{
		final float t;
		final FractionClassificationResult f;
		final ClassificationResult c;
		/**
		 * Store details about the spots that were accepted
		 */
		final byte[] good;
		final long time;

		public RankResult(float t, FractionClassificationResult f, ClassificationResult c, byte[] good, long time)
		{
			this.t = t;
			this.f = f;
			this.c = c;
			this.good = good;
			this.time = time;
		}
	}

	private class RankResults
	{
		final ScoredSpot[] spots;
		/** Store the z-position of the actual spots for later analysis. Size is the number of actual spots */
		@SuppressWarnings("unused")
		final double[] zPosition;

		ArrayList<RankResult> results = new ArrayList<RankResult>();

		public RankResults(ScoredSpot[] spots, double[] zPosition)
		{
			this.spots = spots;
			this.zPosition = zPosition;
		}

	}

	/**
	 * Used to allow multi-threading of the fitting method
	 */
	private class Worker implements Runnable
	{
		volatile boolean finished = false;
		final BlockingQueue<Integer> jobs;
		final ImageStack stack;
		final HashMap<Integer, ArrayList<Coordinate>> actualCoordinates;
		final HashMap<Integer, FilterCandidates> filterCandidates;
		final HashMap<Integer, RankResults> results;
		final int fitting;
		final boolean requireSNR;

		float[] data = null;
		double[] region = null;

		public Worker(BlockingQueue<Integer> jobs, ImageStack stack,
				HashMap<Integer, ArrayList<Coordinate>> actualCoordinates,
				HashMap<Integer, FilterCandidates> filterCandidates)
		{
			this.jobs = jobs;
			this.stack = stack;
			this.actualCoordinates = actualCoordinates;
			this.filterCandidates = filterCandidates;
			this.results = new HashMap<Integer, RankResults>();
			fitting = config.getRelativeFitting();
			requireSNR = true;
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
				while (!finished)
				{
					Integer job = jobs.take();
					if (job == null || job.intValue() < 0 || finished)
						break;
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

			FilterCandidates candidates = filterCandidates.get(frame);

			// Extract the spot intensities
			final ScoredSpot[] spots = candidates.spots;
			final float[] intensity = new float[spots.length];
			for (int i = 0; i < spots.length; i++)
			{
				final Spot spot = spots[i].spot;
				intensity[i] = spot.intensity;
			}

			// Estimate SNR
			double[] snr = null;
			long tSnr = 0;
			if (requireSNR)
			{
				tSnr = System.nanoTime();
				data = ImageConverter.getData(stack.getPixels(frame), stack.getWidth(), stack.getHeight(), null, data);
				final int maxx = stack.getWidth();
				final int maxy = stack.getHeight();
				final ImageExtractor ie = new ImageExtractor(data, maxx, maxy);
				final float noise = FitWorker.estimateNoise(data, maxx, maxy, config.getNoiseMethod());
				snr = new double[spots.length];
				for (int i = 0; i < spots.length; i++)
				{
					final Spot spot = spots[i].spot;
					final Rectangle regionBounds = ie.getBoxRegionBounds(spot.x, spot.y, fitting);
					region = ie.crop(regionBounds, region);
					double sum = 0;
					final int width = regionBounds.width;
					final int height = regionBounds.height;
					final int size = width * height;
					for (int k = size; k-- > 0;)
						sum += region[k];
					// TODO - The number of peaks could use the other candidates in the fit region
					final double b = Gaussian2DFitter.getBackground(region, width, height, 1);
					final double signal = sum - b * size;
					snr[i] = signal / noise;
				}
				tSnr = System.nanoTime() - tSnr;
			}

			Coordinate[] actual = ResultsMatchCalculator.getCoordinates(actualCoordinates, frame);
			final double[] zPosition = new double[actual.length];
			for (int i = 0; i < actual.length; i++)
			{
				PeakResultPoint p = (PeakResultPoint) actual[i];
				zPosition[i] = p.peakResult.error;
			}

			RankResults results = new RankResults(spots, zPosition);
			this.results.put(frame, results);

			long t1 = System.nanoTime();
			FloatHistogram histogram = FloatHistogram.buildHistogram(intensity.clone(), true);
			// Only compact once
			Histogram histogram2 = histogram.compact(compactBins);
			t1 = System.nanoTime() - t1;

			for (AutoThreshold.Method m : methods)
			{
				long t2 = System.nanoTime();
				float t = histogram2.getAutoThreshold(m);
				//// Lower to the next bin down to compensate for histogram compaction
				//if (histogram2 instanceof FloatHistogram)
				//{
				//	FloatHistogram h = (FloatHistogram) histogram2;
				//	for (int i = h.maxBin; i-- > 0;)
				//		if (h.value[i] < t)
				//		{
				//			t = h.value[i];
				//			break;
				//		}
				//}
				t2 = System.nanoTime() - t2;

				// Score
				final byte[] category = new byte[spots.length];
				double tp = 0;
				double fp = 0;
				double tn = 0;
				double fn = 0;
				int itp = 0;
				int ifp = 0;
				int itn = 0;
				int ifn = 0;
				for (int i = 0; i < spots.length; i++)
				{
					if (intensity[i] >= t)
					{
						tp += spots[i].getScore();
						fp += spots[i].antiScore();
						if (spots[i].match)
						{
							category[i] = TP;
							itp++;
						}
						else
						{
							category[i] = FP;
							ifp++;
						}
					}
					else
					{
						fn += spots[i].getScore();
						tn += spots[i].antiScore();
						if (spots[i].match)
						{
							category[i] = FN;
							ifn++;
						}
						else
						{
							category[i] = TN;
							itn++;
						}
					}
				}

				// Store the results using a copy of the original (to preserve the candidates for repeat analysis)
				results.results.add(new RankResult(t, new FractionClassificationResult(tp, fp, tn, fn),
						new ClassificationResult(itp, ifp, itn, ifn), category, t1 + t2));
			}

			for (double l : levels)
			{
				// Get the intensity of the lowest spot
				double t = Double.POSITIVE_INFINITY;

				// Score
				final byte[] category = new byte[spots.length];
				double tp = 0;
				double fp = 0;
				double tn = 0;
				double fn = 0;
				int itp = 0;
				int ifp = 0;
				int itn = 0;
				int ifn = 0;
				for (int i = 0; i < spots.length; i++)
				{
					if (snr[i] >= l)
					{
						tp += spots[i].getScore();
						fp += spots[i].antiScore();
						if (spots[i].match)
						{
							category[i] = TP;
							itp++;
						}
						else
						{
							category[i] = FP;
							ifp++;
						}
						if (t > intensity[i])
							t = intensity[i];
					}
					else
					{
						fn += spots[i].getScore();
						tn += spots[i].antiScore();
						if (spots[i].match)
						{
							category[i] = FN;
							ifn++;
						}
						else
						{
							category[i] = TN;
							itn++;
						}
					}
				}

				// Store the results using a copy of the original (to preserve the candidates for repeat analysis)
				results.results.add(new RankResult((float) t, new FractionClassificationResult(tp, fp, tn, fn),
						new ClassificationResult(itp, ifp, itn, ifn), category, tSnr));
			}

			//// Testing: Correlation between intensity and SNR
			//FastCorrelator c = new FastCorrelator();
			//double[] i2 = new double[spots.length];
			//for (int i = 0; i < spots.length; i++)
			//{
			//	i2[i] = intensity[i];
			//	c.add(Math.round(intensity[i]), Math.round(snr[i]));
			//}
			//System.out.printf("%d = %.2f (%d)\n", frame, c.getCorrelation(), c.getN());
			////Utils.display("I vs SNR", new Plot("I vs SNR", "Intensity", "SNR", i2, snr));
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		extraOptions = Utils.isExtraOptions();

		simulationParameters = CreateData.simulationParameters;
		if (simulationParameters == null)
		{
			IJ.error(TITLE, "No benchmark spot parameters in memory");
			return;
		}
		imp = CreateData.getImage();
		if (imp == null)
		{
			IJ.error(TITLE, "No simulation image");
			return;
		}
		results = MemoryPeakResults.getResults(CreateData.CREATE_DATA_IMAGE_TITLE + " (Create Data)");
		if (results == null)
		{
			IJ.error(TITLE, "No benchmark results in memory");
			return;
		}
		if (BenchmarkSpotFilter.filterResults == null)
		{
			IJ.error(TITLE, "No benchmark spot candidates in memory");
			return;
		}
		if (BenchmarkSpotFilter.simulationId != simulationParameters.id)
		{
			IJ.error(TITLE, "Update the benchmark spot candidates for the latest simulation");
			return;
		}

		if (!showDialog())
			return;

		run();
	}

	private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage(String.format(
				"Rank candidate spots in the benchmark image created by " + CreateData.TITLE +
						" plugin\nand identified by the " + BenchmarkSpotFilter.TITLE +
						" plugin.\nPSF width = %s nm (Square pixel adjustment = %s nm)\n \nConfigure the fitting:",
				Utils.rounded(simulationParameters.s), Utils.rounded(getSa())));

		gd.addSlider("Fraction_positives", 50, 100, fractionPositives);
		gd.addSlider("Fraction_negatives_after_positives", 0, 100, fractionNegativesAfterAllPositives);
		gd.addSlider("Min_negatives_after_positives", 0, 10, negativesAfterAllPositives);
		gd.addCheckbox("Select_methods", selectMethods);
		gd.addNumericField("Compact_bins", compactBins, 0);
		gd.addChoice("Sort", SORT, SORT[sortIndex]);
		gd.addCheckbox("Use_fraction_scores", useFractionScores);

		// Collect options for fitting that may effect ranking
		final double sa = getSa();
		gd.addNumericField("Initial_StdDev", sa / simulationParameters.a, 3);
		gd.addSlider("Fitting_width", 2, 4.5, config.getFitting());
		//		gd.addCheckbox("Include_neighbours", config.isIncludeNeighbours());
		//		gd.addSlider("Neighbour_height", 0.01, 1, config.getNeighbourHeightThreshold());

		// Output options
		gd.addCheckbox("Show_overlay", showOverlay);

		if (extraOptions)
		{
			String[] noiseMethodNames = SettingsManager.getNames((Object[]) Method.values());
			gd.addChoice("Noise_method", noiseMethodNames, noiseMethodNames[config.getNoiseMethod().ordinal()]);
		}

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		fractionPositives = Math.abs(gd.getNextNumber());
		fractionNegativesAfterAllPositives = Math.abs(gd.getNextNumber());
		negativesAfterAllPositives = (int) Math.abs(gd.getNextNumber());
		selectMethods = gd.getNextBoolean();
		compactBins = (int) Math.abs(gd.getNextNumber());
		sortIndex = gd.getNextChoiceIndex();
		useFractionScores = gd.getNextBoolean();

		// Collect options for fitting that may effect ranking
		fitConfig.setInitialPeakStdDev(gd.getNextNumber());
		config.setFitting(gd.getNextNumber());
		//		config.setIncludeNeighbours(gd.getNextBoolean());
		//		config.setNeighbourHeightThreshold(gd.getNextNumber());

		showOverlay = gd.getNextBoolean();

		if (extraOptions)
		{
			config.setNoiseMethod(gd.getNextChoiceIndex());
		}

		if (gd.invalidNumber())
			return false;

		methodNames = thresholdMethodNames.clone();
		if (selectMethods)
		{
			int count = 0, count1 = 0, count2 = 0;
			methods = new AutoThreshold.Method[thresholdMethods.length];
			levels = new double[snrLevels.length];

			gd = new GenericDialog(TITLE);
			gd.addHelp(About.HELP_URL);
			for (int i = 0; i < thresholdMethodNames.length; i++)
				gd.addCheckbox(thresholdMethodNames[i], thresholdMethodOptions[i]);

			gd.showDialog();

			if (gd.wasCanceled())
				return false;

			for (int i = 0, j = 0; i < thresholdMethodNames.length; i++)
			{
				thresholdMethodOptions[i] = gd.getNextBoolean();
				if (thresholdMethodOptions[i])
				{
					methodNames[count++] = thresholdMethodNames[i];
					if (i < thresholdMethods.length)
						methods[count1++] = thresholdMethods[i];
					else
						levels[count2++] = snrLevels[j++];
				}
			}

			methodNames = Arrays.copyOf(methodNames, count);
			methods = Arrays.copyOf(methods, count1);
			levels = Arrays.copyOf(levels, count2);
		}
		else
		{
			// Do them all
			methods = thresholdMethods.clone();
			levels = snrLevels.clone();
		}

		if (methodNames.length == 0)
		{
			IJ.error(TITLE, "No methods selected");
			return false;
		}

		return true;
	}

	/** The total progress. */
	int progress, stepProgress, totalProgress;

	/**
	 * Show progress.
	 */
	private synchronized void showProgress()
	{
		if (++progress % stepProgress == 0)
		{
			if (Utils.showStatus("Frame: " + progress + " / " + totalProgress))
				IJ.showProgress(progress, totalProgress);
		}
	}

	private void run()
	{
		// Extract all the results in memory into a list per frame. This can be cached
		boolean refresh = false;
		if (lastId != simulationParameters.id)
		{
			// Do not get integer coordinates
			// The Coordinate objects will be PeakResultPoint objects that store the original PeakResult
			// from the MemoryPeakResults
			actualCoordinates = ResultsMatchCalculator.getCoordinates(results.getResults(), false);
			lastId = simulationParameters.id;
			refresh = true;
		}

		// Extract all the candidates into a list per frame. This can be cached if the settings have not changed
		if (refresh || lastFilterId != BenchmarkSpotFilter.filterResultsId ||
				lastFractionPositives != fractionPositives ||
				lastFractionNegativesAfterAllPositives != fractionNegativesAfterAllPositives ||
				lastNegativesAfterAllPositives != negativesAfterAllPositives)
		{
			filterCandidates = subsetFilterResults(BenchmarkSpotFilter.filterResults);

			lastFilterId = BenchmarkSpotFilter.filterResultsId;
			lastFractionPositives = fractionPositives;
			lastFractionNegativesAfterAllPositives = fractionNegativesAfterAllPositives;
			lastNegativesAfterAllPositives = negativesAfterAllPositives;
		}

		final ImageStack stack = imp.getImageStack();

		// Create a pool of workers
		final int nThreads = Prefs.getThreads();
		BlockingQueue<Integer> jobs = new ArrayBlockingQueue<Integer>(nThreads * 2);
		List<Worker> workers = new LinkedList<Worker>();
		List<Thread> threads = new LinkedList<Thread>();
		for (int i = 0; i < nThreads; i++)
		{
			Worker worker = new Worker(jobs, stack, actualCoordinates, filterCandidates);
			Thread t = new Thread(worker);
			workers.add(worker);
			threads.add(t);
			t.start();
		}

		// Process the frames
		totalProgress = filterCandidates.size();
		stepProgress = Utils.getProgressInterval(totalProgress);
		progress = 0;
		for (int i : filterCandidates.keySet())
		{
			put(jobs, i);
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

		IJ.showProgress(1);

		if (Utils.isInterrupted())
		{
			IJ.showStatus("Aborted");
			return;
		}

		IJ.showStatus("Collecting results ...");

		rankResultsId++;
		rankResults = new HashMap<Integer, RankResults>();
		for (Worker w : workers)
		{
			rankResults.putAll(w.results);
		}

		summariseResults(rankResults);

		IJ.showStatus("");
	}

	/**
	 * Extract all the filter candidates in order until the desired number of positives have been reached and the number
	 * of negatives matches the configured parameters.
	 * 
	 * @param filterResults
	 * @return The filter candidates
	 */
	private HashMap<Integer, FilterCandidates> subsetFilterResults(HashMap<Integer, FilterResult> filterResults)
	{
		// Convert fractions from percent 
		final double f1 = Math.min(1, fractionPositives / 100.0);
		final double f2 = fractionNegativesAfterAllPositives / 100.0;

		HashMap<Integer, FilterCandidates> subset = new HashMap<Integer, FilterCandidates>();
		fP = fN = 0;
		nP = nN = 0;
		for (Entry<Integer, FilterResult> result : filterResults.entrySet())
		{
			FilterResult r = result.getValue();

			// Determine the number of positives to find. This score may be fractional.
			fP += r.result.getTP();
			fN += r.result.getFP();

			// Q. Is r.result.getTP() not the same as the total of r.spots[i].match
			// A. Not if we used fractional scoring.

			for (ScoredSpot spot : r.spots)
			{
				if (spot.match)
					nP++;
				else
					nN++;
			}

			// Make the target use the fractional score
			final double targetP = r.result.getTP() * f1;

			// Count the number of positive & negatives
			int p = 0, n = 0;
			double np = 0, nn = 0;

			boolean reachedTarget = false;
			int nAfter = 0;

			int count = 0;
			for (ScoredSpot spot : r.spots)
			{
				count++;
				nn += spot.antiScore();
				if (spot.match)
				{
					np += spot.getScore();
					p++;
					if (!reachedTarget)
					{
						reachedTarget = np >= targetP;
					}
				}
				else
				{
					n++;
					if (reachedTarget)
					{
						nAfter++;
					}
				}

				if (reachedTarget)
				{
					// Check if we have reached both the limits
					if (nAfter >= negativesAfterAllPositives && (double) n / (n + p) >= f2)
						break;
				}
			}

			// Debug
			//System.out.printf("Frame %d : %.1f / (%.1f + %.1f). p=%d, n=%d, after=%d, f=%.1f\n", result.getKey().intValue(),
			//		r.result.getTP(), r.result.getTP(), r.result.getFP(), p, n,
			//		nAfter, (double) n / (n + p));

			subset.put(result.getKey(), new FilterCandidates(p, n, np, nn, Arrays.copyOf(r.spots, count)));
		}
		return subset;
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

	private class ScoredResult implements Comparable<ScoredResult>
	{
		int i;
		double score;
		String result;

		public ScoredResult(int i, double score, String result)
		{
			this.i = i;
			this.score = score;
			this.result = result;
		}

		public int compareTo(ScoredResult o)
		{
			if (score < o.score)
				return 1;
			if (score > o.score)
				return -1;
			return 0;
		}
	}

	private void summariseResults(HashMap<Integer, RankResults> rankResults)
	{
		createTable();

		// Summarise the ranking results. 
		StringBuilder sb = new StringBuilder(BenchmarkSpotFilter.resultPrefix);

		// nP and nN is the fractional score of the spot candidates 
		addCount(sb, nP + nN);
		addCount(sb, nP);
		addCount(sb, nN);
		addCount(sb, fP);
		addCount(sb, fN);

		double tp = 0;
		double fp = 0;
		int cTP = 0, cFP = 0;
		for (FilterCandidates result : filterCandidates.values())
		{
			tp += result.np;
			fp += result.nn;
			cTP += result.p;
			cFP += result.n;
		}

		//		// This should be the same
		//		double tp2 = 0;
		//		double fp2 = 0;
		//		int cTP2 = 0, cFP2 = 0;
		//		for (RankResults rr : rankResults.values())
		//		{
		//			for (ScoredSpot spot : rr.spots)
		//			{
		//				if (spot.match)
		//					cTP2++;
		//				else
		//					cFP2++;
		//				tp2 += spot.getScore();
		//				fp2 += spot.antiScore();
		//			}
		//		}
		//		if (tp != tp2 || fp != fp2 || cTP != cTP2 || cFP != cFP2)
		//			System.out.println("Error counting");

		// The fraction of positive and negative candidates that were included
		add(sb, (100.0 * cTP) / nP);
		add(sb, (100.0 * cFP) / nN);

		// Add counts of the the candidates
		add(sb, cTP + cFP);
		add(sb, cTP);
		add(sb, cFP);

		// Add fractional counts of the the candidates
		add(sb, tp + fp);
		add(sb, tp);
		add(sb, fp);

		String resultPrefix = sb.toString();

		// ---
		// TODO
		// Add good label to spot candidates and have the benchmark spot filter respect this before applying the fail count limit.

		// Correlation between intensity and SNR ...

		// --- 

		// Allow using the fitted results from benchmark spot fit. Will it make a difference if we fit the candidates (some will fail
		// if weak).
		// Can this be done by allowing the user to select the input (spot candidates or fitted positions)?

		// Perhaps I need to produce a precision estimate for all simulated spots and then only use those that achieve a certain 
		// precision, i.e. are reasonably in focus. Can this be done? Does the image PSF have a width estimate for the entire stack?

		// Perhaps I should filter, fit and then filter all spots using no fail count. These then become the spots to work with
		// for creating a smart fail count filter. 

		// ---

		// Pre-compute the results and have optional sort
		ArrayList<ScoredResult> list = new ArrayList<ScoredResult>(methodNames.length);

		for (int i = 0; i < methodNames.length; i++)
		{
			tp = 0;
			fp = 0;
			double tn = 0;
			int itp = 0;
			int ifp = 0;
			int itn = 0;
			Statistics s = new Statistics();
			long time = 0;

			for (RankResults rr : rankResults.values())
			{
				RankResult r = rr.results.get(i);
				// Some results will not have a threshold
				if (Float.isFinite(r.t))
					s.add(r.t);
				time += r.time;
				tp += r.f.getTP();
				fp += r.f.getFP();
				tn += r.f.getTN();
				itp += r.c.getTP();
				ifp += r.c.getFP();
				itn += r.c.getTN();
			}

			sb.setLength(0);
			sb.append(resultPrefix);
			add(sb, methodNames[i]);
			if (methodNames[i].startsWith("SNR"))
				sb.append("\t");
			else
				add(sb, compactBins);
			add(sb, s.getMean());
			add(sb, s.getStandardDeviation());
			add(sb, Utils.timeToString(time / 1e6));

			// TP are all accepted candidates that can be matched to a spot
			// FP are all accepted candidates that cannot be matched to a spot
			// TN are all accepted candidates that cannot be matched to a spot
			// FN = The number of missed spots

			// Raw counts of match or no-match
			FractionClassificationResult f1 = new FractionClassificationResult(itp, ifp, itn,
					simulationParameters.molecules - itp);
			double s1 = addScores(sb, f1);

			// Fractional scoring
			FractionClassificationResult f2 = new FractionClassificationResult(tp, fp, tn,
					simulationParameters.molecules - tp);
			double s2 = addScores(sb, f2);

			// Store for sorting
			list.add(new ScoredResult(i, (useFractionScores) ? s2 : s1, sb.toString()));
		}

		if (list.isEmpty())
			return;

		Collections.sort(list);

		for (ScoredResult r : list)
			summaryTable.append(r.result);

		if (showOverlay)
		{
			int bestMethod = list.get(0).i;
			Overlay o = new Overlay();
			for (Entry<Integer, RankResults> entry : rankResults.entrySet())
			{
				int frame = entry.getKey();
				//FilterCandidates candidates = filterCandidates.get(frame);
				RankResults rr = entry.getValue();
				RankResult r = rr.results.get(bestMethod);
				int[] x1 = new int[r.good.length];
				int[] y1 = new int[r.good.length];
				int c1 = 0;
				int[] x2 = new int[r.good.length];
				int[] y2 = new int[r.good.length];
				int c2 = 0;
				int[] x3 = new int[r.good.length];
				int[] y3 = new int[r.good.length];
				int c3 = 0;
				int[] x4 = new int[r.good.length];
				int[] y4 = new int[r.good.length];
				int c4 = 0;
				for (int i = 0; i < x1.length; i++)
				{
					if (r.good[i] == TP)
					{
						x1[c1] = rr.spots[i].spot.x;
						y1[c1] = rr.spots[i].spot.y;
						c1++;
					}
					else if (r.good[i] == FP)
					{
						x2[c2] = rr.spots[i].spot.x;
						y2[c2] = rr.spots[i].spot.y;
						c2++;
					}
					else if (r.good[i] == TN)
					{
						x3[c3] = rr.spots[i].spot.x;
						y3[c3] = rr.spots[i].spot.y;
						c3++;
					}
					else if (r.good[i] == FN)
					{
						x4[c4] = rr.spots[i].spot.x;
						y4[c4] = rr.spots[i].spot.y;
						c4++;
					}
				}
				addToOverlay(o, frame, x1, y1, c1, Color.green);
				addToOverlay(o, frame, x2, y2, c2, Color.red);
				//addToOverlay(o, frame, x3, y3, c3, new Color(153, 255, 153)); // light green
				addToOverlay(o, frame, x4, y4, c4, new Color(255, 153, 153)); // light red
			}

			imp.setOverlay(o);
		}
	}

	private void addToOverlay(Overlay o, int frame, int[] x, int[] y, int c, Color color)
	{
		PointRoi roi = new PointRoi(x, y, c);
		roi.setFillColor(color);
		roi.setStrokeColor(color);
		roi.setPosition(frame);
		roi.setHideLabels(true);
		o.add(roi);
	}

	private static void add(StringBuilder sb, String value)
	{
		sb.append("\t").append(value);
	}

	private static void add(StringBuilder sb, int value)
	{
		sb.append("\t").append(value);
	}

	private static void add(StringBuilder sb, double value)
	{
		add(sb, Utils.rounded(value));
	}

	private static void addCount(StringBuilder sb, double value)
	{
		// Check if the double holds an integer count
		if ((int) value == value)
		{
			sb.append("\t").append((int) value);
		}
		else
		{
			// Otherwise add the counts using at least 2 dp
			if (value > 100)
				sb.append("\t").append(IJ.d2s(value));
			else
				add(sb, Utils.rounded(value));
		}
	}

	private void createTable()
	{
		if (summaryTable == null || !summaryTable.isVisible())
		{
			summaryTable = new TextWindow(TITLE, createHeader(false), "", 1000, 300);
			summaryTable.setVisible(true);
		}
	}

	private String createHeader(boolean extraRecall)
	{
		StringBuilder sb = new StringBuilder(BenchmarkSpotFilter.tablePrefix);
		sb.append("\t");
		sb.append("Spots\t");
		sb.append("nP\t");
		sb.append("nN\t");
		sb.append("fP\t");
		sb.append("fN\t");

		sb.append("% nP\t");
		sb.append("% nN\t");

		sb.append("cTotal\t");
		sb.append("cTP\t");
		sb.append("cFP\t");

		sb.append("cfTotal\t");
		sb.append("cfTP\t");
		sb.append("cfFP\t");

		sb.append("Method\t");
		sb.append("Bins\t");
		sb.append("Threshold Av\t");
		sb.append("Threshold SD\t");
		sb.append("Time\t");

		addScoreColumns(sb, null);
		addScoreColumns(sb, "f ");

		return sb.toString();
	}

	private void addScoreColumns(StringBuilder sb, String prefix)
	{
		SORT = new String[] { "(None)", "tp", "fp", "tn", "fn", "Precision", "Recall", "F0.5", "F1", "F2", "Jaccard",
				"MCC" };
		for (int i = 1; i < SORT.length; i++)
			addScoreColumn(sb, prefix, SORT[i]);
	}

	private double addScores(StringBuilder sb, FractionClassificationResult m)
	{
		double[] scores = new double[SORT.length - 1];
		int i = 0;
		scores[i++] = m.getTP();
		scores[i++] = m.getFP();
		scores[i++] = m.getTN();
		scores[i++] = m.getFN();
		scores[i++] = m.getPrecision();
		scores[i++] = m.getRecall();
		scores[i++] = m.getFScore(0.5);
		scores[i++] = m.getF1Score();
		scores[i++] = m.getFScore(2);
		scores[i++] = m.getJaccard();
		scores[i++] = m.getMCC();
		for (double s : scores)
			add(sb, s);
		return (sortIndex != 0) ? scores[sortIndex - 1] : 0;
	}

	private void addScoreColumn(StringBuilder sb, String prefix, String name)
	{
		if (prefix != null)
			sb.append(prefix);
		sb.append(name);
		sb.append("\t");
	}

	private double getSa()
	{
		final double sa = PSFCalculator.squarePixelAdjustment(simulationParameters.s, simulationParameters.a);
		return sa;
	}
}
