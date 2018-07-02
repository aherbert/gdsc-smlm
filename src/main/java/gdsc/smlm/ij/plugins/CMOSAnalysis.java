/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 * 
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm.ij.plugins;

import java.io.File;
import java.io.FileFilter;
import java.util.Collections;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.JLabel;

import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import org.apache.commons.math3.stat.inference.TestUtils;
import org.apache.commons.math3.util.MathArrays;

import gdsc.core.data.IntegerType;
import gdsc.core.data.SIPrefix;
import gdsc.core.generics.CloseableBlockingQueue;
import gdsc.core.ij.Utils;
import gdsc.core.logging.TrackProgress;
import gdsc.core.math.ArrayMoment;
import gdsc.core.math.IntegerArrayMoment;
import gdsc.core.math.RollingArrayMoment;
import gdsc.core.math.SimpleArrayMoment;
import gdsc.core.utils.DoubleData;
import gdsc.core.utils.PseudoRandomGenerator;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.core.utils.Statistics;
import gdsc.core.utils.StoredData;
import gdsc.core.utils.TextUtils;
import gdsc.core.utils.TurboList;
import gdsc.smlm.ij.SeriesImageSource;
import gdsc.smlm.ij.settings.Constants;
import gdsc.smlm.model.camera.PerPixelCameraModel;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.ExtendedGenericDialog;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.ProgressBar;
import ij.io.FileSaver;
import ij.plugin.PlugIn;
import ij.plugin.WindowOrganiser;

/**
 * Analyse the per pixel offset, variance and gain from a sCMOS camera.
 * <p>
 * See Huang et al (2013) Video-rate nanoscopy using sCMOS camera–specific single-molecule localization algorithms.
 * Nature Methods 10, 653-658 (Supplementary Information).
 */
public class CMOSAnalysis implements PlugIn
{
	private class SimulationWorker implements Runnable
	{
		final RandomGenerator rg;
		final String out;
		final float[] pixelOffset, pixelVariance, pixelGain;
		final int from, to, blockSize, photons;

		public SimulationWorker(long seed, String out, ImageStack stack, int from, int to, int blockSize, int photons)
		{
			rg = new Well19937c(seed);
			pixelOffset = (float[]) stack.getPixels(1);
			pixelVariance = (float[]) stack.getPixels(2);
			pixelGain = (float[]) stack.getPixels(3);
			this.out = out;
			this.from = from;
			this.to = to;
			this.blockSize = blockSize;
			this.photons = photons;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		@Override
		public void run()
		{
			// Avoid the status bar talking to the current image
			WindowManager.setTempCurrentImage(null);

			// Convert variance to SD
			float[] pixelSD = new float[pixelVariance.length];
			for (int i = 0; i < pixelVariance.length; i++)
				pixelSD[i] = (float) Math.sqrt(pixelVariance[i]);

			int size = (int) Math.sqrt(pixelVariance.length);

			// Pre-compute a set of Poisson numbers since this is slow
			int[] poisson = null;
			if (photons != 0)
			{
				// For speed we can precompute a set of random numbers to reuse
				RandomGenerator rg = new PseudoRandomGenerator(pixelVariance.length, this.rg);
				final PoissonDistribution pd = new PoissonDistribution(rg, photons, PoissonDistribution.DEFAULT_EPSILON,
						PoissonDistribution.DEFAULT_MAX_ITERATIONS);
				poisson = new int[pixelVariance.length];
				for (int i = poisson.length; i-- > 0;)
					poisson[i] = pd.sample();
			}

			// Save image in blocks
			ImageStack stack = new ImageStack(size, size);
			int start = from;
			for (int i = from; i < to; i++)
			{
				showProgress();

				// Create image
				final short[] pixels = new short[pixelOffset.length];
				if (photons == 0)
				{
					for (int j = 0; j < pixelOffset.length; j++)
					{
						// Fixed offset per pixel plus a variance
						double p = pixelOffset[j] + rg.nextGaussian() * pixelSD[j];
						pixels[j] = clip16bit(p);
					}
				}
				else
				{
					for (int j = 0; j < pixelOffset.length; j++)
					{
						// Fixed offset per pixel plus a variance plus a 
						// fixed gain multiplied by a Poisson sample of the photons
						double p = pixelOffset[j] + rg.nextGaussian() * pixelSD[j] + (poisson[j] * pixelGain[j]);
						pixels[j] = clip16bit(p);
					}

					// Rotate Poisson numbers.
					// Shuffling what we have is faster than generating new values
					// and we should have enough.
					MathArrays.shuffle(poisson, rg);
				}

				// Save image
				stack.addSlice(null, pixels);
				if (stack.getSize() == blockSize)
				{
					save(stack, start);
					start = i + 1;
					stack = new ImageStack(size, size);
				}
			}
			// This should not happen if we control the to-from range correctly
			if (stack.getSize() != 0)
				save(stack, start);
		}

		final static short MIN_SHORT = 0;
		// Require cast since this is out-of-range for a signed short
		final static short MAX_SHORT = (short) 65335;

		/**
		 * Clip to the range for a 16-bit image.
		 *
		 * @param value
		 *            the value
		 * @return the clipped value
		 */
		private short clip16bit(double value)
		{
			int i = (int) Math.round(value);
			if (i < 0)
				return MIN_SHORT;
			if (i > 65335)
				return MAX_SHORT;
			return (short) i;
		}

		private void save(ImageStack stack, int start)
		{
			ImagePlus imp = new ImagePlus("", stack);
			String path = new File(out, String.format("image%06d.tif", start)).getPath();
			FileSaver fs = new FileSaver(imp);
			fs.saveAsTiffStack(path);
		}
	}

	private class SubDir implements Comparable<SubDir>
	{
		int exposureTime;
		File path;
		String name;

		SubDir(int exposureTime, File path, String name)
		{
			this.exposureTime = exposureTime;
			this.path = path;
			this.name = name;
		}

		@Override
		public int compareTo(SubDir o)
		{
			return Integer.compare(exposureTime, o.exposureTime);
		}
	}

	/**
	 * Used to allow multi-threading of the scoring the filters
	 */
	private class ImageWorker implements Runnable
	{
		volatile boolean finished = false;
		final BlockingQueue<Object> jobs;
		final ArrayMoment moment;
		int bitDepth = 0;

		public ImageWorker(BlockingQueue<Object> jobs, ArrayMoment moment)
		{
			this.jobs = jobs;
			this.moment = moment.newInstance();
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Runnable#run()
		 */
		@Override
		public void run()
		{
			try
			{
				while (true)
				{
					Object pixels = jobs.take();
					if (pixels == null)
						break;
					if (!finished)
						// Only run jobs when not finished. This allows the queue to be emptied.
						run(pixels);
				}
			}
			catch (InterruptedException e)
			{
				if (!finished)
				{
					// This is not expected
					System.out.println(e.toString());
					throw new RuntimeException(e);
				}
			}
			finally
			{
				//Utils.log("Finished");
				finished = true;
			}
		}

		private void run(Object pixels)
		{
			if (Utils.isInterrupted())
			{
				finished = true;
				return;
			}
			showProgress();
			if (bitDepth == 0)
				bitDepth = Utils.getBitDepth(pixels);
			// Most likely first
			if (bitDepth == 16)
				moment.addUnsigned((short[]) pixels);
			else if (bitDepth == 32)
				moment.add((float[]) pixels);
			else if (bitDepth == 8)
				moment.addUnsigned((byte[]) pixels);
			else
				throw new IllegalStateException("Unsupported bit depth");
		}
	}

	private static final String TITLE = "sCMOS Analysis";

	private static String directory = Prefs.get(Constants.sCMOSAnalysisDirectory, "");
	private static String modelDirectory = null;
	private static String modelName = null;
	private static boolean rollingAlgorithm = false;
	private static boolean reuseProcessedData = true;

	// The simulation can default roughly to the values displayed 
	// in the Huang sCMOS paper supplementary figure 1:

	// Offset = Approximately Normal or Poisson. We use Poisson
	// since that is an integer distribution which would be expected
	// for an offset & Poisson approaches the Gaussian at high mean.
	private static double offset = 100;

	// Variance = Exponential (equivalent to chi-squared with k=1, i.e. 
	// sum of the squares of 1 normal distribution). 
	// We want 99.9% @ 400 ADU based on supplementary figure 1.a/1.b 
	// cumul = 1 - e^-lx (with l = 1/mean)
	// => e^-lx = 1 - cumul
	// => -lx = log(1-0.999)
	// => l = -log(0.001) / 400  (since x==400)
	// => 1/l = 57.9
	private static double variance = 57.9; // SD = 7.6

	// Gain = Approximately Normal
	private static double gain = 2.2;
	private static double gainSD = 0.2;

	private static int size = 512;
	private static int frames = 512;

	private boolean extraOptions = false;

	private static int imagejNThreads = Prefs.getThreads();
	private static int lastNThreads = imagejNThreads;

	private int nThreads = 0;
	// The simulated offset, variance and gain
	private ImagePlus simulationImp;
	// The measured offset, variance and gain
	private ImageStack measuredStack;
	// The sub-directories containing the sCMOS images
	private TurboList<SubDir> subDirs;

	/**
	 * Gets the last N threads used in the input dialog.
	 *
	 * @return the last N threads
	 */
	private static int getLastNThreads()
	{
		// See if ImageJ preference were updated
		if (imagejNThreads != Prefs.getThreads())
		{
			lastNThreads = imagejNThreads = Prefs.getThreads();
		}
		// Otherwise use the last user input
		return lastNThreads;
	}

	/**
	 * Gets the threads to use for multi-threaded computation.
	 *
	 * @return the threads
	 */
	private int getThreads()
	{
		if (nThreads == 0)
		{
			nThreads = Prefs.getThreads();
		}
		return nThreads;
	}

	/**
	 * Sets the threads to use for multi-threaded computation.
	 *
	 * @param nThreads
	 *            the new threads
	 */
	private void setThreads(int nThreads)
	{
		this.nThreads = Math.max(1, nThreads);
		// Save user input
		lastNThreads = this.nThreads;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	@Override
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		extraOptions = Utils.isExtraOptions();
		// Avoid the status bar talking to the current image
		WindowManager.setTempCurrentImage(null);

		//@formatter:off
		IJ.log(TextUtils.wrap(
				TITLE + ": Analyse the per-pixel offset, variance and gain of sCMOS images. " + 
				"See Huang et al (2013) Video-rate nanoscopy using sCMOS camera–specific " +
				"single-molecule localization algorithms. Nature Methods 10, 653-658 " +
				"(Supplementary Information).",
				80));
		//@formatter:on

		String dir = Utils.getDirectory(TITLE, directory);
		if (TextUtils.isNullOrEmpty(dir))
			return;
		directory = dir;
		Prefs.set(Constants.sCMOSAnalysisDirectory, dir);

		boolean simulate = "simulate".equals(arg);
		if (simulate || extraOptions)
		{
			if (!showSimulateDialog())
				return;
			simulate();
		}

		if (!showDialog())
			return;

		run();

		if (simulationImp == null)
		{
			// Just in case an old simulation is in the directory
			ImagePlus imp = IJ.openImage(new File(directory, "perPixelSimulation.tif").getPath());
			if (imp != null && imp.getStackSize() == 3 && imp.getWidth() == measuredStack.getWidth() &&
					imp.getHeight() == measuredStack.getHeight())
				simulationImp = imp;
		}

		if (simulationImp != null)
			computeError();
	}

	/** The total progress. */
	int progress, stepProgress, totalProgress;
	ProgressBar progressBar;

	/**
	 * Show progress.
	 */
	private synchronized void showProgress()
	{
		//Utils.log("%d/%d\n", progress, totalProgress);
		if (progress % stepProgress == 0)
		{
			//IJ.showProgress(progress, totalProgress);

			// Use the actual progress bar so we can show progress 
			// when all other IJ commands cannot
			double p = (progress + 1.0) / totalProgress;
			progressBar.show(p, true);
		}
		progress++;
	}

	private void simulate()
	{
		// Create the offset, variance and gain for each pixel
		int n = size * size;
		float[] pixelOffset = new float[n];
		float[] pixelVariance = new float[n];
		float[] pixelGain = new float[n];

		IJ.showStatus("Creating random per-pixel readout");
		long start = System.currentTimeMillis();

		RandomGenerator rg = new Well19937c();
		PoissonDistribution pd = new PoissonDistribution(rg, offset, PoissonDistribution.DEFAULT_EPSILON,
				PoissonDistribution.DEFAULT_MAX_ITERATIONS);
		ExponentialDistribution ed = new ExponentialDistribution(rg, variance,
				ExponentialDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
		totalProgress = n;
		stepProgress = Utils.getProgressInterval(totalProgress);
		for (int i = 0; i < n; i++)
		{
			if (i % n == 0)
				IJ.showProgress(i, n);
			// Q. Should these be clipped to a sensible range?
			pixelOffset[i] = pd.sample();
			pixelVariance[i] = (float) ed.sample();
			pixelGain[i] = (float) (gain + rg.nextGaussian() * gainSD);
		}
		IJ.showProgress(1);

		// Avoid all the file saves from updating the progress bar and status line
		Utils.setShowStatus(false);
		Utils.setShowProgress(false);
		JLabel statusLine = Utils.getStatusLine();
		progressBar = Utils.getProgressBar();

		// Save to the directory as a stack
		ImageStack simulationStack = new ImageStack(size, size);
		simulationStack.addSlice("Offset", pixelOffset);
		simulationStack.addSlice("Variance", pixelVariance);
		simulationStack.addSlice("Gain", pixelGain);
		simulationImp = new ImagePlus("PerPixel", simulationStack);
		// Only the info property is saved to the TIFF file
		simulationImp.setProperty("Info", String.format("Offset=%s; Variance=%s; Gain=%s +/- %s", Utils.rounded(offset),
				Utils.rounded(variance), Utils.rounded(gain), Utils.rounded(gainSD)));
		IJ.save(simulationImp, new File(directory, "perPixelSimulation.tif").getPath());

		// Create thread pool and workers
		ExecutorService executor = Executors.newFixedThreadPool(getThreads());
		TurboList<Future<?>> futures = new TurboList<Future<?>>(nThreads);

		// Simulate the zero exposure input.
		// Simulate 20 - 200 photon images.
		int[] photons = new int[] { 0, 20, 50, 100, 200 };

		totalProgress = photons.length * frames;
		stepProgress = Utils.getProgressInterval(totalProgress);
		progress = 0;
		progressBar.show(0);

		int blockSize = 10; // For saving stacks
		int nPerThread = (int) Math.ceil((double) frames / nThreads);
		// Convert to fit the block size
		nPerThread = (int) Math.ceil((double) nPerThread / blockSize) * blockSize;
		long seed = start;

		for (int p : photons)
		{
			statusLine.setText("Simulating " + TextUtils.pleural(p, "photon"));

			// Create the directory
			File out = new File(directory, String.format("photon%03d", p));
			if (!out.exists())
				out.mkdir();

			for (int from = 0; from < frames;)
			{
				int to = Math.min(from + nPerThread, frames);
				futures.add(executor
						.submit(new SimulationWorker(seed++, out.getPath(), simulationStack, from, to, blockSize, p)));
				from = to;
			}
			// Wait for all to finish
			for (int t = futures.size(); t-- > 0;)
			{
				try
				{
					// The future .get() method will block until completed
					futures.get(t).get();
				}
				catch (Exception e)
				{
					// This should not happen. 
					e.printStackTrace();
				}
			}
			futures.clear();
		}

		Utils.setShowStatus(true);
		Utils.setShowProgress(true);
		IJ.showProgress(1);

		executor.shutdown();

		Utils.log("Simulation time = " + Utils.timeToString(System.currentTimeMillis() - start));
	}

	private boolean showSimulateDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Simulate per-pixel offset, variance and gain of sCMOS images.");

		gd.addNumericField("nThreads", getLastNThreads(), 0);
		gd.addNumericField("Offset (Poisson)", offset, 3);
		gd.addNumericField("Variance (Exponential)", variance, 3);
		gd.addNumericField("Gain (Gaussian)", gain, 3);
		gd.addNumericField("Gain_SD", gainSD, 3);
		gd.addNumericField("Size", size, 0);
		gd.addNumericField("Frames", frames, 0);
		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		setThreads((int) gd.getNextNumber());
		offset = Math.abs(gd.getNextNumber());
		variance = Math.abs(gd.getNextNumber());
		gain = Math.abs(gd.getNextNumber());
		gainSD = Math.abs(gd.getNextNumber());
		size = Math.abs((int) gd.getNextNumber());
		frames = Math.abs((int) gd.getNextNumber());

		// Check arguments
		try
		{
			Parameters.isAboveZero("Offset", offset);
			Parameters.isAboveZero("Variance", variance);
			Parameters.isAboveZero("Gain", gain);
			Parameters.isAboveZero("Gain SD", gainSD);
			Parameters.isAboveZero("Size", size);
			Parameters.isAboveZero("Frames", frames);
		}
		catch (IllegalArgumentException ex)
		{
			Utils.log(TITLE + ": " + ex.getMessage());
			return false;
		}

		return true;
	}

	private boolean showDialog()
	{
		// Determine sub-directories to process
		File dir = new File(directory);
		File[] dirs = dir.listFiles(new FileFilter()
		{
			@Override
			public boolean accept(File pathname)
			{
				return pathname.isDirectory();
			}
		});

		if (dirs.length == 0)
		{
			IJ.error(TITLE, "No sub-directories");
			return false;
		}

		// Get only those with numbers at the end. 
		// These should correspond to exposure times
		subDirs = new TurboList<SubDir>();
		Pattern p = Pattern.compile("([0-9]+)$");
		for (File path : dirs)
		{
			String name = path.getName();
			Matcher m = p.matcher(name);
			if (m.find())
			{
				int t = Integer.parseInt(m.group(1));
				subDirs.add(new SubDir(t, path, name));
			}
		}

		if (subDirs.size() < 2)
		{
			IJ.error(TITLE, "Not enough sub-directories with exposure time suffix");
			return false;
		}

		Collections.sort(subDirs);

		if (subDirs.get(0).exposureTime != 0)
		{
			IJ.error(TITLE, "No sub-directories with exposure time 0");
			return false;
		}

		for (SubDir sd : subDirs)
		{
			Utils.log("Sub-directory: %s. Exposure time = %d", sd.name, sd.exposureTime);
		}

		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		//@formatter:off
		gd.addMessage("Analyse the per-pixel offset, variance and gain of sCMOS images.\n \n" + 
				TextUtils.wrap(
				"See Huang et al (2013) Video-rate nanoscopy using sCMOS camera–specific " +
				"single-molecule localization algorithms. Nature Methods 10, 653-658 " +
				"(Supplementary Information).",
				80));
		//@formatter:on

		gd.addNumericField("nThreads", getLastNThreads(), 0);
		gd.addMessage(TextUtils.wrap(
				"A rolling algorithm can handle any size of data but is slower. Otherwise the camera is assumed to produce a maximum of 16-bit unsigned data.",
				80));
		gd.addCheckbox("Rolling_algorithm", rollingAlgorithm);
		gd.addCheckbox("Re-use_processed_data", reuseProcessedData);
		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		setThreads((int) gd.getNextNumber());
		rollingAlgorithm = gd.getNextBoolean();

		return true;
	}

	private void run()
	{
		long start = System.currentTimeMillis();

		// Avoid all the file saves from updating the progress bar and status line
		Utils.setShowProgress(false);
		Utils.setShowStatus(false);
		JLabel statusLine = Utils.getStatusLine();
		progressBar = Utils.getProgressBar();

		TrackProgress trackProgress = new TrackProgress()
		{
			@Override
			public void progress(double fraction)
			{
				progressBar.show(fraction);
			}

			@Override
			public void progress(long position, long total)
			{
				progressBar.show((double) position / total);
			}

			@Override
			public void incrementProgress(double fraction)
			{

			}

			@Override
			public void log(String format, Object... args)
			{

			}

			@Override
			public void status(String format, Object... args)
			{

			}

			@Override
			public boolean isEnded()
			{
				return false;
			}

			@Override
			public boolean isProgress()
			{
				return true;
			}

			@Override
			public boolean isLog()
			{
				return false;
			}

			@Override
			public boolean isStatus()
			{
				return false;
			}
		};

		// Create thread pool and workers. The system is likely to be IO limited 
		// so reduce the computation threads to allow the reading thread in the 
		// SeriesImageSource to run.
		// If the images are small enough to fit into memory then 3 threads are used, 
		// otherwise it is 1.
		int nThreads = Math.max(1, getThreads() - 3);
		ExecutorService executor = Executors.newFixedThreadPool(nThreads);
		TurboList<Future<?>> futures = new TurboList<Future<?>>(nThreads);
		TurboList<ImageWorker> workers = new TurboList<ImageWorker>(nThreads);

		double[][] data = new double[subDirs.size() * 2][];
		double[] pixelOffset = null, pixelVariance = null;
		Statistics statsOffset = null, statsVariance = null;

		// For each sub-directory compute the mean and variance
		final int nSubDirs = subDirs.size();
		boolean error = false;
		int width = 0, height = 0;
		for (int n = 0; n < nSubDirs; n++)
		{
			SubDir sd = subDirs.getf(n);
			statusLine.setText("Analysing " + sd.name);
			StopWatch sw = StopWatch.createStarted();

			// Option to reuse data
			File file = new File(directory, "perPixel" + sd.name + ".tif");
			boolean found = false;
			if (reuseProcessedData && file.exists())
			{
				ImagePlus imp = IJ.openImage(file.getPath());
				if (imp != null && imp.getStackSize() == 2 && imp.getBitDepth() == 32)
				{
					if (n == 0)
					{
						width = imp.getWidth();
						height = imp.getHeight();
					}
					else
					{
						if (width != imp.getWidth() || height != imp.getHeight())
						{
							error = true;
							IJ.error(TITLE,
									"Image width/height mismatch in image series: " + file.getPath() +
											String.format("\n \nExpected %dx%d, Found %dx%d", width, height,
													imp.getWidth(), imp.getHeight()));
							break;
						}
					}

					ImageStack stack = imp.getImageStack();
					data[2 * n] = SimpleArrayUtils.toDouble((float[]) stack.getPixels(1));
					data[2 * n + 1] = SimpleArrayUtils.toDouble((float[]) stack.getPixels(2));
					found = true;
				}
			}

			if (!found)
			{
				// Open the series
				SeriesImageSource source = new SeriesImageSource(sd.name, sd.path.getPath());
				source.setTrackProgress(trackProgress);
				if (!source.open())
				{
					error = true;
					IJ.error(TITLE, "Failed to open image series: " + sd.path.getPath());
					break;
				}

				if (n == 0)
				{
					width = source.getWidth();
					height = source.getHeight();
				}
				else
				{
					if (width != source.getWidth() || height != source.getHeight())
					{
						error = true;
						IJ.error(TITLE,
								"Image width/height mismatch in image series: " + sd.path.getPath() +
										String.format("\n \nExpected %dx%d, Found %dx%d", width, height,
												source.getWidth(), source.getHeight()));
						break;
					}
				}

				totalProgress = source.getFrames() + 1; // So the bar remains at 99% when workers have finished
				stepProgress = Utils.getProgressInterval(totalProgress);
				progress = 0;
				progressBar.show(0);

				// Open the first frame to get the bit depth.
				// Assume the first pixels are not empty as the source is open.
				Object pixels = source.nextRaw();
				int bitDepth = Utils.getBitDepth(pixels);

				ArrayMoment moment;
				if (rollingAlgorithm)
				{
					moment = new RollingArrayMoment();
				}
				else
				{
					// We assume 16-bit camera at the maximum
					if (bitDepth <= 16 && IntegerArrayMoment.isValid(IntegerType.UNSIGNED_16, source.getFrames()))
						moment = new IntegerArrayMoment();
					else
						moment = new SimpleArrayMoment();
				}

				final CloseableBlockingQueue<Object> jobs = new CloseableBlockingQueue<Object>(nThreads * 2);
				for (int i = 0; i < nThreads; i++)
				{
					final ImageWorker worker = new ImageWorker(jobs, moment);
					workers.add(worker);
					futures.add(executor.submit(worker));
				}

				// Process the raw pixel data
				long lastTime = 0;
				while (pixels != null)
				{
					long time = System.currentTimeMillis();
					if (time - lastTime > 150)
					{
						if (Utils.isInterrupted())
						{
							error = true;
							break;
						}
						lastTime = time;
						statusLine.setText("Analysing " + sd.name + " Frame " + source.getStartFrameNumber());
					}
					put(jobs, pixels);
					pixels = source.nextRaw();
				}
				source.close();

				if (error)
				{
					// Kill the workers
					jobs.close(true);
					for (int t = futures.size(); t-- > 0;)
					{
						try
						{
							workers.get(t).finished = true;
							futures.get(t).cancel(true);
						}
						catch (Exception e)
						{
							// This should not happen. 
							e.printStackTrace();
						}
					}
					break;
				}

				// Finish all the worker threads by passing in a null job
				jobs.close(false);

				// Wait for all to finish
				for (int t = futures.size(); t-- > 0;)
				{
					try
					{
						// The future .get() method will block until completed
						futures.get(t).get();
					}
					catch (Exception e)
					{
						// This should not happen. 
						e.printStackTrace();
					}
				}

				// Create the final aggregate statistics
				for (ImageWorker w : workers)
					moment.add(w.moment);
				data[2 * n] = moment.getFirstMoment();
				data[2 * n + 1] = moment.getVariance();

				// Get the processing speed.
				sw.stop();
				// progress holds the number of calls to showProgress() for 
				// processing a frame (i.e. number of frames)
				double bits = (double) bitDepth * progress * source.getWidth() * source.getHeight();
				double seconds = sw.getNanoTime() / 1e9;
				double bps = bits / seconds;
				SIPrefix prefix = SIPrefix.getPrefix(bps);
				Utils.log("Processed %d frames. Time = %s. Rate = %s %sbits/s", moment.getN(), sw.toString(),
						Utils.rounded(prefix.convert(bps)), prefix.getName());

				// Reset
				futures.clear();
				workers.clear();

				ImageStack stack = new ImageStack(width, height);
				stack.addSlice("Mean", SimpleArrayUtils.toFloat(data[2 * n]));
				stack.addSlice("Variance", SimpleArrayUtils.toFloat(data[2 * n + 1]));
				IJ.save(new ImagePlus("PerPixel", stack), file.getPath());
			}

			Statistics s = new Statistics(data[2 * n]);

			if (n != 0)
			{
				// Compute mean ADU
				Statistics signal = new Statistics();
				double[] mean = data[2 * n];
				for (int i = 0; i < pixelOffset.length; i++)
					signal.add(mean[i] - pixelOffset[i]);
				Utils.log("%s Mean = %s +/- %s. Signal = %s +/- %s ADU", sd.name, Utils.rounded(s.getMean()),
						Utils.rounded(s.getStandardDeviation()), Utils.rounded(signal.getMean()),
						Utils.rounded(signal.getStandardDeviation()));
			}
			else
			{
				pixelOffset = data[0];
				pixelVariance = data[1];
				statsOffset = s;
				statsVariance = new Statistics(pixelVariance);
				Utils.log("%s Offset = %s +/- %s. Variance = %s +/- %s", sd.name, Utils.rounded(s.getMean()),
						Utils.rounded(s.getStandardDeviation()), Utils.rounded(statsVariance.getMean()),
						Utils.rounded(statsVariance.getStandardDeviation()));
			}

			progressBar.show(1);
		}
		progressBar.show(1);

		Utils.setShowStatus(true);
		Utils.setShowProgress(true);
		IJ.showProgress(1);

		if (error)
		{
			executor.shutdownNow();
			statusLine.setText(TITLE + " cancelled");
			return;
		}

		executor.shutdown();

		// Compute the gain
		statusLine.setText("Computing gain");

		double[] pixelGain = new double[pixelOffset.length];
		double[] bibiT = new double[pixelGain.length];
		double[] biaiT = new double[pixelGain.length];

		// Ignore first as this is the 0 exposure image
		for (int n = 1; n < nSubDirs; n++)
		{
			// Use equation 2.5 from the Huang et al paper.
			double[] b = data[2 * n];
			double[] a = data[2 * n + 1];
			for (int i = 0; i < pixelGain.length; i++)
			{
				double bi = b[i] - pixelOffset[i];
				double ai = a[i] - pixelVariance[i];
				bibiT[i] += bi * bi;
				biaiT[i] += bi * ai;
			}
		}
		for (int i = 0; i < pixelGain.length; i++)
			pixelGain[i] = biaiT[i] / bibiT[i];

		//for (int i = 0; i < pixelGain.length; i++)
		//{
		//	// Use equation 2.5 from the Huang et al paper.
		//	double bibiT = 0;
		//	double biaiT = 0;
		//	// Ignore first as this is the 0 exposure image
		//	for (int n = 1; n < nSubDirs; n++)
		//	{
		//		double bi = data[2*n][i] - pixelOffset[i];
		//		double ai = data[2*n+1][i] - pixelVariance[i];
		//		bibiT += bi * bi;
		//		biaiT += bi * ai;
		//	}
		//	pixelGain[i] = biaiT / bibiT;
		//}

		Statistics statsGain = new Statistics(pixelGain);
		Utils.log("Gain Mean = %s +/- %s", Utils.rounded(statsGain.getMean()),
				Utils.rounded(statsGain.getStandardDeviation()));

		// Histogram of offset, variance and gain
		int bins = 2 * Utils.getBinsSturges(pixelGain.length);
		WindowOrganiser wo = new WindowOrganiser();
		showHistogram("Offset (ADU)", pixelOffset, bins, statsOffset, wo);
		showHistogram("Variance (ADU^2)", pixelVariance, bins, statsVariance, wo);
		showHistogram("Gain (ADU/e)", pixelGain, bins, statsGain, wo);
		wo.tile();

		// Save
		float[] bias = SimpleArrayUtils.toFloat(pixelOffset);
		float[] variance = SimpleArrayUtils.toFloat(pixelVariance);
		float[] gain = SimpleArrayUtils.toFloat(pixelGain);
		measuredStack = new ImageStack(width, height);
		measuredStack.addSlice("Offset", bias);
		measuredStack.addSlice("Variance", variance);
		measuredStack.addSlice("Gain", gain);

		ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
		egd.addMessage("Save the sCMOS camera model?");
		if (modelDirectory == null)
		{
			modelDirectory = directory;
			modelName = "sCMOS Camera";
		}
		egd.addStringField("Model_name", modelName, 30);
		egd.addDirectoryField("Model_directory", modelDirectory);
		egd.showDialog();
		if (!egd.wasCanceled())
		{
			modelName = egd.getNextString();
			modelDirectory = egd.getNextString();
			PerPixelCameraModel cameraModel = new PerPixelCameraModel(width, height, bias, gain, variance);
			if (!CameraModelManager.save(cameraModel, new File(directory, modelName).getPath()))
				IJ.error(TITLE, "Failed to save model to file");
		}
		IJ.showStatus(""); // Remove the status from the ij.io.ImageWriter class

		Utils.log("Analysis time = " + Utils.timeToString(System.currentTimeMillis() - start));
	}

	private void showHistogram(String name, double[] values, int bins, Statistics stats, WindowOrganiser wo)
	{
		DoubleData data = new StoredData(values, false);
		double minWidth = 0;
		int removeOutliers = 0;
		int shape = Plot.CIRCLE; // Plot2.BAR; // A bar chart confuses the log plot since it plots lines to zero.
		String label = String.format("Mean = %s +/- %s", Utils.rounded(stats.getMean()),
				Utils.rounded(stats.getStandardDeviation()));
		int id = Utils.showHistogram(TITLE, data, name, minWidth, removeOutliers, bins, shape, label);
		if (Utils.isNewWindow())
			wo.add(id);
		// Redraw using a log scale. This requires a non-zero y-min
		Plot plot = Utils.plot;
		double[] limits = plot.getLimits();
		plot.setLimits(limits[0], limits[1], 1, limits[3]);
		plot.setAxisYLog(true);
		Utils.plot.updateImage();
	}

	private <T> void put(BlockingQueue<T> jobs, T job)
	{
		try
		{
			jobs.put(job);
		}
		catch (InterruptedException e)
		{
			throw new RuntimeException("Unexpected interruption", e);
		}
	}

	private void computeError()
	{
		// Assume the simulation stack and measured stack are not null.
		Utils.log("Comparison to simulation: %s", simulationImp.getInfoProperty());
		ImageStack simulationStack = simulationImp.getImageStack();
		for (int slice = 1; slice <= 3; slice++)
		{
			computeError(slice, simulationStack);
		}
	}

	private void computeError(int slice, ImageStack simulationStack)
	{
		String label = simulationStack.getSliceLabel(slice);
		float[] e = (float[]) simulationStack.getPixels(slice);
		float[] o = (float[]) measuredStack.getPixels(slice);

		// Get the mean error
		Statistics s = new Statistics();
		for (int i = e.length; i-- > 0;)
			s.add(o[i] - e[i]);

		StringBuilder result = new StringBuilder("Error ").append(label);
		result.append(" = ").append(Utils.rounded(s.getMean()));
		result.append(" +/- ").append(Utils.rounded(s.getStandardDeviation()));

		// Do statistical tests
		double[] x = SimpleArrayUtils.toDouble(e), y = SimpleArrayUtils.toDouble(o);

		PearsonsCorrelation c = new PearsonsCorrelation();
		result.append(" : R=").append(Utils.rounded(c.correlation(x, y)));

		// Mann-Whitney U is valid for any distribution, e.g. variance
		MannWhitneyUTest test = new MannWhitneyUTest();
		double p = test.mannWhitneyUTest(x, y);
		result.append(" : Mann-Whitney U p=").append(Utils.rounded(p)).append(' ')
				.append(((p < 0.05) ? "reject" : "accept"));

		if (slice != 2)
		{
			// T-Test is valid for approximately Normal distributions, e.g. offset and gain
			p = TestUtils.tTest(x, y);
			result.append(" : T-Test p=").append(Utils.rounded(p)).append(' ')
					.append(((p < 0.05) ? "reject" : "accept"));
			p = TestUtils.pairedTTest(x, y);
			result.append(" : Paired T-Test p=").append(Utils.rounded(p)).append(' ')
					.append(((p < 0.05) ? "reject" : "accept"));
		}

		Utils.log(result.toString());
	}
}
