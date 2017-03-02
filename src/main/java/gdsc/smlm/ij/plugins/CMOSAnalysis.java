package gdsc.smlm.ij.plugins;

import java.io.File;
import java.io.FileFilter;
import java.util.Collections;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.JLabel;

import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

import gdsc.core.ij.Utils;
import gdsc.core.math.ArrayMoment;
import gdsc.core.math.RollingArrayMoment;
import gdsc.core.math.SimpleArrayMoment;
import gdsc.core.utils.PseudoRandomGenerator;
import gdsc.core.utils.Statistics;
import gdsc.core.utils.TextUtils;
import gdsc.core.utils.TurboList;
import gdsc.smlm.ij.SeriesImageSource;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.ProgressBar;
import ij.io.FileSaver;
import ij.plugin.PlugIn;
import ij.process.ShortProcessor;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

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
		public void run()
		{
			// Avoid the status bar talking to the current image
			WindowManager.setTempCurrentImage(null);

			// Convert variance to SD
			float[] pixelSD = pixelVariance.clone();
			for (int i = 0; i < pixelVariance.length; i++)
				pixelSD[i] = (float) Math.sqrt(pixelVariance[i]);

			int size = (int) Math.sqrt(pixelVariance.length);

			// For speed we can precompute a set of random numbers to reuse
			RandomGenerator rg = new PseudoRandomGenerator(pixelVariance.length, this.rg);

			// Pre-compute a set of Poisson numbers since this is slow
			int nPoisson = pixelVariance.length + 1;
			int[] poisson = null;
			if (photons != 0)
			{
				final PoissonDistribution pd = new PoissonDistribution(rg, photons, PoissonDistribution.DEFAULT_EPSILON,
						PoissonDistribution.DEFAULT_MAX_ITERATIONS);
				poisson = new int[nPoisson];
				for (int i = 0; i < nPoisson; i++)
					poisson[i] = pd.sample();
			}

			// Save image in blocks
			ImageStack stack = new ImageStack(size, size);
			int start = from;
			for (int i = from; i < to; i++)
			{
				showProgress();

				// Create image
				ShortProcessor ip = new ShortProcessor(size, size);
				if (photons == 0)
				{
					for (int j = 0; j < pixelOffset.length; j++)
					{
						double p = pixelOffset[j] + rg.nextGaussian() * pixelSD[j];
						ip.set(j, clip(p));
					}
				}
				else
				{
					for (int j = 0; j < pixelOffset.length; j++)
					{
						double p = pixelOffset[j] + rg.nextGaussian() * pixelSD[j] +
								(poisson[--nPoisson] * pixelGain[j]);
						if (nPoisson == 0) // Rotate Poisson numbers
							nPoisson = poisson.length;
						ip.set(j, clip(p));
					}
				}

				// Save image
				stack.addSlice(null, ip.getPixels());
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

		/**
		 * Clip to the range for a 16-bit image.
		 *
		 * @param value
		 *            the value
		 * @return the clipped value
		 */
		private int clip(double value)
		{
			int i = (int) Math.round(value);
			if (i < 0)
				return 0;
			if (i > 65335)
				return 65335;
			return i;
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

		public int compareTo(SubDir o)
		{
			return exposureTime - o.exposureTime;
		}
	}

	private class ImageJob
	{
		float[] data;

		ImageJob(float[] data)
		{
			this.data = data;
		}
	}

	/**
	 * Used to allow multi-threading of the scoring the filters
	 */
	private class ImageWorker implements Runnable
	{
		volatile boolean finished = false;
		final BlockingQueue<ImageJob> jobs;
		final ArrayMoment moment;

		public ImageWorker(BlockingQueue<ImageJob> jobs, ArrayMoment moment)
		{
			this.jobs = jobs;
			this.moment = moment.newInstance();
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
					ImageJob job = jobs.take();
					if (job == null || job.data == null)
						break;
					if (!finished)
						// Only run jobs when not finished. This allows the queue to be emptied.
						run(job);
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

		private void run(ImageJob job)
		{
			if (Utils.isInterrupted())
			{
				finished = true;
				return;
			}
			showProgress();
			moment.add(job.data);
		}
	}

	private static final String TITLE = "sCMOS Analysis";

	private static String directory = "";
	private static boolean rollingAlgorithm = false;

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
	private ImageStack simulationStack;
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
		if (Utils.isNullOrEmpty(dir))
			return;
		directory = dir;

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

		if (simulationStack != null)
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
		if (progress % stepProgress == 0)
		{
			//IJ.showProgress(progress, totalProgress);

			// Use the actual progress bar so we can show progress 
			// when all other IJ commands cannot
			double p = ((double) progress + 1.0) / totalProgress;
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
			pixelOffset[i] = (float) pd.sample();
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
		simulationStack = new ImageStack(size, size);
		simulationStack.addSlice("Offset", pixelOffset);
		simulationStack.addSlice("Variance", pixelVariance);
		simulationStack.addSlice("Gain", pixelGain);
		IJ.save(new ImagePlus("PerPixel", simulationStack), new File(directory, "perPixelSimulation.tif").getPath());

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
			statusLine.setText("Simulating " + Utils.pleural(p, "photon"));

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
		gd.addCheckbox("Rolling_algorithm", rollingAlgorithm);
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

		// Create thread pool and workers
		ExecutorService executor = Executors.newFixedThreadPool(getThreads());
		TurboList<Future<?>> futures = new TurboList<Future<?>>(nThreads);
		TurboList<ImageWorker> workers = new TurboList<ImageWorker>(nThreads);

		double[][][] data = new double[subDirs.size()][2][];
		double[] pixelOffset = null, pixelVariance = null;

		// For each sub-directory compute the mean and variance
		final int nSubDirs = subDirs.size();
		boolean error = false;
		for (int n = 0; n < nSubDirs; n++)
		{
			SubDir sd = subDirs.getf(n);
			statusLine.setText("Analysing " + sd.name);

			// Open the series
			SeriesImageSource source = new SeriesImageSource(sd.name, sd.path.getPath());
			//source.setLogProgress(true);
			if (!source.open())
			{
				error = true;
				IJ.error(TITLE, "Failed to open image series: " + sd.path.getPath());
				break;
			}

			totalProgress = source.getFrames() + 1; // So the bar remains at 99% when workers have finished
			stepProgress = Utils.getProgressInterval(totalProgress);
			progress = 0;
			progressBar.show(0);

			ArrayMoment moment = (rollingAlgorithm) ? new RollingArrayMoment() : new SimpleArrayMoment();

			final BlockingQueue<ImageJob> jobs = new ArrayBlockingQueue<ImageJob>(nThreads * 2);
			for (int i = 0; i < nThreads; i++)
			{
				final ImageWorker worker = new ImageWorker(jobs, moment);
				workers.add(worker);
				futures.add(executor.submit(worker));
			}

			// Process the data
			for (float[] pixels = source.next(); pixels != null; pixels = source.next())
			{
				put(jobs, new ImageJob(pixels));
			}
			// Finish all the worker threads by passing in a null job
			for (int i = 0; i < nThreads; i++)
			{
				put(jobs, new ImageJob(null));
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

			// Create the final aggregate statistics
			for (ImageWorker w : workers)
				moment.add(w.moment);
			data[n][0] = moment.getFirstMoment();
			data[n][1] = moment.getVariance();

			// Reset
			futures.clear();
			workers.clear();

			Statistics s = new Statistics(data[n][0]);

			if (n != 0)
			{
				// Compute mean ADU
				Statistics signal = new Statistics();
				double[] mean = data[n][0];
				for (int i = 0; i < pixelOffset.length; i++)
					signal.add(mean[i] - pixelOffset[i]);
				Utils.log("%s Mean = %s +/- %s. Signal = %s +/- %s ADU", sd.name, Utils.rounded(s.getMean()),
						Utils.rounded(s.getStandardDeviation()), Utils.rounded(signal.getMean()),
						Utils.rounded(signal.getStandardDeviation()));

				// TODO - optionally save
				ImageStack stack = new ImageStack(source.getWidth(), source.getHeight());
				stack.addSlice("Mean", Utils.toFloat(data[n][0]));
				stack.addSlice("Variance", Utils.toFloat(data[n][1]));
				IJ.save(new ImagePlus("PerPixel", stack), new File(directory, "perPixel" + sd.name + ".tif").getPath());
			}
			else
			{
				pixelOffset = data[0][0];
				pixelVariance = data[0][1];
				Statistics var = new Statistics(pixelVariance);
				Utils.log("%s Offset = %s +/- %s. Variance = %s +/- %s", sd.name, Utils.rounded(s.getMean()),
						Utils.rounded(s.getStandardDeviation()), Utils.rounded(var.getMean()),
						Utils.rounded(var.getStandardDeviation()));
			}
		}

		Utils.setShowStatus(true);
		Utils.setShowProgress(true);
		IJ.showProgress(1);

		executor.shutdown();

		if (error)
			return;

		// Compute the gain
		statusLine.setText("Computing gain");

		double[] pixelGain = new double[pixelOffset.length];

		// Ignore first as this is the 0 exposure image
		for (int i = 0; i < pixelGain.length; i++)
		{
			// Use equation 2.5 from the Huang et al paper.
			double bibiT = 0;
			double biaiT = 0;
			for (int n = 1; n < nSubDirs; n++)
			{
				double bi = data[n][0][i] - pixelOffset[i];
				double ai = data[n][1][i] - pixelVariance[i];
				bibiT += bi * bi;
				biaiT += bi * ai;
			}
			pixelGain[i] = biaiT / bibiT;
		}

		Statistics s = new Statistics(pixelGain);
		Utils.log("Gain Mean = %s +/- %s", Utils.rounded(s.getMean()), Utils.rounded(s.getStandardDeviation()));

		// TODO - Histogram of offset, variance and gain

		// Save
		simulationStack = new ImageStack(size, size);
		simulationStack.addSlice("Offset", Utils.toFloat(pixelOffset));
		simulationStack.addSlice("Variance", Utils.toFloat(pixelVariance));
		simulationStack.addSlice("Gain", Utils.toFloat(pixelGain));
		IJ.save(new ImagePlus("PerPixel", simulationStack), new File(directory, "perPixel.tif").getPath());
		IJ.showStatus(""); // Remove the status from the ij.io.ImageWriter class

		Utils.log("Analysis time = " + Utils.timeToString(System.currentTimeMillis() - start));
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
		// TODO Auto-generated method stub

	}
}
