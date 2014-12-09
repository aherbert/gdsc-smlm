package gdsc.smlm.ij.plugins;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2014 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import gdsc.smlm.fitting.FitConfiguration;
import gdsc.smlm.fitting.FitFunction;
import gdsc.smlm.fitting.FitSolver;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolver;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.ImageConverter;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.utils.Statistics;
import gdsc.smlm.utils.StoredDataStatistics;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.text.TextWindow;

import java.awt.Rectangle;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

/**
 * Fits the benchmark image created by CreateData plugin.
 */
public class BenchmarkFit implements PlugIn
{
	private static final String TITLE = "Benchmark Fit";

	private static int regionSize = 3;
	private static double psfWidth = 1;
	private static boolean showHistograms = false;

	private static TextWindow summaryTable = null;

	private FitConfiguration fitConfig;
	private ImagePlus imp;
	private CreateData.BenchmarkParameters benchmarkParameters;
	private static final String[] NAMES = new String[] { "dB (photons)", "dSignal (photons)", "dAngle (deg)",
			"dX (nm)", "dY (nm)", "dSx (nm)", "dSy (nm)", "Time (ms)", "dN (photons)" };
	private static final int TIME = 7;
	private static final int ACTUAL_SIGNAL = 8;
	private double[] answer = new double[7];
	private double[] lb, ub = null;
	private double[] lc, uc = null;

	private class Worker implements Runnable
	{
		volatile boolean finished = false;
		final BlockingQueue<Integer> jobs;
		final Statistics[] stats = new Statistics[9];
		final ImageStack stack;
		final Rectangle region;
		final double[][] xy;
		final FitConfiguration fitConfig;

		double[] data = null;

		public Worker(BlockingQueue<Integer> jobs, ImageStack stack, Rectangle region, double[][] xy,
				FitConfiguration fitConfig)
		{
			this.jobs = jobs;
			this.stack = stack;
			this.region = region;
			this.fitConfig = (FitConfiguration) fitConfig.clone();
			// Clone the positions array
			this.xy = new double[xy.length][];
			for (int i = 0; i < xy.length; i++)
				this.xy[i] = xy[i].clone();
			for (int i = 0; i < stats.length; i++)
				stats[i] = (showHistograms) ? new StoredDataStatistics() : new Statistics();
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
			// Extract the data
			data = ImageConverter.getDoubleData(stack.getPixels(frame + 1), stack.getWidth(), stack.getHeight(),
					region, data);

			// Get the background and signal estimate
			final double b = getBackground(data, region.width, region.height);
			final double signal = getSignal(data, b);

			// Find centre-of-mass estimate
			getCentreOfMass(data, region.width, region.height, xy[0]);

			double[] initialParams = new double[7];
			initialParams[Gaussian2DFunction.BACKGROUND] = b;
			initialParams[Gaussian2DFunction.SIGNAL] = signal;
			initialParams[Gaussian2DFunction.X_SD] = initialParams[Gaussian2DFunction.Y_SD] = psfWidth;

			// Subtract the bias
			final double bias = benchmarkParameters.bias;
			if (fitConfig.isRemoveBiasBeforeFitting())
			{
				initialParams[0] -= bias;
				for (int i = 0; i < data.length; i++)
					data[i] -= bias;
			}

			double[] error = new double[1];
			double[][] result = new double[xy.length][];
			int c = 0;
			long start = System.nanoTime();
			for (double[] centre : xy)
			{
				// Do fitting
				final double[] params = initialParams.clone();
				params[Gaussian2DFunction.X_POSITION] = centre[0];
				params[Gaussian2DFunction.Y_POSITION] = centre[1];
				fitConfig.initialise(1, region.width, params);
				FunctionSolver solver = fitConfig.getFunctionSolver();
				if (solver.isBounded())
					setBounds(solver);
				if (solver.isConstrained())
					setConstraints(solver);
				if (solver.fit(data.length, data, null, params, null, error, 0) == FitStatus.OK)
				{
					// Subtract the fitted bias from the background
					if (!fitConfig.isRemoveBiasBeforeFitting())
						params[0] -= bias;
					result[c++] = params;
				}
			}
			long time = System.nanoTime() - start;

			// Store the results from each run
			for (int i = 0; i < c; i++)
			{
				final double[] params = result[i];
				for (int j = 0; j < params.length; j++)
				{
					stats[j].add(params[j] - answer[j]);
				}
				stats[7].add(time);
				stats[8].add(params[Gaussian2DFunction.SIGNAL] - benchmarkParameters.p[frame]);
			}
		}
	}

	public void run(String arg)
	{
		if (CreateData.benchmarkParameters == null)
		{
			IJ.error(TITLE, "No benchmark parameters in memory");
			return;
		}
		benchmarkParameters = CreateData.benchmarkParameters;
		imp = WindowManager.getImage(CreateData.CREATE_DATA_IMAGE_TITLE);
		if (imp == null || imp.getStackSize() != benchmarkParameters.frames)
		{
			IJ.error(TITLE, "No benchmark image to match the parameters in memory");
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

		final double sa = PSFCalculator.squarePixelAdjustment(benchmarkParameters.s, benchmarkParameters.a);
		gd.addMessage(String.format(
				"Fits the benchmark image created by CreateData plugin.\nPSF width = %s, adjusted = %s",
				Utils.rounded(benchmarkParameters.s / benchmarkParameters.a), Utils.rounded(sa / benchmarkParameters.a)));

		String filename = SettingsManager.getSettingsFilename();
		GlobalSettings settings = SettingsManager.loadSettings(filename);
		fitConfig = settings.getFitEngineConfiguration().getFitConfiguration();

		gd.addSlider("Region_size", 2, 20, regionSize);
		gd.addNumericField("PSF_width", psfWidth, 3);
		String[] solverNames = SettingsManager.getNames((Object[]) FitSolver.values());
		gd.addChoice("Fit_solver", solverNames, solverNames[fitConfig.getFitSolver().ordinal()]);
		String[] functionNames = SettingsManager.getNames((Object[]) FitFunction.values());
		gd.addChoice("Fit_function", functionNames, functionNames[fitConfig.getFitFunction().ordinal()]);
		gd.addCheckbox("Show_histograms", showHistograms);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		regionSize = (int) Math.abs(gd.getNextNumber());
		psfWidth = Math.abs(gd.getNextNumber());
		fitConfig.setFitSolver(gd.getNextChoiceIndex());
		fitConfig.setFitFunction(gd.getNextChoiceIndex());
		showHistograms = gd.getNextBoolean();

		if (regionSize < 1)
			regionSize = 1;

		if (gd.invalidNumber())
			return false;

		return PeakFit.configureFitSolver(settings, filename, false);
	}

	private void run()
	{
		// Initialise the answer. Convert to units of the image (ADUs and pixels)
		answer[Gaussian2DFunction.BACKGROUND] = benchmarkParameters.b * benchmarkParameters.gain;
		answer[Gaussian2DFunction.SIGNAL] = benchmarkParameters.signal * benchmarkParameters.gain;
		answer[Gaussian2DFunction.X_POSITION] = benchmarkParameters.x;
		answer[Gaussian2DFunction.Y_POSITION] = benchmarkParameters.y;
		answer[Gaussian2DFunction.X_SD] = benchmarkParameters.s / benchmarkParameters.a;
		answer[Gaussian2DFunction.Y_SD] = benchmarkParameters.s / benchmarkParameters.a;

		// Set up the fit region. Always round down since 0.5 is the centre of the pixel.
		int x = (int) benchmarkParameters.x;
		int y = (int) benchmarkParameters.y;
		Rectangle region = new Rectangle(x - regionSize, y - regionSize, 2 * regionSize + 1, 2 * regionSize + 1);
		if (!new Rectangle(0, 0, imp.getWidth(), imp.getHeight()).contains(region))
		{
			IJ.error(TITLE, "Fit region does not fit within the image");
		}

		// Adjust the centre & account for 0.5 pixel offset during fitting
		x -= region.x;
		y -= region.y;
		answer[Gaussian2DFunction.X_POSITION] -= (region.x + 0.5);
		answer[Gaussian2DFunction.Y_POSITION] -= (region.y + 0.5);

		// Configure for fitting
		fitConfig.setBackgroundFitting(true);
		fitConfig.setComputeDeviations(false);

		// Fit using the centre-of-mass estimate and 4 corners of image
		final double[][] xy = new double[][] { { 0, 0 }, { x, y }, { x, y + 1 }, { x + 1, y }, { x + 1, y + 1 } };

		final ImageStack stack = imp.getImageStack();

		// Create a pool of workers
		int nThreads = Prefs.getThreads();
		BlockingQueue<Integer> jobs = new ArrayBlockingQueue<Integer>(nThreads * 3);
		List<Worker> workers = new LinkedList<Worker>();
		List<Thread> threads = new LinkedList<Thread>();
		for (int i = 0; i < nThreads; i++)
		{
			Worker worker = new Worker(jobs, stack, region, xy, fitConfig);
			Thread t = new Thread(worker);
			workers.add(worker);
			threads.add(t);
			t.start();
		}

		// Fit the frames with simulated photons
		final int totalFrames = benchmarkParameters.frames;
		final int step = (totalFrames > 400) ? totalFrames / 200 : 2;
		for (int i = 0; i < totalFrames; i++)
		{
			if (benchmarkParameters.p[i] > 0)
			{
				put(jobs, i);
				if (i % step == 0)
				{
					IJ.showProgress(i, totalFrames);
					IJ.showStatus("Frame: " + i + " / " + totalFrames);
				}
			}
		}
		// Finish all the worker threads by passing in a null job
		for (int i = 0; i < threads.size(); i++)
		{
			put(jobs, -1);
		}
		IJ.showProgress(1);
		IJ.showStatus("Collecting results ...");

		// Collect the results
		Statistics[] stats = new Statistics[NAMES.length];
		for (int i = 0; i < threads.size(); i++)
		{
			try
			{
				threads.get(i).join();
				Statistics[] next = workers.get(i).stats;
				for (int j = 0; j < next.length; j++)
				{
					if (stats[j] == null)
						stats[j] = next[j];
					else
						stats[j].add(next[j]);
				}
			}
			catch (InterruptedException e)
			{
				e.printStackTrace();
			}
		}
		threads.clear();
		workers.clear();

		// Show a table of the results
		summariseResults(stats);

		// TODO Optionally show histograms
		if (showHistograms)
		{

		}
		
		IJ.showStatus("");		
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
	 * Set background using the average value of the edge in the data
	 * 
	 * @param data
	 * @param maxx
	 * @param maxy
	 * @return The background
	 */
	private static double getBackground(double[] data, int maxx, int maxy)
	{
		double background = 0;
		for (int xi = 0; xi < maxx; xi++)
			background += data[xi] + data[maxx * (maxy - 1) + xi];
		for (int yi = 0; yi < maxy; yi++)
			background += data[maxx * yi] + data[maxx * yi + (maxx - 1)];
		background /= 2 * (maxx + maxy);
		return background;
	}

	/**
	 * Sum the intensity above background to estimate the signal
	 * 
	 * @param data
	 * @param b
	 *            background
	 * @return The signal
	 */
	private static double getSignal(double[] data, double b)
	{
		double s = 0;
		for (double d : data)
			s += d;
		// subtract the background per pixel
		return s - b * data.length;
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
	private static void getCentreOfMass(double[] data, int maxx, int maxy, double[] com)
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

	private void setBounds(FunctionSolver solver)
	{
		if (ub == null)
		{
			ub = new double[7];
			lb = new double[7];

			ub[Gaussian2DFunction.BACKGROUND] = 2 * benchmarkParameters.b * benchmarkParameters.gain;
			double signal = benchmarkParameters.signal * benchmarkParameters.gain;
			lb[Gaussian2DFunction.SIGNAL] = signal * 0.5;
			ub[Gaussian2DFunction.SIGNAL] = signal * 2;
			ub[Gaussian2DFunction.X_POSITION] = 2 * regionSize + 1;
			ub[Gaussian2DFunction.Y_POSITION] = 2 * regionSize + 1;
			lb[Gaussian2DFunction.ANGLE] = -Math.PI;
			ub[Gaussian2DFunction.ANGLE] = Math.PI;
			double wf = 1.5;
			double s = benchmarkParameters.s / benchmarkParameters.a;
			lb[Gaussian2DFunction.X_SD] = s / wf;
			ub[Gaussian2DFunction.X_SD] = s * wf;
			lb[Gaussian2DFunction.Y_SD] = s / wf;
			ub[Gaussian2DFunction.Y_SD] = s * wf;
		}
		solver.setBounds(lb, ub);
	}

	private void setConstraints(FunctionSolver solver)
	{
		if (uc == null)
		{
			lc = new double[7];
			uc = new double[7];
			Arrays.fill(lc, Float.NEGATIVE_INFINITY);
			Arrays.fill(uc, Float.POSITIVE_INFINITY);
			lc[Gaussian2DFunction.BACKGROUND] = 0;
			lc[Gaussian2DFunction.SIGNAL] = 0;
		}
		solver.setConstraints(lc, uc);
	}

	private void summariseResults(Statistics[] stats)
	{
		createTable();
		StringBuilder sb = new StringBuilder();
		sb.append(benchmarkParameters.molecules).append("\t");
		sb.append(Utils.rounded(benchmarkParameters.s)).append("\t");
		sb.append(Utils.rounded(benchmarkParameters.a)).append("\t");
		// Report XY in nm from the pixel centre
		sb.append(Utils.rounded(distanceFromCentre(benchmarkParameters.x))).append("\t");
		sb.append(Utils.rounded(distanceFromCentre(benchmarkParameters.y))).append("\t");
		sb.append(Utils.rounded(benchmarkParameters.gain)).append("\t");
		sb.append(Utils.rounded(benchmarkParameters.readNoise)).append("\t");
		sb.append(Utils.rounded(benchmarkParameters.b)).append("\t");
		sb.append(Utils.rounded(benchmarkParameters.precisionN)).append("\t");
		sb.append(Utils.rounded(benchmarkParameters.precisionX)).append("\t");
		sb.append(Utils.rounded(benchmarkParameters.precisionXML)).append("\t");
		int n = 2 * regionSize + 1;
		//sb.append(n).append("x");
		sb.append(n).append("\t");
		sb.append(Utils.rounded(psfWidth * benchmarkParameters.a)).append("\t");
		sb.append(fitConfig.getFitFunction().toString() + ":" + PeakFit.getSolverName(fitConfig));
		if (fitConfig.getFitSolver() == FitSolver.MLE && fitConfig.isModelCamera())
		{
			// Add details of the noise model for the MLE
			sb.append(":EM=").append(fitConfig.isEmCCD());
			sb.append(":G=").append(fitConfig.getGain());
			sb.append(":N=").append(fitConfig.getReadNoise());
		}
		sb.append("\t");
		final double recall = (double) (stats[0].getN() / 5) / benchmarkParameters.molecules;
		sb.append(Utils.rounded(recall));

		// Convert to units of the image (ADUs and pixels)		
		double[] convert = new double[stats.length];
		convert[Gaussian2DFunction.BACKGROUND] = 1 / benchmarkParameters.gain;
		convert[Gaussian2DFunction.SIGNAL] = 1 / benchmarkParameters.gain;
		convert[Gaussian2DFunction.ANGLE] = (fitConfig.isAngleFitting()) ? 180.0 / Math.PI : 0;
		convert[Gaussian2DFunction.X_POSITION] = benchmarkParameters.a;
		convert[Gaussian2DFunction.Y_POSITION] = benchmarkParameters.a;
		convert[Gaussian2DFunction.X_SD] = (fitConfig.isWidth0Fitting()) ? benchmarkParameters.a : 0;
		convert[Gaussian2DFunction.Y_SD] = (fitConfig.isWidth1Fitting()) ? benchmarkParameters.a : 0;
		convert[TIME] = 1e-6;
		convert[ACTUAL_SIGNAL] = 1 / benchmarkParameters.gain;

		for (int i = 0; i < stats.length; i++)
		{
			if (convert[i] != 0)
				sb.append("\t").append(Utils.rounded(stats[i].getMean() * convert[i])).append("\t")
						.append(Utils.rounded(stats[i].getStandardDeviation() * convert[i]));
			else
				sb.append("\t0\t0");
		}
		summaryTable.append(sb.toString());
	}

	private double distanceFromCentre(double x)
	{
		x -= 0.5;
		int i = (int) Math.round(x);
		x = x - i;
		return x * benchmarkParameters.a;
	}

	private void createTable()
	{
		if (summaryTable == null || !summaryTable.isVisible())
		{
			summaryTable = new TextWindow(TITLE, createHeader(), "", 1000, 300);
			summaryTable.setVisible(true);
		}
	}

	private String createHeader()
	{
		StringBuilder sb = new StringBuilder(
				"Molecules\tS (nm)\ta (nm)\tX (nm)\tY (nm)\tGain\tReadNoise (ADUs)\tB (photons)\tLimit N\tLimit X\tLimit X ML\tRegion\tWidth\tMethod\tRecall");
		for (int i = 0; i < NAMES.length; i++)
		{
			sb.append("\t").append(NAMES[i]).append("\t+/-");
		}
		return sb.toString();
	}
}
