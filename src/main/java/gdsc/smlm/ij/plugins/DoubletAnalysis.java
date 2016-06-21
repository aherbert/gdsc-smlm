package gdsc.smlm.ij.plugins;

import gdsc.smlm.engine.DataFilter;
import gdsc.smlm.engine.DataFilterType;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.filters.MaximaSpotFilter;
import gdsc.smlm.filters.Spot;

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
import gdsc.smlm.fitting.Gaussian2DFitter;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.ij.plugins.CreateData.SimulationParameters;
import gdsc.smlm.ij.plugins.ResultsMatchCalculator.PeakResultPoint;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.ImageConverter;
import gdsc.smlm.results.Calibration;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.core.clustering.Cluster;
import gdsc.core.clustering.ClusterPoint;
import gdsc.core.clustering.ClusteringAlgorithm;
import gdsc.core.clustering.ClusteringEngine;
import gdsc.core.ij.Utils;
import gdsc.core.match.BasePoint;
import gdsc.core.match.Coordinate;
import gdsc.core.match.MatchCalculator;
import gdsc.core.match.PointPair;
import gdsc.core.utils.Maths;
import gdsc.core.utils.Statistics;
import gdsc.core.utils.StoredDataStatistics;
import gdsc.core.utils.TextUtils;
import gdsc.core.utils.NoiseEstimator.Method;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot2;
import ij.gui.PlotWindow;
import ij.plugin.PlugIn;
import ij.plugin.WindowOrganiser;
import ij.text.TextWindow;

import java.awt.Rectangle;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.math3.util.FastMath;

/**
 * Fits spots created by CreateData plugin.
 * <p>
 * Performs clustering to determine if spots are either single or doublets. Larger clusters are ignored. Outputs a table
 * of the single and double fit for each spot with metrics. This can be used to determine the best settings for optimum
 * doublet fitting and filtering.
 */
public class DoubletAnalysis implements PlugIn
{
	private static final String TITLE = "Doublet Analysis";

	static FitConfiguration fitConfig;
	private static FitEngineConfiguration config;
	private static Calibration cal;
	private static int lastId = 0;
	static
	{
		cal = new Calibration();
		fitConfig = new FitConfiguration();
		config = new FitEngineConfiguration(fitConfig);
		// Set some default fit settings here ...
		// Ensure all candidates are fitted
		config.setFailuresLimit(-1);
		fitConfig.setFitValidation(true);
		fitConfig.setMinPhotons(0); // Do not allow negative photons 
		fitConfig.setCoordinateShiftFactor(0); // Disable
		fitConfig.setPrecisionThreshold(0);
		fitConfig.setWidthFactor(0);

		fitConfig.setBackgroundFitting(true);
		fitConfig.setMinIterations(0);
		fitConfig.setNoise(0);
		config.setNoiseMethod(Method.QUICK_RESIDUALS_LEAST_MEAN_OF_SQUARES);

		fitConfig.setBackgroundFitting(true);
		fitConfig.setNotSignalFitting(false);
		fitConfig.setComputeDeviations(false);
	}

	private static boolean showHistograms = false;
	private static boolean saveRawData = false;
	private static String rawDataDirectory = "";
	private static int histogramBins = 100;

	private static TextWindow summaryTable = null, analysisTable = null;

	private ImagePlus imp;
	private MemoryPeakResults results;
	private CreateData.SimulationParameters simulationParameters;

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

	/**
	 * Used to allow multi-threading of the fitting method
	 */
	private class Worker implements Runnable
	{
		volatile boolean finished = false;
		final BlockingQueue<Integer> jobs;
		final Statistics[] stats = new Statistics[NAMES.length];
		final ImageStack stack;
		final HashMap<Integer, ArrayList<Coordinate>> actualCoordinates;
		final int fitting;
		final FitConfiguration fitConfig;
		final MaximaSpotFilter spotFilter;
		final double limit;
		final int[] spotHistogram, resultHistogram;

		float[] data = null;
		private double[] lb, ub = null;
		private double[] lc, uc = null;

		public Worker(BlockingQueue<Integer> jobs, ImageStack stack,
				HashMap<Integer, ArrayList<Coordinate>> actualCoordinates, FitConfiguration fitConfig, int maxCount)
		{
			this.jobs = jobs;
			this.stack = stack;
			this.actualCoordinates = actualCoordinates;
			fitting = config.getRelativeFitting();
			this.fitConfig = fitConfig.clone();
			this.spotFilter = config.createSpotFilter(true);
			limit = spotFilter.getSearch() * spotFilter.getSearch();
			spotHistogram = new int[maxCount + 1];
			resultHistogram = new int[spotHistogram.length];

			for (int i = 0; i < stats.length; i++)
				stats[i] = (showHistograms || saveRawData) ? new StoredDataStatistics() : new Statistics();

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

			Coordinate[] actual = ResultsMatchCalculator.getCoordinates(actualCoordinates, frame);

			// Extract the data
			final int width = stack.getWidth();
			final int height = stack.getHeight();
			data = ImageConverter.getData(stack.getPixels(frame), width, height, null, data);

			// Smooth the image and identify spots with a filter
			Spot[] spots = spotFilter.rank(data, width, height);

			// Match the each actual result to the closest filter candidate
			int[] matches = new int[actual.length];
			for (int i = 0; i < actual.length; i++)
			{
				double dmin = limit;
				int match = -1;
				for (int j = 0; j < spots.length; j++)
				{
					double d2 = actual[i].distance2(spots[j].x, spots[j].y);
					if (dmin > d2)
					{
						dmin = d2;
						match = j;
					}
				}
				matches[i] = match;
			}

			// Identify single and doublets (and other)
			int singles = 0, doublets = 0, multiples = 0, total = 0;
			for (int i = 0; i < actual.length; i++)
			{
				if (matches[i] != -1)
				{
					// Count all matches
					int n = 0;
					final int j = matches[i];
					for (int ii = i; ii < matches.length; ii++)
					{
						if (matches[ii] == j)
						{
							n++;
							// Reset to avoid double counting
							matches[ii] = -1;
						}
					}
					switch (n)
					{
						//@formatter:off
						case 1: singles++; break;
						case 2: doublets++; break;
						default: multiples++;
						//@formatter:on
					}
					// Store the number of actual results that match to a spot
					spotHistogram[n]++;
					// Store the number of results that match to a spot with n results
					resultHistogram[n] += n;
					total += n;

					if (n > 1)
					{
						//System.out.printf("Frame %d match %d : %d,%d \n", frame, n, spots[j].x, spots[j].y);
						// TODO - Create an output stack with coloured ROI overlay for each n=1, n=2, n=other
						// to check that the doublets are correctly identified.
					}

					// Fit the candidates (as per the FitWorker logic)
					// (Fit even multiple since this is what the FitWorker will do) 

					// Compute residuals and fit as a doublet

				}
			}
			System.out.printf("Frame %d, singles=%d, doublets=%d, multi=%d\n", frame, singles, doublets, multiples);
			resultHistogram[0] += actual.length - total;

			// Extract the data
			//data = ImageConverter.getDoubleData(stack.getPixels(frame + 1), stack.getWidth(), stack.getHeight(), region,
			//		data);

			// Fit a single

			// Fit a doublet

			//			
			//			
			//
			//			final int size = region.height;
			//
			//			// Get the background and signal estimate
			//			final double b = getBackground(data, size, size);
			//			final double signal = getSignal(data, b);
			//
			//			// Find centre-of-mass estimate
			//			if (comFitting)
			//			{
			//				getCentreOfMass(data, size, size, xy[xy.length - 1]);
			//
			//				double dx = xy[xy.length - 1][0] - answer[Gaussian2DFunction.X_POSITION];
			//				double dy = xy[xy.length - 1][1] - answer[Gaussian2DFunction.Y_POSITION];
			//				if (dx * dx + dy * dy < startOffset * startOffset)
			//					comValid.incrementAndGet();
			//			}
			//
			//			double[] initialParams = new double[7];
			//			initialParams[Gaussian2DFunction.BACKGROUND] = b;
			//			initialParams[Gaussian2DFunction.SIGNAL] = signal;
			//			initialParams[Gaussian2DFunction.X_SD] = initialParams[Gaussian2DFunction.Y_SD] = psfWidth;
			//
			//			// Subtract the bias
			//			final double bias = simulationParameters.bias;
			//			if (fitConfig.isRemoveBiasBeforeFitting())
			//			{
			//				// MLE can handle negative data 
			//				initialParams[Gaussian2DFunction.BACKGROUND] -= bias;
			//				for (int i = 0; i < data.length; i++)
			//					data[i] -= bias;
			//			}
			//
			//			double[][] bounds = null;
			//			double[] error = new double[1];
			//			double[][] result = new double[xy.length][];
			//			long[] time = new long[xy.length];
			//			int c = 0;
			//			int resultPosition = frame;
			//			for (double[] centre : xy)
			//			{
			//				long start = System.nanoTime();
			//
			//				// Do fitting
			//				final double[] params = initialParams.clone();
			//				params[Gaussian2DFunction.X_POSITION] = centre[0];
			//				params[Gaussian2DFunction.Y_POSITION] = centre[1];
			//				fitConfig.initialise(1, size, params);
			//				FunctionSolver solver = fitConfig.getFunctionSolver();
			//				if (solver.isBounded())
			//					bounds = setBounds(solver, initialParams, bounds);
			//				if (solver.isConstrained())
			//					setConstraints(solver);
			//				final FitStatus status = solver.fit(data.length, data, null, params, null, error, 0);
			//				if (isValid(status, params, size))
			//				{
			//					// Subtract the fitted bias from the background
			//					if (!fitConfig.isRemoveBiasBeforeFitting())
			//						params[Gaussian2DFunction.BACKGROUND] -= bias;
			//					result[c] = params;
			//					time[c] = System.nanoTime() - start;
			//					// Store all the results for later analysis
			//					results[resultPosition] = params;
			//					resultsTime[resultPosition] = time[c];
			//					c++;
			//				}
			//				else
			//				{
			//					//System.out.println(status);
			//				}
			//				resultPosition += totalFrames;
			//			}
			//
			//			addResults(stats, answer, simulationParameters.p[frame], sa, time, result, c);

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
				ub = new double[7];
				lb = new double[7];

				// Background could be zero so always have an upper limit
				ub[Gaussian2DFunction.BACKGROUND] = Math.max(0, 2 * simulationParameters.b * simulationParameters.gain);
				if (!fitConfig.isRemoveBiasBeforeFitting())
				{
					lb[Gaussian2DFunction.BACKGROUND] += simulationParameters.bias;
					ub[Gaussian2DFunction.BACKGROUND] += simulationParameters.bias;
				}

				//				double signal = simulationParameters.getSignal() * simulationParameters.gain;
				//				lb[Gaussian2DFunction.SIGNAL] = signal * 0.5;
				//				ub[Gaussian2DFunction.SIGNAL] = signal * 2;
				//				ub[Gaussian2DFunction.X_POSITION] = 2 * regionSize + 1;
				//				ub[Gaussian2DFunction.Y_POSITION] = 2 * regionSize + 1;
				lb[Gaussian2DFunction.ANGLE] = -Math.PI;
				ub[Gaussian2DFunction.ANGLE] = Math.PI;
				double wf = 1.5;
				double s = simulationParameters.s / simulationParameters.a;
				lb[Gaussian2DFunction.X_SD] = s / wf;
				ub[Gaussian2DFunction.X_SD] = s * wf;
				lb[Gaussian2DFunction.Y_SD] = s / wf;
				ub[Gaussian2DFunction.Y_SD] = s * wf;
			}
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
			if (fitConfig.isWidth0Fitting())
			{
				if (params[Gaussian2DFunction.X_SD] < lb[Gaussian2DFunction.X_SD] ||
						params[Gaussian2DFunction.X_SD] > ub[Gaussian2DFunction.X_SD])
				{
					return false;
				}
			}
			if (fitConfig.isWidth1Fitting())
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
			runAnalysis();
		}
		else
		{
			simulationParameters = CreateData.simulationParameters;
			if (simulationParameters == null)
			{
				IJ.error(TITLE, "No benchmark spot parameters in memory");
				return;
			}
			imp = WindowManager.getImage(CreateData.CREATE_DATA_IMAGE_TITLE);
			if (imp == null)
			{
				IJ.error(TITLE, "No benchmark image");
				return;
			}
			results = MemoryPeakResults.getResults(CreateData.CREATE_DATA_IMAGE_TITLE + " (Create Data)");
			if (results == null)
			{
				IJ.error(TITLE, "No benchmark results in memory");
				return;
			}

			if (!showDialog())
				return;

			run();
		}
	}

	private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		final double sa = getSa();
		gd.addMessage(
				String.format("Fits the benchmark image created by CreateData plugin.\nPSF width = %s, adjusted = %s",
						Utils.rounded(simulationParameters.s / simulationParameters.a), Utils.rounded(sa)));

		// For each new benchmark width, reset the PSF width to the square pixel adjustment
		if (lastId != simulationParameters.id)
		{
			double w = sa;
			fitConfig.setInitialPeakStdDev(w);
		}

		// Collect options for fitting
		gd.addNumericField("Initial_StdDev", fitConfig.getInitialPeakStdDev0(), 3);
		String[] filterTypes = SettingsManager.getNames((Object[]) DataFilterType.values());
		gd.addChoice("Spot_filter_type", filterTypes, filterTypes[config.getDataFilterType().ordinal()]);
		String[] filterNames = SettingsManager.getNames((Object[]) DataFilter.values());
		gd.addChoice("Spot_filter", filterNames, filterNames[config.getDataFilter(0).ordinal()]);
		gd.addSlider("Smoothing", 0, 2.5, config.getSmooth(0));
		gd.addSlider("Search_width", 0.5, 2.5, config.getSearch());
		gd.addSlider("Border", 0.5, 2.5, config.getBorder());
		gd.addSlider("Fitting_width", 2, 4.5, config.getFitting());
		String[] solverNames = SettingsManager.getNames((Object[]) FitSolver.values());
		gd.addChoice("Fit_solver", solverNames, solverNames[fitConfig.getFitSolver().ordinal()]);
		String[] functionNames = SettingsManager.getNames((Object[]) FitFunction.values());
		gd.addChoice("Fit_function", functionNames, functionNames[fitConfig.getFitFunction().ordinal()]);

		gd.addCheckbox("Show_histograms", showHistograms);
		gd.addCheckbox("Save_raw_data", saveRawData);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		fitConfig.setInitialPeakStdDev(gd.getNextNumber());
		config.setDataFilterType(gd.getNextChoiceIndex());
		config.setDataFilter(gd.getNextChoiceIndex(), Math.abs(gd.getNextNumber()), 0);
		config.setSearch(gd.getNextNumber());
		config.setBorder(gd.getNextNumber());
		config.setFitting(gd.getNextNumber());
		fitConfig.setFitSolver(gd.getNextChoiceIndex());
		fitConfig.setFitFunction(gd.getNextChoiceIndex());

		showHistograms = gd.getNextBoolean();
		saveRawData = gd.getNextBoolean();

		if (gd.invalidNumber())
			return false;
		GlobalSettings settings = new GlobalSettings();
		settings.setFitEngineConfiguration(config);
		settings.setCalibration(cal);

		// Copy simulation defaults if a new simulation
		if (lastId != simulationParameters.id)
		{
			cal.nmPerPixel = simulationParameters.a;
			cal.gain = simulationParameters.gain;
			cal.amplification = simulationParameters.amplification;
			cal.exposureTime = 100;
			cal.readNoise = simulationParameters.readNoise;
			cal.bias = simulationParameters.bias;
			cal.emCCD = simulationParameters.emCCD;
		}
		if (!PeakFit.configureFitSolver(settings, null, false))
			return false;

		lastId = simulationParameters.id;

		if (true && showHistograms)
		{
			gd = new GenericDialog(TITLE);
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
		final double sa = PSFCalculator.squarePixelAdjustment(simulationParameters.s, simulationParameters.a) /
				simulationParameters.a;
		return sa;
	}

	private void run()
	{
		final ImageStack stack = imp.getImageStack();

		// Get the coordinates per frame
		HashMap<Integer, ArrayList<Coordinate>> actualCoordinates = ResultsMatchCalculator
				.getCoordinates(results.getResults(), false);

		int maxCount = 0;
		for (ArrayList<Coordinate> list : actualCoordinates.values())
			if (maxCount < list.size())
				maxCount = list.size();

		// Create a pool of workers
		int nThreads = Prefs.getThreads();
		BlockingQueue<Integer> jobs = new ArrayBlockingQueue<Integer>(nThreads * 2);
		List<Worker> workers = new LinkedList<Worker>();
		List<Thread> threads = new LinkedList<Thread>();
		for (int i = 0; i < nThreads; i++)
		{
			Worker worker = new Worker(jobs, stack, actualCoordinates, fitConfig, maxCount);
			Thread t = new Thread(worker);
			workers.add(worker);
			threads.add(t);
			t.start();
		}

		// Fit the frames
		final int totalFrames = actualCoordinates.size();
		final int step = Utils.getProgressInterval(totalFrames);
		int progress = 0;
		for (int frame : actualCoordinates.keySet())
		{
			put(jobs, frame);
			if (++progress % step == 0)
			{
				if (Utils.showStatus("Frame: " + progress + " / " + totalFrames))
					IJ.showProgress(progress, totalFrames);
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

		IJ.showProgress(1);
		IJ.showStatus("Collecting results ...");

		// Collect the results

		double[] spotHistogram = new double[maxCount];
		double[] resultHistogram  = new double[maxCount];
		for (int i = 0; i < workers.size(); i++)
		{
			final int[] h1 = workers.get(i).spotHistogram;
			final int[] h2 = workers.get(i).resultHistogram;
			for (int j = 0; j < spotHistogram.length; j++)
			{
				spotHistogram[j] += h1[j];
				resultHistogram[j] += h2[j];
			}
		}
		workers.clear();

		WindowOrganiser o = new WindowOrganiser();
		
		String title = TITLE + " Candidate Histogram";
		Plot2 plot = new Plot2(title, "N results in candidate", "Count");
		double max = Maths.max(spotHistogram);
		plot.setLimits(0, spotHistogram.length, 0, max * 1.05);
		plot.addPoints(Utils.newArray(spotHistogram.length, 0, 1.0), spotHistogram, Plot2.BAR);
		PlotWindow pw = Utils.display(title, plot);
		if (Utils.isNewWindow())
			o.add(pw.getImagePlus().getID());

		title = TITLE + " Assigned Result Histogram";
		plot = new Plot2(title, "N results in assigned spot", "Count");
		max = Maths.max(resultHistogram);
		plot.setLimits(0, resultHistogram.length, 0, max * 1.05);
		plot.addPoints(Utils.newArray(resultHistogram.length, 0, 1.0), resultHistogram, Plot2.BAR);
		pw = Utils.display(title, plot);
		if (Utils.isNewWindow())
			o.add(pw.getImagePlus().getID());
		
		o.tile();
		
		//		
		//		Statistics[] stats = new Statistics[NAMES.length];
		//		for (int i = 0; i < workers.size(); i++)
		//		{
		//			Statistics[] next = workers.get(i).stats;
		//			for (int j = 0; j < next.length; j++)
		//			{
		//				if (stats[j] == null)
		//					stats[j] = next[j];
		//				else
		//					stats[j].add(next[j]);
		//			}
		//		}
		//		workers.clear();
		//
		//		// Show a table of the results
		//		summariseResults(stats);
		//
		//		// Optionally show histograms
		//		if (showHistograms)
		//		{
		//			IJ.showStatus("Calculating histograms ...");
		//
		//			int[] idList = new int[NAMES.length];
		//			int count = 0;
		//			double[] convert = getConversionFactors();
		//
		//			boolean requireRetile = false;
		//			for (int i = 0; i < NAMES.length; i++)
		//			{
		//				if (displayHistograms[i] && convert[i] != 0)
		//				{
		//					// We will have to convert the values...
		//					double[] tmp = ((StoredDataStatistics) stats[i]).getValues();
		//					for (int j = 0; j < tmp.length; j++)
		//						tmp[j] *= convert[i];
		//					StoredDataStatistics tmpStats = new StoredDataStatistics(tmp);
		//					idList[count++] = Utils.showHistogram(TITLE, tmpStats, NAMES[i], 0, 0, histogramBins,
		//							String.format("%s +/- %s", Utils.rounded(tmpStats.getMean()),
		//									Utils.rounded(tmpStats.getStandardDeviation())));
		//					requireRetile = requireRetile || Utils.isNewWindow();
		//				}
		//			}
		//
		//			if (count > 0 && requireRetile)
		//			{
		//				idList = Arrays.copyOf(idList, count);
		//				new WindowOrganiser().tileWindows(idList);
		//			}
		//		}
		//
		//		if (saveRawData)
		//		{
		//			String dir = Utils.getDirectory("Data_directory", rawDataDirectory);
		//			if (dir != null)
		//				saveData(stats, dir);
		//		}

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
		int molecules = 0;
		sb.append(molecules).append("\t");
		sb.append(Utils.rounded(simulationParameters.s)).append("\t");
		sb.append(Utils.rounded(simulationParameters.a)).append("\t");
		sb.append(Utils.rounded(getSa() * simulationParameters.a)).append("\t");
		sb.append(Utils.rounded(simulationParameters.gain)).append("\t");
		sb.append(Utils.rounded(simulationParameters.readNoise)).append("\t");
		sb.append(Utils.rounded(simulationParameters.b)).append("\t");
		sb.append(Utils.rounded(simulationParameters.b2)).append("\t");

		// Compute the noise
		double noise = simulationParameters.b2;
		if (simulationParameters.emCCD)
		{
			// The b2 parameter was computed without application of the EM-CCD noise factor of 2.
			//final double b2 = backgroundVariance + readVariance
			//                = simulationParameters.getBackground() + readVariance
			// This should be applied only to the background variance.
			final double readVariance = noise - simulationParameters.b2;
			noise = simulationParameters.b2 * 2 + readVariance;
		}

		sb.append(Utils.rounded(fitConfig.getInitialPeakStdDev0() * simulationParameters.a)).append("\t");
		sb.append(config.getRelativeFitting()).append("\t");
		sb.append(fitConfig.getFitFunction().toString());
		sb.append(":").append(PeakFit.getSolverName(fitConfig));
		if (fitConfig.getFitSolver() == FitSolver.MLE && fitConfig.isModelCamera())
		{
			sb.append(":Camera\t");

			// Add details of the noise model for the MLE
			sb.append("EM=").append(fitConfig.isEmCCD());
			sb.append(":G=").append(fitConfig.getGain());
			sb.append(":N=").append(fitConfig.getReadNoise());
		}
		else
			sb.append("\t");

		// Now output the actual results ...		
		sb.append("\t");

		summaryTable.append(sb.toString());
	}

	/**
	 * Get the factors to convert the ADUs/pixel units in the statistics array into calibrated photons and nm units. Set
	 * the conversion to zero if the function does not fit the specified statistic.
	 * 
	 * @return The conversion factors
	 */
	private double[] getConversionFactors()
	{
		final double[] convert = new double[NAMES.length];
		convert[Gaussian2DFunction.BACKGROUND] = (fitConfig.isBackgroundFitting()) ? 1 / simulationParameters.gain : 0;
		convert[Gaussian2DFunction.SIGNAL] = (fitConfig.isNotSignalFitting() &&
				fitConfig.getFitFunction() == FitFunction.FIXED) ? 0 : 1 / simulationParameters.gain;
		convert[Gaussian2DFunction.ANGLE] = (fitConfig.isAngleFitting()) ? 180.0 / Math.PI : 0;
		convert[Gaussian2DFunction.X_POSITION] = simulationParameters.a;
		convert[Gaussian2DFunction.Y_POSITION] = simulationParameters.a;
		convert[Gaussian2DFunction.X_SD] = (fitConfig.isWidth0Fitting()) ? simulationParameters.a : 0;
		convert[Gaussian2DFunction.Y_SD] = (fitConfig.isWidth1Fitting()) ? simulationParameters.a : 0;
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
		return x * simulationParameters.a;
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
			sb.append("\t").append(NAMES[i]).append("\t+/-");
		}
		return sb.toString();
	}

	private String createParameterHeader()
	{
		return "Molecules\tN\ts (nm)\ta (nm)\tsa (nm)\tX (nm)\tY (nm)\tZ (nm)\tGain\tReadNoise (ADUs)\tB (photons)\tb2 (photons)\tSNR\tLimit N\tLimit X\tLimit X ML\tRegion\tWidth\tMethod\tOptions";
	}

	private void runAnalysis()
	{
	}
}
