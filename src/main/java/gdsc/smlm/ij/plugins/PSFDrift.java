package gdsc.smlm.ij.plugins;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.fitting.FitConfiguration;
import gdsc.smlm.fitting.FitFunction;
import gdsc.smlm.fitting.FitSolver;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolver;
import gdsc.smlm.fitting.nonlinear.MaximumLikelihoodFitter.SearchMethod;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.PSFOffset;
import gdsc.smlm.ij.settings.PSFSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.model.ImagePSFModel;
import gdsc.smlm.utils.Maths;
import gdsc.smlm.utils.Statistics;
import gdsc.smlm.utils.XmlUtils;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.Plot2;
import ij.gui.PlotWindow;
import ij.plugin.PlugIn;
import ij.plugin.WindowOrganiser;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.Well19937c;

/**
 * Produces an drift curve for a PSF image using fitting.
 * <p>
 * The input images must be a z-stack of a PSF. These can be produced using the PSFCreator plugin.
 */
public class PSFDrift implements PlugIn
{
	private final static String TITLE = "PSF Drift";

	private static String title = "";
	private static boolean useOffset = false;
	private static double scale = 10;
	private static double zDepth = 1000;
	private static int gridSize = 10;
	private static double recallLimit = 0.25;
	private static int regionSize = 5;
	private static boolean backgroundFitting = false;
	private static boolean offsetFitting = true;
	private static double startOffset = 0.5;
	private static boolean comFitting = true;
	private static boolean useSampling = false;
	private static double photons = 1000;
	private static double photonLimit = 0.25;
	private static int positionsToAverage = 5;
	private static double smoothing = 0.1;

	private ImagePlus imp;
	private PSFSettings psfSettings;
	private static FitConfiguration fitConfig;

	private static final double bias = 500;
	static
	{
		// Initialise for fitting
		fitConfig = new FitConfiguration();
		// Set defaults for the MLE method
		fitConfig.setBias(bias);
		fitConfig.setEmCCD(false);
		fitConfig.setGain(1);
		fitConfig.setReadNoise(0);
		fitConfig.setModelCamera(false);
		fitConfig.setSearchMethod(SearchMethod.POWELL);
	}

	private int centrePixel;
	private int total;
	private double[][] results;

	private int[] idList = new int[3];
	private int idCount = 0;

	private class Job
	{
		final int z;
		final double cx, cy;
		final int index;

		public Job(int z, double cx, double cy, int index)
		{
			this.z = z;
			this.cx = cx;
			this.cy = cy;
			this.index = index;
		}

		public Job()
		{
			this(0, 0, 0, -1);
		}

		@Override
		public String toString()
		{
			return String.format("z=%d, cx=%.2f, cy=%.2f", z, cx, cy);
		}
	}

	/**
	 * Used to allow multi-threading of the fitting method
	 */
	private class Worker implements Runnable
	{
		volatile boolean finished = false;
		final ImagePSFModel psf;
		final BlockingQueue<Job> jobs;
		final FitConfiguration fitConfig;
		final double s, a;
		final double[][] xy;
		final int w;
		final int w2;
		final RandomDataGenerator random;

		private double[] lb, ub = null;
		private double[] lc, uc = null;

		public Worker(BlockingQueue<Job> jobs, ImagePSFModel psf, int width, FitConfiguration fitConfig)
		{
			this.jobs = jobs;
			this.psf = psf.copy();
			this.fitConfig = fitConfig.clone();
			s = fitConfig.getInitialPeakStdDev0();
			a = psfSettings.nmPerPixel * scale;
			xy = PSFDrift.getStartPoints(PSFDrift.this);
			w = width;
			w2 = w * w;
			if (useSampling)
				random = new RandomDataGenerator(new Well19937c());
			else
				random = null;

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
					Job job = jobs.take();
					if (job == null || job.index < 0 || finished)
						break;
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

		private void run(Job job)
		{
			if (IJ.escapePressed())
			{
				finished = true;
				return;
			}

			double cx = centrePixel + job.cx;
			double cy = centrePixel + job.cy;

			// Draw the PSF
			double[] data = new double[w2];
			// Fitting is done when there is a bias
			Arrays.fill(data, bias);
			if (useSampling)
			{
				int p = (int) random.nextPoisson(photons);
				psf.sample3D(data, w, w, p, cx, cy, job.z);
			}
			else
				psf.create3D(data, w, w, photons, cx, cy, job.z, false);

			//Utils.display("Data", data, w, w);

			// Fit the PSF. Do this from different start positions.

			// Get the background and signal estimate
			final double b = (backgroundFitting) ? BenchmarkFit.getBackground(data, w, w) : bias;
			final double signal = BenchmarkFit.getSignal(data, b);

			if (comFitting)
			{
				// Get centre-of-mass estimate, then subtract the centre that will be added later 
				BenchmarkFit.getCentreOfMass(data, w, w, xy[xy.length - 1]);
				xy[xy.length - 1][0] -= cx;
				xy[xy.length - 1][1] -= cy;
			}

			double[] initialParams = new double[7];
			initialParams[Gaussian2DFunction.BACKGROUND] = b;
			initialParams[Gaussian2DFunction.SIGNAL] = signal;
			initialParams[Gaussian2DFunction.X_SD] = initialParams[Gaussian2DFunction.Y_SD] = s;

			// Subtract the bias
			if (fitConfig.isRemoveBiasBeforeFitting())
			{
				initialParams[Gaussian2DFunction.BACKGROUND] = Math.max(0,
						initialParams[Gaussian2DFunction.BACKGROUND] - bias);
				for (int i = 0; i < data.length; i++)
					data[i] -= bias;
			}

			double[] error = new double[1];
			int resultPosition = job.index;
			for (double[] centre : xy)
			{
				// Do fitting
				final double[] params = initialParams.clone();
				params[Gaussian2DFunction.X_POSITION] = cx + centre[0];
				params[Gaussian2DFunction.Y_POSITION] = cy + centre[1];
				fitConfig.initialise(1, w, params);
				FunctionSolver solver = fitConfig.getFunctionSolver();
				if (solver.isBounded())
					setBounds(solver);
				if (solver.isConstrained())
					setConstraints(solver);
				final FitStatus status = solver.fit(data.length, data, null, params, null, error, 0);
				// Subtract the fitted bias from the background
				if (!fitConfig.isRemoveBiasBeforeFitting())
					params[Gaussian2DFunction.BACKGROUND] -= bias;
				// Account for 0.5 pixel offset during fitting
				params[Gaussian2DFunction.X_POSITION] += 0.5;
				params[Gaussian2DFunction.Y_POSITION] += 0.5;
				if (isValid(status, params, w))
				{
					// XXX Decide what results are needed for analysis
					// Store all the results for later analysis
					//results[resultPosition] = params;
					// Store only the drift
					results[resultPosition] = new double[] { a * (params[Gaussian2DFunction.X_POSITION] - cx),
							a * (params[Gaussian2DFunction.Y_POSITION] - cy), job.z };
					//System.out.printf("Fit " + job + ". %f,%f\n", results[resultPosition][0],
					//		results[resultPosition][1]);
				}
				else
				{
					//System.out.println("Failed to fit " + job + ". " + status);
				}
				resultPosition += total;
			}
		}

		private void setBounds(FunctionSolver solver)
		{
			solver.setBounds(lb, ub);
		}

		private void createBounds()
		{
			if (ub == null)
			{
				ub = new double[7];
				lb = new double[7];

				// Background could be zero so always have an upper limit
				ub[Gaussian2DFunction.BACKGROUND] = 1;
				lb[Gaussian2DFunction.SIGNAL] = photons * photonLimit;
				ub[Gaussian2DFunction.SIGNAL] = photons * 2;
				ub[Gaussian2DFunction.X_POSITION] = w;
				ub[Gaussian2DFunction.Y_POSITION] = w;
				lb[Gaussian2DFunction.ANGLE] = -Math.PI;
				ub[Gaussian2DFunction.ANGLE] = Math.PI;
				double wf = 1.5;
				lb[Gaussian2DFunction.X_SD] = s / wf;
				ub[Gaussian2DFunction.X_SD] = s * 5;
				lb[Gaussian2DFunction.Y_SD] = s / wf;
				ub[Gaussian2DFunction.Y_SD] = s * 5;
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
				//System.out.println("Failed to fit: " + status);
				return false;
			}

			// Reject fits that are outside the bounds of the data
			if (params[Gaussian2DFunction.X_POSITION] < 0 || params[Gaussian2DFunction.Y_POSITION] < 0 ||
					params[Gaussian2DFunction.X_POSITION] > size || params[Gaussian2DFunction.Y_POSITION] > size)
			{
				//System.out.printf("Failed to fit position: x=%f,y=%f\n", params[Gaussian2DFunction.X_POSITION],
				//		params[Gaussian2DFunction.Y_POSITION]);
				return false;
			}

			// Reject fits that do not correctly estimate the signal
			if (params[Gaussian2DFunction.SIGNAL] < lb[Gaussian2DFunction.SIGNAL] ||
					params[Gaussian2DFunction.SIGNAL] > ub[Gaussian2DFunction.SIGNAL])
			{
				//System.out.printf("Failed to fit signal: %f\n", params[Gaussian2DFunction.SIGNAL]);
				return false;
			}

			// Reject fits that have a background too far from zero
			// TODO - configure this better
			//if (params[Gaussian2DFunction.BACKGROUND] < -10 || params[Gaussian2DFunction.BACKGROUND] > 10)
			//{
			//	return false;
			//}

			// Q. Should we do width bounds checking?
			if (fitConfig.isWidth0Fitting())
			{
				if (params[Gaussian2DFunction.X_SD] < lb[Gaussian2DFunction.X_SD] ||
						params[Gaussian2DFunction.X_SD] > ub[Gaussian2DFunction.X_SD])
				{
					//System.out.printf("Failed to fit x-width: %f\n", params[Gaussian2DFunction.X_SD]);
					return false;
				}
			}
			if (fitConfig.isWidth1Fitting())
			{
				if (params[Gaussian2DFunction.Y_SD] < lb[Gaussian2DFunction.Y_SD] ||
						params[Gaussian2DFunction.Y_SD] > ub[Gaussian2DFunction.Y_SD])
				{
					//System.out.printf("Failed to fit y-width: %f\n", params[Gaussian2DFunction.Y_SD]);
					return false;
				}
			}

			return true;
		}

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		// Build a list of suitable images
		List<String> titles = createImageList();

		if (titles.isEmpty())
		{
			IJ.error(TITLE, "No suitable PSF images");
			return;
		}

		GenericDialog gd = new GenericDialog(TITLE);
		gd.addMessage("Select the input PSF image");
		gd.addChoice("PSF", titles.toArray(new String[titles.size()]), title);
		gd.addCheckbox("Use_offset", useOffset);
		gd.addNumericField("Scale", scale, 2);
		gd.addNumericField("z_depth", zDepth, 2, 6, "nm");
		gd.addNumericField("Grid_size", gridSize, 0);
		gd.addSlider("Recall_limit", 0.01, 1, recallLimit);

		gd.addSlider("Region_size", 2, 20, regionSize);
		gd.addCheckbox("Background_fitting", backgroundFitting);
		String[] solverNames = SettingsManager.getNames((Object[]) FitSolver.values());
		gd.addChoice("Fit_solver", solverNames, solverNames[fitConfig.getFitSolver().ordinal()]);
		String[] functionNames = SettingsManager.getNames((Object[]) FitFunction.values());
		gd.addChoice("Fit_function", functionNames, functionNames[fitConfig.getFitFunction().ordinal()]);
		gd.addCheckbox("Offset_fit", offsetFitting);
		gd.addNumericField("Start_offset", startOffset, 3);
		gd.addCheckbox("Include_CoM_fit", comFitting);
		gd.addCheckbox("Use_sampling", useSampling);
		gd.addNumericField("Photons", photons, 0);
		gd.addSlider("Photon_limit", 0, 1, photonLimit);
		gd.addSlider("Smoothing", 0, 0.5, smoothing);

		gd.showDialog();
		if (gd.wasCanceled())
			return;

		title = gd.getNextChoice();
		useOffset = gd.getNextBoolean();
		scale = gd.getNextNumber();
		zDepth = gd.getNextNumber();
		gridSize = (int) gd.getNextNumber();
		recallLimit = gd.getNextNumber();
		regionSize = (int) Math.abs(gd.getNextNumber());
		backgroundFitting = gd.getNextBoolean();
		fitConfig.setFitSolver(gd.getNextChoiceIndex());
		fitConfig.setFitFunction(gd.getNextChoiceIndex());
		offsetFitting = gd.getNextBoolean();
		startOffset = Math.abs(gd.getNextNumber());
		comFitting = gd.getNextBoolean();
		useSampling = gd.getNextBoolean();
		photons = Math.abs(gd.getNextNumber());
		photonLimit = Math.abs(gd.getNextNumber());
		smoothing = Math.abs(gd.getNextNumber());

		if (!comFitting && !offsetFitting)
		{
			IJ.error(TITLE, "No initial fitting positions");
			return;
		}

		if (regionSize < 1)
			regionSize = 1;

		if (gd.invalidNumber())
			return;

		GlobalSettings settings = new GlobalSettings();
		settings.setFitEngineConfiguration(new FitEngineConfiguration(fitConfig));
		if (!PeakFit.configureFitSolver(settings, null, false, true))
			return;

		imp = WindowManager.getImage(title);
		if (imp == null)
		{
			IJ.error(TITLE, "No PSF image for image: " + title);
			return;
		}
		psfSettings = getPSFSettings(imp);
		if (psfSettings == null)
		{
			IJ.error(TITLE, "No PSF settings for image: " + title);
			return;
		}

		computeDrift();
	}

	private void computeDrift()
	{
		// Create a grid of XY offset positions between 0-1 for PSF insert
		final double[] grid = new double[gridSize];
		for (int i = 0; i < grid.length; i++)
			grid[i] = (double) i / gridSize;

		// Configure fitting region
		final int w = 2 * regionSize + 1;
		centrePixel = w / 2;

		// Check region size using the image PSF
		double newPsfWidth = (double) imp.getWidth() / scale;
		if (Math.ceil(newPsfWidth) > w)
			Utils.log(TITLE + ": Region size %d is larger then the scaled PSF %.1f", w, newPsfWidth);

		// Create robust PSF fitting settings
		final double a = psfSettings.nmPerPixel * scale;
		final double sa = PSFCalculator.squarePixelAdjustment(psfSettings.nmPerPixel *
				(psfSettings.fwhm / PSFCalculator.SD_TO_FWHM_FACTOR), a);
		fitConfig.setInitialPeakStdDev(sa / a);
		fitConfig.setBackgroundFitting(backgroundFitting);
		fitConfig.setNotSignalFitting(false);
		fitConfig.setComputeDeviations(false);
		fitConfig.setFitValidation(false);

		// Create the PSF over the desired z-depth
		int depth = (int) Math.round(zDepth / psfSettings.nmPerSlice);
		int startSlice = psfSettings.zCentre - depth;
		int endSlice = psfSettings.zCentre + depth;
		int nSlices = imp.getStackSize();
		startSlice = (startSlice < 1) ? 1 : (startSlice > nSlices) ? nSlices : startSlice;
		endSlice = (endSlice < 1) ? 1 : (endSlice > nSlices) ? nSlices : endSlice;

		ImagePSFModel psf = createImagePSF(startSlice, endSlice);

		int minz = startSlice - psfSettings.zCentre;
		int maxz = endSlice - psfSettings.zCentre;

		final int nZ = maxz - minz + 1;
		final int gridSize2 = grid.length * grid.length;
		total = nZ * gridSize2;

		// Store all the fitting results
		int nStartPoints = getNumberOfStartPoints();
		results = new double[total * nStartPoints][];

		// TODO - Add ability to iterate this, adjusting the current offset in the PSF
		// each iteration

		// Create a pool of workers
		int nThreads = Prefs.getThreads();
		BlockingQueue<Job> jobs = new ArrayBlockingQueue<Job>(nThreads * 2);
		List<Worker> workers = new LinkedList<Worker>();
		List<Thread> threads = new LinkedList<Thread>();
		for (int i = 0; i < nThreads; i++)
		{
			Worker worker = new Worker(jobs, psf, w, fitConfig);
			Thread t = new Thread(worker);
			workers.add(worker);
			threads.add(t);
			t.start();
		}

		// Fit 
		Utils.showStatus("Fitting ...");
		final int step = (total > 400) ? total / 200 : 2;
		outer: for (int z = minz, i = 0; z <= maxz; z++)
		{
			for (int x = 0; x < grid.length; x++)
				for (int y = 0; y < grid.length; y++, i++)
				{
					if (IJ.escapePressed())
					{
						break outer;
					}
					put(jobs, new Job(z, grid[x], grid[y], i));
					if (i % step == 0)
					{
						IJ.showProgress(i, total);
					}
				}
		}

		// If escaped pressed then do not need to stop the workers, just return
		if (Utils.isInterrupted())
		{
			IJ.showProgress(1);
			return;
		}

		// Finish all the worker threads by passing in a null job
		for (int i = 0; i < threads.size(); i++)
		{
			put(jobs, new Job());
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
		IJ.showStatus("");

		// Plot the average and SE for the drift curve
		// Plot the recall
		double[] zPosition = new double[nZ];
		double[] avX = new double[nZ];
		double[] seX = new double[nZ];
		double[] avY = new double[nZ];
		double[] seY = new double[nZ];
		double[] recall = new double[nZ];
		for (int z = minz, i = 0; z <= maxz; z++, i++)
		{
			Statistics statsX = new Statistics();
			Statistics statsY = new Statistics();
			for (int s = 0; s < nStartPoints; s++)
			{
				int resultPosition = i * gridSize2 + s * total;
				final int endResultPosition = resultPosition + gridSize2;
				while (resultPosition < endResultPosition)
				{
					if (results[resultPosition] != null)
					{
						statsX.add(results[resultPosition][0]);
						statsY.add(results[resultPosition][1]);
					}
					resultPosition++;
				}
			}
			zPosition[i] = z * psfSettings.nmPerSlice;
			avX[i] = statsX.getMean();
			seX[i] = statsX.getStandardError();
			avY[i] = statsY.getMean();
			seY[i] = statsY.getStandardError();
			recall[i] = (double) statsX.getN() / (nStartPoints * gridSize2);
		}

		// Find the range from the z-centre above the recall limit 
		int centre = 0;
		for (int slice = startSlice, i = 0; slice <= endSlice; slice++, i++)
		{
			if (slice == psfSettings.zCentre)
			{
				centre = i;
				break;
			}
		}
		if (recall[centre] < recallLimit)
			return;
		int start = centre, end = centre;
		for (int i = centre; i-- > 0;)
		{
			if (recall[i] < recallLimit)
				break;
			start = i;
		}
		for (int i = centre; ++i < recall.length;)
		{
			if (recall[i] < recallLimit)
				break;
			end = i;
		}

		int iterations = 1;
		LoessInterpolator loess = null;
		if (smoothing > 0)
			loess = new LoessInterpolator(smoothing, iterations);

		double[][] smoothx = displayPlot("Drift X", "X (nm)", zPosition, avX, seX, loess, start, end);
		double[][] smoothy = displayPlot("Drift Y", "Y (nm)", zPosition, avY, seY, loess, start, end);
		displayPlot("Recall", "Recall", zPosition, recall, null, null, start, end);

		WindowOrganiser wo = new WindowOrganiser();
		wo.tileWindows(idList);

		// Ask the user if they would like to store them in the image
		GenericDialog gd = new GenericDialog(TITLE);
		gd.enableYesNoCancel();
		gd.hideCancelButton();
		startSlice = psfSettings.zCentre - (centre - start);
		endSlice = psfSettings.zCentre + (end - centre);
		gd.addMessage(String.format("Save the drift to the PSF?\n \nSlices %d (%s nm) - %d (%s nm) above recall limit",
				startSlice, Utils.rounded(zPosition[start]), endSlice, Utils.rounded(zPosition[end])));
		gd.addMessage("Optionally average the end points to set drift outside the limits.\n(Select zero to ignore)");
		gd.addSlider("Number_of_points", 0, 10, positionsToAverage);
		gd.showDialog();
		if (gd.wasOKed())
		{
			positionsToAverage = Math.abs((int) gd.getNextNumber());
			ArrayList<PSFOffset> offset = new ArrayList<PSFOffset>();
			final double pitch = psfSettings.nmPerPixel;
			int j = 0, jj = 0;
			for (int i = start, slice = startSlice; i <= end; slice++, i++)
			{
				j = findCentre(zPosition[i], smoothx, j);
				if (j == -1)
				{
					Utils.log("Failed to find the offset for depth %.2f", zPosition[i]);
					continue;
				}
				// The offset should store the difference to the centre in pixels so divide by the pixel pitch
				double cx = smoothx[1][j] / pitch;
				double cy = smoothy[1][j] / pitch;
				jj = findOffset(slice, jj);
				if (jj != -1)
				{
					cx += psfSettings.offset[jj].cx;
					cy += psfSettings.offset[jj].cy;
				}
				offset.add(new PSFOffset(slice, cx, cy));
			}
			addMissingOffsets(startSlice, endSlice, nSlices, offset);
			psfSettings.offset = offset.toArray(new PSFOffset[offset.size()]);
			psfSettings.addNote(TITLE,
					String.format("Solver=%s, Region=%d", PeakFit.getSolverName(fitConfig), regionSize));
			imp.setProperty("Info", XmlUtils.toXML(psfSettings));
		}
	}

	private int findCentre(double d, double[][] smoothx, int i)
	{
		while (i < smoothx[0].length)
		{
			if (smoothx[0][i] == d)
				return i;
			i++;
		}
		return -1;
	}

	private int findOffset(int slice, int i)
	{
		if (useOffset)
		{
			while (i < psfSettings.offset.length)
			{
				if (psfSettings.offset[i].slice == slice)
					return i;
				i++;
			}
		}
		return -1;
	}

	private void addMissingOffsets(int startSlice, int endSlice, int nSlices, ArrayList<PSFOffset> offset)
	{
		// Add an offset for the remaining slices 
		if (positionsToAverage > 0)
		{
			double cx = 0, cy = 0;
			int n = 0;
			for (int i = 0; n < positionsToAverage && i < offset.size(); i++)
			{
				cx += offset.get(i).cx;
				cy += offset.get(i).cy;
				n++;
			}
			cx /= n;
			cy /= n;
			double cx2 = 0, cy2 = 0;
			double n2 = 0;
			for (int i = offset.size(); n2 < positionsToAverage && i-- > 0;)
			{
				cx2 += offset.get(i).cx;
				cy2 += offset.get(i).cy;
				n2++;
			}
			cx2 /= n2;
			cy2 /= n2;

			for (int slice = 1; slice < startSlice; slice++)
				offset.add(new PSFOffset(slice, cx, cy));
			for (int slice = endSlice + 1; slice <= nSlices; slice++)
				offset.add(new PSFOffset(slice, cx2, cy2));
			Collections.sort(offset, new Comparator<PSFOffset>()
			{
				public int compare(PSFOffset arg0, PSFOffset arg1)
				{
					return arg0.slice - arg1.slice;
				}
			});
		}
	}

	private double[][] displayPlot(String title, String yLabel, double[] x, double[] y, double[] se,
			LoessInterpolator loess, int start, int end)
	{
		// Extract non NaN numbers
		double[] newX = new double[x.length];
		double[] newY = new double[x.length];
		int c = 0;
		for (int i = 0; i < x.length; i++)
			if (!Double.isNaN(y[i]))
			{
				newX[c] = x[i];
				newY[c] = y[i];
				c++;
			}
		newX = Arrays.copyOf(newX, c);
		newY = Arrays.copyOf(newY, c);

		title = TITLE + " " + title;
		Plot2 plot = new Plot2(title, "z (nm)", yLabel);
		double[] limitsx = Maths.limits(x);
		double[] limitsy = new double[2];
		if (se != null)
		{
			if (c > 0)
			{
				limitsy = new double[] { newY[0] - se[0], newY[0] + se[0] };
				for (int i = 1; i < newY.length; i++)
				{
					limitsy[0] = Maths.min(limitsy[0], newY[i] - se[i]);
					limitsy[1] = Maths.max(limitsy[1], newY[i] + se[i]);
				}
			}
		}
		else
		{
			if (c > 0)
				limitsy = Maths.limits(newY);
		}
		double rangex = Math.max(0.05 * (limitsx[1] - limitsx[0]), 0.1);
		double rangey = Math.max(0.05 * (limitsy[1] - limitsy[0]), 0.1);
		plot.setLimits(limitsx[0] - rangex, limitsx[1] + rangex, limitsy[0] - rangey, limitsy[1] + rangey);

		if (loess == null)
		{
			addPoints(plot, Plot.LINE, newX, newY, x[start], x[end]);
		}
		else
		{
			addPoints(plot, Plot.DOT, newX, newY, x[start], x[end]);
			newY = loess.smooth(newX, newY);
			addPoints(plot, Plot.LINE, newX, newY, x[start], x[end]);
		}
		if (se != null)
		{
			plot.setColor(Color.magenta);
			for (int i = 0; i < x.length; i++)
			{
				if (!Double.isNaN(y[i]))
					plot.drawLine(x[i], y[i] - se[i], x[i], y[i] + se[i]);
			}

			// Draw the start and end lines for the valid range
			plot.setColor(Color.green);
			plot.drawLine(x[start], limitsy[0], x[start], limitsy[1]);
			plot.drawLine(x[end], limitsy[0], x[end], limitsy[1]);
		}
		else
		{
			// draw a line for the recall limit
			plot.setColor(Color.magenta);
			plot.drawLine(limitsx[0] - rangex, recallLimit, limitsx[1] + rangex, recallLimit);
		}
		PlotWindow pw = Utils.display(title, plot);
		if (Utils.isNewWindow())
			idList[idCount++] = pw.getImagePlus().getID();

		return new double[][] { newX, newY };
	}

	private void addPoints(Plot plot, int shape, double[] x, double[] y, double lower, double upper)
	{
		if (x.length == 0)
			return;
		// Split the line into three:
		// 1. All points up to and including lower
		// 2. All points between lower and upper inclusive
		// 3. All point from upper upwards

		// Plot the main curve first 
		addPoints(plot, shape, x, y, lower, upper, Color.blue);
		// Then plot the others
		addPoints(plot, shape, x, y, x[0], lower, Color.red);
		addPoints(plot, shape, x, y, upper, x[x.length - 1], Color.red);
	}

	private void addPoints(Plot plot, int shape, double[] x, double[] y, double lower, double upper, Color color)
	{
		double[] x2 = new double[x.length];
		double[] y2 = new double[y.length];
		int c = 0;
		for (int i = 0; i < x.length; i++)
		{
			if (x[i] >= lower && x[i] <= upper)
			{
				x2[c] = x[i];
				y2[c] = y[i];
				c++;
			}
		}
		if (c == 0)
			return;
		x2 = Arrays.copyOf(x2, c);
		y2 = Arrays.copyOf(y2, c);
		plot.setColor(color);
		plot.addPoints(x2, y2, shape);
	}

	private ImagePSFModel createImagePSF(int lower, int upper)
	{
		int zCentre = psfSettings.zCentre;

		final double unitsPerPixel = 1.0 / scale;
		final double unitsPerSlice = 1; // So we can move from -depth to depth

		// Extract data uses index not slice number as arguments so subtract 1
		double noiseFraction = 1e-3;
		ImagePSFModel model = new ImagePSFModel(CreateData.extractImageStack(imp, lower - 1, upper - 1), zCentre -
				lower, unitsPerPixel, unitsPerSlice, psfSettings.fwhm, noiseFraction);

		// Add the calibrated centres
		if (psfSettings.offset != null && useOffset)
		{
			int sliceOffset = lower;
			for (PSFOffset offset : psfSettings.offset)
			{
				model.setRelativeCentre(offset.slice - sliceOffset, offset.cx, offset.cy);
			}
		}

		return model;
	}

	private void put(BlockingQueue<Job> jobs, Job job)
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
				xy[ii++] = new double[] { 0, 0 };
			}
			else
			{
				// Fit using region surrounding the point. Use -1,-1 : -1:1 : 1,-1 : 1,1 directions at 
				// startOffset pixels total distance 
				final double distance = Math.sqrt(startOffset * startOffset * 0.5);

				for (int x = -1; x <= 1; x += 2)
					for (int y = -1; y <= 1; y += 2)
					{
						xy[ii++] = new double[] { x * distance, y * distance };
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

	public static List<String> createImageList()
	{
		List<String> titles = new LinkedList<String>();
		int[] ids = WindowManager.getIDList();
		if (ids != null)
		{
			for (int id : ids)
			{
				ImagePlus imp = WindowManager.getImage(id);
				if (imp != null)
				{
					// Image must be greyscale
					if (imp.getType() == ImagePlus.GRAY8 || imp.getType() == ImagePlus.GRAY16 ||
							imp.getType() == ImagePlus.GRAY32)
					{
						// Image must be square and a stack of a single channel
						if (imp.getWidth() == imp.getHeight() && imp.getNChannels() == 1)
						{
							// Check if these are PSF images created by the SMLM plugins
							PSFSettings psfSettings = getPSFSettings(imp);
							if (psfSettings != null)
							{
								if (psfSettings.zCentre <= 0)
								{
									Utils.log(TITLE + ": Unknown PSF z-centre setting for image: " + imp.getTitle());
									continue;
								}
								if (psfSettings.nmPerPixel <= 0)
								{
									Utils.log(TITLE + ": Unknown PSF nm/pixel setting for image: " + imp.getTitle());
									continue;
								}
								if (psfSettings.nmPerSlice <= 0)
								{
									Utils.log(TITLE + ": Unknown PSF nm/slice setting for image: " + imp.getTitle());
									continue;
								}
								if (psfSettings.fwhm <= 0)
								{
									Utils.log(TITLE + ": Unknown PSF FWHM setting for image: " + imp.getTitle());
									continue;
								}

								titles.add(imp.getTitle());
							}
						}
					}
				}
			}
		}
		return titles;
	}

	private static PSFSettings getPSFSettings(ImagePlus imp)
	{
		Object info = imp.getProperty("Info");
		if (info != null)
		{
			Object o = XmlUtils.fromXML(info.toString());
			if (o != null && o instanceof PSFSettings)
			{
				return (PSFSettings) o;
			}
		}
		return null;
	}

	public static double[][] getStartPoints(PSFDrift psfDrift)
	{
		return psfDrift.getStartPoints();
	}
}
