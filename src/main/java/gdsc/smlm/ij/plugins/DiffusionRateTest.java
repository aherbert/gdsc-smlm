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

import gdsc.smlm.fitting.JumpDistanceAnalysis;
import gdsc.smlm.fitting.JumpDistanceAnalysis.JumpDistanceFunction;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.ij.settings.CreateDataSettings;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.model.ImageModel;
import gdsc.smlm.model.MoleculeModel;
import gdsc.smlm.model.SphericalDistribution;
import gdsc.smlm.results.Calibration;
import gdsc.smlm.results.ExtendedPeakResult;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.utils.Maths;
import gdsc.smlm.utils.RollingArray;
import gdsc.smlm.utils.Statistics;
import gdsc.smlm.utils.StoredDataStatistics;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.Plot2;
import ij.gui.PlotWindow;
import ij.plugin.LutLoader;
import ij.plugin.PlugIn;
import ij.plugin.WindowOrganiser;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij.text.TextWindow;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction.Parametric;
import org.apache.commons.math3.fitting.CurveFitter;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

/**
 * Move a set of molecules and calculates the diffusion rate. Uses settings from the CreateData
 * plugin so that the diffusion should be equivalent.
 * 
 */
public class DiffusionRateTest implements PlugIn
{
	private static final String TITLE = "Diffusion Rate Test";
	private static TextWindow msdTable = null;

	// Used to allow other plugins to detect if a dataset is simulated
	static double lastSimulatedPrecision = 0;
	static String[] lastSimulatedDataset = new String[2];

	static boolean isSimulated(String name)
	{
		for (String name2 : lastSimulatedDataset)
			if (name.equals(name2))
				return true;
		return false;
	}

	/**
	 * Used to aggregate points into results
	 */
	public class Point
	{
		public int id;
		public double x, y;

		/**
		 * Create a cluster point
		 * 
		 * @param id
		 * @param x
		 * @param y
		 */
		public Point(int id, double x, double y)
		{
			this.id = id;
			this.x = x;
			this.y = y;
		}

		public Point(int id, double[] xyz)
		{
			this(id, xyz[0], xyz[1]);
		}

		public double distance2(Point p)
		{
			final double dx = x - p.x;
			final double dy = y - p.y;
			return dx * dx + dy * dy;
		}

		public double distance2(Point p, double error, RandomGenerator rand)
		{
			final double dx = (x + rand.nextGaussian() * error) - (p.x + rand.nextGaussian() * error);
			final double dy = (y + rand.nextGaussian() * error) - (p.y + rand.nextGaussian() * error);
			return dx * dx + dy * dy;
		}
	}

	private CreateDataSettings settings;
	private static boolean useConfinement = false;
	private static int confinementAttempts = 5;
	private static int fitN = 20;
	private static boolean showDiffusionExample = false;
	private static double magnification = 5;
	private static int aggregateSteps = 10;
	private static int msdAnalysisSteps = 0;
	private static double precision = 0;
	private int myAggregateSteps = 0;
	private int myMsdAnalysisSteps = 0;
	private boolean extraOptions = false;
	private double conversionFactor;
	private double myPrecision = 0;

	private int[] idList = new int[4];
	private int idCount = 0;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		PluginTracker.recordPlugin(this.getClass(), arg);
		
		extraOptions = Utils.isExtraOptions();

		if (!showDialog())
			return;

		lastSimulatedDataset[0] = lastSimulatedDataset[1] = "";
		lastSimulatedPrecision = 0;

		final int totalSteps = (int) Math.ceil(settings.seconds * settings.stepsPerSecond);

		conversionFactor = 1000000.0 / (settings.pixelPitch * settings.pixelPitch);

		// Diffusion rate is um^2/sec. Convert to pixels per simulation frame.
		final double diffusionRateInPixelsPerSecond = settings.diffusionRate * conversionFactor;
		final double diffusionRateInPixelsPerStep = diffusionRateInPixelsPerSecond / settings.stepsPerSecond;
		final double precisionInPixels = myPrecision / settings.pixelPitch;
		final boolean addError = myPrecision != 0;

		Utils.log(TITLE + " : D = %s um^2/sec, Precision = %s nm", Utils.rounded(settings.diffusionRate, 4),
				Utils.rounded(myPrecision, 4));
		Utils.log("Mean-displacement per dimension = %s nm/sec",
				Utils.rounded(1e3 * ImageModel.getRandomMoveDistance(settings.diffusionRate), 4));
		if (extraOptions)
			Utils.log("Step size = %s, precision = %s",
					Utils.rounded(ImageModel.getRandomMoveDistance(diffusionRateInPixelsPerStep)),
					Utils.rounded(precisionInPixels));

		// Convert diffusion co-efficient into the standard deviation for the random walk
		final double diffusionSigma = ImageModel.getRandomMoveDistance(diffusionRateInPixelsPerStep);
		Utils.log("Simulation step-size = %s nm", Utils.rounded(settings.pixelPitch * diffusionSigma, 4));

		// Move the molecules and get the diffusion rate
		IJ.showStatus("Simulating ...");
		final long start = System.nanoTime();
		final long seed = System.currentTimeMillis() + System.identityHashCode(this);
		RandomGenerator[] random = new RandomGenerator[3];
		RandomGenerator[] random2 = new RandomGenerator[3];
		for (int i = 0; i < 3; i++)
		{
			random[i] = new Well19937c(seed + i);
			random2[i] = new Well19937c(seed + i + 3);
		}
		Statistics[] stats2D = new Statistics[totalSteps];
		Statistics[] stats3D = new Statistics[totalSteps];
		StoredDataStatistics jumpDistances = new StoredDataStatistics(totalSteps);
		for (int j = 0; j < totalSteps; j++)
		{
			stats2D[j] = new Statistics();
			stats3D[j] = new Statistics();
		}
		SphericalDistribution dist = new SphericalDistribution(settings.confinementRadius / settings.pixelPitch);
		Statistics asymptote = new Statistics();

		// Save results to memory
		MemoryPeakResults results = new MemoryPeakResults(totalSteps);
		Calibration cal = new Calibration(settings.pixelPitch, 1, 1000.0 / settings.stepsPerSecond);
		results.setCalibration(cal);
		results.setName(TITLE);
		int peak = 0;
		// Store raw coordinates
		ArrayList<Point> points = new ArrayList<Point>(totalSteps);

		for (int i = 0; i < settings.particles; i++)
		{
			if (i % 16 == 0)
			{
				IJ.showProgress(i, settings.particles);
				if (Utils.isInterrupted())
					return;
			}

			peak++; // Increment the frame so that tracing analysis can distinguish traces
			double[] origin = new double[3];
			final int id = i + 1;
			MoleculeModel m = new MoleculeModel(id, origin.clone());
			if (addError)
				origin = addError(origin, precisionInPixels, random);
			if (useConfinement)
			{
				// Note: When using confinement the average displacement should asymptote
				// at the average distance of a point from the centre of a ball. This is 3r/4.
				// See: http://answers.yahoo.com/question/index?qid=20090131162630AAMTUfM
				// The equivalent in 2D is 2r/3. However although we are plotting 2D distance
				// this is a projection of the 3D position onto the plane and so the particles
				// will not be evenly spread (there will be clustering at centre caused by the
				// poles)
				for (int j = 0; j < totalSteps; j++)
				{
					double[] xyz = m.getCoordinates();
					double[] originalXyz = xyz.clone();
					for (int n = confinementAttempts; n-- > 0;)
					{
						if (settings.useGridWalk)
							m.walk(diffusionSigma, random);
						else
							m.move(diffusionSigma, random);
						if (!dist.isWithin(m.getCoordinates()))
						{
							// Reset position
							for (int k = 0; k < 3; k++)
								xyz[k] = originalXyz[k];
						}
						else
						{
							// The move was allowed
							break;
						}
					}

					points.add(new Point(id, xyz));

					if (addError)
						xyz = addError(xyz, precisionInPixels, random2);

					peak = record(xyz, id, peak, stats2D[j], stats3D[j], jumpDistances, origin, results);
				}
				asymptote.add(distance(m.getCoordinates()));
			}
			else
			{
				if (settings.useGridWalk)
				{
					for (int j = 0; j < totalSteps; j++)
					{
						m.walk(diffusionSigma, random);
						double[] xyz = m.getCoordinates();
						points.add(new Point(id, xyz));
						if (addError)
							xyz = addError(xyz, precisionInPixels, random2);
						peak = record(xyz, id, peak, stats2D[j], stats3D[j], jumpDistances, origin, results);
					}
				}
				else
				{
					for (int j = 0; j < totalSteps; j++)
					{
						m.move(diffusionSigma, random);
						double[] xyz = m.getCoordinates();
						points.add(new Point(id, xyz));
						if (addError)
							xyz = addError(xyz, precisionInPixels, random2);
						peak = record(xyz, id, peak, stats2D[j], stats3D[j], jumpDistances, origin, results);
					}
				}
			}

			// Debug: record all the particles so they can be analysed
			// System.out.printf("%f %f %f\n", m.getX(), m.getY(), m.getZ());
		}
		final double time = (System.nanoTime() - start) / 1000000.0;
		IJ.showProgress(1);

		MemoryPeakResults.addResults(results);
		lastSimulatedDataset[0] = results.getName();
		lastSimulatedPrecision = myPrecision;

		// Convert pixels^2/step to um^2/sec
		final double msd = (jumpDistances.getMean() / conversionFactor) /
				(results.getCalibration().exposureTime / 1000);
		Utils.log("Raw data D=%s um^2/s, Precision = %s nm, N=%d, step=%s s, mean=%s um^2, MSD = %s um^2/s",
				Utils.rounded(settings.diffusionRate), Utils.rounded(myPrecision), jumpDistances.getN(),
				Utils.rounded(results.getCalibration().exposureTime / 1000),
				Utils.rounded(jumpDistances.getMean() / conversionFactor), Utils.rounded(msd));

		aggregateIntoFrames(points, addError, precisionInPixels, random2);

		IJ.showStatus("Analysing results ...");

		if (showDiffusionExample)
		{
			showExample(totalSteps, diffusionSigma, random);
		}

		// Plot a graph of mean squared distance
		double[] xValues = new double[stats2D.length];
		double[] yValues2D = new double[stats2D.length];
		double[] yValues3D = new double[stats3D.length];
		double[] upper = new double[stats2D.length];
		double[] lower = new double[stats2D.length];
		final CurveFitter<Parametric> fitter2D = new CurveFitter<Parametric>(new LevenbergMarquardtOptimizer());
		final CurveFitter<Parametric> fitter3D = new CurveFitter<Parametric>(new LevenbergMarquardtOptimizer());
		Statistics gradient2D = new Statistics();
		Statistics gradient3D = new Statistics();
		final int firstN = (useConfinement) ? fitN : totalSteps;
		for (int j = 0; j < totalSteps; j++)
		{
			// Convert steps to seconds
			xValues[j] = (double) (j + 1) / settings.stepsPerSecond;

			// Convert values in pixels^2 to um^2
			final double mean = stats2D[j].getMean() / conversionFactor;
			final double sd = stats2D[j].getStandardDeviation() / conversionFactor;
			yValues2D[j] = mean;
			yValues3D[j] = stats3D[j].getMean() / conversionFactor;
			upper[j] = mean + sd;
			lower[j] = mean - sd;

			if (j < firstN)
			{
				fitter2D.addObservedPoint(xValues[j], yValues2D[j]);
				gradient2D.add(yValues2D[j] / xValues[j]);
				fitter3D.addObservedPoint(xValues[j], yValues3D[j]);
				gradient3D.add(yValues3D[j] / xValues[j]);
			}
		}

		// TODO - Fit using the equation for 2D confined diffusion:
		// MSD = 4s^2 + R^2 (1 - 0.99e^(-1.84^2 Dt / R^2)
		// s = localisation precision
		// R = confinement radius
		// D = 2D diffusion coefficient
		// t = time

		final PolynomialFunction fitted2D, fitted3D;
		if (gradient2D.getMean() > 0)
		{
			// Do linear regression to get diffusion rate
			final double[] init2D = { 0, 1 / gradient2D.getMean() }; // a - b x
			final double[] best2D = fitter2D.fit(new PolynomialFunction.Parametric(), init2D);
			fitted2D = new PolynomialFunction(best2D);

			final double[] init3D = { 0, 1 / gradient3D.getMean() }; // a - b x
			final double[] best3D = fitter3D.fit(new PolynomialFunction.Parametric(), init3D);
			fitted3D = new PolynomialFunction(best3D);

			// For 2D diffusion: d^2 = 4D
			// where: d^2 = mean-square displacement

			double D = best2D[1] / 4.0;
			String msg = "2D Diffusion rate = " + Utils.rounded(D, 4) + " um^2 / sec (" + Utils.timeToString(time) +
					")";
			IJ.showStatus(msg);
			Utils.log(msg);

			D = best3D[1] / 6.0;
			Utils.log("3D Diffusion rate = " + Utils.rounded(D, 4) + " um^2 / sec (" + Utils.timeToString(time) + ")");
		}
		else
		{
			fitted2D = fitted3D = null;
		}

		// Create plot
		String title = TITLE;
		Plot2 plot = new Plot2(title, "Time (seconds)", "Mean-squared Distance (um^2)", xValues, yValues2D);
		double[] limits = Maths.limits(upper);
		limits = Maths.limits(limits, lower);
		limits = Maths.limits(limits, yValues3D);
		plot.setLimits(0, totalSteps / settings.stepsPerSecond, limits[0], limits[1]);
		plot.setColor(Color.blue);
		plot.addPoints(xValues, lower, Plot2.LINE);
		plot.addPoints(xValues, upper, Plot2.LINE);
		plot.setColor(Color.magenta);
		plot.addPoints(xValues, yValues3D, Plot2.LINE);
		if (fitted2D != null)
		{
			plot.setColor(Color.red);
			plot.addPoints(new double[] { xValues[0], xValues[xValues.length - 1] },
					new double[] { fitted2D.value(xValues[0]), fitted2D.value(xValues[xValues.length - 1]) },
					Plot2.LINE);
			plot.setColor(Color.green);
			plot.addPoints(new double[] { xValues[0], xValues[xValues.length - 1] },
					new double[] { fitted3D.value(xValues[0]), fitted3D.value(xValues[xValues.length - 1]) },
					Plot2.LINE);
		}
		plot.setColor(Color.black);

		PlotWindow pw1 = Utils.display(title, plot);
		if (Utils.isNewWindow())
			idList[idCount++] = pw1.getImagePlus().getID();

		// Show the cumulative jump distance plot
		final double factor = conversionFactor;
		double[] values = jumpDistances.getValues();
		for (int i = 0; i < values.length; i++)
			values[i] /= factor;
		title += " Cumulative Jump Distance";
		double[][] jdHistogram = JumpDistanceAnalysis.cumulativeHistogram(values);
		if (settings.useGridWalk)
		{
			// In this case with a large simulation size the the jumps are all 
			// the same distance so the histogram is a single step. Check the plot
			// range will be handled by ImageJ otherwise pad it out a but.
			double[] x = jdHistogram[0];
			double[] y = jdHistogram[1];
			if (x[x.length - 1] - x[0] < 0.01)
			{
				double[] x2 = new double[x.length + 3];
				double[] y2 = new double[y.length + 3];
				System.arraycopy(x, 0, x2, 2, x.length);
				System.arraycopy(y, 0, y2, 2, y.length);
				x2[0] = x[0] - 0.1;
				x2[1] = x[0];
				x2[x2.length - 1] = x[x.length - 1] + 0.1;
				y2[0] = 0;
				y2[1] = 0;
				y2[y2.length - 1] = 1;
				jdHistogram[0] = x2;
				jdHistogram[1] = y2;

				// Add some artificial points to allow the plot to be drawn
				values = Arrays.copyOf(values, values.length + 2);
				values[values.length - 2] = x[0] - 0.1;
				values[values.length - 1] = x[x.length - 1] + 0.1;
			}
		}
		Plot2 jdPlot = new Plot2(title, "Distance (um^2)", "Cumulative Probability", jdHistogram[0], jdHistogram[1]);
		PlotWindow pw2 = Utils.display(title, jdPlot);
		if (Utils.isNewWindow())
			idList[idCount++] = pw2.getImagePlus().getID();

		// Show a histogram of the distribution
		title = TITLE + " Jump";
		jumpDistances = new StoredDataStatistics(values);
		int plotId = Utils.showHistogram(title, jumpDistances, "Distance (um^2)", 0, 0, values.length / 1000);
		if (Utils.isNewWindow())
			idList[idCount++] = plotId;

		// Plot the expected function
		double max = Maths.max(values);
		double[] x = Utils.newArray(1000, 0, max / 1000);
		JumpDistanceAnalysis jd = new JumpDistanceAnalysis();
		double estimatedD = settings.diffusionRate / settings.stepsPerSecond;
		if (precision > 0)
			estimatedD += precision * precision / 1e6;
		JumpDistanceFunction fun = jd.new JumpDistanceFunction(x, estimatedD);
		double[] y = fun.evaluateAll(fun.guess());
		// Scale to have the same area
		if (Utils.xValues.length > 1)
		{
			final double area1 = jumpDistances.getN() * (Utils.xValues[1] - Utils.xValues[0]);
			final double area2 = Maths.sum(y) * (x[1] - x[0]);
			final double scale = area1 / area2;
			for (int i = 0; i < y.length; i++)
				y[i] *= scale;
		}
		jdPlot = Utils.plot;
		jdPlot.setColor(Color.red);
		jdPlot.addPoints(x, y, Plot.LINE);
		Utils.display(WindowManager.getImage(plotId).getTitle(), jdPlot);

		if (idCount > 0)
			new WindowOrganiser().tileWindows(idList);

		if (useConfinement)
			Utils.log("3D asymptote distance = %s nm (expected %.2f)",
					Utils.rounded(asymptote.getMean() * settings.pixelPitch, 4), 3 * settings.confinementRadius / 4);
	}

	/**
	 * Add a random Gaussian XY shift using the specified precision
	 * 
	 * @param xyz
	 * @param precision
	 * @param random
	 * @return The new xyz
	 */
	private double[] addError(double[] xyz, double precision, RandomGenerator[] random)
	{
		final double[] xy = xyz.clone();
		for (int i = 0; i < 2; i++)
		{
			final double shift = random[i].nextGaussian() * precision;
			xy[i] += shift;
		}
		return xy;
	}

	private int record(double[] xyz, int id, int peak, Statistics stats2D, Statistics stats3D,
			StoredDataStatistics jumpDistances, double[] origin, MemoryPeakResults results)
	{
		final double dx = xyz[0] - origin[0];
		final double dy = xyz[1] - origin[1];
		jumpDistances.add(dx * dx + dy * dy);
		for (int i = 0; i < 3; i++)
			origin[i] = xyz[i];

		stats2D.add(squared(xyz));
		stats3D.add(distance2(xyz));

		final float[] params = new float[7];
		params[Gaussian2DFunction.X_POSITION] = (float) xyz[0];
		params[Gaussian2DFunction.Y_POSITION] = (float) xyz[1];
		params[Gaussian2DFunction.SIGNAL] = 10f;
		params[Gaussian2DFunction.X_SD] = params[Gaussian2DFunction.Y_SD] = 1;
		final float noise = 0.1f;
		results.add(new ExtendedPeakResult(peak, (int) params[Gaussian2DFunction.X_POSITION],
				(int) params[Gaussian2DFunction.Y_POSITION], 10, 0, noise, params, null, peak, id));
		return ++peak;
	}

	/**
	 * Get the squared distance from the origin in 2D (using XY coordinates)
	 * 
	 * @param coordinates
	 * @return
	 */
	private double squared(double[] coordinates)
	{
		return coordinates[0] * coordinates[0] + coordinates[1] * coordinates[1];
	}

	/**
	 * Get the distance from the origin in 3D
	 * 
	 * @param coordinates
	 * @return
	 */
	private double distance(double[] coordinates)
	{
		return Math.sqrt(coordinates[0] * coordinates[0] + coordinates[1] * coordinates[1] + coordinates[2] *
				coordinates[2]);
	}

	/**
	 * Get the squared distance from the origin in 3D
	 * 
	 * @param coordinates
	 * @return
	 */
	private double distance2(double[] coordinates)
	{
		return (coordinates[0] * coordinates[0] + coordinates[1] * coordinates[1] + coordinates[2] * coordinates[2]);
	}

	private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);

		GlobalSettings globalSettings = SettingsManager.loadSettings();
		settings = globalSettings.getCreateDataSettings();

		if (settings.stepsPerSecond < 1)
			settings.stepsPerSecond = 1;

		gd.addNumericField("Pixel_pitch (nm)", settings.pixelPitch, 2);
		gd.addNumericField("Seconds", settings.seconds, 1);
		gd.addSlider("Steps_per_second", 1, 15, settings.stepsPerSecond);
		if (extraOptions)
		{
			gd.addSlider("Aggregate_steps", 2, 20, aggregateSteps);
			gd.addNumericField("MSD_analysis_steps", msdAnalysisSteps, 0);
		}
		gd.addNumericField("Particles", settings.particles, 0);
		gd.addNumericField("Diffusion_rate (um^2/sec)", settings.diffusionRate, 2);
		if (extraOptions)
			gd.addNumericField("Precision (nm)", precision, 2);
		gd.addCheckbox("Use_grid_walk", settings.useGridWalk);
		gd.addCheckbox("Use_confinement", useConfinement);
		gd.addSlider("Confinement_attempts", 1, 20, confinementAttempts);
		gd.addSlider("Confinement_radius (nm)", 0, 3000, settings.confinementRadius);
		gd.addSlider("Fit_N", 5, 20, fitN);
		gd.addCheckbox("Show_example", showDiffusionExample);
		gd.addSlider("Magnification", 1, 10, magnification);

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		settings.pixelPitch = Math.abs(gd.getNextNumber());
		settings.seconds = Math.abs(gd.getNextNumber());
		settings.stepsPerSecond = Math.abs(gd.getNextNumber());
		if (extraOptions)
		{
			myAggregateSteps = aggregateSteps = Math.abs((int) gd.getNextNumber());
			myMsdAnalysisSteps = msdAnalysisSteps = Math.abs((int) gd.getNextNumber());
		}
		settings.particles = Math.abs((int) gd.getNextNumber());
		settings.diffusionRate = Math.abs(gd.getNextNumber());
		if (extraOptions)
			myPrecision = precision = Math.abs(gd.getNextNumber());
		settings.useGridWalk = gd.getNextBoolean();
		useConfinement = gd.getNextBoolean();
		confinementAttempts = Math.abs((int) gd.getNextNumber());
		settings.confinementRadius = Math.abs(gd.getNextNumber());
		fitN = Math.abs((int) gd.getNextNumber());
		showDiffusionExample = gd.getNextBoolean();
		magnification = gd.getNextNumber();

		// Save before validation so that the current values are preserved.
		SettingsManager.saveSettings(globalSettings);

		// Check arguments
		try
		{
			Parameters.isAboveZero("Pixel Pitch", settings.pixelPitch);
			Parameters.isAboveZero("Seconds", settings.seconds);
			Parameters.isAboveZero("Steps per second", settings.stepsPerSecond);
			Parameters.isAboveZero("Particles", settings.particles);
			Parameters.isPositive("Diffusion rate", settings.diffusionRate);
			Parameters.isAboveZero("Magnification", magnification);
			Parameters.isAboveZero("Confinement attempts", confinementAttempts);
			Parameters.isAboveZero("Fit N", fitN);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return false;
		}

		if (settings.diffusionRate == 0)
			IJ.error(TITLE, "Warning : Diffusion rate is zero");

		if (gd.invalidNumber())
			return false;

		SettingsManager.saveSettings(globalSettings);

		return true;
	}

	private void showExample(int totalSteps, double diffusionSigma, RandomGenerator[] random)
	{
		MoleculeModel m = new MoleculeModel(0, new double[3]);
		float[] xValues = new float[totalSteps];
		float[] x = new float[totalSteps];
		float[] y = new float[totalSteps];
		for (int j = 0; j < totalSteps; j++)
		{
			if (settings.useGridWalk)
				m.walk(diffusionSigma, random);
			else
				m.move(diffusionSigma, random);
			x[j] = (float) (m.getX());
			y[j] = (float) (m.getY());
			xValues[j] = (float) ((j + 1) / settings.stepsPerSecond);
		}

		// Plot x and y coords on a timeline
		String title = TITLE + " example coordinates";
		Plot2 plot = new Plot2(title, "Time (seconds)", "Distance (um)");
		float[] xUm = convertToUm(x);
		float[] yUm = convertToUm(y);
		float[] limits = Maths.limits(xUm);
		limits = Maths.limits(limits, yUm);
		plot.setLimits(0, totalSteps / settings.stepsPerSecond, limits[0], limits[1]);
		plot.setColor(Color.red);
		plot.addPoints(xValues, xUm, Plot2.LINE);
		plot.setColor(Color.blue);
		plot.addPoints(xValues, yUm, Plot2.LINE);

		Utils.display(title, plot);

		// Scale up and draw 2D position
		for (int j = 0; j < totalSteps; j++)
		{
			x[j] *= magnification;
			y[j] *= magnification;
		}
		float[] limitsx = getLimits(x);
		float[] limitsy = getLimits(y);

		int width = (int) (limitsx[1] - limitsx[0]);
		int height = (int) (limitsy[1] - limitsy[0]);

		// Ensure we draw something, even it is a simple dot at the centre for no diffusion
		if (width == 0)
		{
			width = (int) (32 * magnification);
			limitsx[0] = -width / 2;
		}
		if (height == 0)
		{
			height = (int) (32 * magnification);
			limitsy[0] = -height / 2;
		}

		ImageProcessor ip = new ByteProcessor(width, height);

		// Adjust x and y using the minimum to centre
		x[0] -= limitsx[0];
		y[0] -= limitsy[0];

		for (int j = 1; j < totalSteps; j++)
		{
			// Adjust x and y using the minimum to centre
			x[j] -= limitsx[0];
			y[j] -= limitsy[0];

			// Draw a line
			ip.setColor(32 + (223 * j) / (totalSteps - 1));
			ip.drawLine(round(x[j - 1]), round(y[j - 1]), round(x[j]), round(y[j]));
		}
		// Draw the final position
		ip.putPixel((int) round(x[totalSteps - 1]), (int) round(y[totalSteps - 1]), 255);

		ImagePlus imp = Utils.display(TITLE + " example", ip);

		// Apply the fire lookup table
		WindowManager.setTempCurrentImage(imp);
		LutLoader lut = new LutLoader();
		lut.run("fire");
		WindowManager.setTempCurrentImage(null);

	}

	private float[] convertToUm(float[] x)
	{
		final float factor = (float) (settings.pixelPitch / 1000.0);
		float[] newX = new float[x.length];
		for (int i = 0; i < x.length; i++)
			newX[i] = x[i] * factor;
		return newX;
	}

	private int round(float f)
	{
		return (int) Math.round(f);
	}

	private float[] getLimits(float[] x)
	{
		float[] limits = Maths.limits(x);
		limits[0] = (float) Math.floor(limits[0]);
		limits[1] = (float) Math.ceil(limits[1]);
		return limits;
	}

	private void aggregateIntoFrames(ArrayList<Point> points, boolean addError, double precisionInPixels,
			RandomGenerator[] random)
	{
		if (myAggregateSteps < 1)
			return;

		MemoryPeakResults results = new MemoryPeakResults(points.size() / myAggregateSteps);
		Calibration cal = new Calibration(settings.pixelPitch, 1, myAggregateSteps * 1000.0 / settings.stepsPerSecond);
		results.setCalibration(cal);
		results.setName(TITLE + " Aggregated");
		MemoryPeakResults.addResults(results);
		lastSimulatedDataset[1] = results.getName();
		int id = 0;
		int peak = 1;
		int n = 0;
		double cx = 0, cy = 0;
		// Get the mean square distance
		double sum = 0;
		int count = 0;
		PeakResult last = null;
		for (Point result : points)
		{
			final boolean newId = result.id != id;
			if (n >= myAggregateSteps || newId)
			{
				if (n != 0)
				{
					final float[] params = new float[7];
					double[] xyz = new double[] { cx / n, cy / n };
					if (addError)
						xyz = addError(xyz, precisionInPixels, random);
					params[Gaussian2DFunction.X_POSITION] = (float) xyz[0];
					params[Gaussian2DFunction.Y_POSITION] = (float) xyz[1];
					params[Gaussian2DFunction.SIGNAL] = n;
					params[Gaussian2DFunction.X_SD] = params[Gaussian2DFunction.Y_SD] = 1;
					final float noise = 0.1f;
					PeakResult r = new ExtendedPeakResult(peak, (int) params[Gaussian2DFunction.X_POSITION],
							(int) params[Gaussian2DFunction.Y_POSITION], n, 0, noise, params, null, peak, id);
					results.add(r);
					if (last != null)
					{
						sum += last.distance2(r);
						count++;
					}
					last = r;
					n = 0;
					cx = cy = 0;
					peak++;
				}
				if (newId)
				{
					peak++; // Increment the frame so that tracing analysis can distinguish traces
					last = null;
					id = result.id;
				}
			}
			n++;
			cx += result.x;
			cy += result.y;
		}

		// Final peak
		if (n != 0)
		{
			final float[] params = new float[7];
			double[] xyz = new double[] { cx / n, cy / n };
			if (addError)
				xyz = addError(xyz, precisionInPixels, random);
			params[Gaussian2DFunction.X_POSITION] = (float) xyz[0];
			params[Gaussian2DFunction.Y_POSITION] = (float) xyz[1];
			params[Gaussian2DFunction.SIGNAL] = n;
			params[Gaussian2DFunction.X_SD] = params[Gaussian2DFunction.Y_SD] = 1;
			final float noise = 0.1f;
			PeakResult r = new ExtendedPeakResult(peak, (int) params[Gaussian2DFunction.X_POSITION],
					(int) params[Gaussian2DFunction.Y_POSITION], n, 0, noise, params, null, peak, id);
			results.add(r);
			if (last != null)
			{
				sum += last.distance2(r);
				count++;
			}
		}

		// MSD in pixels^2 / frame
		double msd = sum / count;
		// Convert to um^2/second
		Utils.log("Aggregated data D=%s um^2/s, Precision=%s nm, N=%d, step=%s s, mean=%s um^2, MSD = %s um^2/s",
				Utils.rounded(settings.diffusionRate), Utils.rounded(myPrecision), count,
				Utils.rounded(results.getCalibration().exposureTime / 1000), Utils.rounded(msd / conversionFactor),
				Utils.rounded((msd / conversionFactor) / (results.getCalibration().exposureTime / 1000)));

		msdAnalysis(points);
	}

	/**
	 * Tabulate the observed MSD for different jump distances
	 * 
	 * @param points
	 */
	private void msdAnalysis(ArrayList<Point> points)
	{
		if (myMsdAnalysisSteps == 0)
			return;

		IJ.showStatus("MSD analysis ...");
		IJ.showProgress(1, myMsdAnalysisSteps);

		// This will only be fast if the list is an array
		Point[] list = points.toArray(new Point[points.size()]);

		// Compute the base MSD
		Point origin = new Point(0, 0, 0);
		double sum = origin.distance2(list[0]);
		int count = 1;
		for (int i = 1; i < list.length; i++)
		{
			Point last = list[i - 1];
			Point current = list[i];

			if (last.id == current.id)
			{
				sum += last.distance2(current);
			}
			else
			{
				sum += origin.distance2(current);
			}
			count++;
		}
		createMsdTable((sum / count) * settings.stepsPerSecond / conversionFactor);

		// Create a new set of points that have coordinates that 
		// are the rolling average over the number of aggregate steps
		RollingArray x = new RollingArray(aggregateSteps);
		RollingArray y = new RollingArray(aggregateSteps);

		int id = 0;
		int length = 0;
		for (Point p : points)
		{
			if (p.id != id)
			{
				x.reset();
				y.reset();
			}
			id = p.id;
			x.add(p.x);
			y.add(p.y);
			// Only create a point if the full aggregation size is reached
			if (x.isFull())
			{
				list[length++] = new Point(id, x.getAverage(), y.getAverage());
			}
		}

		// Q - is this useful?
		final double p = myPrecision / settings.pixelPitch;
		final long seed = System.currentTimeMillis() + System.identityHashCode(this);
		RandomGenerator rand = new Well19937c(seed);

		StringBuilder sb = new StringBuilder();
		final int totalSteps = (int) Math.ceil(settings.seconds * settings.stepsPerSecond - aggregateSteps);
		final int limit = Math.min(totalSteps, myMsdAnalysisSteps);
		final int interval = Utils.getProgressInterval(limit);
		for (int step = 1; step <= myMsdAnalysisSteps; step++)
		{
			if (step % interval == 0)
				IJ.showProgress(step, limit);

			sum = 0;
			count = 0;
			for (int i = step; i < length; i++)
			{
				Point last = list[i - step];
				Point current = list[i];

				if (last.id == current.id)
				{
					if (p == 0)
					{
						sum += last.distance2(current);
						count++;
					}
					else
					{
						// This can be varied but the effect on the output with only 1 loop 
						// is the same if enough samples are present
						for (int ii = 1; ii-- > 0;)
						{
							sum += last.distance2(current, p, rand);
							count++;
						}
					}
				}
			}
			if (count == 0)
				break;
			addResult(sb, step, sum, count);

			// Flush to auto-space the columns
			if (step == 9)
			{
				msdTable.append(sb.toString());
				sb.setLength(0);
			}
		}
		msdTable.append(sb.toString());

		IJ.showProgress(1);
	}

	private void createMsdTable(double baseMsd)
	{
		String header = createHeader(baseMsd);
		if (msdTable == null || !msdTable.isVisible())
		{
			msdTable = new TextWindow("MSD Analysis", header, "", 800, 300);
			msdTable.setVisible(true);
		}
		else
		{
			//msdTable.getTextPanel().clear();
		}
	}

	private String prefix;
	private double exposureTime;

	private String createHeader(double baseMsd)
	{
		final double apparentD = baseMsd / 4;
		StringBuilder sb = new StringBuilder();
		sb.append(settings.diffusionRate).append('\t');
		sb.append(myPrecision).append('\t');
		sb.append(Utils.rounded(apparentD)).append('\t');
		sb.append(Utils.rounded(1.0 / settings.stepsPerSecond)).append('\t');
		sb.append(myAggregateSteps).append('\t');
		// Exposure time is the aggregated frame time 
		exposureTime = myAggregateSteps / settings.stepsPerSecond;
		sb.append(Utils.rounded(exposureTime)).append('\t');
		prefix = sb.toString();
		return "D (um^2/s)\tPrecision (nm)\tDsim (um^2/s)\tStep (s)\tResolution\tFrame (s)\tt (s)\tn\tN\tMSD (um^2)\tD (um^2/s)";
	}

	private void addResult(StringBuilder sb, int step, double sum, int count)
	{
		// Exposure time is the aggregated frame time 
		final double msd = (sum / count) / conversionFactor;
		// Jump distance separation is the number of steps
		final double t = step / settings.stepsPerSecond;
		sb.append(prefix);
		sb.append(Utils.rounded(t)).append('\t');
		sb.append(Utils.rounded(t / exposureTime)).append('\t');
		sb.append(count).append('\t');
		// Not rounded to preserve precision 
		sb.append(msd).append('\t');
		sb.append(msd / (4 * t));
		sb.append('\n');
	}
}
