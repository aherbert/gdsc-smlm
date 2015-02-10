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

import gdsc.smlm.ij.settings.CreateDataSettings;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.model.ImageModel;
import gdsc.smlm.model.MoleculeModel;
import gdsc.smlm.model.SphericalDistribution;
import gdsc.smlm.utils.Maths;
import gdsc.smlm.utils.Statistics;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.SuperPlot;
import ij.plugin.LutLoader;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

import java.awt.Color;
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

	private CreateDataSettings settings;
	private static boolean useConfinement = false;
	private static int confinementAttempts = 5;
	private static int fitN = 20;
	private static boolean showDiffusionExample = false;
	private static double magnification = 5;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		if (!showDialog())
			return;

		int totalSteps = settings.seconds * settings.stepsPerSecond;

		final double conversionFactor = 1000000.0 / (settings.pixelPitch * settings.pixelPitch);

		// Diffusion rate is um^2/sec. Convert to pixels per simulation frame.
		final double diffusionRateInPixelsPerSecond = settings.diffusionRate * conversionFactor;
		final double diffusionRateInPixelsPerStep = diffusionRateInPixelsPerSecond / settings.stepsPerSecond;

		Utils.log(TITLE + " : D = %s um^2/sec", Utils.rounded(settings.diffusionRate, 4));
		Utils.log("Mean-displacement per dimension = %s nm/sec",
				Utils.rounded(1e3 * ImageModel.getRandomMoveDistance(settings.diffusionRate), 4));

		// Convert diffusion co-efficient into the standard deviation for the random walk
		final double diffusionSigma = ImageModel.getRandomMoveDistance(diffusionRateInPixelsPerStep);
		Utils.log("Simulation step-size = %s nm", Utils.rounded(settings.pixelPitch * diffusionSigma, 4));

		// Move the molecules and get the diffusion rate
		IJ.showStatus("Simulating ...");
		final long start = System.nanoTime();
		RandomGenerator random = new Well19937c(System.currentTimeMillis() + System.identityHashCode(this));
		Statistics[] stats2D = new Statistics[totalSteps];
		Statistics[] stats3D = new Statistics[totalSteps];
		for (int j = 0; j < totalSteps; j++)
		{
			stats2D[j] = new Statistics();
			stats3D[j] = new Statistics();
		}
		SphericalDistribution dist = new SphericalDistribution(settings.confinementRadius / settings.pixelPitch);
		Statistics asymptote = new Statistics();
		for (int i = 0; i < settings.particles; i++)
		{
			if (i % 16 == 0)
				IJ.showProgress(i, settings.particles);

			MoleculeModel m = new MoleculeModel(i, new double[3]);
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
					double[] originalXyz = Arrays.copyOf(xyz, 3);
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

					stats2D[j].add(squared(m.getCoordinates()));
					stats3D[j].add(distance2(m.getCoordinates()));
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
						stats2D[j].add(squared(m.getCoordinates()));
						stats3D[j].add(distance2(m.getCoordinates()));
					}
				}
				else
				{
					for (int j = 0; j < totalSteps; j++)
					{
						m.move(diffusionSigma, random);
						stats2D[j].add(squared(m.getCoordinates()));
						stats3D[j].add(distance2(m.getCoordinates()));
					}
				}
			}

			// Debug: record all the particles so they can be analysed
			// System.out.printf("%f %f %f\n", m.getX(), m.getY(), m.getZ());
		}
		final double time = (System.nanoTime() - start) / 1000000.0;

		IJ.showStatus("Analysing results ...");
		IJ.showProgress(1);

		if (showDiffusionExample)
		{
			showExample(totalSteps, diffusionSigma, random);
		}

		// SuperPlot a graph of mean squared distance
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

		// Do linear regression to get diffusion rate
		final double[] init2D = { 0, 1 / gradient2D.getMean() }; // a - b x
		final double[] best2D = fitter2D.fit(new PolynomialFunction.Parametric(), init2D);
		final PolynomialFunction fitted2D = new PolynomialFunction(best2D);

		final double[] init3D = { 0, 1 / gradient3D.getMean() }; // a - b x
		final double[] best3D = fitter3D.fit(new PolynomialFunction.Parametric(), init3D);
		final PolynomialFunction fitted3D = new PolynomialFunction(best3D);

		// Create plot
		String title = TITLE;
		SuperPlot plot = new SuperPlot(title, "Time (seconds)", "Mean-squared Distance (um^2)", xValues, yValues2D);
		double[] limits = Maths.limits(upper);
		limits = Maths.limits(limits, lower);
		limits = Maths.limits(limits, yValues3D);
		plot.setLimits(0, totalSteps / settings.stepsPerSecond, limits[0], limits[1]);
		plot.setColor(Color.blue);
		plot.addPoints(xValues, lower, SuperPlot.LINE);
		plot.addPoints(xValues, upper, SuperPlot.LINE);
		plot.setColor(Color.magenta);
		plot.addPoints(xValues, yValues3D, SuperPlot.LINE);
		plot.setColor(Color.red);
		plot.addPoints(new double[] { xValues[0], xValues[xValues.length - 1] },
				new double[] { fitted2D.value(xValues[0]), fitted2D.value(xValues[xValues.length - 1]) }, SuperPlot.LINE);
		plot.setColor(Color.green);
		plot.addPoints(new double[] { xValues[0], xValues[xValues.length - 1] },
				new double[] { fitted3D.value(xValues[0]), fitted3D.value(xValues[xValues.length - 1]) }, SuperPlot.LINE);
		plot.setColor(Color.black);

		Utils.display(title, plot);

		// For 2D diffusion: d^2 = 4D
		// where: d^2 = mean-square displacement

		double D = best2D[1] / 4.0;
		String msg = "2D Diffusion rate = " + Utils.rounded(D, 4) + " um^2 / sec (" + Utils.timeToString(time) + ")";
		IJ.showStatus(msg);
		Utils.log(msg);

		D = best3D[1] / 6.0;
		Utils.log("3D Diffusion rate = " + Utils.rounded(D, 4) + " um^2 / sec (" + Utils.timeToString(time) + ")");

		if (useConfinement)
			Utils.log("3D asymptote distance = %s nm (expected %.2f)",
					Utils.rounded(asymptote.getMean() * settings.pixelPitch, 4), 3 * settings.confinementRadius / 4);
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
		gd.addNumericField("Seconds", settings.seconds, 0);
		gd.addSlider("Steps_per_second", 1, 15, settings.stepsPerSecond);
		gd.addNumericField("Particles", settings.particles, 0);
		gd.addNumericField("Diffusion_rate (um^2/sec)", settings.diffusionRate, 2);
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
		settings.seconds = Math.abs((int) gd.getNextNumber());
		settings.stepsPerSecond = Math.abs((int) gd.getNextNumber());
		settings.particles = Math.abs((int) gd.getNextNumber());
		settings.diffusionRate = Math.abs(gd.getNextNumber());
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

		if (gd.invalidNumber())
			return false;

		SettingsManager.saveSettings(globalSettings);

		return true;
	}

	private void showExample(int totalSteps, double diffusionSigma, RandomGenerator random)
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
			xValues[j] = (float) (j + 1) / settings.stepsPerSecond;
		}

		// SuperPlot x and y coords on a timeline
		String title = TITLE + " example coordinates";
		SuperPlot plot = new SuperPlot(title, "Time (seconds)", "Distance (um)");
		float[] xUm = convertToUm(x);
		float[] yUm = convertToUm(y);
		float[] limits = Maths.limits(xUm);
		limits = Maths.limits(limits, yUm);
		plot.setLimits(0, totalSteps / settings.stepsPerSecond, limits[0], limits[1]);
		plot.setColor(Color.red);
		plot.addPoints(xValues, xUm, SuperPlot.LINE);
		plot.setColor(Color.blue);
		plot.addPoints(xValues, yUm, SuperPlot.LINE);

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
}
