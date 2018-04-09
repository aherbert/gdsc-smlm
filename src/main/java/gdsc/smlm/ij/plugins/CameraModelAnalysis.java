package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Color;
import java.util.Arrays;

import org.apache.commons.math3.distribution.CustomGammaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.random.AbstractRandomGenerator;
import org.apache.commons.math3.random.Well19937c;

import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.Utils;
import gdsc.core.logging.Ticker;
import gdsc.core.threshold.IntHistogram;
import gdsc.core.utils.Maths;
import gdsc.core.utils.PseudoRandomGenerator;
import gdsc.smlm.data.config.GUIProtos.CameraModelAnalysisSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.DialogListener;
import ij.gui.ExtendedGenericDialog;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.Plot2;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.ImageProcessor;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Model the on-chip amplification from an EM-CCD camera, CCD or sCMOS camera.
 */
public class CameraModelAnalysis implements ExtendedPlugInFilter, DialogListener
{
	private static final String TITLE = "Camera Model Analysis";

	private CameraModelAnalysisSettings.Builder settings;

	private PlugInFilterRunner pfr;
	private boolean dirty = true;

	private static String[] MODE = { "CCD", "EM-CCD" };
	private static String[] MODEL = { "Poisson", "Poisson+Gaussian", "Poisson+Poisson", "Poisson+Gamma+Gaussian" };

	private static PseudoRandomGenerator random = new PseudoRandomGenerator(new double[1]);
	private static int seed = 0;

	private PseudoRandomGenerator getRandomGenerator()
	{
		if (random.getLength() < settings.getSamples() || seed != settings.getSeed())
		{
			random = new PseudoRandomGenerator(settings.getSamples(), new Well19937c(settings.getSeed()));
			seed = settings.getSeed();
		}
		return random;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	public int setup(String arg, ImagePlus imp)
	{
		return NO_IMAGE_REQUIRED;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.ExtendedPlugInFilter#showDialog(ij.ImagePlus, java.lang.String,
	 * ij.plugin.filter.PlugInFilterRunner)
	 */
	public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr)
	{
		this.pfr = pfr;

		settings = SettingsManager.readCameraModelAnalysisSettings(0).toBuilder();

		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Simulate on-chip camera applification.");

		gd.addNumericField("Photons", settings.getPhotons(), 2);
		gd.addChoice("Mode", MODE, settings.getMode());
		gd.addNumericField("Gain", settings.getGain(), 2);
		gd.addNumericField("Noise", settings.getNoise(), 2);
		gd.addNumericField("Seed", settings.getSeed(), 0);
		gd.addNumericField("Samples", settings.getSamples(), 0);
		gd.addNumericField("Noise_samples", settings.getNoiseSamples(), 0);
		gd.addChoice("Model", MODEL, settings.getModel());
		gd.addDialogListener(this);
		gd.addPreviewCheckbox(pfr);
		gd.showDialog();

		if (gd.wasCanceled())
			return DONE;

		if (dirty)
			execute();

		SettingsManager.writeSettings(settings);

		return DONE;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.gui.DialogListener#dialogItemChanged(ij.gui.GenericDialog, java.awt.AWTEvent)
	 */
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
	{
		dirty = true;
		settings.setPhotons(gd.getNextNumber());
		settings.setMode(gd.getNextChoiceIndex());
		settings.setGain(gd.getNextNumber());
		settings.setNoise(gd.getNextNumber());
		settings.setSeed((int) gd.getNextNumber());
		settings.setSamples(Math.max(1, (int) gd.getNextNumber()));
		settings.setNoiseSamples(Math.max(1, (int) gd.getNextNumber()));
		settings.setModel(gd.getNextChoiceIndex());
		if (gd.getPreviewCheckbox().getState())
			execute();
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.ExtendedPlugInFilter#setNPasses(int)
	 */
	public void setNPasses(int nPasses)
	{
		// Ignore
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	public void run(ImageProcessor ip)
	{
		execute();
	}

	/**
	 * Execute the analysis
	 */
	private void execute()
	{
		dirty = false;

		IntHistogram h = simulateHistogram();

		// Build cumulative distribution
		double[][] cdf = cumulativeHistogram(h);
		double[] x = cdf[0];
		double[] y = cdf[1];

		String title = TITLE;
		Plot2 plot = new Plot2(title, "Count", "CDF");
		plot.setLimits(x[0], x[x.length - 1], 0, 1.05);
		plot.setColor(Color.blue);
		plot.addPoints(x, y, Plot2.BAR);
		plot.setColor(Color.black);
		plot.addLegend("CDF");
		Utils.display(title, plot);

		// Interpolate to 300 steps faster evaluation?

		// Get likelihood function

		// Create likelihood cumulative distribution

		// Compute Komolgorov distance

		// Plot

	}

	private double[][] cumulativeHistogram(IntHistogram histogram)
	{
		int[] h = histogram.h;
		double[] x = new double[h.length];
		double[] y = new double[x.length];
		double sum = 0;
		for (int i = 0; i < x.length; i++)
		{
			x[i] = histogram.getValue(i);
			sum += h[i];
			y[i] = sum;
		}
		for (int i = 0; i < y.length; i++)
		{
			y[i] /= sum;
		}
		y[y.length - 1] = 1;
		return new double[][] { x, y };
	}

	/**
	 * Simulate the histogram for fitting
	 * 
	 * @param method
	 *            0 - sample from the fitted PDF, 1 - sample from a Poisson-Gamma-Gaussian
	 * @return The histogram
	 */
	private IntHistogram simulateHistogram()
	{
		IJ.showStatus("Simulating histogram ...");
		IntHistogram h;
		switch (settings.getMode())
		{
			case 1:
				h = simulatePoissonGammaGaussian();
				break;

			case 0:
				h = simulatePoissonGaussian();
				break;

			default:
				throw new IllegalStateException();
		}
		IJ.showStatus("");
		return h;
	}

	/**
	 * Randomly generate a histogram from poisson-gamma-gaussian samples
	 * 
	 * @return The histogram
	 */
	private IntHistogram simulatePoissonGammaGaussian()
	{
		// Randomly sample
		PseudoRandomGenerator random = getRandomGenerator();
		PseudoRandomGenerator random2 = random.clone();

		PoissonDistribution poisson = new PoissonDistribution(random, settings.getPhotons(),
				PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);

		CustomGammaDistribution gamma = new CustomGammaDistribution(random.clone(), 1, settings.getGain(),
				GammaDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);

		final double noise = settings.getNoise();
		final int noiseSamples = (noise > 0) ? settings.getNoiseSamples() : 1;
		final int steps = settings.getSamples() * noiseSamples;
		int[] sample = new int[steps];
		Ticker ticker = Ticker.createStarted(new IJTrackProgress(), steps, false);
		for (int n = settings.getSamples(), i = 0; n -- > 0; )
		{
			// Poisson
			double d = poisson.sample();

			// Gamma
			if (d > 0)
			{
				gamma.setShapeUnsafe(d);
				d = gamma.sample();
			}

			// Over-sample the Gaussian
			for (int j = noiseSamples; j-- > 0;)
			{
				// Convert the sample to a count 
				sample[i++] = (int) Math.round(d + noise * random2.nextGaussian());
				ticker.tick();
			}
		}

		ticker.stop();

		int[] limits = Maths.limits(sample);
		int min = limits[0];
		int max = limits[1];

		int[] h = new int[max - min + 1];
		for (int s : sample)
			h[s - min]++;
		return new IntHistogram(h, min);
	}

	/**
	 * Randomly generate a histogram from poisson-gaussian samples
	 * 
	 * @return The histogram
	 */
	private IntHistogram simulatePoissonGaussian()
	{
		// Randomly sample
		PseudoRandomGenerator random = getRandomGenerator();
		PseudoRandomGenerator random2 = random.clone();

		PoissonDistribution poisson = new PoissonDistribution(random, settings.getPhotons(),
				PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);

		// Fixed gain
		final double gain = settings.getGain();
		
		final double noise = settings.getNoise();
		final int noiseSamples = (noise > 0) ? settings.getNoiseSamples() : 1;
		final int steps = settings.getSamples() * noiseSamples;
		int[] sample = new int[steps];
		Ticker ticker = Ticker.createStarted(new IJTrackProgress(), steps, false);
		for (int n = settings.getSamples(), i = 0; n -- > 0; )
		{
			// Poisson
			double d = poisson.sample();

			// Fixed gain
			d *= gain;

			// Over-sample the Gaussian
			for (int j = noiseSamples; j-- > 0;)
			{
				// Convert the sample to a count 
				sample[i++] = (int) Math.round(d + noise * random2.nextGaussian());
				ticker.tick();
			}
		}

		ticker.stop();

		int[] limits = Maths.limits(sample);
		int min = limits[0];
		int max = limits[1];

		int[] h = new int[max - min + 1];
		for (int s : sample)
			h[s - min]++;
		return new IntHistogram(h, min);
	}

	//	@SuppressWarnings("unused")
	//	private void plotPMF()
	//	{
	//		if (!showPMFDialog())
	//			return;
	//
	//		double step = getStepSize(_photons, _gain, _noise);
	//
	//		PDF pdf = pdf(0, step, _photons, _gain, _noise);
	//		double[] pmf = pdf.p;
	//		double[] x = pdf.x;
	//		double yMax = Maths.max(pmf);
	//
	//		// Get the approximation
	//		LikelihoodFunction fun;
	//		double myNoise = _noise;
	//		switch (approximation)
	//		{
	//			case 3:
	//				fun = new PoissonFunction(1.0 / _gain, true);
	//				break;
	//			case 2:
	//				// The mean does not matter (as normalisation is done dynamically for 
	//				// PoissonGaussianFunction.likelihood(double, double) so just use zero
	//				//fun = PoissonGaussianFunction.createWithStandardDeviation(1.0 / _gain, 0, _noise);
	//
	//				// Use adaptive normalisation
	//				fun = PoissonGaussianFunction2.createWithStandardDeviation(1.0 / _gain, _noise);
	//				break;
	//			case 1:
	//				myNoise = 0;
	//			case 0:
	//			default:
	//				PoissonGammaGaussianFunction myFun = new PoissonGammaGaussianFunction(1.0 / _gain, myNoise);
	//				myFun.setMinimumProbability(0);
	//				fun = myFun;
	//		}
	//		double expected = _photons;
	//		if (offset != 0)
	//			expected += offset * expected / 100.0;
	//		expected *= _gain;
	//
	//		// Normalise 
	//		boolean normalise = false;
	//		if (normalise)
	//		{
	//			double sum = Maths.sum(pmf);
	//			for (int i = pmf.length; i-- > 0;)
	//				pmf[i] /= sum;
	//		}
	//
	//		// Get CDF
	//		double sum = 0;
	//		double sum2 = 0;
	//		double prev = 0;
	//		double prev2 = 0;
	//		double[] f = new double[x.length];
	//		double[] cdf1 = new double[pmf.length];
	//		double[] cdf2 = new double[pmf.length];
	//		double step_2 = step / 2;
	//		for (int i = 0; i < cdf1.length; i++)
	//		{
	//			// Trapezoid integration
	//			//sum += (pmf[i] + prev) * step_2;
	//			sum += pmf[i] * step;
	//			cdf1[i] = sum;
	//			prev = pmf[i];
	//			f[i] = fun.likelihood(x[i], expected);
	//			//sum2 += (f[i] + prev2) * step_2;
	//			sum2 += f[i] * step;
	//			cdf2[i] = sum2;
	//			prev2 = f[i];
	//		}
	//
	//		// Truncate x for plotting
	//		int max = 0;
	//		sum = prev = 0;
	//		double p = 1 - tail;
	//		while (sum < p && max < pmf.length)
	//		{
	//			//sum += (pmf[max] + prev) * step_2;
	//			sum += pmf[max] * step;
	//			prev = pmf[max];
	//			if (sum > 0.5 && pmf[max] == 0)
	//				break;
	//			max++;
	//		}
	//
	//		int min = pmf.length;
	//		sum = prev = 0;
	//		p = 1 - head;
	//		while (sum < p && min > 0)
	//		{
	//			min--;
	//			//sum += (pmf[min] + prev) * step_2;
	//			sum += pmf[min] * step;
	//			prev = pmf[min];
	//			if (sum > 0.5 && pmf[min] == 0)
	//				break;
	//		}
	//
	//		//int min = (int) (dummyBias - gaussWidth * _noise);
	//		pmf = Arrays.copyOfRange(pmf, min, max);
	//		x = Arrays.copyOfRange(x, min, max);
	//		f = Arrays.copyOfRange(f, min, max);
	//
	//		if (showApproximation)
	//			yMax = Maths.maxDefault(yMax, f);
	//
	//		String label = String.format("Gain=%s, noise=%s, photons=%s", Utils.rounded(_gain), Utils.rounded(_noise),
	//				Utils.rounded(_photons));
	//
	//		Plot2 plot = new Plot2("PMF", "ADUs", "p");
	//		plot.setLimits(x[0], x[x.length - 1], 0, yMax);
	//		plot.setColor(Color.red);
	//		plot.addPoints(x, pmf, Plot2.LINE);
	//		if (showApproximation)
	//		{
	//			plot.setColor(Color.blue);
	//			plot.addPoints(x, f, Plot2.LINE);
	//		}
	//
	//		plot.setColor(Color.magenta);
	//		plot.drawLine(_photons * _gain, 0, _photons * _gain, yMax);
	//		plot.setColor(Color.black);
	//		plot.addLabel(0, 0, label);
	//		PlotWindow win1 = Utils.display("PMF", plot);
	//
	//		// Plot the difference between the actual and approximation
	//		double[] delta = new double[f.length];
	//		for (int i = 0; i < f.length; i++)
	//		{
	//			if (pmf[i] == 0 && f[i] == 0)
	//				continue;
	//			if (relativeDelta)
	//				delta[i] = DoubleEquality.relativeError(f[i], pmf[i]) * Math.signum(f[i] - pmf[i]);
	//			else
	//				delta[i] = f[i] - pmf[i];
	//		}
	//
	//		Plot2 plot2 = new Plot2("PMF delta", "ADUs", (relativeDelta) ? "Relative delta" : "delta");
	//		double[] limits = Maths.limits(delta);
	//		plot2.setLimits(x[0], x[x.length - 1], limits[0], limits[1]);
	//		plot2.setColor(Color.red);
	//		plot2.addPoints(x, delta, Plot2.LINE);
	//		plot2.setColor(Color.magenta);
	//		plot2.drawLine(_photons * _gain, limits[0], _photons * _gain, limits[1]);
	//		plot2.setColor(Color.black);
	//		plot2.addLabel(0, 0, label + ((offset == 0) ? "" : ", expected = " + Utils.rounded(expected / _gain)));
	//		PlotWindow win2 = Utils.display("PMF delta", plot2);
	//
	//		if (Utils.isNewWindow())
	//		{
	//			Point p2 = win1.getLocation();
	//			p2.y += win1.getHeight();
	//			win2.setLocation(p2);
	//		}
	//
	//		// Plot the CDF of each distribution.
	//		// Compute the Kolmogorov distance as the supremum (maximum) 
	//		// difference between the two cumulative probability distributions.
	//		// https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test
	//		double kolmogorovDistance = 0;
	//		double xd = x[0];
	//		for (int i = 0; i < cdf1.length; i++)
	//		{
	//			double dist = Math.abs(cdf1[i] - cdf2[i]);
	//			if (kolmogorovDistance < dist)
	//			{
	//				kolmogorovDistance = dist;
	//				xd = pdf.x[i];
	//			}
	//		}
	//		cdf1 = Arrays.copyOfRange(cdf1, min, max);
	//		cdf2 = Arrays.copyOfRange(cdf2, min, max);
	//
	//		Plot2 plot3 = new Plot2("CDF", "ADUs", "p");
	//		yMax = 1.05;
	//		plot3.setLimits(x[0], x[x.length - 1], 0, yMax);
	//		plot3.setColor(Color.red);
	//		plot3.addPoints(x, cdf1, Plot2.LINE);
	//		plot3.setColor(Color.blue);
	//		plot3.addPoints(x, cdf2, Plot2.LINE);
	//
	//		plot3.setColor(Color.magenta);
	//		plot3.drawLine(_photons * _gain, 0, _photons * _gain, yMax);
	//		plot3.drawDottedLine(xd, 0, xd, yMax, 2);
	//		plot3.setColor(Color.black);
	//		plot3.addLabel(0, 0, label + ", Kolmogorov distance = " + Utils.rounded(kolmogorovDistance) + " @ " + xd);
	//		plot3.addLegend("CDF\nApprox");
	//		PlotWindow win3 = Utils.display("CDF", plot3);
	//
	//		if (Utils.isNewWindow())
	//		{
	//			Point p2 = win1.getLocation();
	//			p2.x += win1.getWidth();
	//			win3.setLocation(p2);
	//		}
	//	}
}
