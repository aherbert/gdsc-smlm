package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Color;
import java.util.Arrays;

import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.math3.distribution.CustomGammaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.random.Well19937c;

import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.Utils;
import gdsc.core.logging.Ticker;
import gdsc.core.math.Geometry;
import gdsc.core.threshold.IntHistogram;
import gdsc.core.utils.Maths;
import gdsc.core.utils.PseudoRandomGenerator;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.smlm.data.config.GUIProtos.CameraModelAnalysisSettings;
import gdsc.smlm.function.LikelihoodFunction;
import gdsc.smlm.function.PoissonFunction;
import gdsc.smlm.function.PoissonGammaGaussianFunction;
import gdsc.smlm.function.PoissonGaussianFunction2;
import gdsc.smlm.function.PoissonPoissonFunction;
import gdsc.smlm.ij.settings.SettingsManager;
import gnu.trove.list.array.TDoubleArrayList;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.DialogListener;
import ij.gui.ExtendedGenericDialog;
import ij.gui.ExtendedGenericDialog.OptionCollectedEvent;
import ij.gui.ExtendedGenericDialog.OptionCollectedListener;
import ij.gui.ExtendedGenericDialog.OptionListener;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingExtendedGenericDialog;
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
public class CameraModelAnalysis implements ExtendedPlugInFilter, DialogListener, OptionCollectedListener
{
	private static final String TITLE = "Camera Model Analysis";

	private CameraModelAnalysisSettings.Builder settings;

	private boolean extraOptions;
	private boolean dirty = true;
	private CameraModelAnalysisSettings lastSettings = null;
	private ExtendedGenericDialog gd;
	private IntHistogram lastHistogram = null;
	private CameraModelAnalysisSettings lastSimulationSettings = null;

	private static String[] MODE = { "CCD", "EM-CCD" };
	private static String[] MODEL = { "Poisson (Discrete)", "Poisson (Continuous)", "Poisson+Gaussian",
			"Poisson+Poisson", "Poisson+Gamma+Gaussian" };

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
		extraOptions = Utils.isExtraOptions();
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
		settings = SettingsManager.readCameraModelAnalysisSettings(0).toBuilder();

		gd = new NonBlockingExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Simulate on-chip camera applification.");

		gd.addNumericField("Photons", settings.getPhotons(), 2);
		gd.addChoice("Mode", MODE, settings.getMode(), new OptionListener<Integer>()
		{
			public boolean collectOptions(Integer value)
			{
				settings.setMode(value);
				return collectOptions(false);
			}

			public boolean collectOptions()
			{
				return collectOptions(true);
			}

			private boolean collectOptions(boolean silent)
			{
				int mode = settings.getMode();
				ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
				if (mode == 0)
				{
					egd.addNumericField("Gain", settings.getGain(), 2, 6, "Count/electrons");
					egd.addNumericField("Noise", settings.getNoise(), 2, 6, "Count");
				}
				else if (mode == 1)
				{
					egd.addNumericField("Gain", settings.getEmGain(), 2, 6, "Count/electrons");
					egd.addNumericField("Noise", settings.getEmNoise(), 2, 6, "Count");
					egd.addNumericField("EM_samples", settings.getEmSamples(), 0);
				}
				else
					throw new IllegalStateException();
				egd.showDialog(true, gd);
				if (egd.wasCanceled())
					return false;
				if (mode == 0)
				{
					settings.setGain(egd.getNextNumber());
					settings.setNoise(egd.getNextNumber());
				}
				else
				{
					settings.setEmGain(egd.getNextNumber());
					settings.setEmNoise(egd.getNextNumber());
					settings.setEmSamples(Math.max(1, (int) egd.getNextNumber()));
				}
				return true;
			}
		});
		if (extraOptions)
			gd.addNumericField("Seed", settings.getSeed(), 0);
		gd.addNumericField("Samples", settings.getSamples(), 0);
		gd.addNumericField("Noise_samples", settings.getNoiseSamples(), 0);
		gd.addChoice("Model", MODEL, settings.getModel());
		gd.addOptionCollectedListener(this);
		gd.addDialogListener(this);
		gd.addPreviewCheckbox(pfr);
		gd.showDialog();

		SettingsManager.writeSettings(settings);

		if (!gd.wasCanceled() && dirty)
			execute();

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
		if (extraOptions)
			settings.setSeed((int) gd.getNextNumber());
		settings.setSamples(Math.max(1, (int) gd.getNextNumber()));
		settings.setNoiseSamples(Math.max(1, (int) gd.getNextNumber()));
		settings.setModel(gd.getNextChoiceIndex());
		if (gd.getPreviewCheckbox().getState())
		{
			return execute();
		}
		return true;
	}

	public void optionCollected(OptionCollectedEvent e)
	{
		if (gd.getPreviewCheckbox().getState())
		{
			if (!execute())
			{
				gd.getPreviewCheckbox().setState(false);
			}
		}
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
	private boolean execute()
	{
		dirty = false;

		CameraModelAnalysisSettings settings = this.settings.build();
		if (!(getTotalGain(settings) > 0))
		{
			Utils.log(TITLE + "Error: No total gain");
			return false;
		}
		if (!(settings.getPhotons() > 0))
		{
			Utils.log(TITLE + "Error: No photons");
			return false;
		}
		// Avoid repeating the same analysis
		if (settings.equals(lastSettings))
			return true;
		lastSettings = settings;

		IntHistogram h = getHistogram(settings);

		// Build cumulative distribution
		double[][] cdf1 = cumulativeHistogram(h);
		double[] x1 = cdf1[0];
		double[] y1 = cdf1[1];

		// Interpolate to 300 steps faster evaluation?

		// Get likelihood function
		LikelihoodFunction f = getLikelihoodFunction(settings);

		// Create likelihood cumulative distribution
		double[][] cdf2 = cumulativeDistribution(settings, cdf1, f);

		// Compute Komolgorov distance
		double[] distanceAndValue = getDistance(cdf1, cdf2);
		double distance = distanceAndValue[0];
		double value = distanceAndValue[1];
		double area = distanceAndValue[2];

		double[] x2 = cdf2[0];
		double[] y2 = cdf2[1];

		// Fill y1
		int offset = 0;
		while (x2[offset] < x1[0])
			offset++;
		double[] y1b = new double[y2.length];
		System.arraycopy(y1, 0, y1b, offset, y1.length);
		Arrays.fill(y1b, offset + y1.length, y2.length, y1[y1.length - 1]);

		// Plot
		String title = TITLE;
		Plot2 plot = new Plot2(title, "Count", "CDF");
		double max = 1.05 * Maths.maxDefault(1, y2);
		plot.setLimits(x2[0], x2[x2.length - 1], 0, max);
		plot.setColor(Color.blue);
		//plot.addPoints(x2, y1b, Plot2.LINE);
		plot.addPoints(x2, y1b, Plot2.BAR);
		plot.setColor(Color.red);
		//plot.addPoints(x2, y2, Plot2.LINE);
		plot.addPoints(x2, y2, Plot2.BAR);
		plot.setColor(Color.magenta);
		plot.drawLine(value, 0, value, max);
		plot.setColor(Color.black);
		plot.addLegend("CDF\nModel");
		plot.addLabel(0, 0,
				String.format("Distance=%s @ %.0f (Mean=%s)", Utils.rounded(distance), value, Utils.rounded(area)));
		Utils.display(title, plot);

		// Q. Plot PMF as well using the CDF to generate it?
		return true;
	}

	private IntHistogram getHistogram(CameraModelAnalysisSettings settings)
	{
		if (lastHistogram == null || newSimulationSettings(settings, lastSimulationSettings))
		{
			IJ.showStatus("Simulating histogram ...");
			StopWatch sw = StopWatch.createStarted();
			// Ensure this is the same each time by using the same seed
			PseudoRandomGenerator random = getRandomGenerator().clone();
			random.setSeed(settings.getSeed());
			lastHistogram = simulateHistogram(settings, random);
			lastSimulationSettings = settings;
			IJ.showStatus("Simulated in " + sw.toString());
		}
		return lastHistogram;
	}

	private boolean newSimulationSettings(CameraModelAnalysisSettings s1, CameraModelAnalysisSettings s2)
	{
		// Check those settings for the simulation
		if (s1.getPhotons() != s2.getPhotons())
			return true;
		if (s1.getMode() != s2.getMode())
			return true;
		if (s1.getSamples() != s2.getSamples())
			return true;
		if (s1.getNoiseSamples() != s2.getNoiseSamples())
			return true;
		if (s1.getSeed() != s2.getSeed())
			return true;
		if (s1.getMode() == 0)
		{
			if (s1.getGain() != s2.getGain())
				return true;
			if (s1.getNoise() != s2.getNoise())
				return true;
		}
		else
		{
			if (s1.getEmGain() != s2.getEmGain())
				return true;
			if (s1.getEmNoise() != s2.getEmNoise())
				return true;
			if (s1.getEmSamples() != s2.getEmSamples())
				return true;
		}
		return false;
	}

	/**
	 * Simulate the histogram for fitting.
	 *
	 * @param settings
	 *            the settings
	 * @param random
	 *            the random
	 * @return The histogram
	 */
	private static IntHistogram simulateHistogram(CameraModelAnalysisSettings settings, PseudoRandomGenerator random)
	{
		IntHistogram h;
		switch (settings.getMode())
		{
			case 1:
				h = simulatePoissonGammaGaussian(settings, random);
				break;

			case 0:
				h = simulatePoissonGaussian(settings, random);
				break;

			default:
				throw new IllegalStateException();
		}
		return h;
	}

	/**
	 * Randomly generate a histogram from poisson-gamma-gaussian samples.
	 *
	 * @param settings
	 *            the settings
	 * @param random
	 *            the random
	 * @return The histogram
	 */
	private static IntHistogram simulatePoissonGammaGaussian(CameraModelAnalysisSettings settings,
			PseudoRandomGenerator random)
	{
		// Randomly sample
		PseudoRandomGenerator random2 = random.clone();

		PoissonDistribution poisson = new PoissonDistribution(random, settings.getPhotons(),
				PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);

		// Note that applying a separate EM-gain and then the camera gain later (as per Create Data)
		// is the same as applying the total gain in the gamma distribution and no camera gain 
		// later, i.e. the Gamma distribution is just squeezed.
		CustomGammaDistribution gamma = new CustomGammaDistribution(random.clone(), 1, settings.getEmGain(),
				GammaDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);

		final double noise = getEmNoise(settings);
		final int samples = settings.getSamples();
		final int emSamples = settings.getEmSamples();
		final int noiseSamples = (noise > 0) ? settings.getNoiseSamples() : 1;
		int[] sample = new int[samples * emSamples * noiseSamples];
		int count = 0;
		// Too fast for a ticker
		//Ticker ticker = Ticker.createStarted(new IJTrackProgress(), samples, false);
		for (int n = samples; n-- > 0;)
		{
			// Poisson
			double d = poisson.sample();

			if (d != 0)
			{
				// Gamma
				gamma.setShapeUnsafe(d);

				// Over-sample the Gamma
				for (int k = emSamples; k-- > 0;)
				{
					final double d2 = gamma.sample();

					// Over-sample the Gaussian
					for (int j = noiseSamples; j-- > 0;)
					{
						// Convert the sample to a count 
						sample[count++] = (int) Math.round(d2 + noise * random2.nextGaussian());
					}
				}
			}
			else
			{
				// Over-sample the Gaussian
				for (int j = noiseSamples; j-- > 0;)
				{
					// Convert the sample to a count 
					sample[count++] = (int) Math.round(noise * random2.nextGaussian());
				}
			}

			//ticker.tick();
		}

		//ticker.stop();

		return createHistogram(sample, count);
	}

	/**
	 * Randomly generate a histogram from poisson-gaussian samples.
	 *
	 * @param settings
	 *            the settings
	 * @param random
	 *            the random
	 * @return The histogram
	 */
	private static IntHistogram simulatePoissonGaussian(CameraModelAnalysisSettings settings,
			PseudoRandomGenerator random)
	{
		// Randomly sample
		PseudoRandomGenerator random2 = random.clone();

		PoissonDistribution poisson = new PoissonDistribution(random, settings.getPhotons(),
				PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);

		final double gain = settings.getGain();
		final double noise = getNoise(settings);
		final int samples = settings.getSamples();
		final int noiseSamples = (noise > 0) ? settings.getNoiseSamples() : 1;
		int[] sample = new int[samples * noiseSamples];
		int count = 0;
		// Too fast for a ticker
		//Ticker ticker = Ticker.createStarted(new IJTrackProgress(), samples, false);
		for (int n = samples; n-- > 0;)
		{
			// Poisson
			double d = poisson.sample();

			// Fixed camera gain
			d *= gain;

			// Over-sample the Gaussian
			for (int j = noiseSamples; j-- > 0;)
			{
				// Convert the sample to a count 
				sample[count++] = (int) Math.round(d + noise * random2.nextGaussian());
			}
			//ticker.tick();
		}

		//ticker.stop();

		return createHistogram(sample, count);
	}

	private static IntHistogram createHistogram(int[] sample, int count)
	{
		int[] limits = Maths.limits(sample);
		int min = limits[0];
		int max = limits[1];

		int[] h = new int[max - min + 1];
		for (int i = count; i-- > 0;)
			h[sample[i] - min]++;
		return new IntHistogram(h, min);
	}

	private static LikelihoodFunction getLikelihoodFunction(CameraModelAnalysisSettings settings)
	{
		double alpha = 1.0 / getTotalGain(settings);
		double noise = getReadNoise(settings);
		switch (settings.getModel())
		{
			case 0:
				return new PoissonFunction(alpha, false);
			case 1:
				return new PoissonFunction(alpha, true);
			case 2:
				return PoissonGaussianFunction2.createWithStandardDeviation(alpha, noise);
			case 3:
				return PoissonPoissonFunction.createWithStandardDeviation(alpha, noise);
			case 4:
				return new PoissonGammaGaussianFunction(alpha, noise);

			default:
				throw new IllegalStateException();
		}
	}

	private static double getNoise(CameraModelAnalysisSettings settings)
	{
		return settings.getNoise();
	}

	private static double getEmNoise(CameraModelAnalysisSettings settings)
	{
		return settings.getEmNoise();
	}

	private static double getTotalGain(CameraModelAnalysisSettings settings)
	{
		switch (settings.getMode())
		{
			case 0:
				return settings.getGain();
			case 1:
				return settings.getEmGain();
			default:
				throw new IllegalStateException();
		}
	}

	private static double getReadNoise(CameraModelAnalysisSettings settings)
	{
		switch (settings.getMode())
		{
			case 0:
				return getNoise(settings);
			case 1:
				return getEmNoise(settings);
			default:
				throw new IllegalStateException();
		}
	}

	private static double[][] cumulativeHistogram(IntHistogram histogram)
	{
		int[] h = histogram.h;
		double[] x = new double[h.length];
		double[] y = new double[x.length];
		double sum = 0;
		for (int i = 0; i < x.length; i++)
		{
			// The cumulative histogram represents the probability of all values up to this one.
			// However the histogram is discrete so this is the probability of all values up to 
			// but not including the next one.

			x[i] = histogram.getValue(i);
			sum += h[i];
			y[i] = sum;
		}
		for (int i = 0; i < y.length; i++)
		{
			y[i] /= sum;
		}
		y[y.length - 1] = 1; // Ensure total is 1
		return new double[][] { x, y };
	}

	private static double[][] cumulativeDistribution(CameraModelAnalysisSettings settings, double[][] cdf,
			LikelihoodFunction f)
	{
		// Q. How to match this is the discrete cumulative histogram using the continuous 
		// likelihood function:
		// 1. Compute integral up to the value
		// 2. Compute integral up to but not including the next value using trapezoid integration
		// 3. Compute integral up to but not including the next value using flat-top integration
		// Since the function will be used on continuous float data when fitting PSFs the best 
		// match for how it will perform in practice is a continuous (trapezoid) integration.
		// The simplest is a flat-top integration.

		// Compute the probability at each value
		double e = settings.getPhotons() * getTotalGain(settings);
		double[] x = cdf[0];
		double[] y = new double[x.length];
		for (int i = 0; i < x.length; i++)
			y[i] = f.likelihood(x[i], e);
		// Add more until the probability change is marginal
		final double delta = 1e-6;
		double sum = Maths.sum(y);
		TDoubleArrayList list = new TDoubleArrayList(y);
		for (int o = (int) x[x.length - 1] + 1;; o++)
		{
			double p = f.likelihood(o, e);
			if (p == 0)
				break;
			list.add(p);
			sum += p;
			if (p / sum < delta)
				break;
		}
		TDoubleArrayList list2 = new TDoubleArrayList(10);
		for (int o = (int) x[0] - 1;; o--)
		{
			double p = f.likelihood(o, e);
			if (p == 0)
				break;
			list2.add(p);
			sum += p;
			if (p / sum < delta)
				break;
		}
		// Insert at start
		double start = x[0];
		if (!list2.isEmpty())
		{
			start -= list2.size();
			list2.reverse();
			list.insert(0, list2.toArray());
		}

		y = list.toArray();
		x = SimpleArrayUtils.newArray(y.length, start, 1.0);

		// Simple flat-top integration
		sum = 0;
		for (int i = 0; i < y.length; i++)
		{
			sum += y[i];
			y[i] = sum;
		}

		return new double[][] { x, y };
	}

	/**
	 * Compute the Kolmogorov distance as the supremum (maximum)
	 * difference between the two cumulative probability distributions.
	 * https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test
	 * <p>
	 * Also compute the mean distance between the two CDFs over the range of CDF 1.
	 *
	 * @param cdf1
	 *            the cdf 1
	 * @param cdf2
	 *            the cdf 2
	 * @return [distance,value,mean distance]
	 */
	private static double[] getDistance(double[][] cdf1, double[][] cdf2)
	{
		// Find the offset
		int offset = 0;
		double[] x1 = cdf1[0];
		double[] x2 = cdf2[0];
		double[] y1 = cdf1[1];
		double[] y2 = cdf2[1];
		while (x2[offset] < x1[0])
			offset++;
		double distance = 0;
		double value = x1[0];
		double area = 0;
		for (int i = 0; i < x1.length; i++)
		{
			double d = Math.abs(y1[i] - y2[offset++]);
			if (distance < d)
			{
				distance = d;
				value = x1[i];
			}

			// Compute area:

			// This assumes both are discrete distributions
			area += d;

			// Note: This assumes but distributions are continuous between the values
			// and computes the actual area, including intersecting lines. 
			//if (i != 0)
			//{
			//	area += area(y1[i - 1], y1[i], y2[i], y2[i - 1]);
			//}
		}
		return new double[] { distance, value, area / x1.length };
	}

	private static final double[] areaX = { 0, 1, 1, 0 };

	@SuppressWarnings("unused")
	private static double area(double y1, double y2, double y3, double y4)
	{
		// Check if they cross
		if (!((y1 > y4 && y2 > y3) || (y1 < y4 && y2 < y3)))
		{
			double[] intersection = new double[2];
			if (Geometry.getIntersection(0, y1, 1, y2, 1, y3, 0, y4, intersection))
			{
				// Compute area as two triangles
				return Geometry.getArea(0, y1, 0, y4, intersection[0], intersection[1]) +
						Geometry.getArea(1, y2, 1, y3, intersection[0], intersection[1]);
			}
		}

		return Math.abs(Geometry.getArea(areaX, new double[] { y1, y2, y3, y4 }));
	}
}
