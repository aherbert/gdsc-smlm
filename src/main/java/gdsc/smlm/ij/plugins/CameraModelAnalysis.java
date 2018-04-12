package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Color;
import java.util.Arrays;

import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.math3.distribution.CustomGammaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

import gdsc.core.ij.Utils;
import gdsc.core.math.Geometry;
import gdsc.core.threshold.IntHistogram;
import gdsc.core.utils.CachedRandomGenerator;
import gdsc.core.utils.Maths;
import gdsc.core.utils.PseudoRandomGenerator;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.smlm.data.config.GUIProtos.CameraModelAnalysisSettings;
import gdsc.smlm.function.LikelihoodFunction;
import gdsc.smlm.function.PoissonFunction;
import gdsc.smlm.function.PoissonGammaGaussianConvolutionFunction;
import gdsc.smlm.function.PoissonGammaGaussianFunction;
import gdsc.smlm.function.PoissonGaussianConvolutionFunction;
import gdsc.smlm.function.PoissonGaussianFunction2;
import gdsc.smlm.function.PoissonPoissonFunction;
import gdsc.smlm.ij.settings.SettingsManager;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
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
import ij.plugin.WindowOrganiser;
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

	private static String[] MODE = { "CCD", "EM-CCD", "sCMOS" };
	//@formatter:off
	private static String[] MODEL = { 
			"Poisson (Discrete)", 
			"Poisson (Continuous)", 
			"Poisson*Gaussian convolution",
			"Poisson+Gaussian approximation",
			"Poisson+Poisson", 
			"Poisson*Gamma*Gaussian convolution", 
			"Poisson+Gamma+Gaussian approximation", 
			"Poisson*Gamma*Gaussian convolution2", 
			};

	private static abstract class Round
	{
		abstract int round(double d);
	}
	private static class RoundDown extends Round
	{
		int round(double d) { return (int) Math.floor(d); }
	}
	private static class RoundUpDown extends Round
	{
		int round(double d) { return (int) Math.round(d); }
	}
	private static Round ROUND_DOWN = new RoundDown();
	private static Round ROUND_UP_DOWN = new RoundUpDown();
	private static Round getRound(CameraModelAnalysisSettings settings)
	{
		return (settings.getRoundDown()) ? ROUND_DOWN : ROUND_UP_DOWN;
	}
	//@formatter:on

	private static CachedRandomGenerator random = null;
	private static int currentSeed = 0;

	private CachedRandomGenerator getRandomGenerator()
	{
		if (random == null || currentSeed != settings.getSeed())
		{
			currentSeed = settings.getSeed();
			// Ensure some bits are set in the default seed of zero
			long seed = (currentSeed == 0) ? Double.doubleToRawLongBits(Math.PI) : currentSeed;
			random = new CachedRandomGenerator(settings.getSamples(), new Well19937c(seed));
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
				else if (mode == 2)
				{
					egd.addNumericField("Gain", settings.getCmosGain(), 2, 6, "Count/electrons");
					egd.addNumericField("Noise", settings.getCmosNoise(), 2, 6, "Count");
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
				else if (mode == 1)
				{
					settings.setEmGain(egd.getNextNumber());
					settings.setEmNoise(egd.getNextNumber());
					settings.setEmSamples(Math.max(1, (int) egd.getNextNumber()));
				}
				else
				{
					settings.setCmosGain(egd.getNextNumber());
					settings.setCmosNoise(egd.getNextNumber());
				}
				return true;
			}
		});
		if (extraOptions)
			gd.addNumericField("Seed", settings.getSeed(), 0);
		gd.addNumericField("Samples", settings.getSamples(), 0);
		gd.addNumericField("Noise_samples", settings.getNoiseSamples(), 0);
		gd.addCheckbox("Round_down", settings.getRoundDown());
		gd.addChoice("Model", MODEL, settings.getModel());
		gd.addCheckbox("Full_integration", settings.getSimpsonIntegration());
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
		settings.setRoundDown(gd.getNextBoolean());
		settings.setModel(gd.getNextChoiceIndex());
		settings.setSimpsonIntegration(gd.getNextBoolean());
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
			boolean ok = false;
			try
			{
				ok = execute();
			}
			catch (Exception ex)
			{
				// Catch as this is run within a AWT dispatch thread 
				Utils.log(TITLE + "Error: " + ex.getMessage());
			}
			finally
			{
				if (!ok)
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
		if (!(getGain(settings) > 0))
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
		WindowOrganiser wo = new WindowOrganiser();
		String title = TITLE + " CDF";
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
		Utils.display(title, plot, 0, wo);

		// Show the histogram
		title = TITLE + " Histogram";
		plot = new Plot2(title, "Count", "Frequency");
		plot.setLimits(x1[0] - 0.5, x1[x1.length - 1] + 1.5, 0, Maths.max(h.h) * 1.05);
		plot.setColor(Color.blue);
		//plot.addPoints(x2, y1b, Plot2.LINE);
		plot.addPoints(x1, SimpleArrayUtils.toDouble(h.h), Plot2.BAR);
		Utils.display(title, plot, 0, wo);

		wo.tile();

		return true;
	}

	private IntHistogram getHistogram(CameraModelAnalysisSettings settings)
	{
		if (lastHistogram == null || newSimulationSettings(settings, lastSimulationSettings))
		{
			IJ.showStatus("Simulating histogram ...");
			StopWatch sw = StopWatch.createStarted();
			CachedRandomGenerator random = getRandomGenerator();
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
		if (s1.getRoundDown() != s2.getRoundDown())
			return true;
		if (s1.getMode() == 0)
		{
			if (s1.getGain() != s2.getGain())
				return true;
			if (s1.getNoise() != s2.getNoise())
				return true;
		}
		else if (s1.getMode() == 2)
		{
			if (s1.getCmosGain() != s2.getCmosGain())
				return true;
			if (s1.getCmosNoise() != s2.getCmosNoise())
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
	private static IntHistogram simulateHistogram(CameraModelAnalysisSettings settings, CachedRandomGenerator random)
	{
		IntHistogram h;
		switch (settings.getMode())
		{
			// EM-CCD
			case 1:
				h = simulatePoissonGammaGaussian(settings, random);
				break;

			// CCD or sCMOS				
			case 2:
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
			CachedRandomGenerator random)
	{
		int[] poissonSample = simulatePoisson(settings, random);

		// Randomly sample re-using the sequence
		RandomGenerator random2 = random.getPseudoRandomGenerator();

		final double gain = getGain(settings);

		// Note that applying a separate EM-gain and then the camera gain later (as per Create Data)
		// is the same as applying the total gain in the gamma distribution and no camera gain 
		// later, i.e. the Gamma distribution is just squeezed.
		CustomGammaDistribution gamma = new CustomGammaDistribution(random.getPseudoRandomGenerator(), 1, gain,
				GammaDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);

		final double noise = getReadNoise(settings);
		final int samples = settings.getSamples();
		final int emSamples = settings.getEmSamples();
		final int noiseSamples = (noise > 0) ? settings.getNoiseSamples() : 1;
		int[] sample = new int[samples * emSamples * noiseSamples];
		int count = 0;
		Round round = getRound(settings);
		for (int n = poissonSample.length; n-- > 0;)
		{
			if (poissonSample[n] != 0)
			{
				// Gamma
				gamma.setShapeUnsafe(poissonSample[n]);

				// Over-sample the Gamma
				for (int k = emSamples; k-- > 0;)
				{
					final double d2 = gamma.sample();

					// Over-sample the Gaussian
					for (int j = noiseSamples; j-- > 0;)
					{
						// Convert the sample to a count 
						sample[count++] = round.round(d2 + noise * random2.nextGaussian());
					}
				}
			}
			else
			{
				// Still over-sample the Gamma even though it was zero
				for (int k = emSamples; k-- > 0;)
				{
					// Over-sample the Gaussian
					for (int j = noiseSamples; j-- > 0;)
					{
						// Convert the sample to a count 
						sample[count++] = round.round(noise * random2.nextGaussian());
					}
				}
			}
		}

		return createHistogram(sample, count);
	}

	private static int[] simulatePoisson(CameraModelAnalysisSettings settings, CachedRandomGenerator random)
	{
		// Ensure we reuse random numbers if possible
		random.reset();
		PoissonDistribution poisson = new PoissonDistribution(random, settings.getPhotons(),
				PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
		return poisson.sample(settings.getSamples());
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
			CachedRandomGenerator random)
	{
		int[] poissonSample = simulatePoisson(settings, random);

		// Randomly sample re-using the sequence
		PseudoRandomGenerator random2 = random.getPseudoRandomGenerator();

		final double gain = getGain(settings);
		final double noise = getReadNoise(settings);
		final int samples = settings.getSamples();
		final int noiseSamples = (noise > 0) ? settings.getNoiseSamples() : 1;
		int[] sample = new int[samples * noiseSamples];
		int count = 0;
		Round round = getRound(settings);
		for (int n = poissonSample.length; n-- > 0;)
		{
			// Fixed camera gain
			double d = poissonSample[n] * gain;

			// Over-sample the Gaussian
			for (int j = noiseSamples; j-- > 0;)
			{
				// Convert the sample to a count 
				sample[count++] = round.round(d + noise * random2.nextGaussian());
			}
		}

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
		double alpha = 1.0 / getGain(settings);
		double noise = getReadNoise(settings);
		switch (settings.getModel())
		{
			case 0:
				return new PoissonFunction(alpha, false);
			case 1:
				return new PoissonFunction(alpha, true);
			case 2:
				return PoissonGaussianConvolutionFunction.createWithStandardDeviation(alpha, noise);
			case 3:
				return PoissonGaussianFunction2.createWithStandardDeviation(alpha, noise);
			case 4:
				return PoissonPoissonFunction.createWithStandardDeviation(alpha, noise);

			// Add PoissonGammaGaussianConvolution ...
			// Get the Poisson-Gamma from EMGainAnalysis. 
			// Add a convolution with a range as in the PoissonGaussianConvolutionFunction

			case 5:
			case 6:
				PoissonGammaGaussianFunction f = new PoissonGammaGaussianFunction(alpha, noise);
				// The full-integration here does not appear to work !
				f.setUseApproximation(settings.getModel() == 6);
				f.setUseSimpleIntegration(false);
				f.setMinimumProbability(0);
				return f;

			case 7:
				return PoissonGammaGaussianConvolutionFunction.createWithStandardDeviation(alpha, noise);

			default:
				throw new IllegalStateException();
		}
	}

	private static double getGain(CameraModelAnalysisSettings settings)
	{
		switch (settings.getMode())
		{
			case 0:
				return settings.getGain();
			case 1:
				return settings.getEmGain();
			case 2:
				return settings.getCmosGain();
			default:
				throw new IllegalStateException();
		}
	}

	private static double getReadNoise(CameraModelAnalysisSettings settings)
	{
		switch (settings.getMode())
		{
			case 0:
				return settings.getNoise();
			case 1:
				return settings.getEmNoise();
			case 2:
				return settings.getCmosNoise();
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
			LikelihoodFunction fun)
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
		double e = settings.getPhotons() * getGain(settings);
		double[] x = cdf[0];
		double[] y = new double[x.length];
		for (int i = 0; i < x.length; i++)
			y[i] = fun.likelihood(x[i], e);
		// Add more until the probability change is marginal
		final double delta = 1e-6;
		double sum = Maths.sum(y);
		TDoubleArrayList list = new TDoubleArrayList(y);
		for (int o = (int) x[x.length - 1] + 1;; o++)
		{
			double p = fun.likelihood(o, e);
			list.add(p);
			if (p == 0)
				break;
			sum += p;
			if (p / sum < delta)
				break;
		}
		TDoubleArrayList list2 = new TDoubleArrayList(10);
		for (int o = (int) x[0] - 1;; o--)
		{
			double p = fun.likelihood(o, e);
			list2.add(p);
			if (p == 0)
				break;
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

		if (settings.getSimpsonIntegration())
		{
			// Use Simpson's integration with n=4 to get the integral of the probability 
			// over the range of each count.
			int n = 4;
			int n_2 = n / 2;
			double h = 1.0 / n;

			// Compute the extra function points
			double[] f = new double[y.length * n + 1];
			{
				int i = f.length;

				final int mod;
				if (settings.getRoundDown())
				{
					// Do this final value outside the loop as y[i/n] does not exists
					mod = 0;
					i--;
					f[i] = fun.likelihood(start + i * h, e);
				}
				else
				{
					// Used when computing for rounding up/down
					mod = n_2;
					start -= n_2 * h;
				}

				while (i-- > 0)
				{
					if (i % n == mod)
						f[i] = y[i / n];
					else
						f[i] = fun.likelihood(start + i * h, e);
				}
			}

			// Compute indices for the integral
			TIntArrayList tmp = new TIntArrayList(n);
			for (int j = 1; j <= n_2 - 1; j++)
				tmp.add(2 * j);
			int[] i2 = tmp.toArray();
			tmp.resetQuick();
			for (int j = 1; j <= n_2; j++)
				tmp.add(2 * j - 1);
			int[] i4 = tmp.toArray();

			// Compute integral
			for (int i = 0; i < y.length; i++)
			{
				int a = i * n;
				int b = a + n;
				double s = f[a] + f[b];
				for (int j = i2.length; j-- > 0;)
					s += 2 * f[a + i2[j]];
				for (int j = i4.length; j-- > 0;)
					s += 4 * f[a + i4[j]];
				s *= h / 3;
				//System.out.printf("y[%d] = %f => %f\n", i, y[i], s);
				y[i] = s;
			}
		}

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
