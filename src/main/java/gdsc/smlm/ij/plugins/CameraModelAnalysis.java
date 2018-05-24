package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Color;
import java.util.Arrays;

import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.CustomSimpsonIntegrator;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;
import org.apache.commons.math3.distribution.CustomGammaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;
import org.apache.commons.math3.util.FastMath;

import gdsc.core.ij.Utils;
import gdsc.core.math.Geometry;
import gdsc.core.threshold.IntHistogram;
import gdsc.core.utils.CachedRandomGenerator;
import gdsc.core.utils.Maths;
import gdsc.core.utils.PseudoRandomGenerator;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.smlm.data.NamedObject;
import gdsc.smlm.data.config.GUIProtos.CameraModelAnalysisSettings;
import gdsc.smlm.function.InterpolatedPoissonFunction;
import gdsc.smlm.function.LikelihoodFunction;
import gdsc.smlm.function.LogFactorial;
import gdsc.smlm.function.PoissonFunction;
import gdsc.smlm.function.PoissonGammaFunction;
import gdsc.smlm.function.PoissonGammaGaussianConvolutionFunction;
import gdsc.smlm.function.PoissonGammaGaussianFunction;
import gdsc.smlm.function.PoissonGammaGaussianFunction.ConvolutionMode;
import gdsc.smlm.function.PoissonGaussianConvolutionFunction;
import gdsc.smlm.function.PoissonGaussianFunction2;
import gdsc.smlm.function.PoissonPoissonFunction;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.utils.Convolution;
import gdsc.smlm.utils.GaussianKernel;
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
import ij.gui.Plot;
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
	// Add mode to compute the distance over a range of photons at a set gain and variance
	// Compute distance to simulation. Compute distance to another distribution. 

	private static final String TITLE = "Camera Model Analysis";
	private static final KolmogorovSmirnovTest kolmogorovSmirnovTest = new KolmogorovSmirnovTest();

	private CameraModelAnalysisSettings.Builder settings;

	private boolean extraOptions;
	private boolean dirty = true;
	private CameraModelAnalysisSettings lastSettings = null;
	private ExtendedGenericDialog gd;
	private IntHistogram lastHistogram = null;
	private double[][] floatHistogram = null;
	private CameraModelAnalysisSettings lastSimulationSettings = null;

	private static String[] MODE = { "CCD", "EM-CCD", "sCMOS" };
	private static final int MODE_CCD = 0;
	private static final int MODE_EM_CCD = 1;
	private static final int MODE_SCMOS = 2;

	//@formatter:off
	private enum Model implements NamedObject
	{
		///////////////
		// CCD / sCMOS 
		///////////////
		
		POISSON_PMF { public String getName() { return "Poisson PMF"; } },
		POISSON_DISRECTE { public String getName() { return "Poisson (Discrete)"; } },
		POISSON_CONTINUOUS { public String getName() { return "Poisson (Continuous)"; } },
		POISSON_GAUSSIAN_PDF { public String getName() { return "Poisson+Gaussian PDF integration"; } },
		
		// Best for CCD/sCMOS
		POISSON_GAUSSIAN_PMF { public String getName() { return "Poisson+Gaussian PMF integration"; } },
		// Saddle-point approximation.
		// Very good. Relatively worse than POISSON_GAUSSIAN_PMF at very low photons.
		POISSON_GAUSSIAN_APPROX { public String getName() { return "Poisson+Gaussian approximation"; } },
		 // Mixed Poisson distribution (Noise is added as a second Poisson)
		POISSON_POISSON { public String getName() { return "Poisson+Poisson"; } },
		
		///////////////
		// EMCCD 
		///////////////
		// There is no obvious best for EM-CCD:
		// This requires a full range test to determine the best function for which
		// parameters.
		
		// Good when no noise.
		// Under-estimates total probability when gain is low (<15) and photons are low (<2).
		POISSON_GAMMA_PMF { public String getName() { return "Poisson+Gamma PMF"; } },
		
		// Good but relatively worse as the read noise increases.
		// Requires full integration when read noise is low (<1) and photons are low (<5).
		POISSON_GAMMA_GAUSSIAN_APPROX { public String getName() { return "Poisson+Gamma+Gaussian approximation"; } },
		// Good when read noise is >>1.
		// Requires full integration when read noise is low (<1).
		POISSON_GAMMA_GAUSSIAN_PDF_INTEGRATION { public String getName() { return "Poisson+Gamma+Gaussian PDF integration"; } },
		// Good
		// Slow
		POISSON_GAMMA_GAUSSIAN_PMF_INTEGRATION { public String getName() { return "Poisson+Gamma+Gaussian PMF integration"; } },
		// Best for EM-CCD. 
		// Very robust (computes the full convolution
		// of the Gaussian and the Poisson-Gamma plus the delta function PMF contribution).
		POISSON_GAMMA_GAUSSIAN_SIMPSON_INTEGRATION { public String getName() { return "Poisson+Gamma+Gaussian Simpson's integration"; } },
		// Good
		POISSON_GAMMA_GAUSSIAN_LEGENDRE_GAUSS_INTEGRATION { public String getName() { return "Poisson+Gamma+Gaussian Legendre-Gauss integration"; } },
		// Independent implementation of POISSON_GAMMA_GAUSSIAN_PDF_INTEGRATION.
		// Good when read noise is >>1.
		// Requires full integration when read noise is low (<1).
		POISSON_GAMMA_GAUSSIAN_PDF_CONVOLUTION { public String getName() { return "Poisson+Gamma*Gaussian convolution"; } },
		;


		public String getShortName()
		{
			return getName();
		}		
		
		public static Model forNumber(int number)
		{
			Model[] values = Model.values();
			if (number < 0 || number >= values.length)
				number = 0;
			return values[number];
		}
	}
	
	private static String[] MODEL = SettingsManager.getNames((Object[]) Model.values());
	
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
				if (mode == MODE_CCD)
				{
					egd.addNumericField("Gain", settings.getGain(), 2, 6, "Count/electrons");
					egd.addNumericField("Noise", settings.getNoise(), 2, 6, "Count");
				}
				else if (mode == MODE_EM_CCD)
				{
					egd.addNumericField("Gain", settings.getEmGain(), 2, 6, "Count/electrons");
					egd.addNumericField("Noise", settings.getEmNoise(), 2, 6, "Count");
					egd.addNumericField("EM_samples", settings.getEmSamples(), 0);
				}
				else if (mode == MODE_SCMOS)
				{
					egd.addNumericField("Gain", settings.getCmosGain(), 2, 6, "Count/electrons");
					egd.addNumericField("Noise", settings.getCmosNoise(), 2, 6, "Count");
				}
				else
					throw new IllegalStateException();
				egd.showDialog(true, gd);
				if (egd.wasCanceled())
					return false;
				if (mode == MODE_CCD)
				{
					settings.setGain(egd.getNextNumber());
					settings.setNoise(egd.getNextNumber());
				}
				else if (mode == MODE_EM_CCD)
				{
					settings.setEmGain(egd.getNextNumber());
					settings.setEmNoise(egd.getNextNumber());
					settings.setEmSamples(Math.max(1, (int) egd.getNextNumber()));
				}
				else // MODE_SCMOS
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

		// KolmogorovSmirnovTest
		// n is the number of samples used to build the probability distribution.
		int n = (int) Maths.sum(h.h);

		// From KolmogorovSmirnovTest.kolmogorovSmirnovTest(RealDistribution distribution, double[] data, boolean exact):
		// Returns the p-value associated with the null hypothesis that data is a sample from distribution.
		// E.g. If p<0.05 then the null hypothesis is rejected and the data do not match the distribution.
		double p = Double.NaN;
		try
		{
			p = 1d - kolmogorovSmirnovTest.cdf(distance, n);
		}
		catch (Exception e)
		{
			// Ignore
		}

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
		plot.addLabel(0, 0, String.format("Distance=%s @ %.0f (Mean=%s) : p=%s", Utils.rounded(distance), value,
				Utils.rounded(area), Utils.rounded(p)));
		Utils.display(title, plot, Utils.NO_TO_FRONT, wo);

		// Show the histogram
		title = TITLE + " Histogram";
		plot = new Plot2(title, "Count", "Frequency");
		// Update X1 so that the histogram bars are centred over the x value
		for (int i = x1.length; i-- > 0;)
			x1[i] -= 0.5;
		plot.setLimits(x1[0] - 0.5, x1[x1.length - 1] + 1.5, 0, Maths.max(h.h) * 1.05);
		plot.setColor(Color.blue);
		//plot.addPoints(x2, y1b, Plot2.LINE);
		plot.addPoints(x1, SimpleArrayUtils.toDouble(h.h), Plot2.BAR);

		plot.setColor(Color.red);
		double[] x = floatHistogram[0].clone();
		double[] y = floatHistogram[1].clone();
		double scale = n / (Maths.sum(y) * (x[1] - x[0]));
		for (int i = 0; i < y.length; i++)
			y[i] *= scale;
		plot.addPoints(x, y, Plot2.LINE);

		plot.setColor(Color.black);
		plot.addLegend("Sample\nExpected");
		Utils.display(title, plot, Utils.NO_TO_FRONT, wo);

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

			// Convolve 
			//sw.reset();
			floatHistogram = convolveHistogram(settings);
			//IJ.log("Computed histogram ... " + sw.toString());
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
		if (s1.getMode() == MODE_CCD)
		{
			if (s1.getGain() != s2.getGain())
				return true;
			if (s1.getNoise() != s2.getNoise())
				return true;
		}
		else if (s1.getMode() == MODE_SCMOS)
		{
			if (s1.getCmosGain() != s2.getCmosGain())
				return true;
			if (s1.getCmosNoise() != s2.getCmosNoise())
				return true;
		}
		else // MODE_EM_CCD
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

				// Q. The Poisson-Gamma Function does not exactly match this
				// when the mean is around 1. The model from Ulbrich+Isacoff
				// is using a gamma function. This returns a continuous value.
				// Should it be rounded to make it a discrete PMF?
				// Or should the model actually be integrated.
				// Basically the value of the Poisson-Gamma from 0-0.5 should be
				// added to the delta function at c=0 for the count at c=0.
				// This would fix the function. 

				// Should we use the Tubb's model which uses:
				//final double shape = count;
				//final double scale = gain - 1 + 1 / shape;
				//final double electrons = random.nextGamma(shape, scale) - 1;
				//final double output = count + electrons; 

				// The Tubb's model is for additional electrons. So a count of 1
				// can never generate an output of 0. This does not fit the
				// Ulbrich+Isacoff model where c=0 is actually defined, i.e. you can
				// model zero output for EM-gain.

				// Over-sample the Gamma
				for (int k = emSamples; k-- > 0;)
				{
					//final double d2 = gamma.sample();
					// Do rounding to simulate a discrete PMF.
					final double d2 = round.round(gamma.sample());
					//final double d2 = (int)(gamma.sample());

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

	/**
	 * Convolve the histogram. The output is a discrete probability distribution.
	 *
	 * @param settings
	 *            the settings
	 * @return The histogram
	 */
	private static double[][] convolveHistogram(CameraModelAnalysisSettings settings)
	{
		final double LOWER = 1e-6;
		final double UPPER = 1 - LOWER;

		// Find the range of the Poisson
		PoissonDistribution poisson = new PoissonDistribution(random, settings.getPhotons(),
				PoissonDistribution.DEFAULT_EPSILON, PoissonDistribution.DEFAULT_MAX_ITERATIONS);
		int maxn = poisson.inverseCumulativeProbability(UPPER);

		final double gain = getGain(settings);
		final double noise = getReadNoise(settings);

		boolean debug = false;

		// Build the Probabity Mass/Density Function (PDF) of the distribution:
		// either a Poisson (PMF) or Poisson-Gamma (PDF). The PDF is 0 at all 
		// values apart from the step interval.
		// Note: The Poisson-Gamm is computed without the Dirac delta contribution 
		// at c=0. This allows correct convolution with the Gaussian of the dirac delta
		// and the rest of the Poisson-Gamma (so matching the simulation).
		TDoubleArrayList list = new TDoubleArrayList();
		double step;
		String name;

		int upsample = 100;

		// Store the Dirac delta value at c=0. This must be convolved separately.
		double dirac = 0;

		// EM-CCD 
		if (settings.getMode() == MODE_EM_CCD)
		{
			name = "Poisson-Gamma";

			final double m = gain;
			final double p = settings.getPhotons();

			dirac = PoissonGammaFunction.dirac(p);

			// Chose whether to compute a discrete PMF or a PDF using the approximation.
			// Note: The delta function at c=0 is from the PMF of the Poisson. So it is
			// a discrete contribution. This is omitted from the PDF and handled in
			// a separate convolution.
			boolean discrete = false; // noise != 0;
			if (discrete)
			{
				// Note: This is obsolete as the Poisson-Gamma function is continuous. 
				// Sampling it at integer intervals is not valid, especially for low gain.
				// The Poisson-Gamma PDF should be integrated to form a discrete PMF. 

				step = 1.0;

				CustomGammaDistribution gamma = new CustomGammaDistribution(random, settings.getPhotons(), gain,
						GammaDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);

				double upper;
				if (settings.getPhotons() < 20)
					upper = maxn;
				else
					// Approximate reasonable range of Poisson as a Gaussian
					upper = settings.getPhotons() + Math.sqrt(settings.getPhotons());

				gamma.setShapeUnsafe(upper);
				int maxc = (int) gamma.inverseCumulativeProbability(0.999);

				int minn = Math.max(1, poisson.inverseCumulativeProbability(LOWER));

				// See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.

				// Note this is not a convolution of a single Gamma distribution since the shape
				// is modified not the count. So it is a convolution of a distribution made with
				// a gamma of fixed count and variable shape.

				// The count=0 is a special case.
				list.add(PoissonGammaFunction.poissonGammaN(0, p, m));

				long total = (maxn - minn) * (long) maxc;

				if (total < 1000000)
				{
					// Full computation

					//G(c) = sum n {  (1 / n!) p^n e^-p (1 / ((n-1!)m^n)) c^n-1 e^-c/m };
					// Compute as a log
					// - log(n!) + n*log(p)-p -log((n-1)!) - n * log(m) + (n-1) * log(c) -c/m

					// Note: Both methods work

					LogFactorial.increaseTableMaxN(maxn);
					double[] f = new double[maxn + 1];
					final double logm = Math.log(m);
					final double logp = Math.log(p);
					for (int n = minn; n <= maxn; n++)
					{
						f[n] = -LogFactorial.logF(n) + n * logp - p - LogFactorial.logF(n - 1) - n * logm;
					}

					// Use Poisson + Gamma distribution
					//double[] pd = new double[maxn + 1];
					//CustomGammaDistribution[] gd = new CustomGammaDistribution[maxn + 1];
					//for (int n = minn; n <= maxn; n++)
					//{
					//	pd[n] = poisson.probability(n);
					//	gd[n] = new CustomGammaDistribution(null, n, m);
					//}

					//double total = list.getQuick(0);
					//double total2 = total;
					for (int c = 1; c <= maxc; c++)
					{
						double sum = 0;
						final double c_m = c / m;
						final double logc = Math.log(c);
						for (int n = minn; n <= maxn; n++)
						{
							sum += FastMath.exp(f[n] + (n - 1) * logc - c_m);
							//sum2 += pd[n] * gd[n].density(c);
						}
						list.add(sum);
						//total += sum;

						// This should match the approximation
						//double approx = PoissonGammaFunction.poissonGamma(c, p, m);
						//total2 += approx;
						//System.out.printf("c=%d sum=%g approx=%g error=%g\n", c, sum2, approx,
						//		gdsc.core.utils.DoubleEquality.relativeError(sum2, approx));
					}

					//System.out.printf("sum=%g approx=%g error=%g\n", total, total2,
					//		gdsc.core.utils.DoubleEquality.relativeError(total, total2));
				}
				else
				{
					// Approximate
					for (int c = 1; c <= maxc; c++)
						list.add(PoissonGammaFunction.poissonGammaN(c, p, m));
				}
			}
			else
			{
				// This integrates the PDF using the approximation and up-samples together.
				// Compute the sampling interval.
				step = 1.0 / upsample;
				upsample = 1; // Reset

				// Compute the integral of [-step/2:step/2] for each point.
				// Use trapezoid integration.
				double step_2 = step / 2;
				double prev = PoissonGammaFunction.poissonGammaN(0, p, m);
				double next = PoissonGammaFunction.poissonGammaN(step_2, p, m);
				list.add((prev + next) * 0.25);
				double max = 0;
				for (int i = 1;; i++)
				{
					// Each remaining point is modelling a PMF for the range [-step/2:step/2]
					prev = next;
					next = PoissonGammaFunction.poissonGammaN(i * step + step_2, p, m);
					double pp = (prev + next) * 0.5;
					if (max < pp)
						max = pp;
					if (pp / max < 1e-5)
					{
						// Use this if non-zero since it has been calculated
						if (pp != 0)
							list.add(pp);
						break;
					}
					list.add(pp);
				}
			}

			// Ensure the combined sum of PDF and Dirac is 1
			double expected = 1 - dirac;
			// Compute the sum using Simpson's rule:
			// Require an odd number to get an even number (n) of sub-intervals:
			if (list.size() % 2 == 0)
				list.add(0);
			double[] g = list.toArray();
			// Number of sub intervals
			int n = g.length - 1;
			double h = 1; // h = (a-b) / n = sub-interval width 
			double sum2 = 0, sum4 = 0;
			for (int j = 1; j <= n / 2 - 1; j++)
				sum2 += g[2 * j];
			for (int j = 1; j <= n / 2; j++)
				sum4 += g[2 * j - 1];
			double sum = (h / 3) * (g[0] + 2 * sum2 + 4 * sum4 + g[n]);
			// Check
			//System.out.printf("Sum=%g Expected=%g\n", sum * step, expected);
			SimpleArrayUtils.multiply(g, expected / sum);
			list.resetQuick();
			list.add(g);
		}
		else
		{
			name = "Poisson";
			// Apply fixed gain. Just change the step interval of the PMF.
			step = gain;

			for (int n = 0; n <= maxn; n++)
			{
				list.add(poisson.probability(n));
			}
			double p = poisson.probability(list.size());
			if (p != 0)
				list.add(p);
		}

		// Debug
		if (debug)
		{
			String title = name;
			Plot plot = new Plot(title, "x", "y", SimpleArrayUtils.newArray(list.size(), 0, step), list.toArray());
			Utils.display(title, plot);
		}

		double zero = 0;
		double[] g = list.toArray();

		// Sample Gaussian
		if (noise > 0)
		{
			step /= upsample;
			g = list.toArray();

			// Convolve with Gaussian kernel
			double[] kernel = GaussianKernel.makeGaussianKernel(Math.abs(noise) / step, 6, true);

			if (upsample != 1)
			{
				// Use scaled convolution. This is faster that zero filling distribution g.
				g = Convolution.convolve(kernel, g, upsample);
			}
			else
			{
				// The Poisson-Gamma may be stepped at low mean causing wrap artifacts in the FFT. 
				// This is a problem if most of the probability is in the Dirac. 
				// Otherwise it can be ignored.
				if (dirac > 0.01)
					g = Convolution.convolve(kernel, g);
				else
					g = Convolution.convolveFast(kernel, g);
			}
			
			// The convolution will have created a larger array so we must adjust the offset for this
			int radius = kernel.length / 2;
			zero -= radius * step;

			// Add convolution of the dirac delta function.
			if (dirac != 0)
			{
				// We only need to convolve the Gaussian at c=0
				for (int i = 0; i < kernel.length; i++)
					g[i] += kernel[i] * dirac;
			}

			// Debug
			if (debug)
			{
				String title = "Gaussian";
				Plot plot = new Plot(title, "x", "y", SimpleArrayUtils.newArray(kernel.length, -radius * step, step),
						kernel);
				Utils.display(title, plot);

				title = name + "-Gaussian";
				plot = new Plot(title, "x", "y", SimpleArrayUtils.newArray(g.length, zero, step), g);
				Utils.display(title, plot);
			}

			zero = downSampleCDFtoPMF(settings, list, step, zero, g, 1.0);

			g = list.toArray();
			zero = (int) Math.floor(zero);
			step = 1.0;
		}
		else
		{
			// No convolution means we have the Poisson PMF/Poisson-Gamma PDF already
			if (step != 1)
			{
				// Sample to 1.0 pixel step interval.
				if (settings.getMode() == MODE_EM_CCD)
				{
					// Poisson-Gamma PDF
					zero = downSampleCDFtoPMF(settings, list, step, zero, g, 1 - dirac);
					g = list.toArray();
					zero = (int) Math.floor(zero);

					// Add the dirac delta function.
					if (dirac != 0)
					{
						// Note: zero is the start of the x-axis. This value should be -1.
						assert (int) zero == -1;
						// Use as an offset to find the actual zero. 
						g[-(int) zero] += dirac;
					}
				}
				else
				{
					// Poisson PMF

					// Simple non-interpolated expansion.
					// This should be used when there is no Gaussian convolution.
					double[] pd = g;
					list.resetQuick();

					// Account for rounding.
					Round round = getRound(settings);

					int maxc = round.round(pd.length * step + 1);
					g = new double[maxc];
					for (int n = pd.length; n-- > 0;)
					{
						g[round.round(n * step)] += pd[n];
					}

					if (g[0] != 0)
					{
						list.add(0);
						list.add(g);
						g = list.toArray();
						zero--;
					}
				}

				step = 1.0;
			}
			else
			{
				// Add the dirac delta function.
				list.setQuick(0, list.getQuick(0) + dirac);
			}
		}

		return new double[][] { SimpleArrayUtils.newArray(g.length, zero, step), g };
	}

	private static double downSampleCDFtoPMF(CameraModelAnalysisSettings settings, TDoubleArrayList list, double step,
			double zero, double[] g, double sum)
	{
		// Down-sample to 1.0 pixel step interval.

		// Build cumulative distribution.
		double lowerSum = 0;
		for (int i = 0; i < g.length; i++)
		{
			lowerSum += g[i];
			g[i] = lowerSum;
		}
		for (int i = 0; i < g.length; i++)
			g[i] *= sum / lowerSum;
		g[g.length - 1] = sum;

		double offset = (settings.getRoundDown()) ? 0 : -0.5;

		// Note the interpolation of the CDF is good when the step is much smaller than 1.
		// When the step is above 1 then the gain is likely to be very low and thus 
		// unrealistic for modelling. This case is ignored.

		// Subtract the CDF to the upper bounds from the CDF of the lower bound
		// to get the discrete PMF

		//		// Pad the CDF to avoid index-out-of bounds during interpolation
		int padSize = 0;
		//		padSize = (int) Math.ceil(1 / step) + 2;
		//		list.resetQuick();
		//		list.add(new double[padSize]);
		//		list.add(g);
		//		for (int i = padSize; i-- > 0;)
		//			list.add(1);
		//double[] pd = list.toArray();
		double[] pd = g;

		list.resetQuick();

		double[] x = SimpleArrayUtils.newArray(pd.length, zero - padSize * step, step);

		// Q. If the EM-CCD the distribution may have a Dirac delta at c=0 which 
		// could break interpolation using a spline?
		UnivariateInterpolator in =
				//(settings.getMode() == MODE_EM_CCD) ? new LinearInterpolator() : 
				new SplineInterpolator();
		UnivariateFunction f = in.interpolate(x, pd);

		int bound = (int) Math.floor(zero);

		list.add(0);
		zero--;
		double upperSum = 0;
		double min = x[0];
		double max = x[x.length - 1];
		while (upperSum < sum)
		{
			bound++;

			// Find the point at which the CDF should be computed
			lowerSum = upperSum;
			double point = bound + offset;
			if (point < min)
				upperSum = 0;
			else if (point > max)
				upperSum = sum;
			else
				upperSum = f.value(point);
			list.add(upperSum - lowerSum);
		}
		list.add(0);
		return zero;
	}

	private static LikelihoodFunction getLikelihoodFunction(CameraModelAnalysisSettings settings)
	{
		double alpha = 1.0 / getGain(settings);
		double noise = getReadNoise(settings);
		Model model = Model.forNumber(settings.getModel());
		switch (model)
		{
			case POISSON_PMF:
				return new PoissonFunction(alpha);
			case POISSON_DISRECTE:
				return new InterpolatedPoissonFunction(alpha, false);
			case POISSON_CONTINUOUS:
				return new InterpolatedPoissonFunction(alpha, true);
			case POISSON_GAUSSIAN_PDF:
			case POISSON_GAUSSIAN_PMF:
				PoissonGaussianConvolutionFunction f1 = PoissonGaussianConvolutionFunction
						.createWithStandardDeviation(alpha, noise);
				f1.setComputePMF(model == Model.POISSON_GAUSSIAN_PMF);
				return f1;
			case POISSON_GAUSSIAN_APPROX:
				return PoissonGaussianFunction2.createWithStandardDeviation(alpha, noise);
			case POISSON_POISSON:
				return PoissonPoissonFunction.createWithStandardDeviation(alpha, noise);

			case POISSON_GAMMA_GAUSSIAN_PDF_CONVOLUTION:
				return PoissonGammaGaussianConvolutionFunction.createWithStandardDeviation(alpha, noise);

			case POISSON_GAMMA_PMF:
				return PoissonGammaFunction.createWithAlpha(alpha);

			case POISSON_GAMMA_GAUSSIAN_APPROX:
			case POISSON_GAMMA_GAUSSIAN_PDF_INTEGRATION:
			case POISSON_GAMMA_GAUSSIAN_PMF_INTEGRATION:
			case POISSON_GAMMA_GAUSSIAN_SIMPSON_INTEGRATION:
			case POISSON_GAMMA_GAUSSIAN_LEGENDRE_GAUSS_INTEGRATION:
				PoissonGammaGaussianFunction f2 = new PoissonGammaGaussianFunction(alpha, noise);
				f2.setMinimumProbability(0);
				f2.setConvolutionMode(getConvolutionMode(model));
				// The function should return a PMF/PDF depending on how it is used
				f2.setPmfMode(!settings.getSimpsonIntegration());

				return f2;

			default:
				throw new IllegalStateException();
		}
	}

	private static boolean isPoissonGammaLikelihoodFunction(CameraModelAnalysisSettings settings)
	{
		Model model = Model.forNumber(settings.getModel());
		switch (model)
		{
			case POISSON_GAMMA_GAUSSIAN_PDF_CONVOLUTION:
			case POISSON_GAMMA_PMF:
			case POISSON_GAMMA_GAUSSIAN_APPROX:
			case POISSON_GAMMA_GAUSSIAN_PDF_INTEGRATION:
			case POISSON_GAMMA_GAUSSIAN_PMF_INTEGRATION:
			case POISSON_GAMMA_GAUSSIAN_SIMPSON_INTEGRATION:
			case POISSON_GAMMA_GAUSSIAN_LEGENDRE_GAUSS_INTEGRATION:
				return true;

			default:
				return false;
		}
	}

	private static ConvolutionMode getConvolutionMode(Model model)
	{
		switch (model)
		{
			case POISSON_GAMMA_GAUSSIAN_APPROX:
				return ConvolutionMode.APPROXIMATION;
			case POISSON_GAMMA_GAUSSIAN_PDF_INTEGRATION:
				return ConvolutionMode.DISCRETE_PDF;
			case POISSON_GAMMA_GAUSSIAN_PMF_INTEGRATION:
				return ConvolutionMode.DISCRETE_PMF;
			case POISSON_GAMMA_GAUSSIAN_SIMPSON_INTEGRATION:
				return ConvolutionMode.SIMPSON_PDF;
			case POISSON_GAMMA_GAUSSIAN_LEGENDRE_GAUSS_INTEGRATION:
				return ConvolutionMode.LEGENDRE_GAUSS_PDF;

			default:
				throw new IllegalStateException();
		}
	}

	private static double getGain(CameraModelAnalysisSettings settings)
	{
		switch (settings.getMode())
		{
			case MODE_CCD:
				return settings.getGain();
			case MODE_EM_CCD:
				return settings.getEmGain();
			case MODE_SCMOS:
				return settings.getCmosGain();
			default:
				throw new IllegalStateException();
		}
	}

	private static double getReadNoise(CameraModelAnalysisSettings settings)
	{
		switch (settings.getMode())
		{
			case MODE_CCD:
				return settings.getNoise();
			case MODE_EM_CCD:
				return settings.getEmNoise();
			case MODE_SCMOS:
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

	private static class CachingUnivariateFunction implements UnivariateFunction
	{
		final LikelihoodFunction fun;
		final double p;
		final TDoubleArrayList list = new TDoubleArrayList();

		public CachingUnivariateFunction(LikelihoodFunction fun, double p)
		{
			this.fun = fun;
			this.p = p;
		}

		public double value(double x)
		{
			double v = fun.likelihood(x, p);
			list.add(x);
			list.add(v);
			return v;
		}

		public void reset()
		{
			list.resetQuick();
		}
	}

	private static double[][] cumulativeDistribution(CameraModelAnalysisSettings settings, double[][] cdf,
			final LikelihoodFunction fun)
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
		double e = settings.getPhotons();
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
			int c0 = -1;
			double dirac = 0;
			int minc = 0, maxc = 0;
			CachingUnivariateFunction uf = null;

			if (settings.getMode() == MODE_EM_CCD && isPoissonGammaLikelihoodFunction(settings))
			{
				// A spike is expected at c=0 due to the Dirac delta contribution.
				// This breaks integration, especially when noise < 0.1.
				// Fix by integrating around c=0 fully then integrating the rest.
				c0 = Arrays.binarySearch(x, 0);
				double noise = getReadNoise(settings);
				final double p = settings.getPhotons();
				if (noise == 0)
				{
					// Pure Poisson-Gamma. Just subtract the delta, do the simple integration
					// below and add the delta back. Only functions that support noise==0
					// will be allowed so this solution works.
					dirac = PoissonGammaFunction.dirac(p);
					if (c0 != -1)
						y[c0] -= dirac;
				}
				else
				{
					// Fix integration around c=0 using the range of the Gaussian
					minc = (int) Math.max(x[0], Math.floor(-5 * noise));
					maxc = (int) Math.min(x[x.length - 1], Math.ceil(5 * noise));
					uf = new CachingUnivariateFunction(fun, p);
				}
			}

			// Use Simpson's integration with n=4 to get the integral of the probability 
			// over the range of each count.
			int n = 4;
			int n_2 = n / 2;
			double h = 1.0 / n;

			// Note the Poisson-Gamma function cannot be integrated with the 
			// Dirac delta function at c==0

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

			// Fix Poisson-Gamma ...
			if (c0 != -1)
			{
				if (uf != null)
				{
					// Convolved Poisson-Gamma. Fix in the range of the Gaussian around c=0 
					final double relativeAccuracy = 1e-4;
					final double absoluteAccuracy = 1e-8;
					CustomSimpsonIntegrator in = new CustomSimpsonIntegrator(relativeAccuracy, absoluteAccuracy, 3,
							CustomSimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
					double lower = (settings.getRoundDown()) ? 0 : -0.5;
					double upper = lower + 1;
					// Switch from c<=maxc to c<maxc. Avoid double computation at minc==maxc
					if (maxc != minc)
						maxc++;
					maxc++;
					for (int c = minc, i = Arrays.binarySearch(x, minc); c < maxc; c++, i++)
					{
						uf.reset();
						try
						{
							y[i] = in.integrate(2000, uf, c + lower, c + upper);
						}
						catch (TooManyEvaluationsException ex)
						{
							//System.out.printf("Integration failed: c=%g-%g\n", c + lower, c + upper);
							// Q. Is the last sum valid?
							if (in.getLastSum() > 0)
							{
								y[i] = in.getLastSum();
							}
							else
							{
								// Otherwise use all the cached values to compute a sum
								// using the trapezoid rule. This will underestimate the sum.

								// Note: The Simpson integrator will have computed the edge values
								// as the first two values in the cache.
								double[] g = uf.list.toArray();
								double dx = (g[3] - g[1]) / in.getN();
								n = 1 + 2 * ((int) in.getN());
								sum = 0;
								for (int j = 4; j < n; j += 2)
									sum += g[j];
								y[i] = (g[0] + g[2] + 2 * sum) / dx;
							}
						}
					}
				}
				else
				{
					// Pure Poisson-Gamma. Just add back the delta.
					y[c0] += dirac;
				}
			}
		}

		// Simple flat-top integration
		sum = 0;
		for (int i = 0; i < y.length; i++)
		{
			sum += y[i];
			y[i] = sum;
		}

		// Check if sum is approximately 1
		//System.out.printf("gain=%g, sum=%g\n", getGain(settings), sum);		

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
