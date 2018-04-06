package gdsc.smlm.ij.plugins;

import java.awt.Color;
import java.awt.Point;
import java.awt.Rectangle;
import java.util.Arrays;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.distribution.CustomGammaDistribution;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateFunctionMappingAdapter;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CustomPowellOptimizer;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well44497b;
import org.apache.commons.math3.util.FastMath;

import gdsc.core.ij.Utils;
import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Maths;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.core.utils.StoredDataStatistics;
import gdsc.smlm.function.Bessel;
import gdsc.smlm.function.LikelihoodFunction;
import gdsc.smlm.function.PoissonFunction;
import gdsc.smlm.function.PoissonGammaGaussianFunction;
import gdsc.smlm.function.PoissonGaussianFunction2;
import gdsc.smlm.utils.Convolution;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.Plot2;
import ij.gui.PlotWindow;
import ij.gui.Roi;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

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

/**
 * Analysis a white light image from an EM-CCD camera, construct a histogram of pixel intensity and fit the histogram to
 * obtain the bias, EM-gain, read noise and photons per pixel.
 * <p>
 * See Ulbrich & Isacoff (2007) Subunit counting in membrane-bound proteins. Nature Methods 4, 319-321 (Supplementary
 * Information).
 */
public class EMGainAnalysis implements PlugInFilter
{
	private static final String TITLE = "EM-Gain Analysis";
	private final int FLAGS = DOES_8G | DOES_16 | NO_CHANGES | NO_UNDO;
	private final int MINIMUM_PIXELS = 1000000;

	private static double bias = 500, gain = 40, noise = 3;
	private static boolean _simulate = false, showApproximation = false, relativeDelta = false;
	private static String[] APPROXIMATION = { "PoissonGammaGaussian", "PoissonGamma", "PoissonGaussian", "Poisson" };
	private static int approximation = 0;
	private boolean simulate = false, extraOptions = false;
	private static double _photons = 1, _bias = 500, _gain = 40, _noise = 3;
	private static double head = 0.01, tail = 0.025, _offset = 0;
	private static int simulationSize = 20000;
	private static boolean usePDF = false;

	private ImagePlus imp;
	private double offset = 0;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	public int setup(String arg, ImagePlus imp)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		extraOptions = Utils.isExtraOptions();

		if ("pmf".equals(arg))
		{
			plotPMF();
			return DONE;
		}

		if (imp == null && !extraOptions)
		{
			IJ.noImage();
			return DONE;
		}
		this.imp = imp;
		return showDialog();
	}

	public void run(ImageProcessor ip)
	{
		// Calculate the histogram
		final int[] h = (simulate) ? simulateHistogram(usePDF ? 0 : 1) : buildHistogram(imp);

		// We need > 10^7 pixels from flat white-light images under constant exposure ...
		final int size = getSize(h);
		if (imp != null)
		{
			Roi roi = imp.getRoi();
			Rectangle bounds;
			if (roi == null)
				bounds = new Rectangle(0, 0, imp.getWidth(), imp.getHeight());
			else
				bounds = roi.getBounds();
			Utils.log("Analysing %s [x=%d,y=%d,width=%d,height=%d]", imp.getTitle(), bounds.x, bounds.y, bounds.width,
					bounds.height);
		}
		Utils.log("Histogram contains %d pixels", size);
		if (size < MINIMUM_PIXELS)
			Utils.log("WARNING : Recommend at least %d pixels (%sx more)", MINIMUM_PIXELS,
					Utils.rounded((double) MINIMUM_PIXELS / size));

		fit(h);
	}

	/**
	 * Simulate the histogram for fitting
	 * 
	 * @param method
	 *            0 - sample from the fitted PDF, 1 - sample from a Poisson-Gamma-Gaussian
	 * @return The histogram
	 */
	private int[] simulateHistogram(int method)
	{
		IJ.showStatus("Simulating histogram ...");
		int[] h;
		switch (method)
		{
			case 1:
				h = simulateFromPoissonGammaGaussian();
				break;

			case 0:
			default:
				h = simulateFromPDF();
		}
		IJ.showStatus("");
		return h;
	}

	/**
	 * Random sample from the cumulative probability distribution function that is used during fitting
	 * 
	 * @return The histogram
	 */
	private int[] simulateFromPDF()
	{
		final double[] g = pdf(0, _photons, _gain, _noise, (int) _bias);

		// Debug this
		double[] x = SimpleArrayUtils.newArray(g.length, 0, 1.0);
		Utils.display(TITLE + " PDF", new Plot(TITLE + " PDF", "ADU", "P", x, Arrays.copyOf(g, g.length)));

		// Get cumulative probability
		double sum = 0;
		for (int i = 0; i < g.length; i++)
		{
			final double p = g[i];
			g[i] += sum;
			sum += p;
		}

		for (int i = 0; i < g.length; i++)
			g[i] /= sum;
		g[g.length - 1] = 1; // Ensure value of 1 at the end

		// Randomly sample
		RandomGenerator random = new Well44497b(System.currentTimeMillis() + System.identityHashCode(this));
		int[] h = new int[g.length];
		final int steps = simulationSize;
		for (int n = 0; n < steps; n++)
		{
			if (n % 64 == 0)
				IJ.showProgress(n, steps);
			final double p = random.nextDouble();
			for (int i = 0; i < g.length; i++)
				if (p <= g[i])
				{
					h[i]++;
					break;
				}
		}
		return h;
	}

	/**
	 * Randomly generate a histogram from poisson-gamma-gaussian samples
	 * 
	 * @return The histogram
	 */
	private int[] simulateFromPoissonGammaGaussian()
	{
		// Randomly sample
		RandomGenerator random = new Well44497b(System.currentTimeMillis() + System.identityHashCode(this));

		PoissonDistribution poisson = new PoissonDistribution(random, _photons, PoissonDistribution.DEFAULT_EPSILON,
				PoissonDistribution.DEFAULT_MAX_ITERATIONS);

		CustomGammaDistribution gamma = new CustomGammaDistribution(random, _photons, _gain,
				GammaDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);

		final int steps = simulationSize;
		int[] sample = new int[steps];
		for (int n = 0; n < steps; n++)
		{
			if (n % 64 == 0)
				IJ.showProgress(n, steps);

			// Poisson
			double d = poisson.sample();

			// Gamma
			if (d > 0)
			{
				gamma.setShapeUnsafe(d);
				d = gamma.sample();
			}

			// Gaussian
			d += _noise * random.nextGaussian();

			// Convert the sample to a count 
			sample[n] = (int) Math.round(d + _bias);
		}

		int max = Maths.max(sample);
		int[] h = new int[max + 1];
		for (int s : sample)
			h[s]++;
		return h;
	}

	/**
	 * Build a histogram using pixels within the image ROI
	 * 
	 * @param image
	 *            The image
	 * @return The image histogram
	 */
	private static int[] buildHistogram(ImagePlus imp)
	{
		ImageStack stack = imp.getImageStack();
		Roi roi = imp.getRoi();
		int[] data = getHistogram(stack.getProcessor(1), roi);
		for (int n = 2; n <= stack.getSize(); n++)
		{
			int[] tmp = getHistogram(stack.getProcessor(n), roi);
			for (int i = 0; i < tmp.length; i++)
				data[i] += tmp[i];
		}

		// Avoid super-saturated pixels by using 98% of histogram
		long sum = 0;
		for (int i = 0; i < data.length; i++)
		{
			sum += data[i];
		}

		final long sum2 = (long) (sum * 0.99);
		sum = 0;
		for (int i = 0; i < data.length; i++)
		{
			sum += data[i];
			if (sum > sum2)
			{
				for (; i < data.length; i++)
					data[i] = 0;
				break;
			}
		}

		return data;
	}

	private static int[] getHistogram(ImageProcessor ip, Roi roi)
	{
		ip.setRoi(roi);
		int[] h = ip.getHistogram();
		return h;
	}

	/**
	 * Fit the EM-gain distribution (Gaussian * Gamma)
	 * 
	 * @param h
	 *            The distribution
	 */
	private void fit(int[] h)
	{
		final int[] limits = limits(h);
		final double[] x = getX(limits);
		final double[] y = getY(h, limits);

		Plot2 plot = new Plot2(TITLE, "ADU", "Frequency");
		double yMax = Maths.max(y);
		plot.setLimits(limits[0], limits[1], 0, yMax);
		plot.setColor(Color.black);
		plot.addPoints(x, y, Plot2.DOT);
		Utils.display(TITLE, plot);

		// Estimate remaining parameters. 
		// Assuming a gamma_distribution(shape,scale) then mean = shape * scale
		// scale = gain
		// shape = Photons = mean / gain
		double mean = getMean(h) - bias;
		// Note: if the bias is too high then the mean will be negative. Just move the bias.
		while (mean < 0)
		{
			bias -= 1;
			mean += 1;
		}
		double photons = mean / gain;

		if (simulate)
			Utils.log("Simulated bias=%d, gain=%s, noise=%s, photons=%s", (int) _bias, Utils.rounded(_gain),
					Utils.rounded(_noise), Utils.rounded(_photons));

		Utils.log("Estimate bias=%d, gain=%s, noise=%s, photons=%s", (int) bias, Utils.rounded(gain),
				Utils.rounded(noise), Utils.rounded(photons));

		final int max = (int) x[x.length - 1];
		double[] g = pdf(max, photons, gain, noise, (int) bias);

		plot.setColor(Color.blue);
		plot.addPoints(x, g, Plot2.LINE);
		Utils.display(TITLE, plot);

		// Perform a fit
		CustomPowellOptimizer o = new CustomPowellOptimizer(1e-6, 1e-16, 1e-6, 1e-16);
		double[] startPoint = new double[] { photons, gain, noise, bias };
		int maxEval = 3000;

		String[] paramNames = { "Photons", "Gain", "Noise", "Bias" };
		// Set bounds
		double[] lower = new double[] { 0, 0.5 * gain, 0, bias - noise };
		double[] upper = new double[] { 2 * photons, 2 * gain, gain, bias + noise };

		// Restart until converged.
		// TODO - Maybe fix this with a better optimiser. This needs to be tested on real data.
		PointValuePair solution = null;
		for (int iter = 0; iter < 3; iter++)
		{
			IJ.showStatus("Fitting histogram ... Iteration " + iter);

			try
			{
				// Basic Powell optimiser
				MultivariateFunction fun = getFunction(limits, y, max, maxEval);
				PointValuePair optimum = o.optimize(new MaxEval(maxEval), new ObjectiveFunction(fun), GoalType.MINIMIZE,
						new InitialGuess((solution == null) ? startPoint : solution.getPointRef()));
				if (solution == null || optimum.getValue() < solution.getValue())
				{
					double[] point = optimum.getPointRef();
					// Check the bounds
					for (int i = 0; i < point.length; i++)
					{
						if (point[i] < lower[i] || point[i] > upper[i])
						{
							throw new RuntimeException(
									String.format("Fit out of of estimated range: %s %f", paramNames[i], point[i]));
						}
					}
					solution = optimum;
				}
			}
			catch (Exception e)
			{
				IJ.log("Powell error: " + e.getMessage());
				if (e instanceof TooManyEvaluationsException)
				{
					maxEval = (int) (maxEval * 1.5);
				}
			}
			try
			{
				// Bounded Powell optimiser
				MultivariateFunction fun = getFunction(limits, y, max, maxEval);
				MultivariateFunctionMappingAdapter adapter = new MultivariateFunctionMappingAdapter(fun, lower, upper);
				PointValuePair optimum = o.optimize(new MaxEval(maxEval), new ObjectiveFunction(adapter),
						GoalType.MINIMIZE, new InitialGuess(
								adapter.boundedToUnbounded((solution == null) ? startPoint : solution.getPointRef())));
				double[] point = adapter.unboundedToBounded(optimum.getPointRef());
				optimum = new PointValuePair(point, optimum.getValue());

				if (solution == null || optimum.getValue() < solution.getValue())
				{
					solution = optimum;
				}
			}
			catch (Exception e)
			{
				IJ.log("Bounded Powell error: " + e.getMessage());
				if (e instanceof TooManyEvaluationsException)
				{
					maxEval = (int) (maxEval * 1.5);
				}
			}
		}

		IJ.showStatus("");
		IJ.showProgress(1);

		if (solution == null)
		{
			Utils.log("Failed to fit the distribution");
			return;
		}

		double[] point = solution.getPointRef();
		photons = point[0];
		gain = point[1];
		noise = point[2];
		bias = (int) Math.round(point[3]);
		String label = String.format("Fitted bias=%d, gain=%s, noise=%s, photons=%s", (int) bias, Utils.rounded(gain),
				Utils.rounded(noise), Utils.rounded(photons));
		Utils.log(label);

		if (simulate)
		{
			Utils.log("Relative Error bias=%s, gain=%s, noise=%s, photons=%s",
					Utils.rounded(relativeError(bias, _bias)), Utils.rounded(relativeError(gain, _gain)),
					Utils.rounded(relativeError(noise, _noise)), Utils.rounded(relativeError(photons, _photons)));
		}

		// Show the PoissonGammaGaussian approximation
		double[] f = null;
		if (showApproximation)
		{
			f = new double[x.length];
			PoissonGammaGaussianFunction fun = new PoissonGammaGaussianFunction(1.0 / gain, noise);
			final double expected = photons * gain;
			for (int i = 0; i < f.length; i++)
			{
				f[i] = fun.likelihood(x[i] - bias, expected);
				//System.out.printf("x=%d, g=%f, f=%f, error=%f\n", (int) x[i], g[i], f[i],
				//		gdsc.smlm.fitting.utils.DoubleEquality.relativeError(g[i], f[i]));
			}
			yMax = Maths.maxDefault(yMax, f);
		}

		// Replot
		g = pdf(max, photons, gain, noise, (int) bias);
		plot = new Plot2(TITLE, "ADU", "Frequency");
		plot.setLimits(limits[0], limits[1], 0, yMax * 1.05);
		plot.setColor(Color.black);
		plot.addPoints(x, y, Plot2.DOT);
		plot.setColor(Color.red);
		plot.addPoints(x, g, Plot2.LINE);

		plot.addLabel(0, 0, label);

		if (showApproximation)
		{
			plot.setColor(Color.blue);
			plot.addPoints(x, f, Plot2.LINE);
		}

		Utils.display(TITLE, plot);
	}

	private double relativeError(double a, double b)
	{
		//return gdsc.smlm.utils.DoubleEquality.relativeError(a, b);
		final double d = a - b; // Math.abs(a - b);
		return d / b;
	}

	private MultivariateFunction getFunction(final int[] limits, final double[] y, final int max, final int maxEval)
	{
		MultivariateFunction fun = new MultivariateFunction()
		{
			int eval = 0;

			public double value(double[] point)
			{
				IJ.showProgress(++eval, maxEval);
				if (Utils.isInterrupted())
					throw new TooManyEvaluationsException(maxEval);
				// Compute the sum of squares between the two functions
				double photons = point[0];
				double gain = point[1];
				double noise = point[2];
				int bias = (int) Math.round(point[3]);
				//System.out.printf("[%d] = %s\n", eval, Arrays.toString(point));
				final double[] g = pdf(max, photons, gain, noise, bias);
				double ss = 0;
				for (int c = limits[0]; c <= limits[1]; c++)
				{
					final double d = g[c] - y[c];
					ss += d * d;
				}
				return ss;
			}
		};
		return fun;
	}

	/**
	 * Calculate the probability density function for EM-gain. The maximum count to evaluate is calculated dynamically
	 * so that the cumulative probability does not change.
	 * <p>
	 * See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
	 * 
	 * @param p
	 *            The average number of photons per pixel input to the EM-camera
	 * @param m
	 *            The multiplication factor (gain)
	 * @return The PDF
	 */
	private double[] pdfEMGain(final double p, final double m)
	{
		StoredDataStatistics stats = new StoredDataStatistics(100);
		stats.add(FastMath.exp(-p));
		for (int c = 1;; c++)
		{
			final double g = pEMGain(c, p, m);
			stats.add(g);
			final double delta = g / stats.getSum();
			if (delta < 1e-5)
				break;
		}
		return stats.getValues();
	}

	private static final double twoSqrtPi = 2 * Math.sqrt(Math.PI);

	/**
	 * Calculate the probability density function for EM-gain.
	 * <p>
	 * See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
	 * 
	 * @param c
	 *            The count to evaluate
	 * @param p
	 *            The average number of photons per pixel input to the EM-camera
	 * @param m
	 *            The multiplication factor (gain)
	 * @return The PDF
	 */
	private double pEMGain(int c, double p, double m)
	{
		// The default evaluation is:
		//return Math.sqrt(p / (c * m)) * FastMath.exp(-c / m - p) * Bessel.I1(2 * Math.sqrt(c * p / m));

		// Bessel.I1(x) -> Infinity
		// The current implementation of Bessel.I1(x) is Infinity at x==710

		// This has been fixed in the PoissonGammaGaussianFunction class so copy that code.
		final double alpha = 1 / m;

		final double cij = c;
		final double eta = p;

		// Any observed count above zero
		if (cij > 0.0)
		{
			// The observed count converted to photons
			final double nij = alpha * cij;

			// The current implementation of Bessel.I1(x) is Infinity at x==710
			// The limit on eta * nij is therefore (709/2)^2 = 125670.25
			if (eta * nij > 10000)
			{
				// Approximate Bessel function i1(x) when using large x:
				// i1(x) ~ exp(x)/sqrt(2*pi*x)
				// However the entire equation is logged (creating transform),
				// evaluated then raised to e to prevent overflow error on 
				// large exp(x)

				final double transform = 0.5 * Math.log(alpha * eta / cij) - nij - eta + 2 * Math.sqrt(eta * nij) -
						Math.log(twoSqrtPi * Math.pow(eta * nij, 0.25));
				return FastMath.exp(transform);
			}
			else
			{
				// Second part of equation 135
				return Math.sqrt(alpha * eta / cij) * FastMath.exp(-nij - eta) * Bessel.I1(2 * Math.sqrt(eta * nij));
			}
		}
		else if (cij == 0.0)
		{
			return FastMath.exp(-eta);
		}
		else
		{
			return 0;
		}
	}

	/**
	 * Calculate the probability density function for EM-gain.
	 * <p>
	 * See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
	 * 
	 * @param max
	 *            The maximum count to evaluate
	 * @param p
	 *            The average number of photons per pixel input to the EM-camera
	 * @param m
	 *            The multiplication factor (gain)
	 * @return The PDF
	 */
	private double[] pdfEMGain(final int max, final double p, final double m)
	{
		if (max == 0)
			return pdfEMGain(p, m);
		double[] g = new double[max + 1];
		g[0] = FastMath.exp(-p);
		for (int c = 1; c <= max; c++)
		{
			g[c] = pEMGain(c, p, m);
			if (g[c] == 0)
				break;
		}
		return g;
	}

	/**
	 * Calculate the probability density function for EM-gain, convolve with a Gaussian and then add a constant offset.
	 * <p>
	 * See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI.
	 * 
	 * @param max
	 *            The maximum count to evaluate
	 * @param p
	 *            The average number of photons per pixel input to the EM-camera
	 * @param m
	 *            The multiplication factor (gain)
	 * @param s
	 *            The read noise (Gaussian standard deviation)
	 * @param c0
	 *            The constant offset (bias)
	 * @return The PDF
	 */
	private double[] pdf(final int max, final double p, final double m, final double s, int c0)
	{
		double[] g = pdfEMGain(max, p, m);
		double[] gg;

		if (s > 0)
		{
			// Convolve with Gaussian kernel up to 4 times the standard deviation
			final int radius = (int) Math.ceil(Math.abs(s) * 4) + 1;
			double[] kernel = new double[2 * radius + 1];
			final double norm = -0.5 / (s * s);
			for (int i = 0, j = radius, jj = radius; j < kernel.length; i++, j++, jj--)
				kernel[j] = kernel[jj] = FastMath.exp(norm * i * i);
			// Normalise
			double sum = 0;
			for (int j = 0; j < kernel.length; j++)
				sum += kernel[j];
			for (int j = 0; j < kernel.length; j++)
				kernel[j] /= sum;

			gg = Convolution.convolveFast(g, kernel);
			// The convolution will have created a larger array so we must adjust the offset for this
			c0 -= radius;
		}
		else
			gg = g;

		// Pad with constant c0
		double[] g0 = new double[gg.length + c0];
		for (int c = 0; c < gg.length; c++)
			if (c + c0 >= 0)
				g0[c + c0] = gg[c];
		return g0;
	}

	private int[] limits(int[] h)
	{
		int min = 0;
		while (h[min] == 0)
			min++;
		int max = h.length - 1;
		while (h[max] == 0)
			max--;
		return new int[] { min, max };
	}

	private double[] getX(int[] limits)
	{
		final int min = 0; //limits[0];
		final int range = limits[1] - min + 1;
		final double[] x = new double[range];
		for (int i = 0; i < range; i++)
			x[i] = min + i;
		return x;
	}

	private double[] getY(int[] h, int[] limits)
	{
		final int min = 0; // limits[0];
		final int range = limits[1] - min + 1;
		final double[] y = new double[range];
		double sum = 0;
		for (int i = 0; i < range; i++)
		{
			y[i] = h[min + i];
			sum += y[i];
		}
		for (int i = 0; i < range; i++)
			y[i] /= sum;
		return y;
	}

	private int getSize(int[] h)
	{
		int size = 0;
		for (int i = 0; i < h.length; i++)
			size += h[i];
		return size;
	}

	private double getMean(int[] h)
	{
		int size = 0;
		double sum = 0;
		for (int i = 0; i < h.length; i++)
		{
			size += h[i];
			sum += h[i] * i;
		}
		return sum / size;
	}

	private int showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Analyse the white-light histogram of an image stack to determine EM-gain parameters.\n \n" +
				"See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321 (Supplementary Information).");

		if (extraOptions)
		{
			gd.addCheckbox("Simulate", _simulate);
			gd.addNumericField("Bias", _bias, 0);
			gd.addNumericField("Gain", _gain, 2);
			gd.addNumericField("Noise", _noise, 2);
			gd.addNumericField("Photons", _photons, 2);
			gd.addNumericField("Samples", simulationSize, 0);
			gd.addCheckbox("Sample_PDF", usePDF);
		}

		gd.addNumericField("Bias_estimate", bias, 0);
		gd.addNumericField("Gain_estimate", gain, 2);
		gd.addNumericField("Noise_estimate", noise, 2);
		gd.addCheckbox("Show_approximation", showApproximation);
		gd.showDialog();

		if (gd.wasCanceled())
			return DONE;

		if (extraOptions)
		{
			simulate = _simulate = gd.getNextBoolean();
			_bias = gd.getNextNumber();
			_gain = gd.getNextNumber();
			_noise = FastMath.abs(gd.getNextNumber());
			_photons = FastMath.abs(gd.getNextNumber());
			simulationSize = (int) FastMath.abs(gd.getNextNumber());
			usePDF = gd.getNextBoolean();
			if (gd.invalidNumber() || _bias < 0 || _gain < 1 || _photons == 0 || simulationSize == 0)
				return DONE;
		}

		bias = gd.getNextNumber();
		gain = gd.getNextNumber();
		noise = FastMath.abs(gd.getNextNumber());
		showApproximation = gd.getNextBoolean();

		if (gd.invalidNumber() || bias < 0 || gain < 1)
			return DONE;

		return (simulate) ? NO_IMAGE_REQUIRED : FLAGS;
	}

	private void plotPMF()
	{
		if (!showPMFDialog())
			return;

		final int gaussWidth = 5;
		int dummyBias = (int) Math.max(500, gaussWidth * _noise + 1);

		double[] pmf = pdf(0, _photons, _gain, _noise, dummyBias);
		double[] x = SimpleArrayUtils.newArray(pmf.length, 0, 1.0);
		double yMax = Maths.max(pmf);

		// Get the approximation
		LikelihoodFunction fun;
		double myNoise = _noise;
		switch (approximation)
		{
			case 3:
				fun = new PoissonFunction(1.0 / _gain, true);
				break;
			case 2:
				// The mean does not matter (as normalisation is done dynamically for 
				// PoissonGaussianFunction.likelihood(double, double) so just use zero
				//fun = PoissonGaussianFunction.createWithStandardDeviation(1.0 / _gain, 0, _noise);
				
				// Use adaptive normalisation
				fun = PoissonGaussianFunction2.createWithStandardDeviation(1.0 / _gain, _noise);
				break;
			case 1:
				myNoise = 0;
			case 0:
			default:
				PoissonGammaGaussianFunction myFun = new PoissonGammaGaussianFunction(1.0 / _gain, myNoise);
				myFun.setMinimumProbability(0);
				fun = myFun;
		}
		double expected = _photons;
		if (offset != 0)
			expected += offset * expected / 100.0;
		expected *= _gain;
		
		// Get CDF
		double sum = 0;
		double sum2 = 0;
		double[] f = new double[x.length];
		double[] cdf1 = new double[pmf.length];
		double[] cdf2 = new double[pmf.length];
		for (int i = 0; i < cdf1.length; i++)
		{
			sum += pmf[i];
			cdf1[i] = sum;
			// Adjust the x-values to remove the dummy bias
			x[i] -= dummyBias;
			f[i] =  fun.likelihood(x[i], expected);
			sum2 += f[i];
			cdf2[i] = sum2;
		}
		
		// Truncate x for plotting
		int max = 0;
		sum = 0;
		double p = 1 - tail;
		while (sum < p && max < pmf.length)
		{
			sum += pmf[max];
			if (sum > 0.5 && pmf[max] == 0)
				break;
			max++;
		}

		int min = pmf.length;
		sum = 0;
		p = 1 - head;
		while (sum < p && min > 0)
		{
			min--;
			sum += pmf[min];
			if (sum > 0.5 && pmf[min] == 0)
				break;
		}

		//int min = (int) (dummyBias - gaussWidth * _noise);
		pmf = Arrays.copyOfRange(pmf, min, max);
		x = Arrays.copyOfRange(x, min, max);
		f = Arrays.copyOfRange(f, min, max);

		if (showApproximation)
			yMax = Maths.maxDefault(yMax, f);

		String label = String.format("Gain=%s, noise=%s, photons=%s", Utils.rounded(_gain), Utils.rounded(_noise),
				Utils.rounded(_photons));

		Plot2 plot = new Plot2("PMF", "ADUs", "p");
		plot.setLimits(x[0], x[x.length - 1], 0, yMax);
		plot.setColor(Color.red);
		plot.addPoints(x, pmf, Plot2.LINE);
		if (showApproximation)
		{
			plot.setColor(Color.blue);
			plot.addPoints(x, f, Plot2.LINE);
		}

		plot.setColor(Color.magenta);
		plot.drawLine(_photons * _gain, 0, _photons * _gain, yMax);
		plot.setColor(Color.black);
		plot.addLabel(0, 0, label);
		PlotWindow win1 = Utils.display("PMF", plot);

		// Plot the difference between the actual and approximation
		double[] delta = new double[f.length];
		for (int i = 0; i < f.length; i++)
		{
			if (pmf[i] == 0 && f[i] == 0)
				continue;
			if (relativeDelta)
				delta[i] = DoubleEquality.relativeError(f[i], pmf[i]) * Math.signum(f[i] - pmf[i]);
			else
				delta[i] = f[i] - pmf[i];
		}

		Plot2 plot2 = new Plot2("PMF delta", "ADUs", (relativeDelta) ? "Relative delta" : "delta");
		double[] limits = Maths.limits(delta);
		plot2.setLimits(x[0], x[x.length - 1], limits[0], limits[1]);
		plot2.setColor(Color.red);
		plot2.addPoints(x, delta, Plot2.LINE);
		plot2.setColor(Color.magenta);
		plot2.drawLine(_photons * _gain, limits[0], _photons * _gain, limits[1]);
		plot2.setColor(Color.black);
		plot2.addLabel(0, 0, label + ((offset == 0) ? "" : ", expected = " + Utils.rounded(expected / _gain)));
		PlotWindow win2 = Utils.display("PMF delta", plot2);

		if (Utils.isNewWindow())
		{
			Point p2 = win1.getLocation();
			p2.y += win1.getHeight();
			win2.setLocation(p2);
		}

		// Plot the CDF of each distribution.
		// Compute the Kolmogorov distance as the supremum (maximum) 
		// difference between the two cumulative probability distributions.
		// https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test
		double kolmogorovDistance = 0;
		int xd = 0;
		for (int i = 0; i < cdf1.length; i++)
		{
			double dist = Math.abs(cdf1[i] - cdf2[i]);
			if (kolmogorovDistance < dist)
			{
				kolmogorovDistance = dist;
				xd = i;
			}
		}
		xd -= dummyBias;
		cdf1 = Arrays.copyOfRange(cdf1, min, max);
		cdf2 = Arrays.copyOfRange(cdf2, min, max);
		
		Plot2 plot3 = new Plot2("CDF", "ADUs", "p");
		yMax = 1.05;
		plot3.setLimits(x[0], x[x.length - 1], 0, yMax);
		plot3.setColor(Color.red);
		plot3.addPoints(x, cdf1, Plot2.LINE);
		plot3.setColor(Color.blue);
		plot3.addPoints(x, cdf2, Plot2.LINE);

		plot3.setColor(Color.magenta);
		plot3.drawLine(_photons * _gain, 0, _photons * _gain, yMax);
		plot3.drawDottedLine(xd, 0, xd, yMax, 2);
		plot3.setColor(Color.black);
		plot3.addLabel(0, 0, label + ", Kolmogorov distance = " + Utils.rounded(kolmogorovDistance) + " @ " + xd);
		PlotWindow win3 = Utils.display("CDF", plot3);
		
		if (Utils.isNewWindow())
		{
			Point p2 = win1.getLocation();
			p2.x += win1.getWidth();
			win3.setLocation(p2);
		}
	}

	private boolean showPMFDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Plot the probability mass function for EM-gain");

		gd.addNumericField("Gain", _gain, 2);
		gd.addNumericField("Noise", _noise, 2);
		gd.addNumericField("Photons", _photons, 2);
		gd.addChoice("Approx", APPROXIMATION, APPROXIMATION[approximation]);
		gd.addCheckbox("Show_approximation", showApproximation);
		if (extraOptions)
			gd.addNumericField("Approximation_offset (%)", _offset, 2);
		gd.addNumericField("Remove_head", head, 3);
		gd.addNumericField("Remove_tail", tail, 3);
		gd.addCheckbox("Relative_delta", relativeDelta);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		_gain = gd.getNextNumber();
		_noise = FastMath.abs(gd.getNextNumber());
		_photons = FastMath.abs(gd.getNextNumber());
		approximation = gd.getNextChoiceIndex();
		showApproximation = gd.getNextBoolean();
		if (extraOptions)
			offset = _offset = gd.getNextNumber();
		head = FastMath.abs(gd.getNextNumber());
		tail = FastMath.abs(gd.getNextNumber());
		relativeDelta = gd.getNextBoolean();

		if (gd.invalidNumber() || _bias < 0 || _gain < 1 || _photons == 0 || tail > 0.5 || head > 0.5)
			return false;

		return true;
	}
}
