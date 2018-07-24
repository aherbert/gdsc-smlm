/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package uk.ac.sussex.gdsc.smlm.ij.plugins;

import java.awt.Color;
import java.awt.Point;
import java.awt.Rectangle;
import java.util.Arrays;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateFunctionMappingAdapter;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well44497b;
import org.apache.commons.math3.util.FastMath;

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
import uk.ac.sussex.gdsc.core.ij.Utils;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.Maths;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.smlm.function.LikelihoodFunction;
import uk.ac.sussex.gdsc.smlm.function.PoissonFunction;
import uk.ac.sussex.gdsc.smlm.function.PoissonGammaFunction;
import uk.ac.sussex.gdsc.smlm.function.PoissonGammaGaussianFunction;
import uk.ac.sussex.gdsc.smlm.function.PoissonGaussianFunction2;
import uk.ac.sussex.gdsc.smlm.math3.distribution.CustomGammaDistribution;
import uk.ac.sussex.gdsc.smlm.math3.optim.nonlinear.scalar.noderiv.CustomPowellOptimizer;
import uk.ac.sussex.gdsc.smlm.utils.Convolution;

/**
 * Analysis a white light image from an EM-CCD camera, construct a histogram of pixel intensity and fit the histogram to
 * obtain the bias, EM-gain, read noise and photons per pixel.
 * <p>
 * See Ulbrich &amp; Isacoff (2007) Subunit counting in membrane-bound proteins. Nature Methods 4, 319-321 (Supplementary
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
	@Override
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

	@Override
	public void run(ImageProcessor ip)
	{
		// Calculate the histogram
		final int[] h = (simulate) ? simulateHistogram(usePDF ? 0 : 1) : buildHistogram(imp);

		// We need > 10^7 pixels from flat white-light images under constant exposure ...
		final int size = getSize(h);
		if (imp != null)
		{
			final Roi roi = imp.getRoi();
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
		final double step = getStepSize(_photons, _gain, _noise);

		final PDF pdf = pdf(0, step, _photons, _gain, _noise);

		// Debug this
		final double[] g = pdf.p;
		final double[] x = pdf.x;
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
		final RandomGenerator random = new Well44497b(System.currentTimeMillis() + System.identityHashCode(this));
		final int bias = (int) _bias;
		final int[] bins = new int[x.length];
		for (int i = 0; i < x.length; i++)
			bins[i] = bias + (int) x[i];
		final int[] h = new int[bins[bins.length - 1] + 1];
		final int steps = simulationSize;
		for (int n = 0; n < steps; n++)
		{
			if (n % 64 == 0)
				IJ.showProgress(n, steps);
			final double p = random.nextDouble();
			int i = binarySearch(g, p);
			if (i < 0)
				i = -(i + 1);
			h[bins[i]]++;

			//for (int i = 0; i < g.length; i++)
			//	if (p <= g[i])
			//	{
			//		h[i]++;
			//		break;
			//	}
		}
		return h;
	}

	private static int binarySearch(double[] a, double key)
	{
		int low = 0;
		int high = a.length - 1;

		while (low <= high)
		{
			final int mid = (low + high) >>> 1;
			final double midVal = a[mid];

			if (midVal < key)
				low = mid + 1; // Neither val is NaN, thisVal is smaller
			else if (midVal > key)
				high = mid - 1; // Neither val is NaN, thisVal is larger
			else
				return mid; // Key found
		}
		return -(low + 1); // key not found.
	}

	/**
	 * Randomly generate a histogram from poisson-gamma-gaussian samples
	 *
	 * @return The histogram
	 */
	private int[] simulateFromPoissonGammaGaussian()
	{
		// Randomly sample
		final RandomGenerator random = new Well44497b(System.currentTimeMillis() + System.identityHashCode(this));

		final PoissonDistribution poisson = new PoissonDistribution(random, _photons, PoissonDistribution.DEFAULT_EPSILON,
				PoissonDistribution.DEFAULT_MAX_ITERATIONS);

		final CustomGammaDistribution gamma = new CustomGammaDistribution(random, _photons, _gain,
				GammaDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);

		final int steps = simulationSize;
		final int[] sample = new int[steps];
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

		final int max = Maths.max(sample);
		final int[] h = new int[max + 1];
		for (final int s : sample)
			h[s]++;
		return h;
	}

	/**
	 * Build a histogram using pixels within the image ROI
	 *
	 * @param imp
	 *            The image
	 * @return The image histogram
	 */
	private static int[] buildHistogram(ImagePlus imp)
	{
		final ImageStack stack = imp.getImageStack();
		final Roi roi = imp.getRoi();
		final int[] data = getHistogram(stack.getProcessor(1), roi);
		for (int n = 2; n <= stack.getSize(); n++)
		{
			final int[] tmp = getHistogram(stack.getProcessor(n), roi);
			for (int i = 0; i < tmp.length; i++)
				data[i] += tmp[i];
		}

		// Avoid super-saturated pixels by using 98% of histogram
		long sum = 0;
		for (int i = 0; i < data.length; i++)
			sum += data[i];

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
		final int[] h = ip.getHistogram();
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
		plot.addPoints(x, y, Plot.DOT);
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
		plot.addPoints(x, g, Plot.LINE);
		Utils.display(TITLE, plot);

		// Perform a fit
		final CustomPowellOptimizer o = new CustomPowellOptimizer(1e-6, 1e-16, 1e-6, 1e-16);
		final double[] startPoint = new double[] { photons, gain, noise, bias };
		int maxEval = 3000;

		final String[] paramNames = { "Photons", "Gain", "Noise", "Bias" };
		// Set bounds
		final double[] lower = new double[] { 0, 0.5 * gain, 0, bias - noise };
		final double[] upper = new double[] { 2 * photons, 2 * gain, gain, bias + noise };

		// Restart until converged.
		// TODO - Maybe fix this with a better optimiser. This needs to be tested on real data.
		PointValuePair solution = null;
		for (int iter = 0; iter < 3; iter++)
		{
			IJ.showStatus("Fitting histogram ... Iteration " + iter);

			try
			{
				// Basic Powell optimiser
				final MultivariateFunction fun = getFunction(limits, y, max, maxEval);
				final PointValuePair optimum = o.optimize(new MaxEval(maxEval), new ObjectiveFunction(fun), GoalType.MINIMIZE,
						new InitialGuess((solution == null) ? startPoint : solution.getPointRef()));
				if (solution == null || optimum.getValue() < solution.getValue())
				{
					final double[] point = optimum.getPointRef();
					// Check the bounds
					for (int i = 0; i < point.length; i++)
						if (point[i] < lower[i] || point[i] > upper[i])
							throw new RuntimeException(
									String.format("Fit out of of estimated range: %s %f", paramNames[i], point[i]));
					solution = optimum;
				}
			}
			catch (final Exception e)
			{
				IJ.log("Powell error: " + e.getMessage());
				if (e instanceof TooManyEvaluationsException)
					maxEval = (int) (maxEval * 1.5);
			}
			try
			{
				// Bounded Powell optimiser
				final MultivariateFunction fun = getFunction(limits, y, max, maxEval);
				final MultivariateFunctionMappingAdapter adapter = new MultivariateFunctionMappingAdapter(fun, lower, upper);
				PointValuePair optimum = o.optimize(new MaxEval(maxEval), new ObjectiveFunction(adapter),
						GoalType.MINIMIZE, new InitialGuess(
								adapter.boundedToUnbounded((solution == null) ? startPoint : solution.getPointRef())));
				final double[] point = adapter.unboundedToBounded(optimum.getPointRef());
				optimum = new PointValuePair(point, optimum.getValue());

				if (solution == null || optimum.getValue() < solution.getValue())
					solution = optimum;
			}
			catch (final Exception e)
			{
				IJ.log("Bounded Powell error: " + e.getMessage());
				if (e instanceof TooManyEvaluationsException)
					maxEval = (int) (maxEval * 1.5);
			}
		}

		IJ.showStatus("");
		IJ.showProgress(1);

		if (solution == null)
		{
			Utils.log("Failed to fit the distribution");
			return;
		}

		final double[] point = solution.getPointRef();
		photons = point[0];
		gain = point[1];
		noise = point[2];
		bias = (int) Math.round(point[3]);
		final String label = String.format("Fitted bias=%d, gain=%s, noise=%s, photons=%s", (int) bias, Utils.rounded(gain),
				Utils.rounded(noise), Utils.rounded(photons));
		Utils.log(label);

		if (simulate)
			Utils.log("Relative Error bias=%s, gain=%s, noise=%s, photons=%s",
					Utils.rounded(relativeError(bias, _bias)), Utils.rounded(relativeError(gain, _gain)),
					Utils.rounded(relativeError(noise, _noise)), Utils.rounded(relativeError(photons, _photons)));

		// Show the PoissonGammaGaussian approximation
		double[] f = null;
		if (showApproximation)
		{
			f = new double[x.length];
			final PoissonGammaGaussianFunction fun = new PoissonGammaGaussianFunction(1.0 / gain, noise);
			final double expected = photons * gain;
			for (int i = 0; i < f.length; i++)
				f[i] = fun.likelihood(x[i] - bias, expected);
				//System.out.printf("x=%d, g=%f, f=%f, error=%f\n", (int) x[i], g[i], f[i],
				//		uk.ac.sussex.gdsc.smlm.fitting.utils.DoubleEquality.relativeError(g[i], f[i]));
			yMax = Maths.maxDefault(yMax, f);
		}

		// Replot
		g = pdf(max, photons, gain, noise, (int) bias);
		plot = new Plot2(TITLE, "ADU", "Frequency");
		plot.setLimits(limits[0], limits[1], 0, yMax * 1.05);
		plot.setColor(Color.black);
		plot.addPoints(x, y, Plot.DOT);
		plot.setColor(Color.red);
		plot.addPoints(x, g, Plot.LINE);

		plot.addLabel(0, 0, label);

		if (showApproximation)
		{
			plot.setColor(Color.blue);
			plot.addPoints(x, f, Plot.LINE);
		}

		Utils.display(TITLE, plot);
	}

	private static double relativeError(double a, double b)
	{
		//return uk.ac.sussex.gdsc.smlm.utils.DoubleEquality.relativeError(a, b);
		final double d = a - b; // Math.abs(a - b);
		return d / b;
	}

	private static MultivariateFunction getFunction(final int[] limits, final double[] y, final int max,
			final int maxEval)
	{
		final MultivariateFunction fun = new MultivariateFunction()
		{
			int eval = 0;

			@Override
			public double value(double[] point)
			{
				IJ.showProgress(++eval, maxEval);
				if (Utils.isInterrupted())
					throw new TooManyEvaluationsException(maxEval);
				// Compute the sum of squares between the two functions
				final double photons = point[0];
				final double gain = point[1];
				final double noise = point[2];
				final int bias = (int) Math.round(point[3]);
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
	 * @param step
	 *            the step between counts to evaluate
	 * @param p
	 *            The average number of photons per pixel input to the EM-camera
	 * @param m
	 *            The multiplication factor (gain)
	 * @return The PDF
	 */
	private static double[] pdfEMGain(final double step, final double p, final double m)
	{
		final StoredDataStatistics stats = new StoredDataStatistics(100);
		stats.add(FastMath.exp(-p));
		for (int c = 1;; c++)
		{
			final double g = pEMGain(c * step, p, m);
			stats.add(g);
			final double delta = g / stats.getSum();
			if (delta < 1e-5)
				break;
		}
		//System.out.printf("Sum = %f\n", stats.getSum() * step);
		return stats.getValues();
	}

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
	private static double pEMGain(double c, double p, double m)
	{
		return PoissonGammaFunction.poissonGamma(c, p, m);
	}

	/**
	 * Calculate the probability density function for EM-gain.
	 * <p>
	 * See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
	 *
	 * @param max
	 *            The maximum count to evaluate
	 * @param step
	 *            the step between counts to evaluate
	 * @param p
	 *            The average number of photons per pixel input to the EM-camera
	 * @param m
	 *            The multiplication factor (gain)
	 * @return The PDF
	 */
	private static double[] pdfEMGain(final int max, final double step, final double p, final double m)
	{
		if (max == 0)
			return pdfEMGain(step, p, m);
		final double[] g = new double[max + 1];
		g[0] = FastMath.exp(-p);
		for (int c = 1;; c++)
		{
			final double count = c * step;
			g[c] = pEMGain(count, p, m);
			if (g[c] == 0 || count >= max)
				break;
		}
		return g;
	}

	private class PDF
	{
		final double[] x, p;

		PDF(double[] x, double[] p)
		{
			this.x = x;
			this.p = p;
		}
	}

	/**
	 * Calculate the probability density function for EM-gain, convolve with a Gaussian and then add a constant offset.
	 * <p>
	 * See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI.
	 *
	 * @param max
	 *            The maximum count to evaluate
	 * @param step
	 *            the step between counts to evaluate
	 * @param p
	 *            The average number of photons per pixel input to the EM-camera
	 * @param m
	 *            The multiplication factor (gain)
	 * @param s
	 *            The read noise (Gaussian standard deviation)
	 * @return The PDF
	 */
	private PDF pdf(final int max, final double step, final double p, final double m, final double s)
	{
		final double[] g = pdfEMGain(max, step, p, m);
		double[] gg;

		int zero = 0;

		if (s > 0)
		{
			// Convolve with Gaussian kernel up to 4 times the standard deviation
			final int radius = (int) Math.ceil(Math.abs(s) * 4 / step) + 1;
			final double[] kernel = new double[2 * radius + 1];
			final double norm = -0.5 / (s * s);
			for (int i = 0, j = radius, jj = radius; j < kernel.length; i++, j++, jj--)
				kernel[j] = kernel[jj] = FastMath.exp(norm * Maths.pow2(i * step));
			// Normalise
			double sum = 0;
			for (int j = 0; j < kernel.length; j++)
				sum += kernel[j];
			for (int j = 0; j < kernel.length; j++)
				kernel[j] /= sum;

			if (extraOptions)
			{
				// Debug
				String title = "Poisson-Gamma";
				Plot plot = new Plot(title, "x", "y", SimpleArrayUtils.newArray(g.length, 0, step), g);
				Utils.display(title, plot);

				title = "Gaussian";
				plot = new Plot(title, "x", "y", SimpleArrayUtils.newArray(kernel.length, radius * -step, step),
						kernel);
				Utils.display(title, plot);
			}

			gg = Convolution.convolveFast(g, kernel);
			// The convolution will have created a larger array so we must adjust the offset for this
			zero = radius;
		}
		else
			gg = g;

		final double[] x = new double[gg.length];
		for (int i = 0, j = -zero; i < x.length; i++, j++)
			x[i] = j * step;

		return new PDF(x, gg);
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
	private static double[] pdf(final int max, final double p, final double m, final double s, int c0)
	{
		final double[] g = pdfEMGain(max, p, m);
		double[] gg;

		if (s > 0)
		{
			// Convolve with Gaussian kernel up to 4 times the standard deviation
			final int radius = (int) Math.ceil(Math.abs(s) * 4) + 1;
			final double[] kernel = new double[2 * radius + 1];
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
		final double[] g0 = new double[gg.length + c0];
		for (int c = 0; c < gg.length; c++)
			if (c + c0 >= 0)
				g0[c + c0] = gg[c];
		return g0;
	}

	private static int[] limits(int[] h)
	{
		int min = 0;
		while (h[min] == 0)
			min++;
		int max = h.length - 1;
		while (h[max] == 0)
			max--;
		return new int[] { min, max };
	}

	private static double[] getX(int[] limits)
	{
		final int min = 0; //limits[0];
		final int range = limits[1] - min + 1;
		final double[] x = new double[range];
		for (int i = 0; i < range; i++)
			x[i] = min + i;
		return x;
	}

	private static double[] getY(int[] h, int[] limits)
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

	private static int getSize(int[] h)
	{
		int size = 0;
		for (int i = 0; i < h.length; i++)
			size += h[i];
		return size;
	}

	private static double getMean(int[] h)
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
		final GenericDialog gd = new GenericDialog(TITLE);
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

	@SuppressWarnings("unused")
	private void plotPMF()
	{
		if (!showPMFDialog())
			return;

		final double step = getStepSize(_photons, _gain, _noise);

		final PDF pdf = pdf(0, step, _photons, _gain, _noise);
		double[] pmf = pdf.p;
		double[] x = pdf.x;
		double yMax = Maths.max(pmf);

		// Get the approximation
		LikelihoodFunction fun;
		double myNoise = _noise;
		switch (approximation)
		{
			case 3:
				//fun = new InterpolatedPoissonFunction(1.0 / _gain, true);
				fun = new PoissonFunction(1.0 / _gain);
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
				final PoissonGammaGaussianFunction myFun = new PoissonGammaGaussianFunction(1.0 / _gain, myNoise);
				myFun.setMinimumProbability(0);
				fun = myFun;
		}
		double expected = _photons;
		if (offset != 0)
			expected += offset * expected / 100.0;
		//expected *= _gain;

		// Normalise
		final boolean normalise = false;
		if (normalise)
		{
			final double sum = Maths.sum(pmf);
			for (int i = pmf.length; i-- > 0;)
				pmf[i] /= sum;
		}

		// Get CDF
		double sum = 0;
		double sum2 = 0;
		double prev = 0;
		double prev2 = 0;
		double[] f = new double[x.length];
		double[] cdf1 = new double[pmf.length];
		double[] cdf2 = new double[pmf.length];
		final double step_2 = step / 2;
		for (int i = 0; i < cdf1.length; i++)
		{
			// Trapezoid integration
			//sum += (pmf[i] + prev) * step_2;
			sum += pmf[i] * step;
			cdf1[i] = sum;
			prev = pmf[i];
			f[i] = fun.likelihood(x[i], expected);
			//sum2 += (f[i] + prev2) * step_2;
			sum2 += f[i] * step;
			cdf2[i] = sum2;
			prev2 = f[i];
		}

		// Truncate x for plotting
		int max = 0;
		sum = prev = 0;
		double p = 1 - tail;
		while (sum < p && max < pmf.length)
		{
			//sum += (pmf[max] + prev) * step_2;
			sum += pmf[max] * step;
			prev = pmf[max];
			if (sum > 0.5 && pmf[max] == 0)
				break;
			max++;
		}

		int min = pmf.length;
		sum = prev = 0;
		p = 1 - head;
		while (sum < p && min > 0)
		{
			min--;
			//sum += (pmf[min] + prev) * step_2;
			sum += pmf[min] * step;
			prev = pmf[min];
			if (sum > 0.5 && pmf[min] == 0)
				break;
		}

		//int min = (int) (dummyBias - gaussWidth * _noise);
		pmf = Arrays.copyOfRange(pmf, min, max);
		x = Arrays.copyOfRange(x, min, max);
		f = Arrays.copyOfRange(f, min, max);

		if (showApproximation)
			yMax = Maths.maxDefault(yMax, f);

		final String label = String.format("Gain=%s, noise=%s, photons=%s", Utils.rounded(_gain), Utils.rounded(_noise),
				Utils.rounded(_photons));

		final Plot2 plot = new Plot2("PMF", "ADUs", "p");
		plot.setLimits(x[0], x[x.length - 1], 0, yMax);
		plot.setColor(Color.red);
		plot.addPoints(x, pmf, Plot.LINE);
		if (showApproximation)
		{
			plot.setColor(Color.blue);
			plot.addPoints(x, f, Plot.LINE);
		}

		plot.setColor(Color.magenta);
		plot.drawLine(_photons * _gain, 0, _photons * _gain, yMax);
		plot.setColor(Color.black);
		plot.addLabel(0, 0, label);
		final PlotWindow win1 = Utils.display("PMF", plot);

		// Plot the difference between the actual and approximation
		final double[] delta = new double[f.length];
		for (int i = 0; i < f.length; i++)
		{
			if (pmf[i] == 0 && f[i] == 0)
				continue;
			if (relativeDelta)
				delta[i] = DoubleEquality.relativeError(f[i], pmf[i]) * Math.signum(f[i] - pmf[i]);
			else
				delta[i] = f[i] - pmf[i];
		}

		final Plot2 plot2 = new Plot2("PMF delta", "ADUs", (relativeDelta) ? "Relative delta" : "delta");
		final double[] limits = Maths.limits(delta);
		plot2.setLimits(x[0], x[x.length - 1], limits[0], limits[1]);
		plot2.setColor(Color.red);
		plot2.addPoints(x, delta, Plot.LINE);
		plot2.setColor(Color.magenta);
		plot2.drawLine(_photons * _gain, limits[0], _photons * _gain, limits[1]);
		plot2.setColor(Color.black);
		plot2.addLabel(0, 0, label + ((offset == 0) ? "" : ", expected = " + Utils.rounded(expected / _gain)));
		final PlotWindow win2 = Utils.display("PMF delta", plot2);

		if (Utils.isNewWindow())
		{
			final Point p2 = win1.getLocation();
			p2.y += win1.getHeight();
			win2.setLocation(p2);
		}

		// Plot the CDF of each distribution.
		// Compute the Kolmogorov distance as the supremum (maximum)
		// difference between the two cumulative probability distributions.
		// https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test
		double kolmogorovDistance = 0;
		double xd = x[0];
		for (int i = 0; i < cdf1.length; i++)
		{
			final double dist = Math.abs(cdf1[i] - cdf2[i]);
			if (kolmogorovDistance < dist)
			{
				kolmogorovDistance = dist;
				xd = pdf.x[i];
			}
		}
		cdf1 = Arrays.copyOfRange(cdf1, min, max);
		cdf2 = Arrays.copyOfRange(cdf2, min, max);

		final Plot2 plot3 = new Plot2("CDF", "ADUs", "p");
		yMax = 1.05;
		plot3.setLimits(x[0], x[x.length - 1], 0, yMax);
		plot3.setColor(Color.red);
		plot3.addPoints(x, cdf1, Plot.LINE);
		plot3.setColor(Color.blue);
		plot3.addPoints(x, cdf2, Plot.LINE);

		plot3.setColor(Color.magenta);
		plot3.drawLine(_photons * _gain, 0, _photons * _gain, yMax);
		plot3.drawDottedLine(xd, 0, xd, yMax, 2);
		plot3.setColor(Color.black);
		plot3.addLabel(0, 0, label + ", Kolmogorov distance = " + Utils.rounded(kolmogorovDistance) + " @ " + xd);
		plot3.addLegend("CDF\nApprox");
		final PlotWindow win3 = Utils.display("CDF", plot3);

		if (Utils.isNewWindow())
		{
			final Point p2 = win1.getLocation();
			p2.x += win1.getWidth();
			win3.setLocation(p2);
		}
	}

	/**
	 * Gets the step size.
	 * <p>
	 * Note: Currently this just returns 1 as it should be a PMF so only uses discrete values.
	 *
	 * @param photons
	 *            the photons
	 * @param gain
	 *            the gain
	 * @param noise
	 *            the noise
	 * @return the step size
	 */
	private static double getStepSize(double photons, double gain, double noise)
	{
		// Note: This is not valid as the PMF should only accept integer input.
		return 1;

		//		// Determine the best step to plot the PMF.
		//		// Ensure there are enough points on the chart.
		//		double step = 1.0;
		//
		//		int poissonGammaWidth = pdfEMGain(step, photons, gain).length;
		//
		//		while (step > 0.001)
		//		{
		//			int gaussianWidth = (int) Math.ceil(Math.abs(noise) * 4 / step);
		//			if (gaussianWidth * 2 + poissonGammaWidth / step < 300)
		//				step /= 2;
		//			else
		//				break;
		//		}
		//		return step;
	}

	private boolean showPMFDialog()
	{
		final GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Plot the probability mass function for EM-gain");

		gd.addNumericField("Gain", _gain, 2, 6, "Count/electron");
		gd.addNumericField("Noise", _noise, 2, 6, "Count");
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
