package gdsc.smlm.ij.plugins;

import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.utils.Bessel;
import gdsc.smlm.utils.Convolution;
import gdsc.smlm.utils.Maths;
import gdsc.smlm.utils.StoredDataStatistics;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import java.awt.Color;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateFunctionMappingAdapter;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.PowellOptimizer;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well44497b;
import org.apache.commons.math3.util.FastMath;

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
	private static boolean _simulate = false;
	private boolean simulate = false, extraOptions = false;
	private static double _photons = 1, _bias = 500, _gain = 40, _noise = 3;

	private ImagePlus imp;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	public int setup(String arg, ImagePlus imp)
	{
		extraOptions = Utils.isExtraOptions();
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
		final int[] h = (simulate) ? simulateHistogram() : buildHistogram(imp);

		// We need > 10^7 pixels from flat white-light images under constant exposure ...
		final int size = getSize(h);
		Utils.log("Histogram contains %d pixels", size);
		if (size < MINIMUM_PIXELS)
			Utils.log("WARNING : Recommend at least %d pixels (%sx more)", MINIMUM_PIXELS,
					Utils.rounded((double) MINIMUM_PIXELS / size));

		fit(h);
	}

	private int[] simulateHistogram()
	{
		final double[] g = pdf(0, _photons, _gain, _noise, (int) _bias);

		IJ.showStatus("Simulating histogram ...");

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
		final int steps = 20000; // MINIMUM_PIXELS;
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
		IJ.showStatus("");
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
		return data;
	}

	private static int[] getHistogram(ImageProcessor ip, Roi roi)
	{
		ip.setRoi(roi);
		return ip.getHistogram();
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

		Plot plot = new Plot(TITLE, "ADU", "Frequency");
		final double yMax = Maths.max(y);
		plot.setLimits(limits[0], limits[1], 0, yMax);
		plot.setColor(Color.black);
		plot.addPoints(x, y, Plot.DOT);
		Utils.display(TITLE, plot);

		// Estimate remaining parameters. 
		// Assuming a gamma_distribution(shape,scale) then mean = shape * scale
		// scale = gain
		// shape = Photons = mean / gain
		final double mean = getMean(h) - bias;
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
		PowellOptimizer o = new PowellOptimizer(1e-6, 1e-16);
		double[] startPoint = new double[] { photons, gain, noise, bias };
		final int maxEval = 3000;

		// Set bounds
		double[] lower = new double[] { 0, 0.5 * gain, 0, bias - gain };
		double[] upper = new double[] { (limits[1] - limits[0]) / gain, 2 * gain, gain, bias + gain };

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
				PointValuePair optimum = o.optimize(new MaxEval(maxEval), new ObjectiveFunction(fun),
						GoalType.MINIMIZE, new InitialGuess((solution == null) ? startPoint : solution.getPointRef()));
				if (solution == null || optimum.getValue() < solution.getValue())
				{
					solution = optimum;
				}
			}
			catch (Exception e)
			{
			}
			try
			{
				// Bounded Powell optimiser
				MultivariateFunction fun = getFunction(limits, y, max, maxEval);
				MultivariateFunctionMappingAdapter adapter = new MultivariateFunctionMappingAdapter(fun, lower, upper);
				PointValuePair optimum = o.optimize(
						new MaxEval(maxEval),
						new ObjectiveFunction(adapter),
						GoalType.MINIMIZE,
						new InitialGuess(adapter.boundedToUnbounded((solution == null) ? startPoint : solution
								.getPointRef())));
				double[] point = adapter.unboundedToBounded(optimum.getPointRef());
				optimum = new PointValuePair(point, optimum.getValue());

				if (solution == null || optimum.getValue() < solution.getValue())
				{
					solution = optimum;
				}
			}
			catch (Exception e)
			{
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

		// Replot
		g = pdf(max, photons, gain, noise, (int) bias);
		plot = new Plot(TITLE, "ADU", "Frequency");
		plot.setLimits(limits[0], limits[1], 0, yMax);
		plot.setColor(Color.black);
		plot.addPoints(x, y, Plot.DOT);
		plot.setColor(Color.red);
		plot.addPoints(x, g, Plot.LINE);

		plot.addLabel(0, 0, label);

		Utils.display(TITLE, plot);
	}

	private MultivariateFunction getFunction(final int[] limits, final double[] y, final int max, final int maxEval)
	{
		MultivariateFunction fun = new MultivariateFunction()
		{
			int eval = 0;

			@Override
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
		double g = FastMath.exp(-p);
		stats.add(g);
		for (int c = 1;; c++)
		{
			g = Math.sqrt(p / (c * m)) * FastMath.exp(-c / m - p) * Bessel.I1(2 * Math.sqrt(c * p / m));
			stats.add(g);
			final double delta = g / stats.getSum();
			if (delta < 1e-4)
				break;
		}
		return stats.getValues();
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
		double[] g = new double[max];
		g[0] = FastMath.exp(-p);
		for (int c = 1; c < max; c++)
			g[c] = Math.sqrt(p / (c * m)) * FastMath.exp(-c / m - p) * Bessel.I1(2 * Math.sqrt(c * p / m));
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

		double[] gg;

		if (s > 0)
		{
			// Use Fourier Transform when the convolution is large
			if (g.length * kernel.length > 100000)
				gg = Convolution.convolveFFT(g, kernel);
			else
				gg = Convolution.convolve(g, kernel);
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

		gd.addMessage("Analyse the white-light histogram of an image stack to determine EM-gain parameters.\n \n"
				+ "See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321 (Supplementary Information).");

		if (extraOptions)
		{
			gd.addCheckbox("Simulate", _simulate);
			gd.addNumericField("sBias", _bias, 0);
			gd.addNumericField("sGain", _gain, 2);
			gd.addNumericField("sNoise", _noise, 2);
			gd.addNumericField("sPhotons", _photons, 2);
		}

		gd.addNumericField("Bias (estimate)", bias, 0);
		gd.addNumericField("Gain (estimate)", gain, 2);
		gd.addNumericField("Noise (estimate)", noise, 2);
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
			if (gd.invalidNumber() || _bias < 0 || _gain < 1 || _photons == 0)
				return DONE;
		}

		bias = gd.getNextNumber();
		gain = gd.getNextNumber();
		noise = FastMath.abs(gd.getNextNumber());

		if (gd.invalidNumber() || bias < 0 || gain < 1)
			return DONE;

		return (simulate) ? NO_IMAGE_REQUIRED : FLAGS;
	}
}
