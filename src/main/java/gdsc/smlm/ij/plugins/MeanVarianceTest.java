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

import gdsc.smlm.ij.utils.SeriesOpener;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.utils.Maths;
import gdsc.smlm.utils.Statistics;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot2;
import ij.plugin.PlugIn;
import ij.text.TextWindow;

import java.awt.Color;
import java.awt.Frame;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction.Parametric;
import org.apache.commons.math3.fitting.CurveFitter;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.LevenbergMarquardtOptimizer;

/**
 * Opens a folder of images and computes a Mean-Variance Test.
 * <p>
 * Each image must contain a 2-slice stack of white light images. The image filename must contain the exposure time
 * separated by whitespace, e.g. 'MVT 30.tif' for 30 milliseconds.
 * <p>
 * Gain calculations for standard CCD and EM-CCD cameras are based on the paper: Hirsch M, Wareham RJ, Martin-Fernandez
 * ML, Hobson MP, Rolfe DJ (2013) A Stochastic Model for Electron Multiplication Charge-Coupled Devices – From Theory to
 * Practice. PLoS ONE 8(1): e53671. doi:10.1371/journal.pone.0053671
 */
public class MeanVarianceTest implements PlugIn
{
	private static final String TITLE = "Mean Variance Test";
	private static double cameraGain = 0;
	private int exposureCounter = 0;

	private class PairSample
	{
		int slice1, slice2;
		double mean1, mean2, variance;

		public PairSample(int slice1, int slice2, double mean1, double mean2, double variance)
		{
			this.slice1 = slice1;
			this.slice2 = slice2;
			this.mean1 = mean1;
			this.mean2 = mean2;
			this.variance = variance;
		}

		public double getMean()
		{
			return (mean1 + mean2) * 0.5;
		}
	}

	private class ImageSample
	{
		String title;
		float[][] slices;
		double exposure;
		double[] means;
		List<PairSample> samples;

		public ImageSample(ImagePlus imp)
		{
			// Check stack has two slices
			if (imp.getStackSize() < 2)
				throw new IllegalArgumentException("Image must have at least 2-slices: " + imp.getTitle());

			// Count all the valid input images
			exposureCounter++;

			// Extract the exposure time
			exposure = -1;
			String[] tokens = imp.getTitle().split("[ .]");
			for (String token : tokens)
			{
				try
				{
					exposure = Double.parseDouble(token);
					if (exposure >= 0)
						break;
				}
				catch (NumberFormatException e)
				{
					// Ignore
				}
			}

			if (exposure < 0)
			{
				//throw new IllegalArgumentException("Image must have exposure time in the filename: " + imp.getTitle());

				// If no exposure was found: assume exposure 0 for first input image otherwise set an arbitrary exposure
				exposure = (exposureCounter == 1) ? 0 : 9999;
			}

			title = imp.getTitle();

			// Get all the pixels into a float stack. 
			// Look for saturated pixels that will invalidate the test.
			final int size = imp.getStackSize();
			slices = new float[size][];
			final float saturated = getSaturation(imp);
			ImageStack stack = imp.getImageStack();
			for (int slice = 1; slice <= size; slice++)
			{
				final float[] thisSlice = slices[slice - 1] = (float[]) stack.getProcessor(slice).toFloat(0, null)
						.getPixels();
				checkSaturation(slice, thisSlice, saturated);
			}
		}

		private float getSaturation(ImagePlus imp)
		{
			switch (imp.getBitDepth())
			{
				case 8:
				case 24:
					return 255f;
				case 16:
					return 65535f;
				case 32:
					// float images cannot be saturated
					return Float.NaN;
			}
			throw new IllegalArgumentException("Cannot determine saturation level for image: " + imp.getTitle());
		}

		private void checkSaturation(int i, float[] data, float saturated)
		{
			if (saturated == Float.NaN)
				return;
			for (float f : data)
				if (f >= saturated)
					throw new IllegalArgumentException("Image " + title + " has saturated pixels in slice: " + i);
		}

		public void compute()
		{
			final int size = slices.length;
			samples = new ArrayList<PairSample>(((size - 1) * size) / 2);

			// Cache data
			means = new double[size];
			for (int slice1 = 0; slice1 < size; slice1++)
			{
				means[slice1] = new Statistics(slices[slice1]).getMean();
			}

			// Compute mean and variance.
			// See http://www.photometrics.com/resources/whitepapers/mean-variance.php
			for (int slice1 = 0; slice1 < size; slice1++)
			{
				float[] data1 = slices[slice1];
				for (int slice2 = slice1 + 1; slice2 < size; slice2++)
				{
					float[] data2 = slices[slice2];
					Statistics s = new Statistics();
					for (int i = 0; i < data1.length; i++)
						s.add(data1[i] - data2[i]);
					double variance = s.getVariance() / 2.0;
					samples.add(new PairSample(slice1 + 1, slice2 + 1, means[slice1], means[slice2], variance));
				}
				slices[slice1] = null; // Allow garbage collection
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		String inputDirectory = IJ.getDirectory("Select image series ...");
		if (inputDirectory == null)
			return;

		SeriesOpener series = new SeriesOpener(inputDirectory, false);
		series.setVariableSize(true);
		if (series.getNumberOfImages() < 3)
		{
			IJ.error(TITLE, "Not enough images in the selected directory");
			return;
		}
		if (!IJ.showMessageWithCancel(TITLE, String.format("Analyse %d images, first image:\n%s",
				series.getNumberOfImages(), series.getImageList()[0])))
		{
			return;
		}

		IJ.showStatus("Loading images");
		List<ImageSample> images = getImages(series);

		if (images.size() < 3)
		{
			IJ.error(TITLE, "Not enough images for analysis");
			return;
		}
		if (images.get(0).exposure != 0)
		{
			IJ.error(TITLE, "First image in series must have exposure 0 (Bias image)");
			return;
		}

		boolean emMode = (arg != null && arg.contains("em"));
		if (emMode)
		{
			// Ask the user for the camera gain ...
			GenericDialog gd = new GenericDialog(TITLE);
			gd.addMessage("Estimating the EM-gain requires the camera gain without EM readout enabled");
			gd.addNumericField("Camera_gain (ADU/e-)", cameraGain, 4);
			gd.showDialog();
			if (gd.wasCanceled())
				return;
			cameraGain = gd.getNextNumber();
		}

		//IJ.showProgress(currentIndex, finalIndex);
		IJ.showStatus("Computing mean & variance");
		for (int i = 0; i < images.size(); i++)
		{
			IJ.showStatus(String.format("Computing mean & variance %d/%d", i + 1, images.size()));
			IJ.showProgress(i, images.size());
			images.get(i).compute();
		}
		IJ.showProgress(1);

		IJ.showStatus("Computing results");
		TextWindow results = createResultsWindow();

		// Allow user to input multiple bias images
		int start = 0;
		Statistics biasStats = new Statistics();
		Statistics noiseStats = new Statistics();
		while (start < images.size())
		{
			ImageSample sample = images.get(start);
			if (sample.exposure == 0)
			{
				biasStats.add(sample.means);
				for (PairSample pair : sample.samples)
				{
					noiseStats.add(pair.variance);
				}
				start++;
			}
			else
				break;
		}
		final double bias = biasStats.getMean();

		// Get the mean-variance data
		int total = 0;
		for (int i = start; i < images.size(); i++)
			total += images.get(i).samples.size();
		double[] mean = new double[total];
		double[] variance = new double[mean.length];
		Statistics gainStats = new Statistics();
		final CurveFitter<Parametric> fitter = new CurveFitter<Parametric>(new LevenbergMarquardtOptimizer());
		for (int i = start, j = 0; i < images.size(); i++)
		{
			ImageSample sample = images.get(i);
			for (PairSample pair : sample.samples)
			{
				mean[j] = pair.getMean();
				variance[j] = pair.variance;
				// Gain is in ADU / e
				double gain = variance[j] / (mean[j] - bias);
				gainStats.add(gain);
				fitter.addObservedPoint(mean[j], variance[j]);

				if (emMode)
				{
					gain /= (2 * cameraGain);
				}

				StringBuilder sb = new StringBuilder();
				sb.append(sample.title).append("\t");
				sb.append(sample.exposure).append("\t");
				sb.append(pair.slice1).append("\t");
				sb.append(pair.slice2).append("\t");
				sb.append(IJ.d2s(pair.mean1, 2)).append("\t");
				sb.append(IJ.d2s(pair.mean2, 2)).append("\t");
				sb.append(IJ.d2s(mean[j], 2)).append("\t");
				sb.append(IJ.d2s(variance[j], 2)).append("\t");
				sb.append(Utils.rounded(gain, 4));
				results.append(sb.toString());
				j++;
			}
		}

		IJ.showStatus("Computing fit");

		// Sort
		int[] indices = rank(mean);
		mean = reorder(mean, indices);
		variance = reorder(variance, indices);

		// Compute optimal coefficients.
		final double[] init = { 0, 1 / gainStats.getMean() }; // a - b x
		final double[] best = fitter.fit(new PolynomialFunction.Parametric(), init);

		// Construct the polynomial that best fits the data.
		final PolynomialFunction fitted = new PolynomialFunction(best);

		// SuperPlot mean verses variance. Gradient is gain in ADU/e.
		String title = TITLE + " results";
		Plot2 plot = new Plot2(title, "Mean", "Variance");
		double[] xlimits = Maths.limits(mean);
		double[] ylimits = Maths.limits(variance);
		double xrange = (xlimits[1] - xlimits[0]) * 0.05;
		if (xrange == 0)
			xrange = 0.05;
		double yrange = (ylimits[1] - ylimits[0]) * 0.05;
		if (yrange == 0)
			yrange = 0.05;
		plot.setLimits(xlimits[0] - xrange, xlimits[1] + xrange, ylimits[0] - yrange, ylimits[1] + yrange);
		plot.setColor(Color.blue);
		plot.addPoints(mean, variance, Plot2.CROSS);
		plot.setColor(Color.red);
		plot.addPoints(new double[] { mean[0], mean[mean.length - 1] },
				new double[] { fitted.value(mean[0]), fitted.value(mean[mean.length - 1]) }, Plot2.LINE);
		Utils.display(title, plot);

		final double avBiasNoise = Math.sqrt(noiseStats.getMean());

		Utils.log(TITLE);
		Utils.log("  Directory = %s", inputDirectory);
		Utils.log("  Bias = %s +/- %s (ADU)", Utils.rounded(bias, 4), Utils.rounded(avBiasNoise, 4));
		Utils.log("  Variance = %s + %s * mean", Utils.rounded(best[0], 4), Utils.rounded(best[1], 4));
		if (emMode)
		{
			final double emGain = best[1] / (2 * cameraGain);

			// Noise is standard deviation of the bias image divided by the total gain (in ADU/e-)
			final double totalGain = emGain * cameraGain;
			Utils.log("  Read Noise = %s (e-) [%s (ADU)]", Utils.rounded(avBiasNoise / totalGain, 4),
					Utils.rounded(avBiasNoise, 4));

			Utils.log("  Gain = 1 / %s (ADU/e-)", Utils.rounded(1 / cameraGain, 4));
			Utils.log("  EM-Gain = %s", Utils.rounded(emGain, 4));
			Utils.log("  Total Gain = %s (ADU/e-)", Utils.rounded(totalGain, 4));
		}
		else
		{
			// Noise is standard deviation of the bias image divided by the gain (in ADU/e-)
			cameraGain = best[1];
			final double readNoise = avBiasNoise / cameraGain;
			Utils.log("  Read Noise = %s (e-) [%s (ADU)]", Utils.rounded(readNoise, 4),
					Utils.rounded(readNoise * cameraGain, 4));

			Utils.log("  Gain = 1 / %s (ADU/e-)", Utils.rounded(1 / cameraGain, 4));
		}
		IJ.showStatus("");
	}

	private TextWindow createResultsWindow()
	{
		Frame f = WindowManager.getFrame(TITLE);
		if (f instanceof TextWindow)
		{
			return (TextWindow) f;
		}
		return new TextWindow(TITLE, "Image\tExposure\tSlice1\tSlice2\tMean1\tMean2\tMean\tVariance\tGain", "", 800,
				500);
	}

	private List<ImageSample> getImages(SeriesOpener series)
	{
		List<ImageSample> images = new ArrayList<ImageSample>(series.getNumberOfImages());
		ImagePlus imp = series.nextImage();
		while (imp != null)
		{
			try
			{
				images.add(new ImageSample(imp));
			}
			catch (IllegalArgumentException e)
			{
				Utils.log(e.getMessage());
			}
			imp.close();
			imp = series.nextImage();
		}
		// Sort to ensure all 0 exposure images are first, the remaining order is arbitrary
		Collections.sort(images, new Comparator<ImageSample>()
		{
			public int compare(ImageSample o1, ImageSample o2)
			{
				if (o1.exposure < o2.exposure)
					return -1;
				if (o1.exposure > o2.exposure)
					return 1;
				return 0;
			}
		});
		return images;
	}

	/**
	 * Returns a sorted list of indices of the specified double array.
	 * Modified from: http://stackoverflow.com/questions/951848 by N.Vischer.
	 * Copied from ImageJ 1.48 for backwards compatibility
	 */
	public static int[] rank(double[] values)
	{
		int n = values.length;
		final Integer[] indexes = new Integer[n];
		final Double[] data = new Double[n];
		for (int i = 0; i < n; i++)
		{
			indexes[i] = new Integer(i);
			data[i] = new Double(values[i]);
		}
		Arrays.sort(indexes, new Comparator<Integer>()
		{
			public int compare(final Integer o1, final Integer o2)
			{
				return data[o1].compareTo(data[o2]);
			}
		});
		int[] indexes2 = new int[n];
		for (int i = 0; i < n; i++)
			indexes2[i] = indexes[i].intValue();
		return indexes2;
	}

	private double[] reorder(double[] data, int[] indices)
	{
		double[] array = new double[indices.length];
		for (int i = 0; i < array.length; i++)
		{
			array[i] = data[indices[i]];
		}
		return array;
	}
}
