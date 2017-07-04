package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

import gdsc.core.ij.BufferedTextWindow;
import gdsc.core.ij.Utils;
import gdsc.core.threshold.AutoThreshold;
import gdsc.core.utils.Maths;
import gdsc.core.utils.NoiseEstimator;
import gdsc.core.utils.Statistics;
import gdsc.smlm.data.config.FitConfig.NoiseEstimatorMethod;
import gdsc.smlm.data.config.FitConfigHelper;
import gdsc.smlm.engine.DataEstimator;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.ImageConverter;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingExtendedGenericDialog;
import ij.gui.Plot2;
import ij.gui.PlotWindow;
import ij.plugin.WindowOrganiser;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.ImageProcessor;
import ij.text.TextWindow;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Contains methods to find the noise in the provided image data.
 */
public class BackgroundEstimator implements ExtendedPlugInFilter, DialogListener
{
	private static final String TITLE = "Background Estimator";
	private List<double[]> results;
	private final int FLAGS = DOES_8G | DOES_16 | DOES_32 | PARALLELIZE_STACKS | FINAL_PROCESSING | NO_CHANGES;
	private PlugInFilterRunner pfr;
	private ImagePlus imp;
	private NonBlockingExtendedGenericDialog gd;
	private static double percentile;
	private static NoiseEstimatorMethod noiseMethod = NoiseEstimatorMethod.QUICK_RESIDUALS_LEAST_TRIMMED_OF_SQUARES;
	private NoiseEstimator.Method myNoiseMethod = null;
	private static AutoThreshold.Method thresholdMethod = AutoThreshold.Method.DEFAULT;
	private static float fraction;
	private static int histogramSize;
	static
	{
		DataEstimator de = new DataEstimator(new float[0], 0, 0);
		fraction = de.getFraction();
		histogramSize = de.getHistogramSize();
	}

	public int setup(String arg, ImagePlus imp)
	{
		if (arg.equalsIgnoreCase("final"))
		{
			showResults();
			return DONE;
		}
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (imp == null)
		{
			IJ.noImage();
			return DONE;
		}
		this.imp = imp;
		results = Collections.synchronizedList(new ArrayList<double[]>(imp.getStackSize()));
		return FLAGS;
	}

	public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr)
	{
		// If using a stack, provide a preview graph of the noise for two methods
		if (imp.getStackSize() > 1)
		{
			this.pfr = pfr;

			drawPlot();

			gd = new NonBlockingExtendedGenericDialog(TITLE);
			gd.addHelp(About.HELP_URL);

			gd.addSlider("Percential", 0, 100, percentile);

			gd.addChoice("Noise_method", SettingsManager.getNoiseEstimatorMethodNames(), noiseMethod.ordinal());

			// For background based on pixel below a threshold
			String[] thresholdMethods = AutoThreshold.getMethods(true);
			gd.addChoice("Threshold_method", thresholdMethods, thresholdMethods[thresholdMethod.ordinal() - 1]);
			gd.addSlider("Fraction", 0, 0.999, fraction);
			gd.addNumericField("Histogram_size", histogramSize, 0);

			gd.addDialogListener(this);
			gd.addMessage("Click OK to compute table for all slices");
			gd.showDialog();

			if (gd.wasCanceled() || !dialogItemChanged(gd, null))
				return DONE;
		}

		return IJ.setupDialog(imp, FLAGS);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.gui.DialogListener#dialogItemChanged(ij.gui.GenericDialog, java.awt.AWTEvent)
	 */
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
	{
		percentile = gd.getNextNumber();
		noiseMethod = SettingsManager.getNoiseEstimatorMethodValues()[gd.getNextChoiceIndex()];
		myNoiseMethod = FitConfigHelper.convertNoiseEstimatorMethod(noiseMethod);
		thresholdMethod = AutoThreshold.getMethod(gd.getNextChoiceIndex(), true);
		fraction = (float) gd.getNextNumber();
		histogramSize = (int) gd.getNextNumber();
		if (gd.isShowing())
			drawPlot();
		return true;
	}

	/**
	 * Build a plot of the noise estimate from the current frame.
	 * Limit the preview to 100 frames.
	 */
	private void drawPlot()
	{
		IJ.showStatus("Estimating background ...");

		int start = imp.getCurrentSlice();
		int end = FastMath.min(imp.getStackSize(), start + 100);
		int size = end - start + 1;
		double[] xValues = new double[size];
		double[] noise1 = new double[size];
		double[] noise2 = new double[size];
		double[] background = new double[size];
		double[] threshold = new double[size];
		double[] percentile = new double[size];

		ImageStack stack = imp.getImageStack();
		Rectangle bounds = imp.getProcessor().getRoi();
		float[] buffer = null;
		for (int slice = start, i = 0; slice <= end; slice++, i++)
		{
			IJ.showProgress(i, size);
			final ImageProcessor ip = stack.getProcessor(slice);
			buffer = ImageConverter.getData(ip.getPixels(), ip.getWidth(), ip.getHeight(), bounds, buffer);
			final DataEstimator de = new DataEstimator(buffer, bounds.width, bounds.height);
			de.setFraction(fraction);
			de.setHistogramSize(histogramSize);
			de.setThresholdMethod(thresholdMethod);
			xValues[i] = slice;
			try
			{
				noise1[i] = de.getNoise();
				noise2[i] = de.getNoise(myNoiseMethod);
				background[i] = de.getBackground();
				threshold[i] = de.getThreshold();
				percentile[i] = de.getPercentile(BackgroundEstimator.percentile);
			}
			catch (Exception e)
			{
				e.printStackTrace();
				throw new RuntimeException(e);
			}
		}
		IJ.showProgress(1);

		IJ.showStatus("Plotting background ...");

		WindowOrganiser wo = new WindowOrganiser();
		plot(wo, xValues, noise1, noise2, null, "Noise", "Background Noise", "Global Noise", null);
		plot(wo, xValues, background, threshold, percentile, "Background", "Background", "Threshold", "Percentile");
		wo.tile();

		IJ.showStatus("");
	}

	private void plot(WindowOrganiser wo, double[] xValues, double[] data1, double[] data2, double[] data3,
			String title, String title1, String title2, String title3)
	{
		// Get limits
		double[] a = Maths.limits(xValues);
		double[] b = Maths.limits(data1);
		b = Maths.limits(b, data2);
		if (data3 != null)
			b = Maths.limits(b, data3);

		title = imp.getTitle() + " " + title;
		Plot2 plot = new Plot2(title, "Slice", title);
		double range = b[1] - b[0];
		if (range == 0)
			range = 1;
		plot.setLimits(a[0], a[1], b[0] - 0.05 * range, b[1] + 0.05 * range);

		plot.setColor(Color.blue);
		plot.addPoints(xValues, data1, Plot2.LINE);
		plot.draw();
		Statistics stats = new Statistics(data1);
		String label = String.format("%s (Blue) = %s +/- %s", title1, Utils.rounded(stats.getMean()),
				Utils.rounded(stats.getStandardDeviation()));

		plot.setColor(Color.red);
		plot.addPoints(xValues, data2, Plot2.LINE);
		stats = new Statistics(data2);
		label += String.format(", %s (Red) = %s +/- %s", title2, Utils.rounded(stats.getMean()),
				Utils.rounded(stats.getStandardDeviation()));

		if (data3 != null)
		{
			plot.setColor(Color.green);
			plot.addPoints(xValues, data3, Plot2.LINE);
			stats = new Statistics(data3);
			label += String.format(", %s (Green) = %s +/- %s", title3, Utils.rounded(stats.getMean()),
					Utils.rounded(stats.getStandardDeviation()));
		}

		plot.setColor(Color.black);
		plot.addLabel(0, 0, label);

		PlotWindow pw = Utils.display(title, plot);
		if (Utils.isNewWindow())
			wo.add(pw);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	public void run(ImageProcessor ip)
	{
		// Perform all methods and add to the results
		double[] result = new double[8];
		int i = 0;
		result[i++] = (pfr == null) ? 1 : pfr.getSliceNumber();
		Rectangle bounds = ip.getRoi();
		float[] buffer = ImageConverter.getData(ip.getPixels(), ip.getWidth(), ip.getHeight(), bounds, null);
		final DataEstimator de = new DataEstimator(buffer, bounds.width, bounds.height);
		de.setFraction(fraction);
		de.setHistogramSize(histogramSize);
		de.setThresholdMethod(thresholdMethod);
		result[i++] = de.isBackgroundRegion() ? 1 : 0;
		result[i++] = de.getNoise();
		result[i++] = de.getNoise(myNoiseMethod);
		result[i++] = de.getBackground();
		result[i++] = de.getThreshold();
		result[i++] = de.getBackgroundSize();
		result[i++] = de.getPercentile(percentile);
		results.add(result);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.ExtendedPlugInFilter#setNPasses(int)
	 */
	public void setNPasses(int nPasses)
	{
		// Do nothing
	}

	private void showResults()
	{
		Collections.sort(results, new Comparator<double[]>()
		{
			public int compare(double[] o1, double[] o2)
			{
				// Sort on slice number
				return (o1[0] < o2[0]) ? -1 : 1;
			}
		});

		BufferedTextWindow tw = new BufferedTextWindow(
				new TextWindow(imp.getTitle() + " Background", createHeader(), "", 800, 400));
		for (double[] result : results)
			tw.append(createResult(result));
		tw.flush();
	}

	private String createHeader()
	{
		StringBuilder sb = new StringBuilder(imp.getTitle());
		sb.append('\t').append(Utils.rounded(percentile));
		sb.append('\t').append(noiseMethod.toString());
		sb.append('\t').append(thresholdMethod.toString());
		sb.append('\t').append(Utils.rounded(fraction));
		sb.append('\t').append(histogramSize).append('\t');
		prefix = sb.toString();

		sb = new StringBuilder("Image");
		sb.append("\tPercentile");
		sb.append("\tNoiseMethod");
		sb.append("\tThresholdMethod");
		sb.append("\tFraction");
		sb.append("\tHistogramSize");
		sb.append("\tSlice");
		sb.append("\tIsBackground");
		sb.append("\tNoise");
		sb.append("\tGlobalNoise");
		sb.append("\tBackground");
		sb.append("\tThreshold");
		sb.append("\tBackgroundSize");
		sb.append("\tPercentile");
		return sb.toString();
	}

	private String prefix;

	private String createResult(double[] result)
	{
		StringBuilder sb = new StringBuilder(prefix);
		sb.append((int) result[0]);
		for (int i = 1; i < result.length; i++)
		{
			sb.append('\t').append(Utils.rounded(result[i]));
		}
		return sb.toString();
	}
}
