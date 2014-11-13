package gdsc.smlm.ij.plugins;

import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.utils.NoiseEstimator;
import gdsc.smlm.utils.Statistics;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.ImageProcessor;
import ij.text.TextWindow;
import ij.util.Tools;

import java.awt.AWTEvent;
import java.awt.Color;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

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
 * Contains methods to find the noise in the provided image data.
 */
public class Noise implements ExtendedPlugInFilter, DialogListener
{
	private List<double[]> results;
	private final int FLAGS = DOES_8G | DOES_16 | DOES_32 | PARALLELIZE_STACKS | FINAL_PROCESSING | NO_CHANGES;
	private PlugInFilterRunner pfr;
	private ImagePlus imp;
	private GenericDialog gd;
	private static int algorithm = 0;
	private static int algorithm2 = 1;
	private static int lowestPixelsRange = 6;

	public int setup(String arg, ImagePlus imp)
	{
		if (arg.equalsIgnoreCase("final"))
		{
			showResults();
		}
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

			gd = new GenericDialog("Noise Estimator");
			gd.addHelp(About.HELP_URL);

			NoiseEstimator.Method[] methods = NoiseEstimator.Method.values();
			String[] methodNames = new String[methods.length];
			for (int i = 0; i < methods.length; i++)
			{
				methodNames[i] = methods[i].toString();
			}

			gd.addChoice("Method1 (blue)", methodNames, methodNames[algorithm]);
			gd.addChoice("Method2 (red)", methodNames, methodNames[algorithm2]);
			gd.addSlider("Lowest_radius", 1, 15, lowestPixelsRange);

			//gd.addPreviewCheckbox(pfr);
			gd.addDialogListener(this);
			gd.addMessage("Click OK to compute noise table using all methods");
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
		algorithm = gd.getNextChoiceIndex();
		algorithm2 = gd.getNextChoiceIndex();
		lowestPixelsRange = (int) gd.getNextNumber();
		if (gd.invalidNumber() || lowestPixelsRange < 1)
			return false;
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
		NoiseEstimator.Method method1 = NoiseEstimator.Method.values()[algorithm];
		NoiseEstimator.Method method2 = NoiseEstimator.Method.values()[algorithm2];

		IJ.showStatus("Estimating noise ...");

		boolean twoMethods = method1 != method2;
		boolean preserveResiduals = method1.name().contains("Residuals") && method2.name().contains("Residuals") &&
				twoMethods;

		int start = imp.getCurrentSlice();
		int end = FastMath.min(imp.getStackSize(), start + 100);
		int size = end - start + 1;
		double[] xValues = new double[size];
		double[] yValues1 = new double[size];
		double[] yValues2 = (twoMethods) ? new double[size] : null;

		ImageStack stack = imp.getImageStack();
		for (int slice = start, i = 0; slice <= end; slice++, i++)
		{
			IJ.showProgress(i, size);
			ImageProcessor ip = stack.getProcessor(slice);
			NoiseEstimator ne = new NoiseEstimator((float[]) ip.toFloat(0, null).getPixels(), ip.getWidth(),
					ip.getHeight());
			ne.preserveResiduals = preserveResiduals;
			ne.setRange(lowestPixelsRange);
			xValues[i] = slice;
			yValues1[i] = ne.getNoise(method1);
			if (twoMethods)
				yValues2[i] = ne.getNoise(method2);
		}
		IJ.showProgress(1);

		IJ.showStatus("Plotting noise ...");

		// Get limits
		double[] a = Tools.getMinMax(xValues);
		double[] b1 = Tools.getMinMax(yValues1);
		if (twoMethods)
		{
			double[] b2 = Tools.getMinMax(yValues2);
			b1[0] = FastMath.min(b1[0], b2[0]);
			b1[1] = FastMath.max(b1[1], b2[1]);
		}

		String title = imp.getTitle() + " Noise";
		Plot plot = new Plot(title, "Slice", "Noise", xValues, yValues1);
		double range = b1[1] - b1[0];
		if (range == 0)
			range = 1;
		plot.setLimits(a[0], a[1], b1[0]-0.05*range, b1[1]+0.05*range);
		plot.setColor(Color.blue);
		plot.draw();
		String label = String.format("Blue = %s", Utils.rounded(new Statistics(yValues1).getMean()));
		if (twoMethods)
		{
			plot.setColor(Color.red);
			plot.addPoints(xValues, yValues2, Plot.LINE);
			label += String.format(", Red = %s", Utils.rounded(new Statistics(yValues2).getMean()));
		}
		plot.addLabel(0, 0, label);

		Utils.display(title, plot);

		IJ.showStatus("");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	public void run(ImageProcessor ip)
	{
		// Perform all methods and add to the results
		double[] result = new double[NoiseEstimator.Method.values().length + 1];
		int i = 0;
		result[i++] = (pfr == null) ? 1 : pfr.getSliceNumber();
		NoiseEstimator ne = new NoiseEstimator((float[]) ip.toFloat(0, null).getPixels(), ip.getWidth(), ip.getHeight());
		ne.preserveResiduals = true;
		for (NoiseEstimator.Method m : NoiseEstimator.Method.values())
		{
			result[i++] = ne.getNoise(m);
		}
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

		// Slow when there are lots of results ... Could change the output options in the future 
		TextWindow tw = new TextWindow(imp.getTitle() + " Noise", createHeader(), "", 800, 400);
		for (double[] result : results)
		{
			tw.append(createResult(result));
		}

		// TODO - ImageJ plotting is not very good. Change to use a Java plotting library 
		//plotResults();
	}

	private String createHeader()
	{
		StringBuilder sb = new StringBuilder("Slice");
		for (NoiseEstimator.Method m : NoiseEstimator.Method.values())
		{
			sb.append("\t").append(m);
		}
		return sb.toString();
	}

	private String createResult(double[] result)
	{
		StringBuilder sb = new StringBuilder("" + (int) result[0]);
		for (int i = 1; i < result.length; i++)
		{
			sb.append("\t").append(result[i]);
		}
		return sb.toString();
	}
}
