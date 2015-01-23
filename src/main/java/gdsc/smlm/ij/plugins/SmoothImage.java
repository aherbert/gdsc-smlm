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

import gdsc.smlm.filters.AverageFilter;
import gdsc.smlm.ij.settings.Constants;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.AWTEvent;
import java.awt.Rectangle;
import java.util.Arrays;

/**
 * Smooths the selected rectangular ROI using a mean filter.
 */
public class SmoothImage implements ExtendedPlugInFilter, DialogListener
{
	private final static String TITLE = "Smooth Image";
	public static final String[] ALGORITHMS = { "Standard", "Striped", "Rolling block", "Standard internal",
			"Striped internal", "Rolling block internal" };
	public static final int STANDARD = 0;
	public static final int STRIPED = 1;
	public static final int ROLLING_BLOCK = 2;
	public static final int STANDARD_INTERNAL = 3;
	public static final int STRIPED_INTERNAL = 4;
	public static final int ROLLING_BLOCK_INTERNAL = 5;

	private double boxSize = Prefs.get(Constants.boxSize, 1);
	private double boxSize2 = Prefs.get(Constants.boxSize2, 0);
	private int algorithm = (int) Prefs.get(Constants.algorithm, 1);
	private boolean gaussian = Prefs.getBoolean("gdsc.fitting.boxSize", true);

	private int flags = DOES_16 | DOES_8G | DOES_32 | PARALLELIZE_STACKS | FINAL_PROCESSING;
	private boolean is32bit = false;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	public int setup(String arg, ImagePlus imp)
	{
		if (arg.equals("final"))
		{
			imp.resetDisplayRange();
			imp.updateAndDraw();
			return DONE;
		}

		if (imp == null)
		{
			IJ.noImage();
			return DONE;
		}

		is32bit = (imp.getBitDepth() == 32);
		Roi roi = imp.getRoi();
		if (roi != null && roi.getType() != Roi.RECTANGLE)
		{
			IJ.error("Rectangular ROI required");
			return DONE;
		}

		return flags;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.ExtendedPlugInFilter#showDialog(ij.ImagePlus, java.lang.String,
	 * ij.plugin.filter.PlugInFilterRunner)
	 */
	public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr)
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		boxSize = Prefs.get(Constants.boxSize, 1);
		boxSize2 = Prefs.get(Constants.boxSize2, 0);
		algorithm = (int) Prefs.get(Constants.algorithm, 1);

		gd.addMessage("Smooth image:\n" + "- Within a 2n+1 box\n");
		gd.addSlider("Box_size", 0.1, 5, boxSize);
		if (is32bit)
			gd.addSlider("Box_size2", 0, 15, boxSize2);
		gd.addChoice("Algorithm", ALGORITHMS, ALGORITHMS[algorithm]);
		gd.addCheckbox("Gaussian 3x3", gaussian);

		gd.addPreviewCheckbox(pfr);
		gd.addDialogListener(this);
		gd.showDialog();

		if (gd.wasCanceled() || !dialogItemChanged(gd, null))
			return DONE;

		return IJ.setupDialog(imp, flags);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.gui.DialogListener#dialogItemChanged(ij.gui.GenericDialog, java.awt.AWTEvent)
	 */
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
	{
		boxSize = gd.getNextNumber();
		if (is32bit)
			boxSize2 = gd.getNextNumber();
		if (gd.invalidNumber())
			return false;
		algorithm = gd.getNextChoiceIndex();
		gaussian = gd.getNextBoolean();

		// Check arguments
		try
		{
			Parameters.isAboveZero("Box size", boxSize);
			Parameters.isPositive("Box size2", boxSize2);
		}
		catch (IllegalArgumentException ex)
		{
			IJ.error(TITLE, ex.getMessage());
			return false;
		}

		Prefs.set(Constants.boxSize, boxSize);
		if (is32bit)
			Prefs.set(Constants.boxSize2, boxSize2);
		Prefs.set(Constants.algorithm, algorithm);
		Prefs.set("gdsc.fitting.boxSize", gaussian);

		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	public void run(ImageProcessor ip)
	{
		Rectangle bounds = ip.getRoi();

		// Crop to the ROI
		FloatProcessor fp = ip.crop().toFloat(0, null);

		float[] data = (float[]) fp.getPixels();
		float[] data2 = null;
		int width = fp.getWidth();
		int height = fp.getHeight();

		if (is32bit && boxSize2 > boxSize)
		{
			data2 = Arrays.copyOf(data, data.length);
			smooth(data2, width, height, boxSize2);
			//new ij.ImagePlus("smoothImage2", new ij.process.FloatProcessor(width, height, data2, null)).show();
		}

		smooth(data, width, height, boxSize);
		//new ij.ImagePlus("smoothImage", new ij.process.FloatProcessor(width, height, Arrays.copyOf(data, data.length),
		//		null)).show();

		if (data2 != null)
		{
			for (int i = 0; i < data.length; i++)
				data[i] -= data2[i];
		}

		//new ImagePlus("smoothImageDiff", new FloatProcessor(width, height, data, null)).show();

		for (int index = data.length; index-- > 0;)
		{
			ip.setf(bounds.x + index % width, bounds.y + index / width, data[index]);
		}
		ip.resetMinAndMax();
	}

	private void smooth(float[] data, final int width, final int height, final double boxSize)
	{
		AverageFilter av = new AverageFilter();

		if (boxSize <= 1)
		{
			if (gaussian)
				av.blockGaussian3x3(data, width, height);
			else
				av.blockAverage3x3(data, width, height, (float) boxSize);
		}
		else
		{
			int iBoxSize = (int) boxSize;
			switch (algorithm)
			{
				case ROLLING_BLOCK_INTERNAL:
					av.rollingBlockAverageInternal(data, width, height, iBoxSize);
					break;
				case STRIPED_INTERNAL:
					av.stripedBlockAverageInternal(data, width, height, iBoxSize);
					break;
				case STANDARD_INTERNAL:
					av.blockAverageInternal(data, width, height, iBoxSize);
					break;
				case ROLLING_BLOCK:
					av.rollingBlockAverage(data, width, height, iBoxSize);
					break;
				case STRIPED:
					av.stripedBlockAverage(data, width, height, iBoxSize);
					break;
				case STANDARD:
				default:
					av.blockAverage(data, width, height, iBoxSize);
					break;
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
		// Nothing to do		
	}
}
