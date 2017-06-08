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

import gdsc.smlm.engine.DataFilter;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.filters.DataProcessor;
import gdsc.smlm.filters.DifferenceSpotFilter;
import gdsc.smlm.filters.MaximaSpotFilter;
import gdsc.smlm.filters.SingleSpotFilter;
import gdsc.smlm.ij.settings.SettingsManager;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.AWTEvent;
import java.awt.Rectangle;

/**
 * Smooths the selected rectangular ROI using a mean filter.
 */
public class SmoothImage implements ExtendedPlugInFilter, DialogListener
{
	private final static String TITLE = "Smooth Image";
	private static final String[] filterNames;
	private static final DataFilter[] filters;
	static
	{
		filters = DataFilter.values();
		filterNames = SettingsManager.dataFilterNames;
	}

	private static int filter1 = 0;
	private static double smooth1 = 1;
	private static boolean differenceFilter = false;
	private static int filter2 = 0;
	private static double smooth2 = 3;

	private int flags = DOES_16 | DOES_8G | DOES_32 | PARALLELIZE_STACKS | FINAL_PROCESSING;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	public int setup(String arg, ImagePlus imp)
	{
		if (arg.equals("final"))
		{
			//imp.resetDisplayRange();
			imp.updateAndDraw();

			return DONE;
		}

		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (imp == null)
		{
			IJ.noImage();
			return DONE;
		}

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
		// Note: We cannot use a NonBlockinnericDialog as scrolling through the image
		// throws away the snap shot. The pixel data for the previous slice is then fixed
		// with the preview. So we can only support a single slice.
		
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Smooth image:");
		gd.addChoice("Spot_filter", filterNames, filterNames[filter1]);
		gd.addSlider("Smoothing", 0, 4.5, smooth1);
		gd.addCheckbox("Difference_filter", differenceFilter);
		gd.addChoice("Spot_filter2", filterNames, filterNames[filter2]);
		gd.addSlider("Smoothing2", 1.5, 6, smooth2);

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
		filter1 = gd.getNextChoiceIndex();
		smooth1 = gd.getNextNumber();
		if (differenceFilter = gd.getNextBoolean())
		{
			filter2 = gd.getNextChoiceIndex();
			smooth2 = gd.getNextNumber();
		}
		return !gd.invalidNumber();
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

		MaximaSpotFilter filter = createSpotFilter();
		int width = fp.getWidth();
		int height = fp.getHeight();
		data = filter.preprocessData(data, width, height);

		//System.out.println(filter.getDescription());

		fp = new FloatProcessor(width, height, data);
		ip.insert(fp, bounds.x, bounds.y);
		//ip.resetMinAndMax();
		ip.setMinAndMax(fp.getMin(), fp.getMax());
	}

	private MaximaSpotFilter createSpotFilter()
	{
		final int search = 1;
		final int border = 0;
		DataProcessor processor0 = FitEngineConfiguration.createDataProcessor(border, filters[filter1], smooth1);
		if (differenceFilter)
		{
			DataProcessor processor1 = FitEngineConfiguration.createDataProcessor(border, filters[filter2], smooth2);
			return new DifferenceSpotFilter(search, border, processor0, processor1);
		}
		return new SingleSpotFilter(search, border, processor0);
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
