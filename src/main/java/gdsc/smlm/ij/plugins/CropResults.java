package gdsc.smlm.ij.plugins;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;

import gdsc.core.utils.TurboList;
import gdsc.smlm.data.config.SMLMSettings.DistanceUnit;
import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.procedures.XYRResultProcedure;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.ExtendedGenericDialog;
import ij.plugin.PlugIn;

/**
 * Filters PeakFit results that are stored in memory using various fit criteria.
 */
public class CropResults implements PlugIn
{
	private static final String TITLE = "Crop Results";
	private static String inputOption = "";
	private static double border = 0;
	private static double x = 0, y = 0, width = 0, height = 0;
	private static boolean selectRegion, overwrite, useRoi;
	private static String roiImage = "";
	private boolean myUseRoi;

	private MemoryPeakResults results;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (MemoryPeakResults.isMemoryEmpty())
		{
			IJ.error(TITLE, "There are no fitting results in memory");
			return;
		}

		// Show a dialog allowing the results set to be filtered
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addMessage("Select a dataset to crop");
		ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		inputOption = ResultsManager.getInputSource(gd);
		results = ResultsManager.loadInputResults(inputOption, false);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			IJ.showStatus("");
			return;
		}

		if (!showDialog())
			return;

		cropResults();
	}

	private boolean showDialog()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		// TODO - add option to crop using the bounding rectangle of an ROI on an open image.
		// See PC-PALM Molecules.

		// Build a list of all images with a region ROI
		TurboList<String> titles = new TurboList<String>(WindowManager.getWindowCount());
		if (WindowManager.getWindowCount() > 0)
		{
			for (int imageID : WindowManager.getIDList())
			{
				ImagePlus imp = WindowManager.getImage(imageID);
				if (imp != null && imp.getRoi() != null && imp.getRoi().isArea())
					titles.add(imp.getTitle());
			}
		}

		Rectangle bounds = results.getBounds(true);

		gd.addMessage(String.format("x=%d,y=%d,w=%d,h=%d", bounds.x, bounds.y, bounds.width, bounds.height));
		gd.addNumericField("Border", border, 2);
		gd.addCheckbox("Select_region", selectRegion);
		gd.addNumericField("X", x, 2);
		gd.addNumericField("Y", y, 2);
		gd.addNumericField("Width", width, 2);
		gd.addNumericField("Height", height, 2);
		if (!titles.isEmpty())
		{
			gd.addCheckbox("Use_ROI", useRoi);
			String[] items = titles.toArray(new String[titles.size()]);
			gd.addChoice("Image", items, roiImage);
		}
		gd.addCheckbox("Overwrite", overwrite);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		border = Math.max(0, gd.getNextNumber());
		selectRegion = gd.getNextBoolean();
		x = gd.getNextNumber();
		y = gd.getNextNumber();
		width = Math.max(0, gd.getNextNumber());
		height = Math.max(0, gd.getNextNumber());
		if (!titles.isEmpty())
		{
			myUseRoi = useRoi = gd.getNextBoolean();
			roiImage = gd.getNextChoice();
		}
		overwrite = gd.getNextBoolean();

		return true;
	}

	/**
	 * Apply the filters to the data
	 */
	private void cropResults()
	{
		final MemoryPeakResults newResults = new MemoryPeakResults();

		// These bounds are integer. But this is because the results are meant to come from an image.
		Rectangle intergerBounds = results.getBounds(true);

		// The crop bounds can be floating point...

		// Border
		double xx = intergerBounds.x + border;
		double yy = intergerBounds.y + border;
		double w = Math.max(0, intergerBounds.width - 2 * border);
		double h = Math.max(0, intergerBounds.height - 2 * border);
		Rectangle2D pixelBounds = new Rectangle2D.Double(xx, yy, w, h);

		// Bounding box
		if (selectRegion)
		{
			Rectangle2D boxBounds = new Rectangle2D.Double(x, y, width, height);
			pixelBounds = pixelBounds.createIntersection(boxBounds);
		}

		// If an ROI was chosen from an image, scale the roi to the bounds of this dataset
		// and create another intersection
		if (myUseRoi)
		{
			ImagePlus imp = WindowManager.getImage(roiImage);
			if (imp != null && imp.getRoi() != null)
			{
				Rectangle roi = imp.getRoi().getBounds();
				int roiImageWidth = imp.getWidth();
				int roiImageHeight = imp.getHeight();

				double xscale = (double) roiImageWidth / intergerBounds.width;
				double yscale = (double) roiImageHeight / intergerBounds.height;

				Rectangle2D roiBounds = new Rectangle2D.Double(roi.x / xscale, roi.y / yscale, roi.width / xscale,
						roi.height / yscale);
				pixelBounds = pixelBounds.createIntersection(roiBounds);
			}
		}

		final Rectangle2D bounds = pixelBounds;
		
		if (bounds.getWidth() > 0 && bounds.getHeight() > 0)
		{
			results.forEach(DistanceUnit.PIXEL, new XYRResultProcedure()
			{
				public void executeXYR(float x, float y, PeakResult result)
				{
					if (bounds.contains(result.getXPosition(), result.getYPosition()))
						newResults.add(result);
				}
			});
		}

		newResults.copySettings(results);
		newResults.setBounds(new Rectangle((int) Math.floor(bounds.getX()), (int) Math.floor(bounds.getY()),
				(int) Math.ceil(bounds.getWidth()), (int) Math.ceil(bounds.getHeight())));
		if (!overwrite)
		{
			newResults.setName(results.getName() + " Cropped");
		}
		MemoryPeakResults.addResults(newResults);

		IJ.showStatus(newResults.size() + " Cropped localisations");
	}
}
