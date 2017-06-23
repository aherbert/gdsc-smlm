package gdsc.smlm.ij.plugins;

import gdsc.smlm.data.config.UnitConfig.DistanceUnit;

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

import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.utils.ObjectAnalyzer;
import gdsc.core.ij.Utils;
import gdsc.smlm.results.Counter;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.procedures.XYRResultProcedure;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.ExtendedGenericDialog;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

/**
 * Splits PeakFit results into separate datasets using an input mask of objects.
 */
public class SplitResults implements PlugIn
{
	private static final String TITLE = "Split Results";
	private static String inputOption = "";
	private static String objectMask = "";
	private static boolean showObjectMask = false;
	private static boolean nonMaskDataset = false;

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
		String[] items = Utils.getImageList(Utils.GREY_8_16);
		if (items.length == 0)
		{
			IJ.error(TITLE, "There are no suitable mask images");
			return;
		}

		// Show a dialog allowing the results set to be filtered
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addMessage("Select a dataset to split");
		ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);
		gd.addChoice("Object_mask", items, objectMask);
		gd.addCheckbox("Show_object_mask", showObjectMask);
		gd.addCheckbox("Non_mask_dataset", nonMaskDataset);
		gd.showDialog();
		if (gd.wasCanceled())
			return;

		inputOption = ResultsManager.getInputSource(gd);
		objectMask = gd.getNextChoice();
		showObjectMask = gd.getNextBoolean();
		nonMaskDataset = gd.getNextBoolean();

		MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, false);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			return;
		}

		ImagePlus imp = WindowManager.getImage(objectMask);
		if (imp == null)
		{
			IJ.error(TITLE, "No object mask could be found");
			return;
		}

		splitResults(results, imp.getProcessor());
	}

	private void splitResults(MemoryPeakResults results, ImageProcessor ip)
	{
		IJ.showStatus("Splitting " + Utils.pleural(results.size(), "result"));

		// Create an object mask
		ObjectAnalyzer objectAnalyzer = new ObjectAnalyzer(ip, false);

		final int maxx = ip.getWidth();
		final int maxy = ip.getHeight();

		final float scaleX = (float) results.getBounds().width / maxx;
		final float scaleY = (float) results.getBounds().height / maxy;

		// Create a results set for each object
		final int maxObject = objectAnalyzer.getMaxObject();
		final MemoryPeakResults[] resultsSet = new MemoryPeakResults[maxObject + 1];
		for (int object = 0; object <= maxObject; object++)
		{
			MemoryPeakResults newResults = new MemoryPeakResults();
			newResults.copySettings(results);
			newResults.setName(results.getName() + " " + object);
			resultsSet[object] = newResults;
		}

		final int[] mask = objectAnalyzer.getObjectMask();

		if (showObjectMask)
		{
			ImageProcessor objectIp = (maxObject <= 255) ? new ByteProcessor(maxx, maxy) : new ShortProcessor(maxx,
					maxy);
			for (int i = 0; i < mask.length; i++)
				objectIp.set(i, mask[i]);
			ImagePlus imp = Utils.display(objectMask + " Objects", objectIp);
			imp.setDisplayRange(0, maxObject);
			imp.updateAndDraw();
		}

		// Process the results mapping them to their objects
		final Counter i = new Counter();
		final int size = results.size();
		final int step = Utils.getProgressInterval(size);
		results.forEach(DistanceUnit.PIXEL, new XYRResultProcedure()
		{
			public void executeXYR(float xx, float yy, PeakResult result)
			{
				if (i.incrementAndGet() % step == 0)
					IJ.showProgress(i.getCount(), size);

				// Map to the mask objects
				final int object;
				int x = (int) (xx / scaleX);
				int y = (int) (yy / scaleY);
				if (x < 0 || x >= maxx || y < 0 || y >= maxy)
				{
					object = 0;
				}
				else
				{
					final int index = y * maxx + x;
					// Q. Is this bounds check needed given the above check shows that x,y 
					// is within the bounds of the image processor?
					if (index < 0 || index >= mask.length)
						object = 0;
					else
						object = mask[index];
				}
				resultsSet[object].add(result);
			}
		});
		IJ.showProgress(1);

		// Add the new results sets to memory
		i.reset();
		for (int object = (nonMaskDataset) ? 0 : 1; object <= maxObject; object++)
		{
			if (resultsSet[object].isNotEmpty())
			{
				MemoryPeakResults.addResults(resultsSet[object]);
				i.increment();
			}
		}

		IJ.showStatus("Split " + Utils.pleural(results.size(), "result") + " into " + Utils.pleural(i.getCount(), "set"));
	}
}
