package gdsc.smlm.ij.plugins;

/*----------------------------------------------------------------------------- 
 * GDSC Plugins for ImageJ
 * 
 * Copyright (C) 2011 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;
import gdsc.smlm.results.Trace;
import gdsc.smlm.results.TraceManager;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;

/**
 * Compares the coordinates in sets of traced results and computes the match statistics.
 */
public class DrawTraces implements PlugIn
{
	private static String TITLE = "Draw Traces";

	private static String inputOption = "";
	private static String title = "";

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		if (MemoryPeakResults.countMemorySize() == 0)
		{
			IJ.error(TITLE, "No localisations in memory");
			return;
		}

		if (!showDialog())
			return;

		// Load the results
		MemoryPeakResults results = ResultsManager.loadInputResults(inputOption, false);
		if (results == null || results.size() == 0)
		{
			IJ.error(TITLE, "No results could be loaded");
			return;
		}

		// Get the traces
		Trace[] traces = TraceManager.convert(results);
		if (traces == null || traces.length == 0)
		{
			IJ.error(TITLE, "No traces could be loaded");
			return;
		}

		// Filter traces to a min size
		int count = 0;
		for (int i = 0; i < traces.length; i++)
		{
			if (traces[i].size() > 1)
			{
				Utils.log("Trace ID %d = %d", traces[i].getId(), traces[i].size());
				traces[count++] = traces[i];
			}
		}

		if (count == 0)
		{
			IJ.error(TITLE, "No traces achieved the minimum size");
			return;
		}

		Utils.log("Traces %d / %d (%d)", count, traces.length, results.size());
		
		Rectangle bounds = results.getBounds(true);
		ImagePlus imp = WindowManager.getImage(title);
		if (imp == null)
		{
			// Create a default image using 100 pixels as the longest edge
			double maxD = (bounds.width > bounds.height) ? bounds.width : bounds.height;
			int w, h;
			if (maxD == 0)
			{
				w = h = 100;
			}
			else
			{
				w = (int) (100.0 * bounds.width / maxD);
				h = (int) (100.0 * bounds.height / maxD);
			}
			imp = new ImagePlus(TITLE, new ByteProcessor(w, h));
			imp.show();
		}

		final float xScale = (float) (imp.getWidth() / bounds.getWidth());
		final float yScale = (float) (imp.getHeight() / bounds.getHeight());

		// Draw the traces as ROIs on an overlay
		Overlay o = new Overlay();
		for (int i = 0; i < count; i++)
		{
			Trace trace = traces[i];
			int nPoints = trace.size();
			Utils.log("Trace %d  = %d", i+1, nPoints);
			float[] xPoints = new float[nPoints];
			float[] yPoints = new float[nPoints];
			int j = 0;
			for (PeakResult result : trace.getPoints())
			{
				xPoints[j] = (result.getXPosition() - bounds.x) * xScale;
				yPoints[j] = (result.getYPosition() - bounds.y) * yScale;
				j++;
			}
			// TODO - Also try Roi.FREELINE
			PolygonRoi roi = new PolygonRoi(xPoints, yPoints, nPoints, Roi.POLYLINE);
			// TODO - Get a colour using the trace number
			Color c = Color.YELLOW;
			roi.setStrokeColor(c);
			o.add(roi);
		}

		imp.setOverlay(o);

		IJ.showStatus("Finished " + Utils.pleural(count, "trace"));
	}

	private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);

		ArrayList<String> titles = new ArrayList<String>(WindowManager.getImageCount());
		titles.add("[None]");
		int[] idList = WindowManager.getIDList();
		if (idList != null)
			for (int id : idList)
			{
				ImagePlus imp = WindowManager.getImage(id);
				if (imp != null)
					titles.add(imp.getTitle());
			}

		gd.addMessage("Draw the traces on an image");
		ResultsManager.addInput(gd, "Results", inputOption, InputSource.MEMORY_CLUSTERED);
		gd.addChoice("Image", titles.toArray(new String[0]), title);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption = ResultsManager.getInputSource(gd);
		title = gd.getNextChoice();

		return true;
	}
}
