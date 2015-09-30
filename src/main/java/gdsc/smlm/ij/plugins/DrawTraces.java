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
import gdsc.smlm.utils.Sort;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.LUT;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;

/**
 * Compares the coordinates in sets of traced results and computes the match statistics.
 */
public class DrawTraces implements PlugIn
{
	private static final String TITLE = "Draw Traces";
	private static final String[] sorts = new String[] { "None", "ID", "Time", "Size", "Length" };
	private static final String[] luts = new String[] { "red-hot", "ice", "rainbow", "fire", "red-yellow" };

	private static String inputOption = "";
	private static String title = "";
	private static int minSize = 2;
	private static int sort = 0;
	private static boolean splineFit = false;
	private static int lut = 0;

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
			if (traces[i].size() >= minSize)
			{
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
			imp = Utils.display(TITLE, new ByteProcessor(w, h));
		}

		final float xScale = (float) (imp.getWidth() / bounds.getWidth());
		final float yScale = (float) (imp.getHeight() / bounds.getHeight());

		// Create ROIs and store data to sort them
		Roi[] rois = new Roi[count];
		int[] indices = Utils.newArray(count, 0, 1);
		double[] values = new double[count];
		for (int i = 0; i < count; i++)
		{
			Trace trace = traces[i];
			int nPoints = trace.size();
			float[] xPoints = new float[nPoints];
			float[] yPoints = new float[nPoints];
			int j = 0;
			for (PeakResult result : trace.getPoints())
			{
				xPoints[j] = (result.getXPosition() - bounds.x) * xScale;
				yPoints[j] = (result.getYPosition() - bounds.y) * yScale;
				j++;
			}
			PolygonRoi roi = new PolygonRoi(xPoints, yPoints, nPoints, Roi.POLYLINE);
			if (splineFit)
				roi.fitSpline();
			rois[i] = roi;
			switch (sort)
			{
				case 0:
				default:
					break;
				case 1: // Sort by ID
					values[i] = traces[i].getId();
					break;
				case 2: // Sort by time
					values[i] = traces[i].getHead().peak;
					break;
				case 3: // Sort by size descending
					values[i] = -traces[i].size();
					break;
				case 4: // Sort by length descending
					values[i] = -roi.getLength();
					break;
			}
		}

		if (sort > 0)
			Sort.sort(indices, values);

		// Draw the traces as ROIs on an overlay
		Overlay o = new Overlay();
		LUT lut = createLUT();
		double scale = 256.0 / count;
		for (int i = 0; i < count; i++)
		{
			// Get a colour using the trace number
			Color c = new Color(lut.getRGB((int) (i * scale)));
			Roi roi = rois[indices[i]];
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
		gd.addSlider("Min_size", 2, 15, minSize);
		gd.addChoice("Sort", sorts, sorts[sort]);
		gd.addCheckbox("Spline_fit", splineFit);
		gd.addChoice("LUT", luts, luts[lut]);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption = ResultsManager.getInputSource(gd);
		title = gd.getNextChoice();
		minSize = (int) Math.abs(gd.getNextNumber());
		sort = gd.getNextChoiceIndex();
		splineFit = gd.getNextBoolean();
		lut = gd.getNextChoiceIndex();

		return true;
	}

	/**
	 * Build a custom LUT that helps show the classes
	 * 
	 * @return
	 */
	private LUT createLUT()
	{
		byte[] reds = new byte[256];
		byte[] greens = new byte[256];
		byte[] blues = new byte[256];
		int nColors;
		switch (lut)
		{
			case 4: // red-yellow
				nColors = setColours(reds, greens, blues, Color.red, Color.yellow);
				break;
			case 3:
				nColors = fire(reds, greens, blues);
				break;
			case 2:
				nColors = rainbow(reds, greens, blues);
				break;
			case 1:
				nColors = ice(reds, greens, blues);
				break;
			case 0: // red-hot
			default: 
				nColors = setColours(reds, greens, blues, Color.red, Color.yellow, Color.WHITE);
				break;
		}
		if (nColors < 256)
			interpolate(reds, greens, blues, nColors);
		return new LUT(reds, greens, blues);
	}

	private int rainbow(byte[] reds, byte[] greens, byte[] blues)
	{
		// Using HSV vary the Hue from 300 (magenta) to Red (0)
		int n = 0;
		for (int h = 300; h >= 0; h -= 2)
		{
			Color c = Color.getHSBColor(h / 360.0f, 1, 1);
			reds[n] = (byte) c.getRed();
			greens[n] = (byte) c.getGreen();
			blues[n] = (byte) c.getBlue();
			n++;
		}
		return n;
	}

	private int setColours(byte[] reds, byte[] greens, byte[] blues, Color... colours)
	{
		for (int i = 0; i < colours.length; i++)
		{
			reds[i] = (byte) colours[i].getRed();
			greens[i] = (byte) colours[i].getGreen();
			blues[i] = (byte) colours[i].getBlue();
		}
		return colours.length;
	}

	/**
	 * Copied from ij.plugin.LutLoader
	 * 
	 * @param reds
	 * @param greens
	 * @param blues
	 * @return
	 */
	private int ice(byte[] reds, byte[] greens, byte[] blues)
	{
		int[] r = { 0, 0, 0, 0, 0, 0, 19, 29, 50, 48, 79, 112, 134, 158, 186, 201, 217, 229, 242, 250, 250, 250, 250,
				251, 250, 250, 250, 250, 251, 251, 243, 230 };
		int[] g = { 156, 165, 176, 184, 190, 196, 193, 184, 171, 162, 146, 125, 107, 93, 81, 87, 92, 97, 95, 93, 93,
				90, 85, 69, 64, 54, 47, 35, 19, 0, 4, 0 };
		int[] b = { 140, 147, 158, 166, 170, 176, 209, 220, 234, 225, 236, 246, 250, 251, 250, 250, 245, 230, 230, 222,
				202, 180, 163, 142, 123, 114, 106, 94, 84, 64, 26, 27 };
		for (int i = 0; i < r.length; i++)
		{
			reds[i] = (byte) r[i];
			greens[i] = (byte) g[i];
			blues[i] = (byte) b[i];
		}
		return r.length;
	}

	/**
	 * Adapted from ij.plugin.LutLoader to remove the dark colours
	 * 
	 * @param reds
	 * @param greens
	 * @param blues
	 * @return
	 */
	private int fire(byte[] reds, byte[] greens, byte[] blues)
	{
		int[] r = { //0, 0, 1, 25, 49, 
		73, 98, 122, 146, 162, 173, 184, 195, 207, 217, 229, 240, 252, 255, 255, 255, 255, 255, 255, 255, 255, 255,
				255, 255, 255, 255, 255 };
		int[] g = { //0, 0, 0, 0, 0, 
		0, 0, 0, 0, 0, 0, 0, 0, 14, 35, 57, 79, 101, 117, 133, 147, 161, 175, 190, 205, 219, 234, 248, 255, 255, 255,
				255 };
		int[] b = { //0, 61, 96, 130, 165, 
		192, 220, 227, 210, 181, 151, 122, 93, 64, 35, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 35, 98, 160, 223, 255 };
		for (int i = 0; i < r.length; i++)
		{
			reds[i] = (byte) r[i];
			greens[i] = (byte) g[i];
			blues[i] = (byte) b[i];
		}
		return r.length;
	}

	/**
	 * Copied from ij.plugin.LutLoader.
	 * 
	 * @param reds
	 * @param greens
	 * @param blues
	 * @param nColors
	 */
	private void interpolate(byte[] reds, byte[] greens, byte[] blues, int nColors)
	{
		byte[] r = new byte[nColors];
		byte[] g = new byte[nColors];
		byte[] b = new byte[nColors];
		System.arraycopy(reds, 0, r, 0, nColors);
		System.arraycopy(greens, 0, g, 0, nColors);
		System.arraycopy(blues, 0, b, 0, nColors);
		double scale = nColors / 256.0;
		int i1, i2;
		double fraction;
		for (int i = 0; i < 256; i++)
		{
			i1 = (int) (i * scale);
			i2 = i1 + 1;
			if (i2 == nColors)
				i2 = nColors - 1;
			fraction = i * scale - i1;
			//IJ.write(i+" "+i1+" "+i2+" "+fraction);
			reds[i] = (byte) ((1.0 - fraction) * (r[i1] & 255) + fraction * (r[i2] & 255));
			greens[i] = (byte) ((1.0 - fraction) * (g[i1] & 255) + fraction * (g[i2] & 255));
			blues[i] = (byte) ((1.0 - fraction) * (b[i1] & 255) + fraction * (b[i2] & 255));
		}
	}
}
