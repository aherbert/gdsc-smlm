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
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.ImageWindow;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.FloatPolygon;
import ij.process.LUT;

import java.awt.Color;
import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;

/**
 * Compares the coordinates in sets of traced results and computes the match statistics.
 */
public class DrawClusters implements PlugIn
{
	private static final String TITLE = "Draw Clusters";
	private static final String[] sorts = new String[] { "None", "ID", "Time", "Size", "Length", "MSD", "Mean/Frame" };
	private static final String[] luts = new String[] { "Red-Hot", "Ice", "Rainbow", "Fire", "Red-Yellow", "Red",
			"Green", "Blue", "Cyan", "Magenta", "Yellow" };

	private static String inputOption = "";
	private static String title = "";
	private static int imageSize = 20;
	private static boolean expandToSingles = false;
	private static int minSize = 2;
	private static int maxSize = 0;
	private static boolean drawLines = true;
	private static int sort = 0;
	private static boolean splineFit = false;
	private static boolean useStackPosition = false;
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
		int maxFrame = 0;
		int count = 0;
		final int myMaxSize = (maxSize < minSize) ? Integer.MAX_VALUE : maxSize;
		final boolean myDrawLines = (myMaxSize < 2) ? false : drawLines;
		for (int i = 0; i < traces.length; i++)
		{
			if (expandToSingles)
				traces[i].expandToSingles();
			if (traces[i].size() >= minSize && traces[i].size() <= myMaxSize)
			{
				traces[count++] = traces[i];
				traces[i].sort();
				if (maxFrame < traces[i].getTail().peak)
					maxFrame = traces[i].getTail().peak;
			}
		}

		if (count == 0)
		{
			IJ.error(TITLE, "No traces achieved the size limits");
			return;
		}

		String msg = String.format(TITLE + ": %d / %s (%s)", count, Utils.pleural(traces.length, "trace"),
				Utils.pleural(results.size(), "localisation"));
		IJ.showStatus(msg);
		//Utils.log(msg);

		Rectangle bounds = results.getBounds(true);
		ImagePlus imp = WindowManager.getImage(title);
		boolean isUseStackPosition = useStackPosition;
		if (imp == null)
		{
			// Create a default image using 100 pixels as the longest edge
			double maxD = (bounds.width > bounds.height) ? bounds.width : bounds.height;
			int w, h;
			if (maxD == 0)
			{
				w = h = imageSize;
			}
			else
			{
				w = (int) (imageSize * bounds.width / maxD);
				h = (int) (imageSize * bounds.height / maxD);
			}
			ByteProcessor bp = new ByteProcessor(w, h);
			if (isUseStackPosition)
			{
				ImageStack stack = new ImageStack(w, h, maxFrame);
				for (int i = 1; i <= maxFrame; i++)
					stack.setPixels(bp.getPixels(), i); // Do not clone as the image is empty
				imp = Utils.display(TITLE, stack);
			}
			else
				imp = Utils.display(TITLE, bp);

			// Enlarge
			ImageWindow iw = imp.getWindow();
			for (int i = 9; i-- > 0 && iw.getWidth() < 500 && iw.getHeight() < 500;)
			{
				iw.getCanvas().zoomIn(imp.getWidth() / 2, imp.getHeight() / 2);
			}
		}
		else
		{
			// Check if the image has enough frames for all the traces
			if (maxFrame > imp.getNFrames())
				isUseStackPosition = false;
		}

		final float xScale = (float) (imp.getWidth() / bounds.getWidth());
		final float yScale = (float) (imp.getHeight() / bounds.getHeight());

		// Create ROIs and store data to sort them
		Roi[] rois = new Roi[count];
		int[][] frames = (isUseStackPosition) ? new int[count][] : null;
		int[] indices = Utils.newArray(count, 0, 1);
		double[] values = new double[count];
		for (int i = 0; i < count; i++)
		{
			Trace trace = traces[i];
			int nPoints = trace.size();
			float[] xPoints = new float[nPoints];
			float[] yPoints = new float[nPoints];
			int j = 0;
			if (isUseStackPosition)
				frames[i] = new int[nPoints];
			for (PeakResult result : trace.getPoints())
			{
				xPoints[j] = (result.getXPosition() - bounds.x) * xScale;
				yPoints[j] = (result.getYPosition() - bounds.y) * yScale;
				if (isUseStackPosition)
					frames[i][j] = result.peak;
				j++;
			}
			Roi roi;
			if (myDrawLines)
			{
				roi = new PolygonRoi(xPoints, yPoints, nPoints, Roi.POLYLINE);
				if (splineFit)
					((PolygonRoi) roi).fitSpline();
			}
			else
			{
				roi = new PointRoi(xPoints, yPoints, nPoints);
				((PointRoi) roi).setHideLabels(true);
			}

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
				case 5: // Mean Square Displacement
					values[i] = -traces[i].getMSD();
					break;
				case 6: // Mean / Frame
					values[i] = -traces[i].getMeanPerFrame();
					break;
			}
		}

		if (sort > 0)
			Sort.sort(indices, values);

		// Draw the traces as ROIs on an overlay
		Overlay o = new Overlay();
		LUT lut = createLUT();
		double scale = 256.0 / count;
		if (isUseStackPosition)
		{
			// Add the tracks on the frames containing the results
			final boolean isHyperStack = imp.isDisplayedHyperStack();
			for (int i = 0; i < count; i++)
			{
				final int index = indices[i];
				final Color c = new Color(lut.getRGB((int) (i * scale)));
				final PolygonRoi roi = (PolygonRoi) rois[index];
				roi.setFillColor(c);
				roi.setStrokeColor(c);
				final FloatPolygon fp = roi.getNonSplineFloatCoordinates();
				final Rectangle2D.Double pos = roi.getFloatBounds();
				// For each frame in the track, add the ROI track and a point ROI for the current position
				for (int j = 0; j < frames[index].length; j++)
				{
					addToOverlay(o, (Roi) roi.clone(), isHyperStack, frames[index][j]);
					PointRoi pointRoi = new PointRoi(pos.x + fp.xpoints[j], pos.y + fp.ypoints[j]);
					pointRoi.setFillColor(c);
					pointRoi.setStrokeColor(Color.black);
					addToOverlay(o, pointRoi, isHyperStack, frames[index][j]);
				}
			}
		}
		else
		{
			// Add the tracks as a single overlay
			for (int i = 0; i < count; i++)
			{
				final Roi roi = rois[indices[i]];
				roi.setStrokeColor(new Color(lut.getRGB((int) (i * scale))));
				o.add(roi);
			}
		}
		imp.setOverlay(o);

		IJ.showStatus(msg);
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

		gd.addMessage("Draw the clusters on an image");
		ResultsManager.addInput(gd, "Input", inputOption, InputSource.MEMORY_CLUSTERED);
		gd.addChoice("Image", titles.toArray(new String[0]), title);
		gd.addNumericField("Image_size", imageSize, 0);
		gd.addCheckbox("Expand_to_singles", expandToSingles);
		gd.addSlider("Min_size", 1, 15, minSize);
		gd.addSlider("Max_size", 0, 20, maxSize);
		gd.addCheckbox("Traces (draw lines)", drawLines);
		gd.addChoice("Sort", sorts, sorts[sort]);
		gd.addCheckbox("Spline_fit (traces only)", splineFit);
		gd.addCheckbox("Use_stack_position", useStackPosition);
		gd.addChoice("LUT", luts, luts[lut]);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		inputOption = ResultsManager.getInputSource(gd);
		title = gd.getNextChoice();
		imageSize = (int) Math.abs(gd.getNextNumber());
		if (imageSize < 1)
			imageSize = 1;
		expandToSingles = gd.getNextBoolean();
		minSize = (int) Math.abs(gd.getNextNumber());
		maxSize = (int) Math.abs(gd.getNextNumber());
		drawLines = gd.getNextBoolean();
		sort = gd.getNextChoiceIndex();
		splineFit = gd.getNextBoolean();
		useStackPosition = gd.getNextBoolean();
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
			case 0: // red-hot
			default:
				nColors = setColours(reds, greens, blues, Color.red, Color.yellow, Color.WHITE);
				break;
			case 1:
				nColors = ice(reds, greens, blues);
				break;
			case 2:
				nColors = rainbow(reds, greens, blues);
				break;
			case 3:
				nColors = fire(reds, greens, blues);
				break;
			case 4: // red-yellow
				nColors = setColours(reds, greens, blues, Color.red, Color.yellow);
				break;
			case 5:
				nColors = setColours(reds, greens, blues, Color.red);
				break;
			case 6:
				nColors = setColours(reds, greens, blues, Color.green);
				break;
			case 7:
				nColors = setColours(reds, greens, blues, Color.blue);
				break;
			case 8:
				nColors = setColours(reds, greens, blues, Color.cyan);
				break;
			case 9:
				nColors = setColours(reds, greens, blues, Color.magenta);
				break;
			case 10:
				nColors = setColours(reds, greens, blues, Color.yellow);
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
		int n = 0;
		if (colours.length == 1)
		{
			reds[n] = (byte) (colours[0].getRed() / 2);
			greens[n] = (byte) (colours[0].getGreen() / 2);
			blues[n] = (byte) (colours[0].getBlue() / 2);
			n++;
		}

		for (Color colour : colours)
		{
			reds[n] = (byte) colour.getRed();
			greens[n] = (byte) colour.getGreen();
			blues[n] = (byte) colour.getBlue();
			n++;
		}
		return n;
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

	private void addToOverlay(Overlay o, Roi roi, boolean isHyperStack, int frame)
	{
		if (isHyperStack)
			roi.setPosition(0, 0, frame);
		else
			roi.setPosition(frame);
		o.add(roi);
	}
}
