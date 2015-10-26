package gdsc.smlm.ij.utils;

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

import gdsc.smlm.utils.Sort;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.text.TextPanel;

import java.awt.Color;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

/**
 * Attaches to a text panel and listens for mouse events. Upon double click it obtains the coordinates from a provider
 * and draws a point ROI on the named image. Supports drawing multi-point ROI when multiple lines in the table are
 * selected.
 * <p>
 * The provider should provide [slice,x,y] coordinates for the image ROI.
 */
public class ImageROIPainter implements MouseListener
{
	private TextPanel textPanel;
	private String title;
	private CoordinateProvider coordProvider;

	/**
	 * @param textPanel
	 *            The text panel to listen to for mouse events
	 * @param title
	 *            The title of the image to add the ROI to
	 * @param coordProvider
	 *            Provides coordinates from the lines selected in the text panel
	 */
	public ImageROIPainter(TextPanel textPanel, String title, CoordinateProvider coordProvider)
	{
		// Check if the image is displayed
		this.title = title;
		this.textPanel = textPanel;
		this.coordProvider = coordProvider;
		textPanel.addMouseListener(this);
	}

	public void mouseClicked(MouseEvent e)
	{
		// Show the result that was double clicked in the result table
		if (e.getClickCount() > 1)
		{
			ImagePlus imp = WindowManager.getImage(title);
			if (imp == null)
				return;

			int index = textPanel.getSelectionStart();
			if (index == -1)
				return;

			double[] position = coordProvider.getCoordinates(textPanel.getLine(index));

			if (position == null || position.length < 3)
				return;

			int slice = (int) position[0];
			double x = position[1];
			double y = position[2];

			addRoi(imp, slice, new PointRoi(x, y));

			Utils.adjustSourceRect(imp, 0, (int) x, (int) y);
		}
	}

	public void mousePressed(MouseEvent e)
	{
		// If a multiple-line selection is made then show all the points
		int index = textPanel.getSelectionStart();
		if (index == -1)
			return;
		int index2 = textPanel.getSelectionEnd();
		if (index == index2)
			return;
		ImagePlus imp = WindowManager.getImage(title);
		if (imp == null)
			return;

		// Show all
		int points = 0;
		float[] x = new float[index2 - index + 1];
		float[] y = new float[x.length];
		int[] slice = new int[x.length];
		while (index <= index2)
		{
			double[] position = coordProvider.getCoordinates(textPanel.getLine(index));

			if (position == null || position.length < 3)
				continue;

			slice[points] = (int) position[0];
			x[points] = (float) position[1];
			y[points] = (float) position[2];
			points++;
			index++;
		}

		if (points == 0)
			return;

		// Simple code to add the ROI onto a single slice: addRoi(imp, slice[0], new PointRoi(x, y, points));

		// Add the ROI to each relevant slice

		// Sort the slices
		int[] indices = new int[points];
		for (int i = 0; i < points; i++)
			indices[i] = i;

		Sort.sort(indices, slice);

		Overlay o = new Overlay();

		// Create an ROI for each slice
		int start = 0;
		for (int i = 0; i < points; i++)
		{
			if (slice[indices[i]] != slice[indices[start]])
			{
				appendRoi(x, y, slice, indices, o, start, i);
				start = i;
			}
		}
		appendRoi(x, y, slice, indices, o, start, points);

		// Choose the first slice and add the final overlay
		imp.setSlice(slice[indices[start]]);
		if (imp.getWindow() != null)
			imp.getWindow().toFront();
		o.setStrokeColor(Color.green);
		imp.setOverlay(o);
	}

	/**
	 * Adds a new ROI to the overlay using the coordinates from start to end (non-inclusive)
	 * 
	 * @param x
	 * @param y
	 * @param slice
	 * @param indices
	 * @param o
	 * @param start
	 * @param end
	 */
	private void appendRoi(float[] x, float[] y, int[] slice, int[] indices, Overlay o, int start, int end)
	{
		int p = end - start;
		float[] x2 = new float[p];
		float[] y2 = new float[p];
		for (int j = start, ii = 0; j < end; j++, ii++)
		{
			x2[ii] = x[indices[start]];
			y2[ii] = y[indices[start]];
		}
		PointRoi roi = new PointRoi(x2, y2, p);
		roi.setPosition(slice[indices[start]]);
		o.add(roi);
	}

	public static void addRoi(ImagePlus imp, int slice, PointRoi roi)
	{
		if (imp != null && slice > 0 && slice <= imp.getStackSize())
		{
			imp.setSlice(slice);
			if (imp.getWindow() != null)
				imp.getWindow().toFront();

			if (roi != null)
			{
				//imp.setRoi(roi);

				if (imp.getStackSize() > 1)
					roi.setPosition(slice);
				Overlay o = new Overlay(roi);
				o.setStrokeColor(Color.green);
				imp.setOverlay(o);
			}
			else
			{
				imp.setOverlay(null);
			}
		}
	}

	public void mouseReleased(MouseEvent e)
	{
		// Ignore
	}

	public void mouseEntered(MouseEvent e)
	{
		// Ignore
	}

	public void mouseExited(MouseEvent e)
	{
		// Ignore
	}

	/**
	 * @return the title of the image
	 */
	public String getTitle()
	{
		return title;
	}

	/**
	 * @param title
	 *            the title of the image
	 */
	public void setTitle(String title)
	{
		this.title = title;
	}

	/**
	 * @return the coordProvider
	 */
	public CoordinateProvider getCoordProvider()
	{
		return coordProvider;
	}

	/**
	 * @param coordProvider
	 *            the coordProvider to set
	 */
	public void setCoordProvider(CoordinateProvider coordProvider)
	{
		this.coordProvider = coordProvider;
	}
}
