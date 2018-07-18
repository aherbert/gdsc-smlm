/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm.ij.utils;

import java.awt.Color;

import gdsc.core.ij.Utils;
import gdsc.core.utils.Sort;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.text.TextPanel;

/**
 * Attaches to a text panel and listens for mouse events. Upon double click it obtains the coordinates from a provider
 * and draws a point ROI on the named image. Supports drawing multi-point ROI when multiple lines in the table are
 * selected.
 * <p>
 * The provider should provide [slice,x,y] coordinates for the image ROI.
 */
public class ImageROIPainter extends TextPanelMouseListener
{
	private final TextPanel textPanel;
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
		super(textPanel);

		// Check if the image is displayed
		this.title = title;
		this.textPanel = textPanel;
		this.coordProvider = coordProvider;
		textPanel.addMouseListener(this);
	}

	/**
	 * Trigger the ROI painter using the selected index from the text panel.
	 *
	 * @param selectedIndex
	 *            the selected index
	 * @see gdsc.smlm.ij.utils.TextPanelMouseListener#selected(int)
	 */
	@Override
	public void selected(int selectedIndex)
	{
		if (selectedIndex < 0 || selectedIndex >= textPanel.getLineCount())
			return;

		final ImagePlus imp = WindowManager.getImage(title);
		if (imp == null)
			return;

		final double[] position = coordProvider.getCoordinates(textPanel.getLine(selectedIndex));

		if (position == null || position.length < 3)
			return;

		final int slice = (int) position[0];
		final double x = position[1];
		final double y = position[2];

		addRoi(imp, slice, new PointRoi(x, y));

		Utils.adjustSourceRect(imp, 0, (int) x, (int) y);
	}

	/**
	 * Trigger the ROI painter using the selection from the text panel.
	 *
	 * @param selectionStart
	 *            the selection start
	 * @param selectionEnd
	 *            the selection end
	 * @see gdsc.smlm.ij.utils.TextPanelMouseListener#selected(int, int)
	 */
	@Override
	public void selected(int selectionStart, int selectionEnd)
	{
		if (selectionStart < 0 || selectionStart >= textPanel.getLineCount())
			return;
		if (selectionEnd < selectionStart || selectionEnd >= textPanel.getLineCount())
			return;
		final ImagePlus imp = WindowManager.getImage(title);
		if (imp == null)
			return;

		// Show all
		int points = 0;
		final float[] x = new float[selectionEnd - selectionStart + 1];
		final float[] y = new float[x.length];
		final int[] slice = new int[x.length];
		while (selectionStart <= selectionEnd)
		{
			final double[] position = coordProvider.getCoordinates(textPanel.getLine(selectionStart));

			if (position == null || position.length < 3)
				continue;

			slice[points] = (int) position[0];
			x[points] = (float) position[1];
			y[points] = (float) position[2];
			points++;
			selectionStart++;
		}

		if (points == 0)
			return;

		// Simple code to add the ROI onto a single slice: addRoi(imp, slice[0], new PointRoi(x, y, points));

		// Add the ROI to each relevant slice

		// Sort the slices
		final int[] indices = new int[points];
		for (int i = 0; i < points; i++)
			indices[i] = i;

		Sort.sort(indices, slice);

		final Overlay o = new Overlay();

		// Create an ROI for each slice
		int start = 0;
		for (int i = 0; i < points; i++)
			if (slice[indices[i]] != slice[indices[start]])
			{
				appendRoi(x, y, slice, indices, o, start, i);
				start = i;
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
	 * Adds a new ROI to the overlay using the coordinates from start to end (non-inclusive).
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param slice
	 *            the slice
	 * @param indices
	 *            the indices
	 * @param o
	 *            the o
	 * @param start
	 *            the start
	 * @param end
	 *            the end
	 */
	private static void appendRoi(float[] x, float[] y, int[] slice, int[] indices, Overlay o, int start, int end)
	{
		final int p = end - start;
		final float[] x2 = new float[p];
		final float[] y2 = new float[p];
		for (int j = start, ii = 0; j < end; j++, ii++)
		{
			x2[ii] = x[indices[j]];
			y2[ii] = y[indices[j]];
		}
		final PointRoi roi = new PointRoi(x2, y2, p);
		roi.setPosition(slice[indices[start]]);
		o.add(roi);
	}

	/**
	 * Adds the roi to the specified slice in the image using an overlay.
	 *
	 * @param imp
	 *            the image
	 * @param slice
	 *            the slice
	 * @param roi
	 *            the roi
	 */
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
				final Overlay o = new Overlay(roi);
				o.setStrokeColor(Color.green);
				imp.setOverlay(o);
			}
			else
				imp.setOverlay(null);
		}
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
