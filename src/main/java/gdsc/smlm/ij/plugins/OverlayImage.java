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
package gdsc.smlm.ij.plugins;

import java.awt.Rectangle;
import java.util.Arrays;

import ij.IJ;
import ij.ImagePlus;
import ij.Undo;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.ImageRoi;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.plugin.PlugIn;

/**
 * This plugin is extracted from ij.plugins.OverlayCommands to allow an image to be added with a transparent
 * background.
 */
public class OverlayImage implements PlugIn
{
	private static int opacity = 100;
	private static boolean transparent = true, replace = true;
	private static String title = "";

	@Override
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		addImage();
	}

	/**
	 * Adapted from ij.plugins.OverlayCommands#addImage(boolean) with the additional option for setting the zero pixels
	 * to transparent.
	 */
	void addImage()
	{
		ImagePlus imp = IJ.getImage();
		int[] wList = WindowManager.getIDList();
		if (wList == null || wList.length < 2)
		{
			IJ.error("Add Image...", "The command requires at least two open images.");
			return;
		}
		String[] titles = new String[wList.length];
		int count = 0;
		for (int i = 0; i < wList.length; i++)
		{
			ImagePlus imp2 = WindowManager.getImage(wList[i]);
			if (imp2 != null && imp2 != imp && imp.getWidth() >= imp2.getWidth() && imp.getHeight() >= imp2.getHeight())
				titles[count++] = imp2.getTitle();
		}
		if (count < 1)
		{
			IJ.error("Add Image...", "The command requires at least one valid overlay image.");
			return;
		}
		titles = Arrays.copyOf(titles, count);

		int x = 0, y = 0;
		Roi roi = imp.getRoi();
		if (roi != null && roi.isArea())
		{
			Rectangle r = roi.getBounds();
			x = r.x;
			y = r.y;
		}

		GenericDialog gd = new GenericDialog("Add Image...");
		gd.addChoice("Image to add:", titles, title);
		gd.addNumericField("X location:", x, 0);
		gd.addNumericField("Y location:", y, 0);
		gd.addNumericField("Opacity (0-100%):", opacity, 0);
		gd.addCheckbox("Transparent background", transparent);
		gd.addCheckbox("Replace overlay", replace);

		gd.showDialog();
		if (gd.wasCanceled())
			return;

		title = gd.getNextChoice();
		x = (int) gd.getNextNumber();
		y = (int) gd.getNextNumber();
		opacity = (int) gd.getNextNumber();
		transparent = gd.getNextBoolean();
		replace = gd.getNextBoolean();

		ImagePlus overlay = WindowManager.getImage(title);
		if (overlay == imp)
		{
			IJ.error("Add Image...", "Image to be added cannot be the same as\n\"" + imp.getTitle() + "\".");
			return;
		}
		if (overlay.getWidth() > imp.getWidth() && overlay.getHeight() > imp.getHeight())
		{
			IJ.error("Add Image...", "Image to be added cannnot be larger than\n\"" + imp.getTitle() + "\".");
			return;
		}

		ImageRoi roi2 = new ImageRoi(x, y, overlay.getProcessor());
		roi2.setZeroTransparent(transparent);
		roi2.setName(overlay.getShortTitle());
		if (opacity != 100)
			roi2.setOpacity(opacity / 100.0);

		Overlay overlayList = imp.getOverlay();
		if (overlayList == null || replace)
			overlayList = new Overlay();
		overlayList.add(roi2);
		imp.setOverlay(overlayList);
		Undo.setup(Undo.OVERLAY_ADDITION, imp);
	}
}
