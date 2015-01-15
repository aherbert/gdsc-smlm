package gdsc.smlm.ij.plugins;

import ij.IJ;
import ij.ImagePlus;
import ij.Undo;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.ImageRoi;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.plugin.PlugIn;

import java.awt.Rectangle;

/**
 * This plugin is extracted from ij.plugins.OverlayCommands to allow an image to be added with a transparent
 * background.
 */
public class OverlayImage implements PlugIn
{
	private static int opacity = 100;
	private static boolean transparent = true;

	public void run(String arg)
	{
		addImage(false);
	}

	/**
	 * Adapted from ij.plugins.OverlayCommands#addImage(boolean) with the additional option for setting the zero pixels
	 * to transparent.
	 * 
	 * @param createImageRoi
	 */
	void addImage(boolean createImageRoi)
	{
		ImagePlus imp = IJ.getImage();
		int[] wList = WindowManager.getIDList();
		if (wList == null || wList.length < 2)
		{
			IJ.error("Add Image...", "The command requires at least two open images.");
			return;
		}
		String[] titles = new String[wList.length];
		for (int i = 0; i < wList.length; i++)
		{
			ImagePlus imp2 = WindowManager.getImage(wList[i]);
			titles[i] = imp2 != null ? imp2.getTitle() : "";
		}
		int x = 0, y = 0;
		Roi roi = imp.getRoi();
		if (roi != null && roi.isArea())
		{
			Rectangle r = roi.getBounds();
			x = r.x;
			y = r.y;
		}
		int index = 0;
		if (wList.length == 2)
		{
			ImagePlus i1 = WindowManager.getImage(wList[0]);
			ImagePlus i2 = WindowManager.getImage(wList[1]);
			if (i2.getWidth() < i1.getWidth() && i2.getHeight() < i1.getHeight())
				index = 1;
		}
		else if (imp.getID() == wList[0])
			index = 1;

		String title = createImageRoi ? "Create Image ROI" : "Add Image...";
		GenericDialog gd = new GenericDialog(title);
		if (createImageRoi)
			gd.addChoice("Image:", titles, titles[index]);
		else
		{
			gd.addChoice("Image to add:", titles, titles[index]);
			gd.addNumericField("X location:", x, 0);
			gd.addNumericField("Y location:", y, 0);
			gd.addCheckbox("Transparent background", transparent);
		}
		gd.addNumericField("Opacity (0-100%):", opacity, 0);
		//gd.addCheckbox("Create image selection", createImageRoi);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		index = gd.getNextChoiceIndex();
		if (!createImageRoi)
		{
			x = (int) gd.getNextNumber();
			y = (int) gd.getNextNumber();
			transparent = gd.getNextBoolean();
		}
		opacity = (int) gd.getNextNumber();
		//createImageRoi = gd.getNextBoolean();
		ImagePlus overlay = WindowManager.getImage(wList[index]);
		if (wList.length == 2)
		{
			ImagePlus i1 = WindowManager.getImage(wList[0]);
			ImagePlus i2 = WindowManager.getImage(wList[1]);
			if (i2.getWidth() < i1.getWidth() && i2.getHeight() < i1.getHeight())
			{
				imp = i1;
				overlay = i2;
			}
		}
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
		if (createImageRoi && x == 0 && y == 0)
		{
			x = imp.getWidth() / 2 - overlay.getWidth() / 2;
			y = imp.getHeight() / 2 - overlay.getHeight() / 2;
		}
		ImageRoi roi2 = new ImageRoi(x, y, overlay.getProcessor());
		roi2.setZeroTransparent(transparent);
		roi = roi2;
		roi.setName(overlay.getShortTitle());
		if (opacity != 100)
			((ImageRoi) roi).setOpacity(opacity / 100.0);
		if (createImageRoi)
			imp.setRoi(roi);
		else
		{
			Overlay overlayList = imp.getOverlay();
			if (overlayList == null)
				overlayList = new Overlay();
			overlayList.add(roi);
			imp.setOverlay(overlayList);
			Undo.setup(Undo.OVERLAY_ADDITION, imp);
		}
	}

}
