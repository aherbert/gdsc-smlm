package gdsc.smlm.ij.plugins;

import gdsc.core.ij.Utils;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

/**
 * This plugin creates a mask image stack using an XY and XZ mask image
 */
public class DepthMask implements PlugIn
{
	private static final String TITLE = "Depth Mask";

	private static String titleXY = "";
	private static String titleXZ = "";
	private static String titleYZ = "";

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (!showDialog())
			return;

		createMask();
	}

	private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Create a mask stack using XY, XZ and YZ mask images");

		String[] maskList = Utils.getImageList(Utils.SINGLE);
		gd.addChoice("Mask_XY", maskList, titleXY);
		gd.addChoice("Mask_XZ", maskList, titleXZ);
		gd.addChoice("Mask_YZ", maskList, titleYZ);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		titleXY = gd.getNextChoice();
		titleXZ = gd.getNextChoice();
		titleYZ = gd.getNextChoice();

		return true;
	}

	private void createMask()
	{
		ImagePlus impXY = WindowManager.getImage(titleXY);
		ImagePlus impXZ = WindowManager.getImage(titleXZ);
		ImagePlus impYZ = WindowManager.getImage(titleYZ);
		if (impXY == null)
		{
			IJ.error(TITLE, "No XY mask");
			return;
		}
		if (impXZ == null)
		{
			IJ.error(TITLE, "No XZ mask");
			return;
		}
		if (impYZ == null)
		{
			IJ.error(TITLE, "No YZ mask");
			return;
		}
		if (impXY.getWidth() != impXZ.getWidth())
		{
			IJ.error(TITLE, "XY mask width does not match XZ mask width");
			return;
		}
		if (impXY.getHeight() != impYZ.getWidth())
		{
			IJ.error(TITLE, "XY mask height does not match YZ mask width");
			return;
		}
		if (impXZ.getHeight() != impYZ.getHeight())
		{
			IJ.error(TITLE, "XZ mask height does not match YZ mask height");
			return;
		}

		final int maxx = impXY.getWidth();
		final int maxy = impXY.getHeight();
		final int maxz = impXZ.getHeight();
		ImageStack stack = new ImageStack(maxx, maxy, maxz);
		byte[] maskXY = getMask(impXY);
		byte[] maskXZ = getMask(impXZ);
		byte[] maskYZ = getMask(impYZ);
		for (int z = 0; z < maxz; z++)
		{
			byte[] mask = maskXY.clone();

			//// Simple method
			//for (int y = 0, i = 0; y < maxy; y++, i++)
			//	for (int x = 0; x < maxx; x++, i++)
			//	{
			//		if (maskXZ[z * maxx + x] == 0)
			//			mask[i] = 0;
			//		else if (maskYZ[z * maxy + y] == 0)
			//			mask[i] = 0;
			//	}
			
			for (int x = 0, i = maxx * z; x < maxx; x++, i++)
			{
				if (maskXZ[i] == 0)
				{
					// Blank all the (x,y) for this X
					for (int y = 0, xy = x; y < maxy; y++, xy += maxx)
						mask[xy] = 0;
				}
			}

			for (int y = 0, i = maxy * z; y < maxy; y++, i++)
			{
				if (maskYZ[i] == 0)
				{
					// Blank all the (x,y) for this Y
					for (int x = 0, xy = y * maxx; x < maxx; x++, xy++)
						mask[xy] = 0;
				}
			}

			stack.setPixels(mask, z + 1);
		}
		Utils.display(TITLE, stack);
	}

	private byte[] getMask(ImagePlus impXY)
	{
		final byte[] mask = (byte[]) impXY.getProcessor().convertToByte(false).getPixels();
		// Make binary
		final byte ON = (byte) 255;
		for (int i = 0; i < mask.length; i++)
			if (mask[i] != 0)
				mask[i] = ON;
		return mask;
	}
}
