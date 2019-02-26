package gdsc.smlm.ij.plugins;

import uk.ac.sussex.gdsc.core.ij.ImageJUtils; import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils; import uk.ac.sussex.gdsc.core.utils.TextUtils; import uk.ac.sussex.gdsc.core.utils.MathUtils;
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

		gd.addMessage("Create a mask stack using XY and XZ mask images");

		String[] maskList = ImageJUtils.getImageList(ImageJUtils.SINGLE);
		gd.addChoice("Mask_XY", maskList, titleXY);
		gd.addChoice("Mask_XZ", maskList, titleXZ);

		gd.showDialog();

		if (gd.wasCanceled())
			return false;

		titleXY = gd.getNextChoice();
		titleXZ = gd.getNextChoice();

		return true;
	}

	private void createMask()
	{
		ImagePlus impXY = WindowManager.getImage(titleXY);
		ImagePlus impXZ = WindowManager.getImage(titleXZ);
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
		if (impXY.getWidth() != impXZ.getWidth())
		{
			IJ.error(TITLE, "XY mask width does not match XZ mask width");
			return;
		}

		final int maxx = impXY.getWidth();
		final int maxy = impXY.getHeight();
		final int maxz = impXZ.getHeight();
		ImageStack stack = new ImageStack(maxx, maxy, maxz);
		byte[] maskXY = getMask(impXY);
		byte[] maskXZ = getMask(impXZ);
		byte[] strip = new byte[maxx];
		for (int z = 0, p = 0; z < maxz; z++, p += maxx)
		{
			byte[] mask = maskXY.clone();
			System.arraycopy(maskXZ, p, strip, 0, maxx);
			for (int y = 0, i = 0; y < maxy; y++)
			{
				for (int x = 0; x < maxx; x++, i++)
				{
					if (strip[x] == 0)
						mask[i] = 0;
				}
			}
			stack.setPixels(mask, z + 1);
		}
		ImageJUtils.display(TITLE, stack);
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
