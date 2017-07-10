package gdsc.smlm.ij.plugins;

import gdsc.smlm.data.config.PSFProtos.ImagePSF;

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

import gdsc.smlm.ij.settings.ImagePSFHelper;
import gdsc.smlm.ij.utils.ImageConverter;
import gdsc.smlm.utils.XmlUtils;
import gdsc.core.ij.Utils;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.Blitter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.util.LinkedList;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

/**
 * Produces an average PSF image from multiple PSF images.
 * <p>
 * The input images must be a z-stack of a PSF. These can be produced using the PSFCreator plugin.
 */
public class PSFCombiner implements PlugIn
{
	private final static String TITLE = "PSF Combiner";
	private static int zDepth = 100;

	private List<String> titles = new LinkedList<String>();
	private List<PSF> input = new LinkedList<PSF>();

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		// Build a list of suitable images
		titles = createImageList();

		if (titles.isEmpty())
		{
			IJ.error(TITLE, "No suitable PSF images");
			return;
		}

		try
		{
			while (selectNextImage())
				;
		}
		catch (Exception e)
		{
			IJ.error(TITLE, e.getMessage());
			return;
		}

		if (input.isEmpty())
		{
			return;
		}

		if (input.size() < 2)
		{
			IJ.error(TITLE, "Require at least 2 PSF images to combine");
			return;
		}

		GenericDialog gd = new GenericDialog(TITLE);
		gd.addMessage("Set the maximum z-depth +/- from the PSF centre");
		gd.addSlider("Z-depth", 20, 200, zDepth);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		zDepth = Math.abs((int) gd.getNextNumber());

		for (PSF psf : input)
			psf.crop(zDepth);

		combineImages();
	}

	public static List<String> createImageList()
	{
		List<String> titles = new LinkedList<String>();
		int[] ids = WindowManager.getIDList();
		if (ids != null)
		{
			for (int id : ids)
			{
				ImagePlus imp = WindowManager.getImage(id);
				if (imp != null)
				{
					// Image must be greyscale
					if (imp.getType() == ImagePlus.GRAY8 || imp.getType() == ImagePlus.GRAY16 ||
							imp.getType() == ImagePlus.GRAY32)
					{
						// Image must be square and a stack of a single channel
						if (imp.getWidth() == imp.getHeight() && imp.getNChannels() == 1)
						{
							// Check if these are PSF images created by the SMLM plugins
							if (containsPSF(imp))
								titles.add(imp.getTitle());
						}
					}
				}
			}
		}
		return titles;
	}

	private static boolean containsPSF(ImagePlus imp)
	{
		Object info = imp.getProperty("Info");
		if (info != null)
		{
			return ImagePSFHelper.fromString(info.toString()) != null;
		}
		return false;
	}

	private boolean selectNextImage()
	{
		// Show a dialog allowing the user to select an input image
		if (titles.isEmpty())
			return false;

		GenericDialog gd = new GenericDialog(TITLE);
		gd.addMessage("Select the next input PSF image.\n(Each PSF must have the nm/pixel scale)");
		int n = (input.size() + 1);

		// If in macro mode then we must just use the String input field to allow the macro
		// IJ to return the field values from the macro arguments. Using a Choice input
		// will always return a field value.

		if (IJ.isMacro())
			gd.addStringField("PSF_" + n, "");
		else
			gd.addChoice("PSF_" + n, titles.toArray(new String[titles.size()]), "");

		gd.addMessage("Cancel to finish");
		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		String title;
		if (IJ.isMacro())
			title = gd.getNextString();
		else
			title = gd.getNextChoice();

		// Check the image exists. If not then exit. This is mainly relevant for Macro mode since
		// the
		// loop will continue otherwise since the titles list is not empty.
		ImagePlus imp = WindowManager.getImage(title);
		if (imp == null)
			return false;

		titles.remove(title);
		input.add(new PSF(title));
		return true;
	}

	private void combineImages()
	{
		double nmPerPixel = getNmPerPixel();
		if (nmPerPixel <= 0)
			return;
		double nmPerSlice = getNmPerSlice();
		if (nmPerPixel <= 0)
			return;

		// Find the lowest start point
		int min = 0;
		for (PSF psf : input)
		{
			if (min > psf.start)
				min = psf.start;
		}

		// Shift all stacks and find the dimensions
		final int shift = -min;
		int max = 0, size = 0;
		int totalImages = 0;
		for (PSF psf : input)
		{
			psf.start += shift;
			totalImages += psf.psfSettings.getImageCount();
			if (max < psf.getEnd())
				max = psf.getEnd();
			if (size < psf.getSize())
				size = psf.getSize();
		}

		// Create a stack to hold all the images
		ImageStack stack = new ImageStack(size, size, max);
		for (int n = 1; n <= max; n++)
			stack.setPixels(new float[size * size], n);

		// Insert all the PSFs
		IJ.showStatus("Creating combined image ...");
		int imageNo = 0;
		double fraction = 1.0 / input.size();
		for (PSF psf : input)
		{
			double progress = imageNo * fraction;
			ImageStack psfStack = psf.psfStack;
			final int offsetXY = (psf.getSize() - size) / 2;
			final int offsetZ = psf.start;
			final int w = psf.getSize();
			final double weight = (1.0 * psf.psfSettings.getImageCount()) / totalImages;
			final double increment = fraction / psfStack.getSize();
			for (int n = 1; n <= psfStack.getSize(); n++)
			{
				IJ.showProgress(progress += increment);

				// Get the data and adjust using the weight
				float[] psfData = ImageConverter.getData(psfStack.getProcessor(n));
				for (int i = 0; i < psfData.length; i++)
					psfData[i] *= weight;

				// Insert into the combined PSF
				ImageProcessor ip = stack.getProcessor(n + offsetZ);
				ip.copyBits(new FloatProcessor(w, w, psfData, null), offsetXY, offsetXY, Blitter.ADD);
			}
			imageNo++;
		}

		// IJ.showStatus("Normalising ...");
		// PSFCreator.normalise(stack, 1 + shift);

		IJ.showProgress(1);
		IJ.showStatus("");

		ImagePlus imp = Utils.display("Combined PSF", stack);
		imp.setSlice(1 + shift);
		imp.resetDisplayRange();
		imp.updateAndDraw();

		final double fwhm = getFWHM();
		imp.setProperty("Info", ImagePSFHelper
				.toString(ImagePSFHelper.create(imp.getSlice(), nmPerPixel, nmPerSlice, totalImages, fwhm)));

		Utils.log("%s : z-centre = %d, nm/Pixel = %s, nm/Slice = %s, %d images, FWHM = %s\n", imp.getTitle(),
				imp.getSlice(), Utils.rounded(nmPerPixel), Utils.rounded(nmPerSlice), totalImages, Utils.rounded(fwhm));
	}

	private double getNmPerPixel()
	{
		final double nmPerPixel = input.get(0).psfSettings.getPixelSize();
		for (PSF psf : input)
			if (psf.psfSettings.getPixelSize() != nmPerPixel)
			{
				IJ.error(TITLE, "Different pixel size resolutions for the input PSFs");
				return -1;
			}
		return nmPerPixel;
	}

	private double getNmPerSlice()
	{
		final double nmPerSlice = input.get(0).psfSettings.getPixelDepth();
		for (PSF psf : input)
			if (psf.psfSettings.getPixelDepth() != nmPerSlice)
			{
				IJ.error(TITLE, "Different pixel depth resolutions for the input PSFs");
				return -1;
			}
		return nmPerSlice;
	}

	private double getFWHM()
	{
		double fwhm = 0;
		for (PSF psf : input)
		{
			fwhm += psf.psfSettings.getFwhm();
		}
		return fwhm / input.size();
	}

	private class PSF
	{
		private ImagePlus imp;
		ImagePSF psfSettings;
		int start;
		ImageStack psfStack;

		public PSF(String title)
		{
			imp = WindowManager.getImage(title);
			if (imp == null)
				throw new RuntimeException("No image with title: " + title);

			Object o = XmlUtils.fromXML(imp.getProperty("Info").toString());
			if (!(o != null && o instanceof ImagePSF))
				throw new RuntimeException("Unknown PSF settings for image: " + title);
			this.psfSettings = (ImagePSF) o;
			int zCentre = psfSettings.getCentreImage();
			if (zCentre < 1 || zCentre > imp.getStackSize())
				throw new RuntimeException("z-centre must be within the stack size: " + imp.getStackSize());
			start = 1 - zCentre;
			psfStack = imp.getImageStack();
		}

		/**
		 * Remove frames above and below the centre to the specified depth
		 * 
		 * @param zDepth
		 */
		public void crop(int zDepth)
		{
			int minZ = FastMath.max(1, psfSettings.getCentreImage() - zDepth);
			int maxZ = FastMath.min(psfStack.getSize(), psfSettings.getCentreImage() + zDepth);
			psfStack = psfStack.crop(0, 0, minZ, psfStack.getWidth(), psfStack.getHeight(), maxZ - minZ + 1);

			// Update range
			start += minZ - 1;
		}

		public int getEnd()
		{
			return start + psfStack.getSize();
		}

		public int getSize()
		{
			return imp.getWidth();
		}
	}
}
