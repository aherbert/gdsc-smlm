package gdsc.smlm.ij.plugins;

import java.awt.AWTEvent;

import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.Utils;
import gdsc.core.logging.Ticker;
import gdsc.smlm.filters.KernelFilter;
import gdsc.smlm.filters.ZeroKernelFilter;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Convolve an image with another image.
 */
public class ConvolveFilter implements ExtendedPlugInFilter, DialogListener
{
	private static final String TITLE = "Convolve Filter";
	private final int FLAGS = DOES_8G | DOES_16 | DOES_32 | KEEP_PREVIEW | PARALLELIZE_STACKS | CONVERT_TO_FLOAT;

	private static String title = "";
	private static int border = 0;
	private static boolean zero = false;

	// Ensure not null
	private Ticker ticker = Ticker.INSTANCE;

	private int lastId = 0;
	private boolean lastZero;
	private KernelFilter kf = null;
	private ImagePlus imp;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	public int setup(String arg, ImagePlus imp)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (imp == null)
		{
			IJ.noImage();
			return DONE;
		}
		return FLAGS;
	}

	public void run(ImageProcessor ip)
	{
		ticker.tick();
		kf.convolve((float[]) ip.getPixels(), ip.getWidth(), ip.getHeight(), border);
	}

	public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr)
	{
		// Get available kernels
		String[] names = Utils.getImageList(Utils.GREY_SCALE | Utils.SINGLE);
		if (names.length == 0)
		{
			IJ.error(TITLE, "No suitable kernel images");
			return DONE;
		}

		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Convolve an image using another image as the convolution kernel");

		gd.addChoice("Kernel_image", names, title);
		gd.addSlider("Border", 0, 10, border);
		gd.addCheckbox("Zero_outside_image", zero);

		gd.addDialogListener(this);
		gd.addPreviewCheckbox(pfr);

		gd.showDialog();

		if (gd.wasCanceled() || !dialogItemChanged(gd, null))
			return DONE;

		return IJ.setupDialog(imp, FLAGS);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.gui.DialogListener#dialogItemChanged(ij.gui.GenericDialog, java.awt.AWTEvent)
	 */
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
	{
		title = gd.getNextChoice();
		border = (int) gd.getNextNumber();
		zero = gd.getNextBoolean();

		imp = WindowManager.getImage(title);
		if (imp == null)
			return false;

		return true;
	}

	public void setNPasses(int nPasses)
	{
		// Create the kernel from the image
		if (kf == null || imp.getID() != lastId || zero != lastZero)
		{
			FloatProcessor fp = imp.getProcessor().toFloat(0, null);
			lastId = imp.getID();
			lastZero = zero;
			int kw = fp.getWidth();
			int kh = fp.getHeight();
			// Ensure odd size(to avoid exceptions)
			if ((kw & 1) != 1)
				kw++;
			if ((kh & 1) != 1)
				kh++;
			if (kw != fp.getWidth() || kh != fp.getHeight())
			{
				FloatProcessor fp2 = new FloatProcessor(kw, kh);
				fp2.insert(fp, 0, 0);
				fp = fp2;
			}
			float[] kernel = (float[]) fp.getPixels();
			kf = (zero) ? new ZeroKernelFilter(kernel, kw, kh) : new KernelFilter(kernel, kw, kh);
		}

		ticker = Ticker.create(new IJTrackProgress(), nPasses, true);
		ticker.start();
	}
}
