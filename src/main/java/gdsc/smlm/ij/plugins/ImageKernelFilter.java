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

import java.awt.AWTEvent;

import gdsc.core.ij.IJTrackProgress;
import gdsc.core.ij.Utils;
import gdsc.core.logging.Ticker;
import gdsc.smlm.filters.FHTFilter;
import gdsc.smlm.filters.FHTFilter.Operation;
import gdsc.smlm.filters.KernelFilter;
import gdsc.smlm.filters.ZeroKernelFilter;
import gdsc.smlm.ij.settings.SettingsManager;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.ExtendedGenericDialog;
import ij.gui.GenericDialog;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/**
 * Convolve an image with a kernel from another image.
 */
public class ImageKernelFilter implements ExtendedPlugInFilter, DialogListener
{
	private static final String TITLE = "Image Kernel Filter";
	private final int FLAGS = DOES_8G | DOES_16 | DOES_32 | KEEP_PREVIEW | PARALLELIZE_STACKS | CONVERT_TO_FLOAT |
			FINAL_PROCESSING;

	private static final String[] METHODS = { "Spatial domain", "FHT" };
	private static final int METHOD_SPATIAL = 0;
	private static final int METHOD_FHT = 1;
	private static final String[] FILTERS;
	static
	{
		FILTERS = SettingsManager.getNames((Object[]) Operation.values());
	}

	private static String title = "";
	private static int method = METHOD_FHT;
	private static int filter = Operation.CORRELATION.ordinal();
	private static int border = 0;
	private static boolean zero = false;

	// Ensure not null
	private Ticker ticker = Ticker.INSTANCE;

	private int lastId = 0;
	private int lastMethod = -1;
	private int lastFilter = -1;
	private boolean lastZero;
	private KernelFilter kf = null;
	private FHTFilter ff = null;
	private ImagePlus dataImp, kernelImp;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	@Override
	public int setup(String arg, ImagePlus imp)
	{
		if ("final".equals(arg))
		{
			imp.getProcessor().resetMinAndMax();
			imp.updateAndDraw();
			return DONE;
		}

		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (imp == null)
		{
			IJ.noImage();
			return DONE;
		}
		return FLAGS;
	}

	@Override
	public void run(ImageProcessor ip)
	{
		float[] data = (float[]) ip.getPixels();
		int w = ip.getWidth();
		int h = ip.getHeight();
		if (method == METHOD_SPATIAL)
		{
			kf.convolve(data, w, h, border);
		}
		else
		{
			// Use a clone for thread safety
			FHTFilter f = (ticker.getTotal() > 1l) ? ff.clone() : ff;
			f.filter(data, w, h, border);
		}
		if (ticker.getTotal() == 1l)
			ip.resetMinAndMax();
		ticker.tick();
	}

	@Override
	public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr)
	{
		// Get available kernels
		String[] names = Utils.getImageList(Utils.GREY_SCALE | Utils.SINGLE);
		if (names.length == 0)
		{
			IJ.error(TITLE, "No suitable kernel images");
			return DONE;
		}

		this.dataImp = imp;

		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Convolve an image using another image as the convolution kernel");

		gd.addChoice("Kernel_image", names, title);
		gd.addChoice("Method", METHODS, method);
		gd.addChoice("Filter", FILTERS, filter);
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
	@Override
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e)
	{
		title = gd.getNextChoice();
		method = gd.getNextChoiceIndex();
		filter = gd.getNextChoiceIndex();
		border = (int) gd.getNextNumber();
		zero = gd.getNextBoolean();

		kernelImp = WindowManager.getImage(title);
		if (kernelImp == null)
			return false;

		return true;
	}

	@Override
	public void setNPasses(int nPasses)
	{
		// Create the kernel from the image
		boolean build = kernelImp.getID() != lastId || method != lastMethod || filter != lastFilter;
		build = build || (method == METHOD_SPATIAL && kf == null);
		build = build || (method == METHOD_FHT && ff == null);
		if (build)
		{
			Operation operation = Operation.forOrdinal(filter);
			FloatProcessor fp = kernelImp.getProcessor().toFloat(0, null);
			if (method == METHOD_SPATIAL)
			{
				if (kf == null || kernelImp.getID() != lastId || zero != lastZero)
				{
					fp = KernelFilter.pad(fp);
					int kw = fp.getWidth();
					int kh = fp.getHeight();
					float[] kernel = (float[]) fp.getPixels();
					kf = (zero) ? new ZeroKernelFilter(kernel, kw, kh) : new KernelFilter(kernel, kw, kh);
				}
				switch (operation)
				{
					case CONVOLUTION:
						kf.setConvolution(true);
						break;
					case CORRELATION:
						kf.setConvolution(false);
						break;
					case DECONVOLUTION:
					default:
						Utils.log("Unsupported operation (%s), default to correlation", operation.getName());
						kf.setConvolution(false);
						break;
				}
			}
			else
			{
				if (ff == null || kernelImp.getID() != lastId)
				{
					int kw = fp.getWidth();
					int kh = fp.getHeight();
					float[] kernel = (float[]) fp.getPixels();
					ff = new FHTFilter(kernel, kw, kh);
					ff.initialiseKernel(dataImp.getWidth(), dataImp.getHeight());
				}
				ff.setOperation(operation);
			}
			lastId = kernelImp.getID();
			lastMethod = method;
			lastFilter = filter;
			lastZero = zero;
		}

		ticker = Ticker.create(new IJTrackProgress(), nPasses, nPasses != 1);
		ticker.start();
	}
}
