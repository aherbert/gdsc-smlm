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

import org.apache.commons.math3.util.FastMath;

import gdsc.core.ij.Utils;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.ZProjector;
import ij.plugin.filter.GaussianBlur;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

/**
 * Produces a background intensity image and a mask from a sample image.
 * <p>
 * The input image should be representative of the super-resolution imaging conditions and so will produce suitable
 * input for the Create Data plugin to create realistic images.
 */
public class ImageBackground implements PlugInFilter
{
	private final static String TITLE = "Image Background";

	private static float bias = 500;
	private static double sigma = 2;

	private int flags = DOES_16 | DOES_8G | DOES_32 | NO_CHANGES;
	private ImagePlus imp;

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	@Override
	public int setup(String arg, ImagePlus imp)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (imp == null)
		{
			IJ.noImage();
			return DONE;
		}

		this.imp = imp;

		return showDialog();
	}

	private int showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Creates a background and mask image from a sample input stack\nusing a median projection");

		gd.addNumericField("Bias", bias, 0);
		gd.addSlider("Blur", 0, 20, sigma);

		gd.showDialog();

		if (gd.wasCanceled())
			return DONE;

		bias = (float) gd.getNextNumber();
		sigma = gd.getNextNumber();

		// Check arguments
		try
		{
			Parameters.isPositive("Bias", bias);
		}
		catch (IllegalArgumentException e)
		{
			IJ.error(TITLE, e.getMessage());
			return DONE;
		}

		return flags;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	@Override
	public void run(ImageProcessor ip)
	{
		ImageProcessor median = getProjection();
		//Utils.display("Median", median);

		ImageProcessor background = applyBlur(median);
		subtractBias(background);

		Utils.display("Background", background);

		// Q. Is there a better way to do the thresholding for foreground pixels. 
		// Ideally we want to outline cell shapes. 
		ImageProcessor mask = median.convertToByte(true);
		mask.autoThreshold();

		Utils.display("Mask", mask);
	}

	private ImageProcessor getProjection()
	{
		// Get median intensity projection
		ZProjector p = new ZProjector(imp);
		p.setMethod(ZProjector.MEDIAN_METHOD);
		p.doProjection();
		ImageProcessor median = p.getProjection().getProcessor();
		return median;
	}

	private ImageProcessor applyBlur(ImageProcessor median)
	{
		ImageProcessor blur = median;
		if (sigma > 0)
		{
			blur = median.duplicate();
			GaussianBlur gb = new GaussianBlur();
			gb.blurGaussian(blur, sigma, sigma, 0.0002);
		}
		return blur;
	}

	private void subtractBias(ImageProcessor background)
	{
		float[] data = (float[]) background.getPixels();
		for (int i = 0; i < data.length; i++)
			data[i] = FastMath.max(0f, data[i] - bias);
		background.resetMinAndMax();
	}
}
