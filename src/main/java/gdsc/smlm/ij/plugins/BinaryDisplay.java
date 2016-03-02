package gdsc.smlm.ij.plugins;

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

import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/**
 * Set the display values of an image to render it as a binary mask, all non-zero values are white.
 */
public class BinaryDisplay implements PlugInFilter
{
	private ImagePlus imp;

	/* (non-Javadoc)
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	public int setup(String arg, ImagePlus imp)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);
		
		if (imp == null)
			return DONE;
		
		if (arg.equals("reset"))
		{
			ImageProcessor ip = imp.getProcessor();
			ip.reset();
			imp.setProcessor(ip);
			imp.resetDisplayRange();
			imp.updateAndDraw();
			return DONE;
		}
		
		this.imp = imp;
		return DOES_ALL;
	}

	/* (non-Javadoc)
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	public void run(ImageProcessor ip)
	{
//		float min = Float.POSITIVE_INFINITY;
//		for (int i=0; i<ip.getPixelCount(); i++)
//		{
//			final float value = ip.getf(i);
//			if (value == 0)
//				continue;
//			if (value < min)
//				min = value; 
//		}
//		ip.setMinAndMax(0, min);
//		imp.updateAndDraw();

		FloatProcessor fp = new FloatProcessor(ip.getWidth(), ip.getHeight());
		float[] data = (float[])fp.getPixels();
		for (int i=0; i<ip.getPixelCount(); i++)
		{
			final float value = ip.getf(i);
			if (value == 0)
				continue;
			data[i] = 1;
		}
		
		ip.snapshot();
		ip.setPixels(0, fp);
		ip.setMinAndMax(0, 1);
		imp.updateAndDraw();
	}
}
