package gdsc.smlm.results;

import java.awt.Rectangle;

import com.thoughtworks.xstream.XStream;

import gdsc.smlm.ij.utils.ImageConverter;
import gdsc.smlm.results.ImageSource;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.io.FileInfo;
import ij.process.ImageProcessor;

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

/**
 * Represent a named results source. Does not support data provision.
 */
public class NullSource extends ImageSource
{
	/**
	 * Create a new image source
	 * 
	 * @param name
	 */
	public NullSource(String name)
	{
		super(name);
	}

	@Override
	public boolean open()
	{
		return false;
	}

	@Override
	public float[] next(Rectangle bounds)
	{
		return null;
	}

	@Override
	public float[] get(int frame, Rectangle bounds)
	{
		return null;
	}
}
