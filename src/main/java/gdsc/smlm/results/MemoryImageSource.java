package gdsc.smlm.results;

import java.awt.Rectangle;

import com.thoughtworks.xstream.XStream;

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
public class MemoryImageSource extends ImageSource
{
	private int counter;
	private float[][] data;
	private boolean freeMemoryOnClose;

	/**
	 * Create a new image source
	 */
	public MemoryImageSource(int width, int height, float[]... data)
	{
		super("Memory");
		if (width < 1)
			throw new IllegalArgumentException("Width must be above zero");
		if (height < 1)
			throw new IllegalArgumentException("Height must be above zero");
		if (data == null)
			throw new IllegalArgumentException("Image data must not be null");
		final int length = width * height;
		for (float[] f : data)
		{
			if (f == null)
				throw new IllegalArgumentException("Image data must not be null");
			if (f.length != length)
				throw new IllegalArgumentException("Image data length does not match width*height");
		}
		this.width = width;
		this.height = height;
		this.data = data;
		this.frames = data.length;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#openSource()
	 */
	@Override
	protected boolean openSource()
	{
		counter = 0;
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#close()
	 */
	@Override
	public void close()
	{
		// Free the memory
		if (freeMemoryOnClose)
			data = new float[0][0];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#nextFrame(java.awt.Rectangle)
	 */
	@Override
	protected float[] nextFrame(Rectangle bounds)
	{
		if (counter < data.length)
			return getFrame(++counter, bounds);
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#getFrame(int, java.awt.Rectangle)
	 */
	@Override
	protected float[] getFrame(int frame, Rectangle bounds)
	{
		if (frame > 0 && frame <= data.length)
		{
			if (checkBounds(bounds))
			{
				final float[] pixels = data[frame - 1];
				float[] pixels2 = new float[bounds.width * bounds.height];
				for (int ys = bounds.y; ys < bounds.y + bounds.height; ys++)
				{
					int offset1 = (ys - bounds.y) * bounds.width;
					int offset2 = ys * width + bounds.x;
					for (int xs = 0; xs < bounds.width; xs++)
						pixels2[offset1++] = pixels[offset2++];
				}
				return pixels2;
			}
			else
			{
				return data[frame - 1];
			}
		}
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#isValid(int)
	 */
	@Override
	public boolean isValid(int frame)
	{
		return (frame > 0 && frame <= data.length);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#init(com.thoughtworks.xstream.XStream)
	 */
	@Override
	public void init(XStream xs)
	{
		super.init(xs);
		xs.omitField(MemoryImageSource.class, "counter");
	}

	/**
	 * @return Set to true if freeing the memory on calling {@link #close()}
	 */
	public boolean isFreeMemoryOnClose()
	{
		return freeMemoryOnClose;
	}

	/**
	 * Set this to true to free the memory when the {@link #close()} method is called. No subsequent calls to
	 * {@link #open()} will be valid.
	 * 
	 * @param freeMemoryOnClose
	 *            Set to true to free the memory on calling {@link #close()}
	 */
	public void setFreeMemoryOnClose(boolean freeMemoryOnClose)
	{
		this.freeMemoryOnClose = freeMemoryOnClose;
	}
}
