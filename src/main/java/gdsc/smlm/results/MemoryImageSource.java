package gdsc.smlm.results;

import com.thoughtworks.xstream.annotations.XStreamOmitField;

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
 * Represent a named results source providing data from memory.
 */
public class MemoryImageSource extends ImageSource
{
	@XStreamOmitField
	private int counter;
	private float[][] data;
	private boolean freeMemoryOnClose;

	/**
	 * Create a new image source
	 */
	public MemoryImageSource(int width, int height, float[]... data)
	{
		this(0, 0, width, height, data);
	}

	/**
	 * Create a new image source
	 */
	public MemoryImageSource(int xOrigin, int yOrigin, int width, int height, float[]... data)
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
		this.xOrigin = xOrigin;
		this.yOrigin = yOrigin;
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
	protected float[] nextRawFrame()
	{
		if (counter < data.length)
			return getRawFrame(++counter);
		return null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.ImageSource#getFrame(int, java.awt.Rectangle)
	 */
	@Override
	protected float[] getRawFrame(int frame)
	{
		if (frame > 0 && frame <= data.length)
		{
			return data[frame - 1];
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
