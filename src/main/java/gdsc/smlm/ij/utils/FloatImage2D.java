package gdsc.smlm.ij.utils;

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
 * Store a 2D image in a single float array. Forms a base for 2D DHT transform using the JTransforms library.
 */
public class FloatImage2D extends Image2D
{
	protected float[] data;

	/**
	 * Instantiates a new 2D image.
	 *
	 * @param nc
	 *            the number of columns
	 * @param nr
	 *            the number of rows
	 * @throws IllegalArgumentException
	 *             If the combined dimensions is too large for an array
	 */
	public FloatImage2D(int nc, int nr) throws IllegalArgumentException
	{
		super(nc, nr);
	}

	/**
	 * Instantiates a new 2D image.
	 *
	 * @param image
	 *            the image
	 * @throws IllegalArgumentException
	 *             If the combined dimensions is too large for an array
	 */
	public FloatImage2D(ImageProcessor image) throws IllegalArgumentException
	{
		super(image);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.utils.Image2D#createData(int)
	 */
	@Override
	protected void createData(int size)
	{
		data = new float[size];
	}

	/**
	 * Instantiates a new 2D image.
	 *
	 * @param nc
	 *            the number of columns
	 * @param nr
	 *            the number of rows
	 * @param data
	 *            the data
	 * @throws IllegalArgumentException
	 *             If the data is not the correct length
	 */
	public FloatImage2D(int nc, int nr, float[] data) throws IllegalArgumentException
	{
		// Avoid constructor that calls createData(int)
		super(nc, nr, false);
		if (data == null || data.length != checkSize(nc, nr, true))
			throw new IllegalArgumentException("Data is not correct length");
		this.data = data;
	}

	/**
	 * Instantiates a new 2D image.
	 *
	 * @param nc
	 *            the number of columns
	 * @param nr
	 *            the number of rows
	 * @param data
	 *            the data
	 * @param dummy
	 *            the dummy flag
	 */
	protected FloatImage2D(int nc, int nr, float[] data, boolean dummy)
	{
		// No checks as this is used internally		
		super(nc, nr, false);
		this.data = data;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.utils.Image2D#copy()
	 */
	public FloatImage2D copy()
	{
		return new FloatImage2D(nc, nr, data.clone(), false);
	}

	/**
	 * Gets the data.
	 *
	 * @return the data
	 */
	public float[] getData()
	{
		return data;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.utils.Image2D#getDataLength()
	 */
	@Override
	public int getDataLength()
	{
		return data.length;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.utils.Image2D#crop(int, int, int, int)
	 */
	@Override
	public FloatImage2D crop(int x, int y, int w, int h) throws IllegalArgumentException
	{
		return crop(x, y, w, h, null);
	}

	/**
	 * Crop a sub-region of the data. The target dimensions must be positive.
	 *
	 * @param x
	 *            the x index
	 * @param y
	 *            the y index
	 * @param w
	 *            the width
	 * @param h
	 *            the height
	 * @param region
	 *            the cropped data (will be reused if the correct size)
	 * @return the cropped data
	 * @throws IllegalArgumentException
	 *             if the region is not within the data
	 */
	public FloatImage2D crop(int x, int y, int w, int h, float[] region) throws IllegalArgumentException
	{
		// Check the region range
		if (x < 0 || w < 1 || (long) x + w > nc || y < 0 || h < 1 || (long) y + h > nr)
			throw new IllegalArgumentException("Region not within the data");
		int size = h * w;
		if (region == null || region.length != size)
			region = new float[size];
		int base = y * nc + x;
		for (int r = 0, i = 0; r < h; r++)
		{
			System.arraycopy(data, base, region, i, w);
			base += nc;
			i += w;
		}
		return new FloatImage2D(w, h, region, false);
	}

	/**
	 * Crop a sub-region of the data. The target dimensions must be positive.
	 *
	 * @param image
	 *            the image
	 * @param x
	 *            the x index
	 * @param y
	 *            the y index
	 * @param w
	 *            the width
	 * @param h
	 *            the height
	 * @param region
	 *            the cropped data (will be reused if the correct size)
	 * @return the cropped data
	 * @throws IllegalArgumentException
	 *             if the region is not within the data
	 */
	public static FloatImage2D crop(ImageProcessor image, int x, int y, int w, int h, float[] region)
			throws IllegalArgumentException
	{
		int nc = image.getWidth();
		int nr = image.getHeight();

		// Check the region range
		if (x < 0 || w < 1 || (long) x + w > nc || y < 0 || h < 1 || (long) y + h > nr)
			throw new IllegalArgumentException("Region not within the data");
		int size = checkSize(w, h, true);
		if (region == null || region.length != size)
			region = new float[size];

		if (image.getBitDepth() != 32)
		{
			// Handle non-float data
			int base = y * nc + x;
			for (int r = 0, i = 0; r < h; r++)
			{
				for (int c = 0; c < w; c++)
				{
					region[i++] = image.getf(base + c);
				}
				base += nc;
			}
		}
		else
		{
			float[] data = (float[]) image.getPixels();
			int base = y * nc + x;
			for (int r = 0, i = 0; r < h; r++)
			{
				System.arraycopy(data, base, region, i, w);
				base += nc;
				i += w;
			}
		}
		return new FloatImage2D(w, h, region, false);
	}

	@Override
	public void insert(int x, int y, Image2D image) throws IllegalArgumentException
	{
		if (image instanceof FloatImage2D)
		{
			insert(x, y, (FloatImage2D) image);
		}
		else
		{
			super.insert(x, y, image);
		}
	}

	/**
	 * Insert a sub-region.
	 *
	 * @param x
	 *            the x position
	 * @param y
	 *            the y position
	 * @param image
	 *            the image
	 * @throws IllegalArgumentException
	 *             if the region is not within the data
	 */
	public void insert(int x, int y, FloatImage2D image) throws IllegalArgumentException
	{
		// Check the region range
		int w = image.getWidth();
		int h = image.getHeight();
		if (w < 1 || h < 1)
			return;
		if (x < 0 || (long) x + w > nc || y < 0 || (long) y + h > nr)
			throw new IllegalArgumentException("Region not within the data");
		float[] region = image.data;
		int base = y * nc + x;
		for (int r = 0, i = 0; r < h; r++)
		{
			System.arraycopy(region, i, data, base, w);
			base += nc;
			i += w;
		}
	}

	@Override
	protected void copyTo(int i, float[] buffer, int j, int size)
	{
		System.arraycopy(data, i, buffer, j, size);
	}

	@Override
	protected void copyFrom(float[] buffer, int j, int size, int i)
	{
		System.arraycopy(buffer, j, data, i, size);
	}

	@Override
	public double get(int i)
	{
		return data[i];
	}

	@Override
	public void set(int i, double value)
	{
		data[i] = (float) value;
	}

	@Override
	public float getf(int i)
	{
		return data[i];
	}

	@Override
	public void setf(int i, float value)
	{
		data[i] = value;
	}
}