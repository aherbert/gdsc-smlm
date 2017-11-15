package gdsc.smlm.ij.utils;

import ij.ImageStack;
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
 * Store a 3D image in a single double array
 */
public class DoubleImage3D extends Image3D
{
	protected double[] data;

	/**
	 * Instantiates a new 3D image.
	 *
	 * @param nc
	 *            the number of columns
	 * @param nr
	 *            the number of rows
	 * @param ns
	 *            the number of slices
	 * @throws IllegalArgumentException
	 *             If the combined dimensions is too large for an array
	 */
	public DoubleImage3D(int nc, int nr, int ns) throws IllegalArgumentException
	{
		super(nc, nr, ns);
	}

	/**
	 * Instantiates a new 3D image
	 *
	 * @param stack
	 *            the stack
	 * @throws IllegalArgumentException
	 *             If the combined dimensions is too large for an array
	 */
	public DoubleImage3D(ImageStack stack) throws IllegalArgumentException
	{
		super(stack);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.utils.Image3D#createData(int)
	 */
	@Override
	protected void createData(int size)
	{
		data = new double[size];
	}

	/**
	 * Instantiates a new 3D image.
	 *
	 * @param stack
	 *            the stack
	 * @param data
	 *            the data
	 */
	private DoubleImage3D(ImageStack stack, double[] data)
	{
		super(stack.getWidth(), stack.getHeight(), stack.getSize(), stack.getWidth() * stack.getHeight());

		// This is used internally so the data is the correct length
		this.data = data;

		for (int s = 1, i = 0; s <= ns; s++)
		{
			ImageProcessor ip = stack.getProcessor(s);
			for (int j = 0; i < nr_by_nc; j++)
				data[i++] = ip.getf(j);
		}
	}

	/**
	 * Instantiates a new 3D image.
	 *
	 * @param nc
	 *            the number of columns
	 * @param nr
	 *            the number of rows
	 * @param ns
	 *            the number of slices
	 * @param data
	 *            the data
	 * @throws IllegalArgumentException
	 *             If the data is not the correct length
	 */
	public DoubleImage3D(int nc, int nr, int ns, double[] data) throws IllegalArgumentException
	{
		// Avoid constructor that calls createData(int)
		super(nc, nr, ns, nr * nc);
		if (data == null || data.length != checkSize(nc, nr, ns, true))
			throw new IllegalArgumentException("Data is not correct length");
		this.data = data;
	}

	/**
	 * Instantiates a new 3D image.
	 *
	 * @param nc
	 *            the number of columns
	 * @param nr
	 *            the number of rows
	 * @param ns
	 *            the number of slices
	 * @param nr_by_nc
	 *            the number of rows multiplied by the number of columns
	 * @param data
	 *            the data
	 */
	protected DoubleImage3D(int nc, int nr, int ns, int nr_by_nc, double[] data)
	{
		// No checks as this is used internally		
		super(nc, nr, ns, nr_by_nc);
		this.data = data;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.utils.Image3D#copy()
	 */
	public DoubleImage3D copy()
	{
		return new DoubleImage3D(nc, nr, ns, nr_by_nc, data.clone());
	}

	/**
	 * Gets the data.
	 *
	 * @return the data
	 */
	public double[] getData()
	{
		return data;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.utils.Image3D#getDataLength()
	 */
	@Override
	public int getDataLength()
	{
		return data.length;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.utils.Image3D#crop(int, int, int, int, int, int)
	 */
	@Override
	public DoubleImage3D crop(int x, int y, int z, int w, int h, int d) throws IllegalArgumentException
	{
		return crop(x, y, z, w, h, d, null);
	}

	/**
	 * Crop a sub-region of the data. The target dimensions must be positive.
	 *
	 * @param x
	 *            the x index
	 * @param y
	 *            the y index
	 * @param z
	 *            the z index
	 * @param w
	 *            the width
	 * @param h
	 *            the height
	 * @param d
	 *            the depth
	 * @param region
	 *            the cropped data (will be reused if the correct size)
	 * @return the cropped data
	 * @throws IllegalArgumentException
	 *             if the region is not within the data
	 */
	public DoubleImage3D crop(int x, int y, int z, int w, int h, int d, double[] region) throws IllegalArgumentException
	{
		// Check the region range
		if (x < 0 || w < 1 || (long) x + w > nc || y < 0 || h < 1 || (long) y + h > nr || z < 0 || d < 1 ||
				(long) z + d > ns)
			throw new IllegalArgumentException("Region not within the data");
		int size = d * h * w;
		if (region == null || region.length != size)
			region = new double[size];
		for (int s = 0, i = 0; s < d; s++, z++)
		{
			int base = z * nr_by_nc + y * nc + x;
			for (int r = 0; r < h; r++)
			{
				System.arraycopy(data, base, region, i, w);
				base += nc;
				i += w;
			}
		}
		return new DoubleImage3D(w, h, d, w * h, region);
	}

	/**
	 * Crop a sub-region of the data. The target dimensions must be positive.
	 *
	 * @param stack
	 *            the stack
	 * @param x
	 *            the x index
	 * @param y
	 *            the y index
	 * @param z
	 *            the z index
	 * @param w
	 *            the width
	 * @param h
	 *            the height
	 * @param d
	 *            the depth
	 * @param region
	 *            the cropped data (will be reused if the correct size)
	 * @return the cropped data
	 * @throws IllegalArgumentException
	 *             if the region is not within the data
	 */
	public static DoubleImage3D crop(ImageStack stack, int x, int y, int z, int w, int h, int d, double[] region)
			throws IllegalArgumentException
	{
		int nc = stack.getWidth();
		int nr = stack.getHeight();
		int ns = stack.getSize();

		// Check the region range
		if (x < 0 || w < 1 || (long) x + w > nc || y < 0 || h < 1 || (long) y + h > nr || z < 0 || d < 1 ||
				(long) z + d > ns)
			throw new IllegalArgumentException("Region not within the data");
		int size = checkSize(w, h, d, true);
		if (region == null || region.length != size)
			region = new double[size];
		for (int s = 0, i = 0; s < d; s++, z++)
		{
			ImageProcessor ip = stack.getProcessor(1 + z);
			int base = y * nc + x;
			for (int r = 0; r < h; r++)
			{
				for (int c = 0; c < w; c++)
				{
					region[i++] = ip.getf(base + c);
				}
				base += nc;
			}
		}
		return new DoubleImage3D(w, h, d, w * h, region);
	}

	@Override
	public void insert(int x, int y, int z, Image3D image) throws IllegalArgumentException
	{
		if (image instanceof DoubleImage3D)
		{
			insert(x, y, z, (DoubleImage3D) image);
		}
		else
		{
			super.insert(x, y, z, image);
		}
	}

	/**
	 * Insert a sub-region.
	 *
	 * @param x
	 *            the x position
	 * @param y
	 *            the y position
	 * @param z
	 *            the z position
	 * @param image
	 *            the image
	 * @throws IllegalArgumentException
	 *             if the region is not within the data
	 */
	public void insert(int x, int y, int z, DoubleImage3D image) throws IllegalArgumentException
	{
		// Check the region range
		int w = image.getWidth();
		int h = image.getHeight();
		int d = image.getSize();
		if (w < 1 || h < 1 || d < 1)
			return;
		if (x < 0 || (long) x + w > nc || y < 0 || (long) y + h > nr || z < 0 || (long) z + d > ns)
			throw new IllegalArgumentException("Region not within the data");
		double[] region = image.data;
		for (int s = 0, i = 0; s < d; s++, z++)
		{
			int base = z * nr_by_nc + y * nc + x;
			for (int r = 0; r < h; r++)
			{
				System.arraycopy(region, i, data, base, w);
				base += nc;
				i += w;
			}
		}
	}

	@Override
	protected void copyTo(int i, float[] buffer, int j, int size)
	{
		while (size-- > 0)
			buffer[j++] = (float) data[i++];
	}

	@Override
	protected void copyFrom(float[] buffer, int j, int size, int i)
	{
		while (size-- > 0)
			data[i++] = buffer[j++];
	}

	@Override
	public double get(int i)
	{
		return data[i];
	}

	@Override
	public void set(int i, double value)
	{
		data[i] = value;
	}

	@Override
	public float getf(int i)
	{
		return (float) data[i];
	}

	@Override
	public void setf(int i, float value)
	{
		data[i] = value;
	}
}