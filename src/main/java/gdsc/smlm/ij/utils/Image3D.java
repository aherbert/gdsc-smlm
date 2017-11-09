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
 * Store a 3D image in a single array
 */
public class Image3D
{
	/**
	 * The largest array size for which a regular 1D Java array is used to store the
	 * data (2^30)
	 */
	public static final int maxSizeOf32bitArray = 1073741824;

	/** The number of slices (max z). */
	public final int ns;
	/** The number of rows (max y). */
	public final int nr;
	/** The number of columns (max x). */
	public final int nc;

	/** The number of rows multiplied by the number of columns */
	public final int nr_by_nc;

	protected final float[] data;

	/**
	 * Instantiates a new 3D image
	 *
	 * @param stack
	 *            the stack
	 * @throws IllegalArgumentException
	 *             If the combined dimensions is too large for an array
	 */
	public Image3D(ImageStack stack) throws IllegalArgumentException
	{
		ns = stack.getSize();
		nr = stack.getHeight();
		nc = stack.getWidth();

		long size = (long) ns * nr * nc;
		if (size > maxSizeOf32bitArray)
			throw new IllegalArgumentException("3D data too large");

		data = new float[(int) size];

		nr_by_nc = nr * nc;
		if (stack.getBitDepth() == 32)
		{
			for (int s = 0; s < ns; s++)
			{
				System.arraycopy(stack.getPixels(s + 1), 0, data, s * nr_by_nc, nr_by_nc);
			}
		}
		else
		{
			for (int s = 1, i = 0; s <= ns; s++)
			{
				ImageProcessor ip = stack.getProcessor(s);
				for (int j = 0; i < nr_by_nc; j++)
					data[i++] = ip.getf(j);
			}
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
	 *             If any dimension is less than 2, or if the data is not the correct length
	 */
	public Image3D(int nc, int nr, int ns, float[] data) throws IllegalArgumentException
	{
		long size = (long) ns * nr * nc;
		if (data == null || data.length != size)
			throw new IllegalArgumentException("Data is not correct length");
		this.ns = ns;
		this.nr = nr;
		this.nc = nc;
		nr_by_nc = nr * nc;
		this.data = data;
	}

	/**
	 * Instantiates a new 3D image
	 *
	 * @param ns
	 *            the number of slices
	 * @param nr
	 *            the number of rows
	 * @param nc
	 *            the number of columns
	 * @param nr_by_nc
	 *            the number of rows multiplied by the number of columns
	 * @param data
	 *            the data
	 */
	protected Image3D(int ns, int nr, int nc, int nr_by_nc, float[] data)
	{
		// No checks as this is used internally		
		this.ns = ns;
		this.nr = nr;
		this.nc = nc;
		this.nr_by_nc = nr_by_nc;
		this.data = data;
	}

	/**
	 * Return a copy of the 3D image.
	 *
	 * @return the copy
	 */
	public Image3D copy()
	{
		return new Image3D(ns, nr, nc, nr_by_nc, data.clone());
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

	/**
	 * Convert to an image stack.
	 *
	 * @return the image stack
	 */
	public ImageStack getImageStack()
	{
		ImageStack stack = new ImageStack(nc, nr);
		for (int s = 0; s < ns; s++)
		{
			float[] pixels = new float[nr_by_nc];
			System.arraycopy(data, s * nr_by_nc, pixels, 0, nr_by_nc);
			stack.addSlice(null, pixels);
		}
		return stack;
	}

	/**
	 * Gets the xyz components of the index.
	 *
	 * @param i
	 *            the index
	 * @return the xyz components
	 * @throws IllegalArgumentException
	 *             if the index is not within the data
	 */
	public int[] getXYZ(int i) throws IllegalArgumentException
	{
		if (i < 0 || i >= data.length)
			throw new IllegalArgumentException("Index in not in the correct range: 0 <= i < " + data.length);
		int[] xyz = new int[3];
		xyz[2] = i / nr_by_nc;
		int j = i % nr_by_nc;
		xyz[1] = j / nc;
		xyz[0] = j % nc;
		return xyz;
	}

	/**
	 * Gets the xyz components of the index.
	 *
	 * @param i
	 *            the index
	 * @param xyz
	 *            the xyz components (must be an array of at least length 3)
	 * @throws IllegalArgumentException
	 *             if the index is not within the data
	 */
	public void getXYZ(int i, int[] xyz) throws IllegalArgumentException
	{
		if (i < 0 || i >= data.length)
			throw new IllegalArgumentException("Index in not in the correct range: 0 <= i < " + data.length);
		xyz[2] = i / nr_by_nc;
		int j = i % nr_by_nc;
		xyz[1] = j / nc;
		xyz[0] = j % nc;
	}

	/**
	 * Crop a sub-region of the data.
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
	public Image3D crop(int x, int y, int z, int w, int h, int d, float[] region)
	{
		// Check the region range
		if (x < 0 || x + w >= nc || y < 0 || y + h >= nr || z < 0 || z + d >= ns)
			throw new IllegalArgumentException("Region not within the data");
		int size = d * h * w;
		if (region == null || region.length != size)
			region = new float[size];
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
		return new Image3D(w, h, d, w * d, region);
	}

	/**
	 * Crop a sub-region of the data.
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
	 * @return the cropped data
	 * @throws IllegalArgumentException
	 *             if the region is not within the data
	 */
	public ImageStack cropToStack(int x, int y, int z, int w, int h, int d)
	{
		// Check the region range
		if (x < 0 || x + w >= nc || y < 0 || y + h >= nr || z < 0 || z + d >= ns)
			throw new IllegalArgumentException("Region not within the data");
		int size = w * h;
		ImageStack stack = new ImageStack(w, h);
		for (int s = 0; s < d; s++, z++)
		{
			int base = z * nr_by_nc + y * nc + x;
			float[] region = new float[size];
			for (int r = 0, i = 0; r < h; r++)
			{
				System.arraycopy(data, base, region, i, w);
				base += nc;
				i += w;
			}
			stack.addSlice(null, region);
		}
		return stack;
	}
}