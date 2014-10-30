package gdsc.smlm.utils;

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

import java.awt.Rectangle;

/**
 * Contains methods for extracting parts of an image
 */
public class ImageExtractor
{
	private float[] data;
	private int width;
	private int height;

	/**
	 * Constructor
	 * 
	 * @param data
	 *            The image data
	 * @param width
	 *            The image width
	 * @param height
	 *            The image height
	 */
	public ImageExtractor(float[] data, int width, int height)
	{
		this.data = data;
		this.width = width;
		this.height = height;
	}

	/**
	 * Extract a region from the image.
	 * 
	 * @param regionBounds
	 *            The region to extract
	 * @return The image region (with dimensions specified in the dimensions array)
	 */
	public float[] crop(Rectangle regionBounds)
	{
		return crop(regionBounds, (float[]) null);
	}

	/**
	 * Extract a region from the image. The output array can be truncated
	 * using the {@link #truncate(float[], int)} method.
	 * 
	 * @param regionBounds
	 *            The region to extract
	 * @param region
	 *            A reusable buffer for the region
	 * @return The image region (with dimensions specified in the dimensions array)
	 */
	public float[] crop(Rectangle regionBounds, float[] region)
	{
		region = allocate(region, regionBounds.width * regionBounds.height);

		int offset1 = 0;
		for (int ys = regionBounds.y; ys < regionBounds.y + regionBounds.height; ys++)
		{
			int offset2 = ys * width + regionBounds.x;
			for (int xs = 0; xs < regionBounds.width; xs++)
				region[offset1++] = data[offset2++];
		}

		return region;
	}

	private static float[] allocate(float[] buffer, int size)
	{
		if (buffer == null || buffer.length < size)
			buffer = new float[size];
		return buffer;
	}

	/**
	 * Truncate the input data to the given length. Does nothing if the data is shorter or null.
	 * 
	 * @param data
	 * @param length
	 * @return The truncated data
	 */
	public static float[] truncate(float[] data, int length)
	{
		if (data != null && data.length > length)
		{
			float[] newData = new float[length];
			for (int i = length; i-- > 0;)
			{
				newData[i] = data[i];
			}
			return newData;
		}
		return data;
	}

	/**
	 * Extract a region from the image.
	 * 
	 * @param regionBounds
	 *            The region to extract
	 * @return The image region (with dimensions specified in the dimensions array)
	 */
	public double[] cropToDouble(Rectangle regionBounds)
	{
		return crop(regionBounds, (double[]) null);
	}
	
	/**
	 * Extract a region from the image. The output array can be truncated
	 * using the {@link #truncate(double[], int)} method.
	 * 
	 * @param regionBounds
	 *            The region to extract
	 * @param region
	 *            A reusable buffer for the region
	 * @return The image region (with dimensions specified in the dimensions array)
	 */
	public double[] crop(Rectangle regionBounds, double[] region)
	{
		region = allocate(region, regionBounds.width * regionBounds.height);

		int offset1 = 0;
		for (int ys = regionBounds.y; ys < regionBounds.y + regionBounds.height; ys++)
		{
			int offset2 = ys * width + regionBounds.x;
			for (int xs = 0; xs < regionBounds.width; xs++)
				region[offset1++] = data[offset2++];
		}

		return region;
	}

	private static double[] allocate(double[] buffer, int size)
	{
		if (buffer == null || buffer.length < size)
			buffer = new double[size];
		return buffer;
	}

	/**
	 * Truncate the input data to the given length. Does nothing if the data is shorter or null.
	 * 
	 * @param data
	 * @param length
	 * @return The truncated data
	 */
	public static double[] truncate(double[] data, int length)
	{
		if (data != null && data.length > length)
		{
			double[] newData = new double[length];
			for (int i = length; i-- > 0;)
			{
				newData[i] = data[i];
			}
			return newData;
		}
		return data;
	}

	/**
	 * Calculate a square region of size 2n+1 around the given coordinates. Respects the image boundaries and
	 * so may return a non-square region.
	 * 
	 * @param x
	 * @param y
	 * @param n
	 * @return The region
	 */
	public Rectangle getBoxRegionBounds(int x, int y, int n)
	{
		Rectangle r1 = new Rectangle(x - n, y - n, 2 * n + 1, 2 * n + 1);
		return r1.intersection(new Rectangle(0, 0, width, height));
	}
}
