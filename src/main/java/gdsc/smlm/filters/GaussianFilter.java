package gdsc.smlm.filters;

import org.apache.commons.math3.util.FastMath;

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
 * Computes a Gaussian convolution in the spatial domain for each point within the array.
 * <p>
 * Note: Due to lack of small dimension checking the routines will fail if maxx or maxy are less than 2. All routines
 * are OK for 3x3 images and larger.
 */
public class GaussianFilter implements Cloneable
{
	private float[] floatDataBuffer = null;
	
	private float floatGetAverage(float sum, final float divisor)
	{
		return sum * divisor;
	}

	/**
	 * Compute the block average within a 2n+1 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 */
	public void blockAverageInternal(float[] data, final int maxx, final int maxy, final int n)
	{
		if (n == 1)
			blockAverage3x3Internal(data, maxx, maxy);
		else
			blockAverageNxNInternal(data, maxx, maxy, n);
	}

	/**
	 * Compute the block average within a 2n+1 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 */
	public void blockAverageNxNInternal(float[] data, final int maxx, final int maxy, final int n)
	{
		float[] newData = floatBuffer(data.length);

		// Boundary control
		final int xwidth = FastMath.min(n, maxx - 1);
		final int ywidth = FastMath.min(n, maxy - 1);

		int[] offset = new int[(2 * xwidth + 1) * (2 * ywidth + 1) - 1];
		for (int y = -ywidth, d = 0; y <= ywidth; y++)
			for (int x = -xwidth; x <= xwidth; x++)
				if (x != 0 || y != 0)
				{
					offset[d] = maxx * y + x;
					d++;
				}

		final float divisor = (float) (1.0 / ((2 * n + 1) * (2 * n + 1)));

		for (int y = n; y < maxy - n; y++)
		{
			int index = y * maxx + n;
			for (int x = n; x < maxx - n; x++, index++)
			{
				float sum = data[index];

				// Sweep neighbourhood - 
				// No check for boundaries as this should be an internal sweep.
				for (int offset_d : offset)
				{
					sum += data[index + offset_d];
				}

				newData[index] = floatGetAverage(sum, divisor);
			}
		}

		// Copy back
		for (int y = n; y < maxy - n; y++)
		{
			int index = y * maxx + n;
			for (int x = n; x < maxx - n; x++, index++)
			{
				data[index] = newData[index];
			}
		}
	}

	/**
	 * Compute the block average within a 3x3 size block around each point.
	 * Only pixels with a full block are processed. Pixels within border regions
	 * are unchanged.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 */
	public void blockAverage3x3Internal(float[] data, final int maxx, final int maxy)
	{
		float[] newData = floatBuffer(data.length);

		final float divisor = (float) (1.0 / 9);

		for (int y = 1; y < maxy - 1; y++)
		{
			int index0 = (y - 1) * maxx + 1;
			int index1 = y * maxx + 1;
			int index2 = (y + 1) * maxx + 1;
			for (int x = 1; x < maxx - 1; x++)
			{
				float sum = data[index0 - 1] + data[index0] + data[index0 + 1] + data[index1 - 1] + data[index1] +
						data[index1 + 1] + data[index2 - 1] + data[index2] + data[index2 + 1];
				newData[index1] = floatGetAverage(sum, divisor);
				index0++;
				index1++;
				index2++;
			}
		}

		// Copy back
		for (int y = 1; y < maxy - 1; y++)
		{
			int index = y * maxx + 1;
			for (int x = 1; x < maxx - 1; x++, index++)
			{
				data[index] = newData[index];
			}
		}
	}

	private float[] floatBuffer(int size)
	{
		if (floatDataBuffer == null || floatDataBuffer.length < size)
		{
			floatDataBuffer = new float[size];
		}
		return floatDataBuffer;
	}

	/**
	 * Compute the block average within a 2n+1 size block around each point.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 */
	public void blockAverage(float[] data, final int maxx, final int maxy, final int n)
	{
		if (n == 1)
			blockAverage3x3(data, maxx, maxy);
		else
			blockAverageNxN(data, maxx, maxy, n);
	}

	/**
	 * Compute the block average within a 2n+1 size block around each point.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param n
	 *            The block size
	 */
	public void blockAverageNxN(float[] data, final int maxx, final int maxy, final int n)
	{
		float[] newData = floatBuffer(data.length);

		// Boundary control
		final int xwidth = FastMath.min(n, maxx - 1);
		final int ywidth = FastMath.min(n, maxy - 1);
		final int xlimit = maxx - xwidth;
		final int ylimit = maxy - ywidth;

		int[] offset = new int[(2 * xwidth + 1) * (2 * ywidth + 1) - 1];
		int[] xoffset = new int[offset.length];
		int[] yoffset = new int[offset.length];
		for (int y = -ywidth, d = 0; y <= ywidth; y++)
			for (int x = -xwidth; x <= xwidth; x++)
				if (x != 0 || y != 0)
				{
					offset[d] = maxx * y + x;
					xoffset[d] = x;
					yoffset[d] = y;
					d++;
				}

		final float divisor = (float) (1.0 / ((2 * n + 1) * (2 * n + 1)));

		int index = 0;
		for (int y = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++, index++)
			{
				float sum = data[index];

				// Flag to indicate this pixels has a complete (2n+1) neighbourhood 
				boolean isInnerXY = (y >= ywidth && y < ylimit) && (x >= xwidth && x < xlimit);

				// Sweep neighbourhood
				if (isInnerXY)
				{
					for (int d = offset.length; d-- > 0;)
					{
						sum += data[index + offset[d]];
					}
				}
				else
				{
					for (int d = offset.length; d-- > 0;)
					{
						// Get the pixel with boundary checking
						int yy = y + yoffset[d];
						int xx = x + xoffset[d];
						if (xx <= 0)
							xx = 0;
						else if (xx >= maxx)
							xx = maxx - 1;
						if (yy <= 0)
							yy = 0;
						else if (yy >= maxy)
							yy = maxy - 1;
						sum += data[xx + yy * maxx];
					}
				}
				newData[index] = floatGetAverage(sum, divisor);
			}
		}

		// Copy back
		for (index = data.length; index-- > 0;)
		{
			data[index] = newData[index];
		}
	}

	/**
	 * Compute the block average within a 3x3 size block around each point.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 */
	public void blockAverage3x3(float[] data, final int maxx, final int maxy)
	{
		float[] newData = floatBuffer(data.length);

		// Boundary control
		final int xwidth = 1;
		final int ywidth = 1;
		final int xlimit = maxx - xwidth;
		final int ylimit = maxy - ywidth;

		int[] offset = new int[8];
		int[] xoffset = new int[offset.length];
		int[] yoffset = new int[offset.length];
		for (int y = -ywidth, d = 0; y <= ywidth; y++)
			for (int x = -xwidth; x <= xwidth; x++)
				if (x != 0 || y != 0)
				{
					offset[d] = maxx * y + x;
					xoffset[d] = x;
					yoffset[d] = y;
					d++;
				}

		final float divisor = (float) (1.0 / 9);

		for (int y = 0; y < maxy; y++)
		{
			int index0 = (y - 1) * maxx;
			int index1 = y * maxx;
			int index2 = (y + 1) * maxx;
			for (int x = 0; x < maxx; x++)
			{
				// Flag to indicate this pixels has a complete (2n+1) neighbourhood 
				boolean isInnerXY = (y >= ywidth && y < ylimit) && (x >= xwidth && x < xlimit);

				// Sweep neighbourhood
				if (isInnerXY)
				{
					float sum = data[index0 - 1] + data[index0] + data[index0 + 1] + data[index1 - 1] + data[index1] +
							data[index1 + 1] + data[index2 - 1] + data[index2] + data[index2 + 1];
					newData[index1] = floatGetAverage(sum, divisor);
				}
				else
				{
					float sum = data[index1];

					for (int d = offset.length; d-- > 0;)
					{
						// Get the pixel with boundary checking
						int yy = y + yoffset[d];
						int xx = x + xoffset[d];
						if (xx <= 0)
							xx = 0;
						else if (xx >= maxx)
							xx = maxx - 1;
						if (yy <= 0)
							yy = 0;
						else if (yy >= maxy)
							yy = maxy - 1;
						sum += data[xx + yy * maxx];
					}
					newData[index1] = floatGetAverage(sum, divisor);
				}
				index0++;
				index1++;
				index2++;
			}
		}

		// Copy back
		for (int index = data.length; index-- > 0;)
		{
			data[index] = newData[index];
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public Object clone()
	{
		try
		{
			GaussianFilter o = (GaussianFilter) super.clone();
			o.floatDataBuffer = null;
			return o;
		}
		catch (CloneNotSupportedException e)
		{
			// Ignore
		}
		return null;
	}
}