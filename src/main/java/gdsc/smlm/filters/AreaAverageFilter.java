package gdsc.smlm.filters;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2015 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Computes the area average for each point within the array.
 * <p>
 * The algorithm computes the average in the specified area. If the width is an integer the area is defined using an
 * entire block. If the width is non-integer it is rounded down and up to the adjacent integers. The sum of each block
 * is computed. The larger block is subtracted from the smaller block to give the edge sum. The average is then computed
 * using the sum of the smaller block and a proportion of the edge sum.
 * <p>
 * Note: Due to lack of small dimension checking the routines will fail if maxx or maxy are less than 2. All routines
 * are OK for 3x3 images and larger.
 */
public class AreaAverageFilter implements Cloneable
{
	private SumFilter filter = new SumFilter();
	private AverageFilter avFilter = new AverageFilter();

	private boolean simpleInterpolation = true;

	/**
	 * Compute the block average within a 2w+1 size block around each point.
	 * Pixels within border regions (width = ceil(w)) are unchanged.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param w
	 *            The width
	 */
	public void areaAverageUsingSumsInternal(float[] data, final int maxx, final int maxy, final double w)
	{
		//if (w < 1)
		//{
		//	// For small widths then use a dedicated filter
		//	avFilter.blockAverage3x3Internal(data, maxx, maxy, (float) w);
		//	return;
		//}

		if (w <= 0)
			return;

		final int n = (int) w;
		final int n1 = (n == w) ? n : n + 1;
		final int blockSize = 2 * n1 + 1;
		if (maxx < blockSize || maxy < blockSize)
			return;

		if (n == n1)
		{
			// There is no edge
			avFilter.rollingBlockAverageInternal(data, maxx, maxy, n);
			return;
		}

		// Calculate the sum in the n and n+1 regions
		final float[] sum1;
		if (n == 0)
		{
			sum1 = data;
		}
		else
		{
			sum1 = data.clone();
			filter.rollingBlockSum(sum1, maxx, maxy, n);
		}
		final float[] sum2 = data.clone();
		filter.rollingBlockSumInternal(sum2, maxx, maxy, n1);

		// Get the average by adding the inner sum to the weighted edge pixels.  
		final double area = (2 * w + 1) * (2 * w + 1);

		final float edgeWeight;

		if (simpleInterpolation)
		{
			edgeWeight = (float) (w - n);
		}
		else
		{
			// Use the area to produce the weighting
			final int inner = (2 * n + 1) * (2 * n + 1);
			final int outer = (2 * n1 + 1) * (2 * n1 + 1);
			edgeWeight = (float) ((area - inner) / (outer - inner));
		}

		final float norm = (float) (1.0 / area);

		for (int y = n1; y < maxy - n1; y++)
		{
			int index = y * maxx + n1;
			for (int x = n1; x < maxx - n1; x++, index++)
			{
				data[index] = norm * (sum1[index] + edgeWeight * (sum2[index] - sum1[index]));
			}
		}
	}

	/**
	 * Compute the block average within a 2w+1 size block around each point.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param w
	 *            The width
	 */
	public void areaAverageUsingSums(float[] data, final int maxx, final int maxy, final double w)
	{
		//if (w < 1)
		//{
		//	// For small widths then use a dedicated filter
		//	avFilter.blockAverage3x3(data, maxx, maxy, (float) w);
		//	return;
		//}

		if (w <= 0)
			return;

		final int n = (int) w;
		final int n1 = (n == w) ? n : n + 1;

		if (n == n1 || (maxx < n1 && maxy < n1))
		{
			// There is no edge
			avFilter.rollingBlockAverage(data, maxx, maxy, n);
			return;
		}

		// Calculate the sum in the n and n+1 regions
		final float[] sum1;
		if (n == 0)
		{
			sum1 = data;
		}
		else
		{
			sum1 = data.clone();
			filter.rollingBlockSum(sum1, maxx, maxy, n);
		}
		final float[] sum2 = data.clone();
		filter.rollingBlockSum(sum2, maxx, maxy, n1);

		// Get the average by adding the inner sum to the weighted edge pixels.  
		final double area = (2 * w + 1) * (2 * w + 1);

		final float edgeWeight;

		if (simpleInterpolation)
		{
			edgeWeight = (float) (w - n);
		}
		else
		{
			// Use the area to produce the weighting
			final int inner = (2 * n + 1) * (2 * n + 1);
			final int outer = (2 * n1 + 1) * (2 * n1 + 1);
			edgeWeight = (float) ((area - inner) / (outer - inner));
		}

		final float norm = (float) (1.0 / area);

		for (int index = 0; index < sum1.length; index++)
		{
			data[index] = norm * (sum1[index] + edgeWeight * (sum2[index] - sum1[index]));
		}
	}

	/**
	 * Compute the block average within a 2w+1 size block around each point.
	 * Pixels within border regions (width = ceil(w)) are unchanged.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param w
	 *            The width
	 */
	public void areaAverageInternal(float[] data, final int maxx, final int maxy, final double w)
	{
		//if (w < 1)
		//{
		//	// For small widths then use a dedicated filter
		//	avFilter.blockAverage3x3Internal(data, maxx, maxy, (float) w);
		//	return;
		//}

		if (w <= 0)
			return;

		final int n = (int) w;
		final int n1 = (n == w) ? n : n + 1;
		final int blockSize = 2 * n1 + 1;
		if (maxx < blockSize || maxy < blockSize)
			return;

		if (n == n1)
		{
			// There is no edge
			avFilter.rollingBlockAverageInternal(data, maxx, maxy, n);
			return;
		}

		// Calculate the sum in the n and n+1 regions
		final float[] av1;
		if (n == 0)
		{
			av1 = data;
		}
		else
		{
			av1 = data.clone();
			avFilter.rollingBlockAverage(av1, maxx, maxy, n);
		}
		final float[] av2 = data.clone();
		avFilter.rollingBlockAverageInternal(av2, maxx, maxy, n1);

		// Get the average by weighting the two
		final float outerWeight;

		if (simpleInterpolation)
		{
			outerWeight = (float) (w - n);
		}
		else
		{
			// Use the area to produce the weighting
			final int inner = (2 * n + 1) * (2 * n + 1);
			final int outer = (2 * n1 + 1) * (2 * n1 + 1);
			final double area = (2 * w + 1) * (2 * w + 1);
			outerWeight = (float) ((area - inner) / (outer - inner));
		}

		final float innerWeight = 1 - outerWeight;

		for (int y = n1; y < maxy - n1; y++)
		{
			int index = y * maxx + n1;
			for (int x = n1; x < maxx - n1; x++, index++)
			{
				data[index] = av1[index] * innerWeight + av2[index] * outerWeight;
			}
		}
	}

	/**
	 * Compute the block average within a 2w+1 size block around each point.
	 * <p>
	 * Note: the input data is destructively modified
	 * 
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 * @param w
	 *            The width
	 */
	public void areaAverage(float[] data, final int maxx, final int maxy, final double w)
	{
		//if (w < 1)
		//{
		//	// For small widths then use a dedicated filter
		//	avFilter.blockAverage3x3(data, maxx, maxy, (float) w);
		//	return;
		//}

		if (w <= 0)
			return;

		final int n = (int) w;
		final int n1 = (n == w) ? n : n + 1;

		if (n == n1 || (maxx < n1 && maxy < n1))
		{
			// There is no edge
			avFilter.rollingBlockAverage(data, maxx, maxy, n);
			return;
		}

		// Calculate the sum in the n and n+1 regions
		final float[] av1;
		if (n == 0)
		{
			av1 = data;
		}
		else
		{
			av1 = data.clone();
			avFilter.rollingBlockAverage(av1, maxx, maxy, n);
		}
		final float[] av2 = data.clone();
		avFilter.rollingBlockAverage(av2, maxx, maxy, n1);

		// Get the average by weighting the two
		final float outerWeight;

		if (simpleInterpolation)
		{
			outerWeight = (float) (w - n);
		}
		else
		{
			// Use the area to produce the weighting
			final int inner = (2 * n + 1) * (2 * n + 1);
			final int outer = (2 * n1 + 1) * (2 * n1 + 1);
			final double area = (2 * w + 1) * (2 * w + 1);
			outerWeight = (float) ((area - inner) / (outer - inner));
		}

		final float innerWeight = 1 - outerWeight;

		for (int index = 0; index < av1.length; index++)
		{
			data[index] = av1[index] * innerWeight + av2[index] * outerWeight;
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
			AreaAverageFilter o = (AreaAverageFilter) super.clone();
			o.filter = (SumFilter) filter.clone();
			o.avFilter = (AverageFilter) avFilter.clone();
			return o;
		}
		catch (CloneNotSupportedException e)
		{
			// Ignore
		}
		return null;
	}

	/**
	 * @return true for simple interpolation
	 */
	public boolean isSimpleInterpolation()
	{
		return simpleInterpolation;
	}

	/**
	 * The average for block size n and n+1 is linear interpolated. Set this to true to use a weight of (w-n) for the
	 * outer average. Set to false to use a weight based on the area of the edge pixels in the (2w+1) region.
	 * 
	 * @param simpleInterpolation
	 *            true for simple interpolation
	 */
	public void setSimpleInterpolation(boolean simpleInterpolation)
	{
		this.simpleInterpolation = simpleInterpolation;
	}
}