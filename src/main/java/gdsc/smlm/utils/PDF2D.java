package gdsc.smlm.utils;

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
 * Provides sampling from a 2D histogram
 * <p>
 * Adapted from The GNU Scientific library (http://www.gnu.org/software/gsl/)
 */
public class PDF2D
{
	private final double[] sum;
	public final int nx, ny, n;
	/**
	 * The cumulative sum of the original input data
	 */
	public final double cumulative;

	/**
	 * Default constructor. Assumes the x-range and y-range increment from zero in integers.
	 * 
	 * @param data The data (packed in XY order, i = nx*y + x) 
	 * @param nx The X-dimension size
	 * @param ny The y-dimension size
	 * @throws InvalidArgumentException
	 *             if the dimensions are not above zero
	 * @throws InvalidArgumentException
	 *             if the input data length is not at least nx * ny
	 * @throws InvalidArgumentException
	 *             if the input data contains negatives
	 */
	public PDF2D(double[] data, int nx, int ny)
	{
		if (nx < 1 || ny < 1)
			throw new IllegalArgumentException("Dimensions must be above zero");
		this.nx = nx;
		this.ny = ny;
		n = nx * ny;

		if (data == null || data.length < n)
			throw new IllegalArgumentException("Input data must be at least equal to nx * ny");
		this.sum = new double[n + 1];

		double mean = 0, sum = 0;
		double c = 0;

		for (int i = 0; i < n; i++)
		{
			if (data[i] < 0)
				throw new IllegalArgumentException("Histogram bins must be non-negative");
			mean += (data[i] - mean) / ((double) (i + 1));
			c += data[i];
		}

		cumulative = c;
		
		this.sum[0] = 0;

		for (int i = 0; i < n; i++)
		{
			sum += (data[i] / mean) / n;
			this.sum[i + 1] = sum;
		}
	}

	/**
	 * Sample from the histogram using two uniform random numbers
	 * 
	 * @param r1
	 * @param r2
	 * @param point
	 *            The output coordinates buffer
	 * @return true if a sample was produced
	 */
	public boolean sample(double r1, double r2, double[] point)
	{
		if (point == null || point.length < 2)
			return false;

		/*
		 * Wrap the exclusive top of the bin down to the inclusive bottom of
		 * the bin. Since this is a single point it should not affect the
		 * distribution.
		 */

		if (r2 >= 1.0 || r2 < 0)
		{
			r2 = 0.0;
		}
		if (r1 >= 1.0 || r1 < 0)
		{
			r1 = 0.0;
		}

		int k = find(r1);

		if (k == -1)
			return false;

		int y = k / nx;
		int x = k - (y * nx);
		double delta = (r1 - sum[k]) / (sum[k + 1] - sum[k]);

		// Assume the x-range and y-range increment from zero in integers.
		// We could extend this class to support non-uniform ranges as per the GSL library
		//point[0]= xrange[x] + delta * (xrange[x + 1] - xrange[x]);
		//point[0] = yrange[y] + r2 * (yrange[y + 1] - yrange[y]);
		point[0] = x + delta;
		point[1] = y + r2;

		return true;
	}

	private int find(double x)
	{
		if (x >= sum[n])
		{
			return -1;
		}

		/* perform binary search */

		int upper = n;
		int lower = 0;

		while (upper - lower > 1)
		{
			final int mid = (upper + lower) / 2;

			if (x >= sum[mid])
			{
				lower = mid;
			}
			else
			{
				upper = mid;
			}
		}

		/* sanity check the result */

		if (x < sum[lower] || x >= sum[lower + 1])
		{
			return -1;
		}

		return lower;
	}

	/**
	 * Return the cumulative probability for the given coordinates
	 * @param x
	 * @param y
	 * @return p
	 */
	double get(int x, int y)
	{
		if (x < 0 || x >= nx || y < 0 || y > ny)
			return 0;
		return sum[y * nx + x];
	}
}
