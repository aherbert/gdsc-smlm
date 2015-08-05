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
 * Provides sampling from a 1D histogram
 * <p>
 * Adapted from The GNU Scientific library (http://www.gnu.org/software/gsl/)
 */
public class PDF
{
	private final double[] sum;

	/**
	 * The cumulative sum of the original input data
	 */
	public final double cumulative;

	/**
	 * Default constructor. Assumes the range increments from zero in integers.
	 * 
	 * @param data
	 *            The data
	 * @throws InvalidArgumentException
	 *             if the input data length is not at least 1
	 * @throws InvalidArgumentException
	 *             if the input data contains negatives
	 */
	public PDF(double[] data)
	{
		if (data == null || data.length < 1)
			throw new IllegalArgumentException("Input data must be at least 1");

		this.sum = new double[data.length + 1];

		double mean = 0, sum = 0;
		double c = 0;

		for (int i = 0; i < data.length; i++)
		{
			if (data[i] < 0)
				throw new IllegalArgumentException("Histogram bins must be non-negative");
			mean += (data[i] - mean) / ((double) (i + 1));
			c += data[i];
		}

		cumulative = c;

		this.sum[0] = 0;

		for (int i = 0; i < data.length; i++)
		{
			sum += (data[i] / mean) / data.length;
			this.sum[i + 1] = sum;
		}
	}

	/**
	 * Sample from the PDF using a uniform random number
	 * 
	 * @param r1
	 * @return the sample (or -1 on error)
	 */
	public double sample(double r1)
	{
		/*
		 * Wrap the exclusive top of the bin down to the inclusive bottom of
		 * the bin. Since this is a single point it should not affect the
		 * distribution.
		 */

		if (r1 >= 1.0 || r1 < 0)
		{
			r1 = 0.0;
		}

		int k = find(r1);

		if (k == -1)
			return -1;

		double delta = (r1 - sum[k]) / (sum[k + 1] - sum[k]);

		// Assume the x-range and y-range increment from zero in integers.
		// We could extend this class to support non-uniform ranges as per the GSL library:
		// x = xrange[x] + delta * (xrange[x + 1] - xrange[x]);

		return k + delta;
	}

	private int find(double x)
	{
		if (x >= sum[sum.length - 1])
		{
			return -1;
		}

		/* perform binary search */

		int upper = sum.length - 1;
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
	 * 
	 * @param x
	 * @return p
	 */
	double get(int x)
	{
		if (x < 0 || x >= sum[sum.length - 1])
			return 0;
		return sum[x];
	}
}
