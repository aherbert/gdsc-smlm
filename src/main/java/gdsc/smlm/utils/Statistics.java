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

/**
 * Simple class to calculate the mean and standard deviation of data
 */
public class Statistics
{
	protected int n = 0;
	protected double s = 0, ss = 0;

	public Statistics()
	{
	}

	public Statistics(float[] data)
	{
		add(data);
	}

	public Statistics(double[] data)
	{
		add(data);
	}

	public Statistics(int[] data)
	{
		add(data);
	}

	/**
	 * Add the data
	 * 
	 * @param data
	 */
	public void add(float[] data)
	{
		if (data == null)
			return;
		for (final float value : data)
		{
			n++;
			s += value;
			ss += value * value;
		}
	}

	/**
	 * Add the data
	 * 
	 * @param data
	 */
	public void add(double[] data)
	{
		if (data == null)
			return;
		for (final double value : data)
		{
			n++;
			s += value;
			ss += value * value;
		}
	}

	/**
	 * Add the data
	 * 
	 * @param data
	 */
	public void add(int[] data)
	{
		if (data == null)
			return;
		for (final double value : data)
		{
			n++;
			s += value;
			ss += value * value;
		}
	}

	/**
	 * Add the value
	 * 
	 * @param value
	 */
	public void add(final double value)
	{
		n++;
		s += value;
		ss += value * value;
	}

	/**
	 * Add the data. Synchronized for thread safety.
	 * 
	 * @param data
	 */
	synchronized public void safeAdd(float[] data)
	{
		if (data == null)
			return;
		for (final float value : data)
		{
			n++;
			s += value;
			ss += value * value;
		}
	}

	/**
	 * Add the data. Synchronized for thread safety.
	 * 
	 * @param data
	 */
	synchronized public void safeAdd(double[] data)
	{
		if (data == null)
			return;
		for (final double value : data)
		{
			n++;
			s += value;
			ss += value * value;
		}
	}

	/**
	 * Add the value. Synchronized for thread safety.
	 * 
	 * @param value
	 */
	synchronized public void safeAdd(final double value)
	{
		n++;
		s += value;
		ss += value * value;
	}

	/**
	 * @return The number of data points
	 */
	public int getN()
	{
		return n;
	}

	/**
	 * @return The sum of the data points
	 */
	public double getSum()
	{
		return s;
	}

	/**
	 * @return The sum of squares of the data points
	 */
	public double getSumOfSquares()
	{
		return ss;
	}

	/**
	 * @return The mean of the data points
	 */
	public double getMean()
	{
		return s / n;
	}

	/**
	 * @return The unbiased standard deviation of the data points
	 */
	public double getStandardDeviation()
	{
		final double d = n;
		double stdDev = ss - ((double) s * s) / d;
		if (stdDev > 0.0)
			stdDev = Math.sqrt(stdDev / (d - 1.0));
		else
			stdDev = 0.0;
		return stdDev;
	}

	/**
	 * @return The unbiased variance of the data points
	 */
	public double getVariance()
	{
		final double d = n;
		double variance = ss - ((double) s * s) / d;
		if (variance > 0.0)
			variance = variance / (d - 1.0);
		else
			variance = 0.0;
		return variance;
	}

	/**
	 * The standard error is the standard deviation of the sample-mean's estimate of a population mean.
	 * <p>
	 * Uses the unbiased standard deviation divided by the square root of the sample size.
	 * 
	 * @return The standard error
	 */
	public double getStandardError()
	{
		if (n > 0)
		{
			return getStandardDeviation() / Math.sqrt(n);
		}
		return 0;
	}

	/**
	 * Add the statistics to the data
	 * 
	 * @param statistics
	 */
	public void add(Statistics statistics)
	{
		n += statistics.n;
		s += statistics.s;
		ss += statistics.ss;
	}

	/**
	 * Add the statistics to the data. Synchronized for thread safety.
	 * 
	 * @param statistics
	 */
	synchronized public void safeAdd(Statistics statistics)
	{
		n += statistics.n;
		s += statistics.s;
		ss += statistics.ss;
	}
}