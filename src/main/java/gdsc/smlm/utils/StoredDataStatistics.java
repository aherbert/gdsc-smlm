package gdsc.smlm.utils;

import java.util.Arrays;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

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
 * Calculate the mean and standard deviation of data. Stores the data for later retrieval.
 */
public class StoredDataStatistics extends Statistics
{
	private double[] values = new double[0];
	private DescriptiveStatistics stats = null;

	public StoredDataStatistics()
	{
	}

	public StoredDataStatistics(int capacity)
	{
		values = new double[capacity];
	}

	public StoredDataStatistics(float[] data)
	{
		add(data);
	}

	public StoredDataStatistics(double[] data)
	{
		add(data);
	}

	public StoredDataStatistics(int[] data)
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
		checkCapacity(data.length);
		for (final float value : data)
		{
			values[n++] = value;
			s += value;
			ss += value * value;
		}
	}

	/**
	 * Ensure that the specified number of elements can be added to the array.
	 * <p>
	 * This is not synchronized. However any class using the safeAdd() methods in different threads should be using the
	 * same synchronized method to add data thus this method will be within synchronized code.
	 * 
	 * @param length
	 */
	private void checkCapacity(int length)
	{
		stats = null;
		final int minCapacity = n + length;
		final int oldCapacity = values.length;
		if (minCapacity > oldCapacity)
		{
			int newCapacity = (oldCapacity * 3) / 2 + 1;
			if (newCapacity < minCapacity)
				newCapacity = minCapacity;
			double[] newValues = new double[newCapacity];
			System.arraycopy(values, 0, newValues, 0, n);
			values = newValues;
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
		checkCapacity(data.length);
		for (final double value : data)
		{
			values[n++] = value;
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
		checkCapacity(data.length);
		for (final int value : data)
		{
			values[n++] = value;
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
		checkCapacity(1);
		values[n++] = value;
		s += value;
		ss += value * value;
	}

	/**
	 * Add the data. Synchronized for thread safety. (Multiple threads must all use the same safeAdd method to ensure
	 * thread safety.)
	 * 
	 * @param data
	 */
	synchronized public void safeAdd(float[] data)
	{
		if (data == null)
			return;
		checkCapacity(data.length);
		for (final float value : data)
		{
			values[n++] = value;
			s += value;
			ss += value * value;
		}
	}

	/**
	 * Add the data. Synchronized for thread safety. (Multiple threads must all use the same safeAdd method to ensure
	 * thread safety.)
	 * 
	 * @param data
	 */
	synchronized public void safeAdd(double[] data)
	{
		if (data == null)
			return;
		checkCapacity(data.length);
		for (final double value : data)
		{
			values[n++] = value;
			s += value;
			ss += value * value;
		}
	}

	/**
	 * Add the value. Synchronized for thread safety. (Multiple threads must all use the same safeAdd method to ensure
	 * thread safety.)
	 * 
	 * @param value
	 */
	synchronized public void safeAdd(final double value)
	{
		checkCapacity(1);
		values[n++] = value;
		s += value;
		ss += value * value;
	}

	/**
	 * @return A copy of the values added
	 */
	public double[] getValues()
	{
		return Arrays.copyOf(values, n);
	}

	/**
	 * @return A copy of the values added
	 */
	public float[] getFloatValues()
	{
		float[] data = new float[n];
		for (int i = 0; i < n; i++)
			data[i] = (float) values[i];
		return data;
	}

	/**
	 * @return object used to compute descriptive statistics. The object is cached
	 * @see {@link org.apache.commons.math3.stat.descriptive.DescriptiveStatistics }
	 */
	public DescriptiveStatistics getStatistics()
	{
		if (stats == null)
			stats = new DescriptiveStatistics(values);
		return stats;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.utils.Statistics#add(gdsc.smlm.utils.Statistics)
	 */
	@Override
	public void add(Statistics statistics)
	{
		if (statistics instanceof StoredDataStatistics)
		{
			StoredDataStatistics extra = (StoredDataStatistics) statistics;
			if (extra.n > 0)
			{
				checkCapacity(extra.n);
				System.arraycopy(extra.values, 0, values, n, extra.n);
			}
		}
		super.add(statistics);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.utils.Statistics#safeAdd(gdsc.smlm.utils.Statistics)
	 */
	@Override
	synchronized public void safeAdd(Statistics statistics)
	{
		this.add(statistics);
	}

	/**
	 * @return The median
	 */
	public double getMedian()
	{
		// Check for negatives
		for (double d : values)
		{
			if (d < 0)
			{
				if (n == 0)
					return Double.NaN;
				if (n == 1)
					return values[0];

				double[] data = getValues();
				Arrays.sort(data);
				return (data[(data.length - 1) / 2] + data[data.length / 2]) * 0.5;
			}
		}

		// This does not work when the array contains negative data due to the 
		// implementation of the library using partially sorted data
		return getStatistics().getPercentile(50);
	}
}