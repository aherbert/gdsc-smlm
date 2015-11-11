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
 * Provide a rolling array
 */
public class RollingArray
{
	private final double[] data;
	private final int size;
	private int index, count;
	private double sum;

	/**
	 * Create a rolling average
	 * 
	 * @param size
	 */
	public RollingArray(int size)
	{
		this.size = size;
		this.data = new double[size];
	}

	public void reset()
	{
		sum = 0;
		index = 0;
		count = 0;
	}

	public void add(double d)
	{
		// Add to the total
		sum += d;
		// If at capacity
		if (isFull())
		{
			// Subtract the item to be replaced
			sum -= data[index];
		}
		else
		{
			// Otherwise increase the count
			count++;
		}
		// Replace the item
		data[index++] = d;
		// Wrap the index
		if (index == size)
			index = 0;
	}

	public int getCount()
	{
		return count;
	}

	public int getSize()
	{
		return size;
	}

	public double getSum()
	{
		return sum;
	}

	public double getAverage()
	{
		return sum / count;
	}

	public boolean isFull()
	{
		return count == size;
	}
}
