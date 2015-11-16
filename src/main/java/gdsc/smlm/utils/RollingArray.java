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

	/**
	 * Remove all the numbers from the array 
	 */
	public void reset()
	{
		sum = 0;
		index = 0;
		count = 0;
	}

	/**
	 * Add a number to the array
	 * @param d The number
	 */
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

	/**
	 * @return The count of numbers stored in the array
	 */
	public int getCount()
	{
		return count;
	}

	/**
	 * @return The capacity of the array
	 */
	public int getSize()
	{
		return size;
	}

	/**
	 * @return The sum using the rolling sum of the numbers (may accumulate errors)
	 */
	public double getSum()
	{
		return sum;
	}

	/**
	 * @return The recomputed sum using the current set of numbers
	 */
	public double getSum2()
	{
		double s = 0;
		// If full 'count' will be the length of the data array
		for (int i = 0; i < count; i++)
			s += data[i];

		// Debug
		//System.out.printf("Error = %g\n", DoubleEquality.relativeError(s, sum));

		// Reset the sum
		return sum = s;
	}

	/**
	 * @return The average using the rolling sum of the numbers
	 */
	public double getAverage()
	{
		return sum / count;
	}

	/**
	 * @return The average using a recomputed sum of the current numbers
	 */
	public double getAverage2()
	{
		return getSum2() / count;
	}

	/**
	 * @return True if full
	 */
	public boolean isFull()
	{
		return count == size;
	}
}
