package gdsc.smlm.search;

import java.util.Arrays;

import gdsc.core.utils.Maths;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Specify the dimensions for a search
 */
public class SearchDimension implements Cloneable
{
	public final double min;
	public final double max;
	public final double minIncrement;
	public final int nIncrement;
	public final boolean active;

	private double centre;
	private double increment;
	private double reduceFactor = 0.5;

	/**
	 * Instantiates a new inactive search dimension. The centre can be set to any value, the default is zero.
	 */
	public SearchDimension()
	{
		this(0);
	}

	/**
	 * Instantiates a new inactive search dimension. The centre can be set to any value.
	 *
	 * @param centre
	 *            the centre
	 */
	public SearchDimension(double centre)
	{
		this(0, 0, 0, 0);
		setCentre(centre);
	}

	/**
	 * Instantiates a new search dimension.
	 *
	 * @param min
	 *            the minimum of the range
	 * @param max
	 *            the maximum of the range
	 * @param minIncrement
	 *            the min increment to use around the centre
	 * @param nIncrement
	 *            the number of increments to use around the centre
	 */
	public SearchDimension(double min, double max, double minIncrement, int nIncrement)
	{
		this(min, max, minIncrement, nIncrement, min, max);
	}

	/**
	 * Instantiates a new search dimension.
	 *
	 * @param min
	 *            the minimum of the range
	 * @param max
	 *            the maximum of the range
	 * @param minIncrement
	 *            the min increment to use around the centre
	 * @param nIncrement
	 *            the number of increments to use around the centre
	 * @param lower
	 *            the current lower bound of the range
	 * @param upper
	 *            the current upper bound of the range
	 */
	public SearchDimension(double min, double max, double minIncrement, int nIncrement, double lower, double upper)
	{
		if (isInvalid(min))
			throw new IllegalArgumentException("Min is not a valid number: " + min);
		if (isInvalid(max))
			throw new IllegalArgumentException("Max is not a valid number: " + max);
		if (isInvalid(lower))
			throw new IllegalArgumentException("Lower is not a valid number: " + lower);
		if (isInvalid(upper))
			throw new IllegalArgumentException("Upper is not a valid number: " + upper);
		if (isInvalid(minIncrement))
			throw new IllegalArgumentException("Min increment is not a valid number: " + minIncrement);
		if (max < min)
			throw new IllegalArgumentException("Max is less than min");
		this.active = min < max;
		if (active && nIncrement < 1)
			throw new IllegalArgumentException("Steps must be more than 0: " + nIncrement);
		if (minIncrement < 0)
			throw new IllegalArgumentException("Min increment is negative: " + minIncrement);
		if (upper < lower)
			throw new IllegalArgumentException("Upper is less than lower");
		if (upper < min || upper > max)
			throw new IllegalArgumentException("Upper is outside min/max range");
		if (lower < min || lower > max)
			throw new IllegalArgumentException("Lower is outside min/max range");

		// We round to the min increment so that the values returned should be identical if the centre is moved by a factor of the increment.
		this.minIncrement = minIncrement;
		this.min = round(min);
		this.max = round(max);
		this.nIncrement = nIncrement;

		setCentre((upper + lower) / 2);
		setIncrement((upper - lower) / (2 * nIncrement));
	}

	private static boolean isInvalid(double d)
	{
		return Double.isNaN(d) || Double.isInfinite(d);
	}

	private double round(double d)
	{
		if (minIncrement != 0)
			return Maths.round(d, minIncrement);
		return d;
	}

	/**
	 * Sets the centre.
	 *
	 * @param centre
	 *            the new centre of the range in the dimension
	 */
	public void setCentre(double centre)
	{
		if (active && (centre < min || centre > max))
			throw new IllegalArgumentException("Centre is outside min/max range");
		this.centre = round(centre);
	}

	/**
	 * Gets the centre of the range in the dimension
	 *
	 * @return the centre of the range in the dimension
	 */
	public double getCentre()
	{
		return centre;
	}

	/**
	 * Gets the increment.
	 *
	 * @return the increment
	 */
	public double getIncrement()
	{
		return increment;
	}

	/**
	 * Checks if is at min increment.
	 *
	 * @return true, if is at min increment
	 */
	public boolean isAtMinIncrement()
	{
		return !active || increment == minIncrement;
	}

	/**
	 * Sets the increment.
	 *
	 * @param increment
	 *            the new increment
	 */
	public void setIncrement(double increment)
	{
		if (increment < minIncrement)
			increment = minIncrement;
		this.increment = round(increment);
	}

	/**
	 * Gets the current lower bound of the range
	 *
	 * @return the current lower bound of the range
	 */
	public double getLower()
	{
		final double v = centre - nIncrement * increment;
		if (v < min)
			return min;
		return round(v);
	}

	/**
	 * Gets the current upper bound of the range
	 *
	 * @return the current upper bound of the range
	 */
	public double getUpper()
	{
		final double v = centre + nIncrement * increment;
		if (v > max)
			return max;
		return round(v);
	}

	/**
	 * Checks if the value is at the bounds of the current dimension range.
	 *
	 * @param v
	 *            the value
	 * @return true, if is at bounds
	 */
	public boolean isAtBounds(double v)
	{
		if (v == getLower())
			return true;
		if (v == getUpper())
			return true;
		return false;
	}

	/**
	 * Get the values of the dimension around the current centre using the configured increments.
	 * Note: If the values are outside the min/max range then the number of values may be reduced.
	 *
	 * @return the values
	 */
	public double[] values()
	{
		if (!active)
			return new double[] { centre };

		final double[] values = new double[getMaxLength()];
		int size = 0;
		for (int i = 1; i <= nIncrement; i++)
		{
			double value = round(centre - i * increment);
			if (value < min)
				value = min;
			values[size++] = value;
		}
		values[size++] = centre; // Already rounded
		for (int i = 1; i <= nIncrement; i++)
		{
			double value = round(centre + i * increment);
			if (value > max)
				value = max;
			values[size++] = value;
		}

		// Check for duplicates if at the limits
		if (values[0] == min || values[size - 1] == max)
		{
			double last = values[0];
			size = 1;
			for (int i = 1; i < values.length; i++)
			{
				if (last == values[i])
					continue;
				values[size++] = last = values[i];
			}
		}

		return (size == values.length) ? values : Arrays.copyOf(values, size);
	}

	/**
	 * Gets the max length of the values array
	 *
	 * @return the max length
	 */
	public int getMaxLength()
	{
		return 2 * nIncrement + 1;
	}

	/**
	 * @return the reduceFactor
	 */
	public double getReduceFactor()
	{
		return reduceFactor;
	}

	/**
	 * @param reduceFactor
	 *            the reduceFactor to set
	 */
	public void setReduceFactor(double reduceFactor)
	{
		if (reduceFactor <= 0 || reduceFactor >= 1)
			throw new IllegalArgumentException("Reduce factor must be between 0 and 1 (exclusive)");
		this.reduceFactor = reduceFactor;
	}

	/**
	 * Reduce the size of the increment by multiplying by the reduce factor
	 */
	public void reduce()
	{
		setIncrement(increment * reduceFactor);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	public SearchDimension clone()
	{
		try
		{
			return (SearchDimension) super.clone();
		}
		catch (CloneNotSupportedException e)
		{
			return null;
		}

	}
}
