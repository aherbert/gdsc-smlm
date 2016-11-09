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
	private double[] values;
	private boolean pad = false;

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
		this(0, 0, 0, 1);
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
		
		// Rounding changes the range so bring the upper and lower back within
		lower = Math.min(this.max, Math.max(lower, this.min));
		upper = Math.min(this.max, Math.max(upper, this.min));

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
		values = null;
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
		values = null;
	}

	/**
	 * Gets the current lower bound of the range
	 *
	 * @return the current lower bound of the range
	 */
	public double getLower()
	{
		return values()[0];
	}

	/**
	 * Gets the current upper bound of the range
	 *
	 * @return the current upper bound of the range
	 */
	public double getUpper()
	{
		return values()[0];
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
		values();
		if (v == values[0])
			return true;
		if (v == values[values.length - 1])
			return true;
		return false;
	}

	/**
	 * Get the values of the dimension around the current centre using the configured increments.
	 * Note: If the values are outside the min/max range then by default the number of values may be reduced.
	 * <p>
	 * If the pad setting is enabled then the number of values should remain constant as the values are padded in the
	 * opposite direction.
	 *
	 * @return the values
	 */
	public double[] values()
	{
		if (values != null)
			return values;

		if (!active)
			return values = new double[] { centre };

		values = new double[getMaxLength()];
		int size = 0;

		double value = round(centre - nIncrement * increment);
		if (value < min)
		{
			values[size++] = min;

			// Avoid further values below the min
			for (int i = nIncrement - 1; i >= 1; i--)
			{
				value = round(centre - i * increment);
				if (value < min)
					continue;
				values[size++] = value;
			}

			if (centre != min)
				values[size++] = centre;
		}
		else
		{
			// Not at the limit
			for (int i = nIncrement; i >= 1; i--)
			{
				values[size++] = round(centre - i * increment);
			}
			
			values[size++] = centre; // Already rounded and within range
		}
		
		for (int i = 1; i <= nIncrement; i++)
		{
			value = round(centre + i * increment);
			if (value > max)
			{
				if (centre != max)
					values[size++] = max;
				// Avoid further values outside the range
				break;
			}
			values[size++] = value;
		}

		//		double[] check = values.clone();
		//		Arrays.sort(check);
		//		for (int i=0; i<check.length; i++)
		//			if (check[i] != values[i])
		//				throw new RuntimeException("Not sorted");

		// Check for duplicates if at the limits
		if (size != values.length)
		{
			// Option to pad in the opposite direction
			if (pad)
			{
				if (values[0] == min)
				{
					if (values[size - 1] != max)
					{
						// Pad up
						for (int i = nIncrement + 1; size < values.length; i++)
						{
							value = round(centre + i * increment);
							if (value > max)
							{
								values[size++] = max;
								break;
							}
							values[size++] = value;
						}
					}
				}
				else
				{
					// Pad down
					for (int i = nIncrement + 1; size < values.length; i++)
					{
						value = round(centre - i * increment);
						if (value < min)
						{
							values[size++] = min;
							break;
						}
						values[size++] = value;
					}
					// Simple option is to sort.
					// A better option is to copy the values to the correct place.
					Arrays.sort(values, 0, size);
				}

				// In case we could not pad enough
				if (size != values.length)
					values = Arrays.copyOf(values, size);
			}
			else
			{
				// No padding so truncate 
				values = Arrays.copyOf(values, size);
			}
		}

		return values;
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

	/**
	 * Checks if padding the values in the opposite direction when the range overlaps the min/max.
	 *
	 * @return true, if padding the values
	 */
	public boolean isPad()
	{
		return pad;
	}

	/**
	 * Set to true if padding the values in the opposite direction when the range overlaps the min/max
	 *
	 * @param pad
	 *            true, if padding the values in the opposite direction when the range overlaps the min/max
	 */
	public void setPad(boolean pad)
	{
		this.pad = pad;
		values = null;
	}
}
