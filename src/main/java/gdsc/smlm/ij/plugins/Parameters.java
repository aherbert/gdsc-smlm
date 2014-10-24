package gdsc.smlm.ij.plugins;

import java.util.EnumSet;

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
 * Utility methods for checking parameters
 */
public class Parameters
{
	/**
	 * Define the requirements for the parameter
	 */
	public enum Requirement
	{
		/**
		 * Value is greater than zero
		 */
		ABOVE_ZERO,
		/**
		 * Value is zero or above
		 */
		POSITIVE
	}
	
	/**
	 * Check if the named parameter value greater than zero
	 * @param name
	 * @param value
	 * @throws IllegalArgumentException
	 */
	public static void isAboveZero(String name, double value)
	{
		if (value <= 0)
			throw new IllegalArgumentException(name + " should be above zero");
	}
	
	/**
	 * Check if the named parameter value greater than the given limit
	 * @param name
	 * @param value
	 * @param limit
	 * @throws IllegalArgumentException
	 */
	public static void isAbove(String name, double value,  double limit)
	{
		if (value <= limit)
			throw new IllegalArgumentException(name + " should be > " + limit);
	}
	
	/**
	 * Check if the named parameter value greater than or equal to the given limit
	 * @param name
	 * @param value
	 * @param limit
	 * @throws IllegalArgumentException
	 */
	public static void isEqualOrAbove(String name, double value,  double limit)
	{
		if (value < limit)
			throw new IllegalArgumentException(name + " should be >= " + limit);
	}

	
	/**
	 * Check if the named parameter value is less than the given limit
	 * @param name
	 * @param value
	 * @param limit
	 * @throws IllegalArgumentException
	 */
	public static void isBelow(String name, double value,  double limit)
	{
		if (value >= limit)
			throw new IllegalArgumentException(name + " should be < " + limit);
	}
	
	/**
	 * Check if the named parameter value is less then or equal to the given limit
	 * @param name
	 * @param value
	 * @param limit
	 * @throws IllegalArgumentException
	 */
	public static void isEqualOrBelow(String name, double value,  double limit)
	{
		if (value > limit)
			throw new IllegalArgumentException(name + " should be <= " + limit);
	}
	
	/**
	 * Check if the named parameter value is zero or greater
	 * @param name
	 * @param value
	 * @throws IllegalArgumentException
	 */
	public static void isPositive(String name, double value)
	{
		if (value < 0)
			throw new IllegalArgumentException(name + " should be positive");
	}
	
	/**
	 * Check if the named parameter meets the requirements
	 * @param name
	 * @param requirements
	 */
	public static void isValid(String name, double value, EnumSet<Requirement> requirements)
	{
		if (requirements.contains(Requirement.ABOVE_ZERO))
			isAboveZero(name, value);
		if (requirements.contains(Requirement.POSITIVE))
			isPositive(name, value);
	}
}
