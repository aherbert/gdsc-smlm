package gdsc.smlm.data.units;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Create a Rounder interface implementation
 */
public class RounderFactory
{
	/**
	 * Creates the rounder. If the precision is less than 1 then an instance will be created that does not perform
	 * rounding.
	 *
	 * @param precision
	 *            the precision
	 * @return the rounder
	 */
	public static Rounder create(int precision)
	{
		return (precision > 0) ? new MathContextRounder(precision) : new NonRounder();
	}
}