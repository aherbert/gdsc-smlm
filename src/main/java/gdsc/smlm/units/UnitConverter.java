package gdsc.smlm.units;

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
 * Define conversion of a unit
 */
public interface UnitConverter<T extends Unit>
{
	/**
	 * Convert the value.
	 *
	 * @param value the value
	 * @return the new value
	 */
	public double convert(double value);
	
	/**
	 * Specify the source unit to be converted from
	 *
	 * @return the source unit
	 */
	public T from();
	
	/**
	 * Specify the destination unit to be converted to
	 *
	 * @return the destination unit
	 */
	public T to();
}