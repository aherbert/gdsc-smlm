package gdsc.smlm.data.utils;

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
 * Interface for rounding
 */
public interface Rounder
{
	/**
	 * Round the value.
	 *
	 * @param value
	 *            the value
	 * @return the rounded value
	 */
	public double round(double value);

	/**
	 * Round the value to a string.
	 *
	 * @param value
	 *            the value
	 * @return the rounded string value
	 */
	public String toString(double value);
	
	/**
	 * Round the value.
	 *
	 * @param value
	 *            the value
	 * @return the rounded value
	 */
	public float round(float value);

	/**
	 * Round the value to a string.
	 *
	 * @param value
	 *            the value
	 * @return the rounded string value
	 */
	public String toString(float value);
}