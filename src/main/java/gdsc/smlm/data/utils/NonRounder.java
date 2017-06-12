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
 * Class to the Rounder interface that does not perform rounding
 */
public class NonRounder implements Rounder
{
	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.data.Rounder#round(double)
	 */
	public double round(double value)
	{
		return value;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.data.Rounder#toString(double)
	 */
	public String toString(double value)
	{
		return Double.toString(value);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.data.Rounder#round(float)
	 */
	public float round(float value)
	{
		return value;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.data.Rounder#toString(float)
	 */
	public String toString(float value)
	{
		return Float.toString(value);
	}
}