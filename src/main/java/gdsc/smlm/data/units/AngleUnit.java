/*
 * 
 */
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
 * Unit for measuring angle
 */
public enum AngleUnit implements Unit
{
	/** Radian units */
	RADIAN
	{
		public String getName()
		{
			return "radian";
		}

		public String getShortName()
		{
			return "rad";
		}
	},

	/** Degree units */
	DEGREE
	{
		public String getName()
		{
			return "degree";
		}

		public String getShortName()
		{
			return "°";
		}
	},;

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Enum#toString()
	 */
	public String toString()
	{
		return getName();
	}
}