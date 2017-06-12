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
 * Unit for measuring distance
 */
public enum DistanceUnit implements Unit
{
	/** Pixel units */
	PIXEL
	{
		public String getName()
		{
			return "pixel";
		}

		public String getShortName()
		{
			return "px";
		}
	},

	/** Micrometer units */
	UM
	{
		public String getName()
		{
			return "micrometer";
		}

		public String getShortName()
		{
			return "um";
		}
	},

	/** Nanometer units */
	NM
	{
		public String getName()
		{
			return "nanometer";
		}

		public String getShortName()
		{
			return "nm";
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