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
 * Unit for measuring intensity
 */
public enum IntensityUnit implements Unit
{
	/** Camera count units */
	COUNT
	{
		public String getName()
		{
			return "count";
		}

		public String getShortName()
		{
			return "count"; // Nothing suitable
		}
	},

	/** Photon units */
	PHOTON
	{
		public String getName()
		{
			return "photon";
		}

		public String getShortName()
		{
			return "photon"; // Nothing suitable
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