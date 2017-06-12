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
 * Unit for measuring time
 */
public enum TimeUnit implements Unit
{
	/** Frame units */
	FRAME
	{
		public String getName()
		{
			return "frame";
		}

		public String getShortName()
		{
			return "t"; // Follow XYZCT conversion for 5D image stacks
		}
	},

	/** Second units */
	SECOND
	{
		public String getName()
		{
			return "second";
		}

		public String getShortName()
		{
			return "s";
		}
	},

	/** Millesecond units */
	MILLISECOND
	{
		public String getName()
		{
			return "millisecond";
		}

		public String getShortName()
		{
			return "ms";
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