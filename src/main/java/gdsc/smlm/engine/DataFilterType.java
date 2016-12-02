package gdsc.smlm.engine;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2015 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Define the type of filter used for identifying candidate peaks
 */
public enum DataFilterType
{
	//@formatter:off
	/**
	 * Use a single filter
	 */
	SINGLE{ public String getName() { return "Single"; }},
	/**
	 * Use a difference filter (the second subtracted from the first).
	 */
	DIFFERENCE{ public String getName() { return "Difference"; }},
	/**
	 * Use a jury of multiple filters
	 */
	JURY{ public String getName() { return "Jury"; }};
	//@formatter:on

	@Override
	public String toString()
	{
		return getName();
	}

	/**
	 * Gets the name.
	 *
	 * @return the name
	 */
	abstract public String getName();
}