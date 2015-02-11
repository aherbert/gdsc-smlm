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
	/**
	 * Use a single filter
	 */
	SINGLE("Single"),
	/**
	 * Use a difference filter (the second subtracted from the first).
	 */
	DIFFERENCE("Difference"),
	/**
	 * Use a jury of multiple filters
	 */
	JURY("Jury");

	private String name;

	private DataFilterType(String name)
	{
		this.name = name;
	}

	@Override
	public String toString()
	{
		return name;
	}
}