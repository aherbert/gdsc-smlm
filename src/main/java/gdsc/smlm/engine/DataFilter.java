package gdsc.smlm.engine;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Define the criteria used to filter the input data before identifying candidate peaks
 */
public enum DataFilter
{
	/**
	 * Use a mean within a specified box area
	 */
	MEAN("Mean"),
	/**
	 * Use a Gaussian with a specified radius
	 */
	GAUSSIAN("Gaussian"),
	/**
	 * Use a median within a specified box area
	 */
	MEDIAN("Median");

	private String name;

	private DataFilter(String name)
	{
		this.name = name;
	}

	@Override
	public String toString()
	{
		return name;
	}
}