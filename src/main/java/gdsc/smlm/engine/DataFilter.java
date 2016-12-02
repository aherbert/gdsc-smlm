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
	//@formatter:off
	/**
	 * Use a mean within a specified area
	 */
	MEAN{ public String getName() { return "Mean"; }},
	/**
	 * Use a mean within a specified box area. The box has integer size.
	 */
	BLOCK_MEAN{ public String getName() { return "Block mean"; }},
	/**
	 * Use a mean within a specified circle area
	 */
	CIRCULAR_MEAN{ public String getName() { return "Circular mean"; }},
	/**
	 * Use a Gaussian with a specified radius
	 */
	GAUSSIAN{ public String getName() { return "Gaussian"; }},
	/**
	 * Use a median within a specified box area. The box has integer size.
	 */
	MEDIAN{ public String getName() { return "Median"; }};
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