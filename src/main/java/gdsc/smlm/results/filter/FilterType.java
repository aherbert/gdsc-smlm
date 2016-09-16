package gdsc.smlm.results.filter;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Define the type of filter
 */
public enum FilterType
{
	/**
	 * A basic filter
	 */
	STANDARD("Standard"),
	/**
	 * A direct filter. This can perform filtering on PreprocessedPeakResult objects.
	 */
	DIRECT("Direct");

	private String name;

	private FilterType(String name)
	{
		this.name = name;
	}

	@Override
	public String toString()
	{
		return name;
	}
}