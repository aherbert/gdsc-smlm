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
 * Contains a set of components of the multi filter.
 */
public class MultiFilterComponentSet1 extends MultiFilterComponentSet
{
	private final MultiFilterComponent component0;

	public MultiFilterComponentSet1(MultiFilterComponent[] components)
	{
		this.component0 = components[0];
	}

	public int validate(final PreprocessedPeakResult peak)
	{
		//@formatter:off
		if (component0.fail(peak)) return component0.getType();
		//@formatter:oon
		return 0;
	}
}