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
public class MultiFilterComponentSetDefault extends MultiFilterComponentSet
{
	private final MultiFilterComponent[] components;

	public MultiFilterComponentSetDefault(MultiFilterComponent[] components)
	{
		this.components = components;
	}

	public int validate(final PreprocessedPeakResult peak)
	{
		for (int i = 0; i < components.length; i++)
			if (components[i].fail(peak))
				return components[i].getType();
		return 0;
	}
}