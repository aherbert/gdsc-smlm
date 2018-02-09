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
	private MultiFilterComponent[] components;

	public MultiFilterComponentSetDefault(MultiFilterComponent[] components)
	{
		this.components = components;
	}

	@Override
	public int getValidationFlags()
	{
		int flags = 0;
		for (int i = 0; i < components.length; i++)
			flags |= components[i].getType();
		return flags;
	}

	@Override
	public int validate(final PreprocessedPeakResult peak)
	{
		for (int i = 0; i < components.length; i++)
			if (components[i].fail(peak))
				return components[i].getType();
		return 0;
	}

	@Override
	void replace0(MultiFilterComponent c)
	{
		if (components.length > 0)
			components[0] = c;
	}
	
	@Override
	public MultiFilterComponentSet clone()
	{
		// Copy the array
		MultiFilterComponent[] c = new MultiFilterComponent[components.length];
		if (c.length > 0)
			System.arraycopy(components, 0, c, 0, c.length);
		return new MultiFilterComponentSetDefault(c);
	}
}