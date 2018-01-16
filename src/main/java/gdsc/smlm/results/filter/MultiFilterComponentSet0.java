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
public class MultiFilterComponentSet0 extends MultiFilterComponentSet
{
	public MultiFilterComponentSet0(MultiFilterComponent[] components)
	{
	}

	@Override
	public int getValidationFlags()
	{
		return 0;
	}
	
	@Override
	public int validate(final PreprocessedPeakResult peak)
	{
		return 0;
	}
	
	@Override
	void replace0(MultiFilterComponent c)
	{
	}
}