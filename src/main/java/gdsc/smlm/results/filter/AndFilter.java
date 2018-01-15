package gdsc.smlm.results.filter;

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

import gdsc.smlm.results.PeakResult;

/**
 * Filter results using the combination of two filters. Results must pass both filters
 */
public class AndFilter extends CombinedFilter
{
	public AndFilter(Filter filter1, Filter filter2)
	{
		super(filter1, filter2);
	}

	@Override
	protected String getOperator()
	{
		return "&&";
	}

	@Override
	public boolean accept(PeakResult peak)
	{
		return accept1(peak) && accept2(peak);
	}

	public int validate(final PreprocessedPeakResult peak)
	{
		if (accept1(peak) && accept2(peak))
			return 0;
		// We get here if filter 1 failed; or filter 1 passed but filter 2 failed.
		return (result1 == 0) ? result2 : result1;
	}

	@Override
	protected Filter createFilter(Filter f1, Filter f2)
	{
		return new AndFilter(f1, f2);
	}

	@Override
	public Filter clone()
	{
		return new AndFilter(filter1.clone(), filter2.clone());
	}
}