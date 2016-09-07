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

import gdsc.smlm.results.ClassifiedPeakResult;

/**
 * Filter results using the combination of two filters. Results can pass either filter
 */
public class OrFilter extends CombinedFilter
{
	public OrFilter(Filter filter1, Filter filter2)
	{
		super(filter1, filter2);
	}

	@Override
	protected String getOperator()
	{
		return "||";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#accept(gdsc.smlm.results.PeakResult)
	 */
	@Override
	public boolean accept(ClassifiedPeakResult peak)
	{
		return filter1.accept(peak) || filter2.accept(peak);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getDescription()
	 */
	@Override
	public String getDescription()
	{
		return "Filter results using the combination of two filters. Results can pass either filter.";
	}
	
	@Override
	protected Filter createFilter(Filter f1, Filter f2)
	{
		return new OrFilter(f1, f2);
	}
}