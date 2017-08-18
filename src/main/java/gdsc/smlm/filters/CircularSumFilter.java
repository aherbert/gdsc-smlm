package gdsc.smlm.filters;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Computes the sum using a circular mask.
 * <p>
 * Adapted from ij.plugin.filter.RankFilters
 */
public class CircularSumFilter extends CircularFilter
{
	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.CircularFilter#getValue(double, int)
	 */
	@Override
	protected float getValue(double sum, int nPoints)
	{
		return (float) sum;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public CircularSumFilter clone()
	{
		CircularSumFilter o = (CircularSumFilter) super.clone();
		return o;
	}
}