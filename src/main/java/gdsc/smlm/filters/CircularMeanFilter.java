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
 * Computes the mean using a circular mask.
 * <p>
 * Adapted from ij.plugin.filter.RankFilters
 */
public class CircularMeanFilter extends CircularFilter
{
	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.CircularFilter#computeWeightedNormaliser(double)
	 */
	@Override
	protected Normaliser computeWeightedNormaliser(double radius)
	{
		float[] nPoints = weights.clone();
		CircularSumFilter sum = new CircularSumFilter();
		sum.convolve(nPoints, weightWidth, weightHeight, radius);
		return new PerPixelNormaliser(nPoints);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.CircularFilter#computeNormaliser(int)
	 */
	@Override
	protected Normaliser computeNormaliser(int nPoints)
	{
		return new FixedNormaliser(nPoints);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public CircularMeanFilter clone()
	{
		CircularMeanFilter o = (CircularMeanFilter) super.clone();
		return o;
	}
}