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
	 * @see gdsc.smlm.filters.CircularFilter#getValue(double, float)
	 */
	@Override
	protected float getValue(double sum, float nPoints)
	{
		return (float) (sum / nPoints);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.CircularFilter#computeWeightedNPoints(int, int, double)
	 */
	@Override
	protected float[] computeWeightedNPoints(int maxx, int maxy, double radius)
	{
		float[] nPoints = weights.clone();
		CircularSumFilter sum = new CircularSumFilter();
		sum.convolve(nPoints, maxx, maxy, radius);
		return nPoints;
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