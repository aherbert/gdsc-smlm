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
 * Computes the sum using a square block mask.
 * <p>
 * Adapted from ij.plugin.filter.RankFilters
 */
public class BlockSumFilter extends BlockFilter
{
	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.BlockFilter#getValue(float, float)
	 */
	@Override
	protected float getValue(float sum, float divisor)
	{
		return sum;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.BlockFilter#computeWeightedDivisor(int, int, float)
	 */
	@Override
	protected float[] computeWeightedDivisor(int maxx, int maxy, float n)
	{
		// The concept of a divisor for a weighted sum filter is invalid.
		
		// To avoid array index exceptions we create a empty divisor.
		//float[] divisor = new float[weights.length];
		//Arrays.fill(divisor, 1.0f);
		//return divisor;

		// Since the divisor will not be used in getValue() just return the weights
		return weights;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public BlockSumFilter clone()
	{
		BlockSumFilter o = (BlockSumFilter) super.clone();
		return o;
	}
}