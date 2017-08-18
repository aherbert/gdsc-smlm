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
 * Computes the mean using a square block mask.
 * <p>
 * Adapted from ij.plugin.filter.RankFilters
 */
public class BlockMeanFilter extends BlockFilter
{
	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.BlockFilter#getValue(float, float)
	 */
	@Override
	protected float getValue(float sum, float divisor)
	{
		return sum / divisor;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.BlockFilter#computeWeightedDivisor(int, int, float)
	 */
	@Override
	protected float[] computeWeightedDivisor(final int maxx, final int maxy, final float n)
	{
		float[] divisor = weights.clone();
		
		// Use a sum filter to get the sum of the weights in the region
		BlockSumFilter sum = new BlockSumFilter();
		if ((int) n == n)
		{
			sum.rollingBlockFilter(divisor, maxx, maxy, (int) n);
		}
		else
		{
			sum.stripedBlockFilter(divisor, maxx, maxy, n);
		}
		return divisor;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public BlockMeanFilter clone()
	{
		BlockMeanFilter o = (BlockMeanFilter) super.clone();
		return o;
	}
}