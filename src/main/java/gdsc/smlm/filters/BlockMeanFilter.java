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
 */
public class BlockMeanFilter extends BlockFilter
{
	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.BlockFilter#computeNormaliser(float)
	 */
	@Override
	protected Normaliser computeWeightedNormaliser(float n)
	{
		float[] divisor = weights.clone();

		// Use a sum filter to get the sum of the weights in the region
		BlockSumFilter sum = new BlockSumFilter();
		if ((int) n == n)
		{
			sum.rollingBlockFilter(divisor, weightWidth, weightHeight, (int) n);
			//sum.blockFilter(divisor, weightWidth, weightHeight, (int) n);
		}
		else
		{
			sum.stripedBlockFilter(divisor, weightWidth, weightHeight, n);
			//sum.blockFilter(divisor, weightWidth, weightHeight, n);
		}
		return new PerPixelNormaliser(divisor);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.BlockFilter#computeNormaliser(float)
	 */
	@Override
	protected Normaliser computeNormaliser(float n)
	{
		return new FixedNormaliser(pow2(2 * n + 1));
	}

	private float pow2(float f)
	{
		return f * f;
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