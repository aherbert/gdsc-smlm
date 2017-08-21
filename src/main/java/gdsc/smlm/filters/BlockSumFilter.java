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
	private static NonNormaliser normaliser = new NonNormaliser();

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.BlockFilter#computeNormaliser(float)
	 */
	@Override
	protected Normaliser computeWeightedNormaliser(float n)
	{
		return normaliser;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.filters.BlockFilter#computeNormaliser(float)
	 */
	@Override
	protected Normaliser computeNormaliser(float n)
	{
		return normaliser;
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