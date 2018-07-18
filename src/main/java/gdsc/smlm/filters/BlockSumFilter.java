/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm.filters;

/**
 * Computes the sum using a square block mask.
 */
public class BlockSumFilter extends BlockFilter
{
	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.filters.BlockFilter#computeNormaliser(float)
	 */
	@Override
	protected Normaliser computeWeightedNormaliser(float n)
	{
		final float[] divisor = weights.clone();

		// Use a mean filter to get the mean of the weights in the region
		final BlockMeanFilter sum = new BlockMeanFilter();
		if ((int) n == n)
			sum.rollingBlockFilter(divisor, weightWidth, weightHeight, (int) n);
			//sum.blockFilter(divisor, weightWidth, weightHeight, (int) n);
		else
			sum.stripedBlockFilter(divisor, weightWidth, weightHeight, n);
			//sum.blockFilter(divisor, weightWidth, weightHeight, n);
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
		return NonNormaliser.INSTANCE;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see java.lang.Object#clone()
	 */
	@Override
	public BlockSumFilter clone()
	{
		final BlockSumFilter o = (BlockSumFilter) super.clone();
		return o;
	}
}
