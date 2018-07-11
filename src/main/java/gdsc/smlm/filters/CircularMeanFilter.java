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
	@Override
	public CircularMeanFilter clone()
	{
		CircularMeanFilter o = (CircularMeanFilter) super.clone();
		return o;
	}
}
