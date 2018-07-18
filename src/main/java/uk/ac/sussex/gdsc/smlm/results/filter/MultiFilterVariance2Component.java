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
package uk.ac.sussex.gdsc.smlm.results.filter;

/**
 * Filter results using Precision with local background
 */
public class MultiFilterVariance2Component extends MultiFilterComponent
{
	private final double variance;

	/**
	 * Instantiates a new multi filter variance 2 component.
	 *
	 * @param precision
	 *            the precision
	 */
	public MultiFilterVariance2Component(double precision)
	{
		this.variance = Filter.getDUpperSquaredLimit(precision);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.filter.MultiFilterComponent#fail(uk.ac.sussex.gdsc.smlm.results.filter.PreprocessedPeakResult)
	 */
	@Override
	public boolean fail(final PreprocessedPeakResult peak)
	{
		return (peak.getLocationVariance2() > variance);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see uk.ac.sussex.gdsc.smlm.results.filter.MultiFilterComponent#getType()
	 */
	@Override
	public int getType()
	{
		return IDirectFilter.V_LOCATION_VARIANCE2;
	}
}
