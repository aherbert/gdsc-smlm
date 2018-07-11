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
package gdsc.smlm.results.filter;

/**
 * Filter results using z-depth
 */
public class MultiFilterZComponent extends MultiFilterComponent
{
	final float minZ, maxZ;

	public MultiFilterZComponent(double minZ, double maxZ)
	{
		this.minZ = (float) minZ;
		this.maxZ = (float) maxZ;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.MultiFilterComponent#fail(gdsc.smlm.results.filter.PreprocessedPeakResult)
	 */
	@Override
	public boolean fail(final PreprocessedPeakResult peak)
	{
		final float z = peak.getZ();
		return (z > maxZ || z < minZ);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.MultiFilterComponent#getType()
	 */
	@Override
	public int getType()
	{
		return IDirectFilter.V_Z;
	}
}
