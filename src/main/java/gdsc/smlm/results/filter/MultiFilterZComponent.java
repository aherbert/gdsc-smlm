package gdsc.smlm.results.filter;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

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
	public int getType()
	{
		return IDirectFilter.V_Z;
	}
}