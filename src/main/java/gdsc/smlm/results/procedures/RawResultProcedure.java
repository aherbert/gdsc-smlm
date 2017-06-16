package gdsc.smlm.results.procedures;

import gdsc.core.data.DataException;
import gdsc.smlm.results.MemoryPeakResults;

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
 * Contains functionality to obtain the standard data for results.
 */
//@formatter:off
public class RawResultProcedure extends AbstractResultProcedure implements 
        BIXYZResultProcedure 
//@formatter:on
{
	/** The background. */
	public float[] background;

	/** The intensity. */
	public float[] intensity;

	/** The x. */
	public float[] x;

	/** The y. */
	public float[] y;

	/** The z. */
	public float[] z;

	/**
	 * Instantiates a new standard result procedure.
	 *
	 * @param results            the results
	 */
	public RawResultProcedure(MemoryPeakResults results)
	{
		super(results);
	}

	/**
	 * Gets the BIXYZ data in the configured units.
	 * 
	 * @throws DataException
	 *             if conversion to the required units is not possible
	 */
	public void getBIXYZ() throws DataException
	{
		i = 0;
		this.background = allocate(this.background);
		this.intensity = allocate(this.intensity);
		this.x = allocate(this.x);
		this.y = allocate(this.y);
		this.z = allocate(this.z);
		results.forEachNative((BIXYZResultProcedure) this);
	}

	public void executeBIXYZ(float background, float intensity, float x, float y, float z)
	{
		this.background[i] = background;
		this.intensity[i] = intensity;
		this.x[i] = x;
		this.y[i] = y;
		this.z[i] = z;
		i++;
	}
}