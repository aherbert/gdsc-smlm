package gdsc.smlm.results.procedures;

import gdsc.core.data.DataException;
import gdsc.smlm.data.config.SMLMSettings.IntensityUnit;
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
 * Contains functionality to obtain the standard calibrated data for results.
 */
public class HeightResultProcedure extends UnitResultProcedure implements HResultProcedure
{
	/** The height. */
	public float[] h;

	/**
	 * Instantiates a new width result procedure.
	 *
	 * @param results
	 *            the results
	 * @param intensityUnit
	 *            the intensity unit
	 */
	public HeightResultProcedure(MemoryPeakResults results, IntensityUnit intensityUnit)
	{
		super(results, intensityUnit);
	}

	/**
	 * Instantiates a new width result procedure.
	 *
	 * @param results
	 *            the results
	 */
	public HeightResultProcedure(MemoryPeakResults results)
	{
		super(results);
	}

	/**
	 * Gets the height data in the configured units.
	 * 
	 * @throws DataException
	 *             if conversion to the required units is not possible
	 */
	public void getH() throws DataException
	{
		i = 0;
		this.h = allocate(this.h);
		results.forEach(getIntensityUnit(), this);
	}

	public void executeH(float a)
	{
		this.h[i++] = a;
	}
}