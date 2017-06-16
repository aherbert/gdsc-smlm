package gdsc.smlm.results.procedures;

import gdsc.core.data.DataException;
import gdsc.smlm.data.config.SMLMSettings.DistanceUnit;
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
//@formatter:off
public class WidthResultProcedure extends UnitResultProcedure implements 
	WResultProcedure,
	WxWyResultProcedure
//@formatter:on
{
	/** The x width. */
	public float[] wx;

	/** The y width. */
	public float[] wy;

	/**
	 * Instantiates a new width result procedure.
	 *
	 * @param results
	 *            the results
	 * @param distanceUnit
	 *            the distance unit
	 */
	public WidthResultProcedure(MemoryPeakResults results, DistanceUnit distanceUnit)
	{
		super(results, distanceUnit);
	}

	/**
	 * Instantiates a new width result procedure.
	 *
	 * @param results
	 *            the results
	 */
	public WidthResultProcedure(MemoryPeakResults results)
	{
		super(results);
	}

	/**
	 * Gets the W data in the configured units.
	 * 
	 * @throws DataException
	 *             if conversion to the required units is not possible
	 */
	public void getW() throws DataException
	{
		i = 0;
		this.wx = allocate(this.wx);
		results.forEach(getDistanceUnit(), (WResultProcedure) this);
	}

	public void executeW(float w)
	{
		this.wx[i++] = w;
	}

	/**
	 * Gets the WxWy data in the configured units.
	 * 
	 * @throws DataException
	 *             if conversion to the required units is not possible
	 */
	public void getWxWy() throws DataException
	{
		i = 0;
		this.wx = allocate(this.wx);
		this.wy = allocate(this.wy);
		results.forEach(getDistanceUnit(), (WxWyResultProcedure) this);
	}

	public void executeWxWy(float wx, float wy)
	{
		this.wx[i] = wx;
		this.wy[i] = wy;
		i++;
	}
}