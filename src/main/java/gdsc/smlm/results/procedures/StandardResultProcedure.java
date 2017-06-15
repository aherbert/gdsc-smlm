package gdsc.smlm.results.procedures;

import gdsc.smlm.data.config.SMLMSettings.DistanceUnit;
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
//@formatter:off
public class StandardResultProcedure extends AbstractResultProcedure implements 
        BIXYResultProcedure, 
        BIXYZResultProcedure, 
        IResultProcedure, 
        IXYResultProcedure, 
        IXYZResultProcedure,
		TXYResultProcedure, 
		XYResultProcedure, 
		XYZResultProcedure
//@formatter:on
{
	/** The frame. */
	public int[] frame;
	
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
	 * @param results
	 *            the results
	 * @param distanceUnit
	 *            the distance unit
	 * @param intensityUnit
	 *            the intensity unit
	 */
	public StandardResultProcedure(MemoryPeakResults results, DistanceUnit distanceUnit, IntensityUnit intensityUnit)
	{
		super(results, distanceUnit, intensityUnit);
	}

	/**
	 * Instantiates a new standard result procedure.
	 *
	 * @param results
	 *            the results
	 * @param distanceUnit
	 *            the distance unit
	 */
	public StandardResultProcedure(MemoryPeakResults results, DistanceUnit distanceUnit)
	{
		super(results, distanceUnit);
	}

	/**
	 * Instantiates a new standard result procedure.
	 *
	 * @param results
	 *            the results
	 * @param intensityUnit
	 *            the intensity unit
	 */
	public StandardResultProcedure(MemoryPeakResults results, IntensityUnit intensityUnit)
	{
		super(results, intensityUnit);
	}

	/**
	 * Instantiates a new standard result procedure.
	 *
	 * @param results
	 *            the results
	 */
	public StandardResultProcedure(MemoryPeakResults results)
	{
		super(results);
	}

	/**
	 * Gets the BIXY data in the configured units.
	 * 
	 * @throws ConversionException
	 *             if conversion to the required units is not possible
	 */
	public void getBIXY()
	{
		i = 0;
		this.background = allocate(this.background);
		this.intensity = allocate(this.intensity);
		this.x = allocate(this.x);
		this.y = allocate(this.y);
		results.forEach((BIXYResultProcedure) this, getIntensityUnit(), getDistanceUnit());
	}

	public void executeBIXY(float background, float intensity, float x, float y)
	{
		this.background[i] = background;
		this.intensity[i] = intensity;
		this.x[i] = x;
		this.y[i] = y;
		i++;
	}
	
	/**
	 * Gets the BIXYZ data in the configured units.
	 * 
	 * @throws ConversionException
	 *             if conversion to the required units is not possible
	 */
	public void getBIXYZ()
	{
		i = 0;
		this.background = allocate(this.background);
		this.intensity = allocate(this.intensity);
		this.x = allocate(this.x);
		this.y = allocate(this.y);
		this.z = allocate(this.z);
		results.forEach((BIXYZResultProcedure) this, getIntensityUnit(), getDistanceUnit());
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

	/**
	 * Gets the I data in the configured units.
	 * 
	 * @throws ConversionException
	 *             if conversion to the required units is not possible
	 */
	public void getI()
	{
		i = 0;
		this.intensity = allocate(this.intensity);
		results.forEach((IResultProcedure) this, getIntensityUnit());
	}

	public void executeI(float intensity)
	{
		this.intensity[i] = intensity;
		i++;
	}

	/**
	 * Gets the IXY data in the configured units.
	 * 
	 * @throws ConversionException
	 *             if conversion to the required units is not possible
	 */
	public void getIXY()
	{
		i = 0;
		this.intensity = allocate(this.intensity);
		this.x = allocate(this.x);
		this.y = allocate(this.y);
		results.forEach((IXYResultProcedure) this, getIntensityUnit(), getDistanceUnit());
	}

	public void executeIXY(float intensity, float x, float y)
	{
		this.intensity[i] = intensity;
		this.x[i] = x;
		this.y[i] = y;
		i++;
	}

	/**
	 * Gets the IXYZ data in the configured units.
	 * 
	 * @throws ConversionException
	 *             if conversion to the required units is not possible
	 */
	public void getIXYZ()
	{
		i = 0;
		this.intensity = allocate(this.intensity);
		this.x = allocate(this.x);
		this.y = allocate(this.y);
		this.z = allocate(this.z);
		results.forEach((IXYZResultProcedure) this, getIntensityUnit(), getDistanceUnit());
	}

	public void executeIXYZ(float intensity, float x, float y, float z)
	{
		this.intensity[i] = intensity;
		this.x[i] = x;
		this.y[i] = y;
		this.z[i] = z;
		i++;
	}

	/**
	 * Gets the TXY data in the configured units.
	 * 
	 * @throws ConversionException
	 *             if conversion to the required units is not possible
	 */
	public void getTXY()
	{
		i = 0;
		this.frame = allocate(this.frame);
		this.x = allocate(this.x);
		this.y = allocate(this.y);
		results.forEach((TXYResultProcedure) this, getDistanceUnit());
	}

	public void executeTXY(int frame, float x, float y)
	{
		this.frame[i]=frame;
		this.x[i] = x;
		this.y[i] = y;
		i++;
	}

	/**
	 * Gets the XY data in the configured units.
	 * 
	 * @throws ConversionException
	 *             if conversion to the required units is not possible
	 */
	public void getXY()
	{
		i = 0;
		this.x = allocate(this.x);
		this.y = allocate(this.y);
		results.forEach((XYResultProcedure) this, getDistanceUnit());
	}

	public void executeXY(float x, float y)
	{
		this.x[i] = x;
		this.y[i] = y;
		i++;
	}

	/**
	 * Gets the XYZ data in the configured units.
	 * 
	 * @throws ConversionException
	 *             if conversion to the required units is not possible
	 */
	public void getXYZ()
	{
		i = 0;
		this.x = allocate(this.x);
		this.y = allocate(this.y);
		this.z = allocate(this.z);
		results.forEach((XYZResultProcedure) this, getDistanceUnit());
	}

	public void executeXYZ(float x, float y, float z)
	{
		this.x[i] = x;
		this.y[i] = y;
		this.z[i] = z;
		i++;
	}
}