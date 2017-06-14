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
 * Contains core functionality to for result procedures.
 */
public abstract class AbstractResultProcedure
{
	/** The results. */
	final MemoryPeakResults results;

	/** The distance unit. */
	private DistanceUnit distanceUnit;

	/** The intensity unit. */
	private IntensityUnit intensityUnit;

	/** The counter for procedures */
	protected int i;

	/**
	 * Instantiates a new abstract result procedure.
	 *
	 * @param results
	 *            the results
	 * @param distanceUnit
	 *            the distance unit
	 * @param intensityUnit
	 *            the intensity unit
	 */
	public AbstractResultProcedure(MemoryPeakResults results, DistanceUnit distanceUnit, IntensityUnit intensityUnit)
	{
		if (results == null)
			throw new IllegalArgumentException("results must not be null");
		this.results = results;
		this.setDistanceUnit(distanceUnit);
		this.setIntensityUnit(intensityUnit);
	}

	/**
	 * Instantiates a new abstract result procedure.
	 *
	 * @param results
	 *            the results
	 * @param distanceUnit
	 *            the distance unit
	 */
	public AbstractResultProcedure(MemoryPeakResults results, DistanceUnit distanceUnit)
	{
		this(results, distanceUnit, IntensityUnit.PHOTON);
	}

	/**
	 * Instantiates a new abstract result procedure.
	 *
	 * @param results
	 *            the results
	 * @param intensityUnit
	 *            the intensity unit
	 */
	public AbstractResultProcedure(MemoryPeakResults results, IntensityUnit intensityUnit)
	{
		this(results, DistanceUnit.PIXEL, intensityUnit);
	}

	/**
	 * Instantiates a new abstract result procedure.
	 *
	 * @param results
	 *            the results
	 */
	public AbstractResultProcedure(MemoryPeakResults results)
	{
		this(results, DistanceUnit.PIXEL, IntensityUnit.PHOTON);
	}

	/**
	 * Gets the distance unit.
	 *
	 * @return the distance unit
	 */
	public DistanceUnit getDistanceUnit()
	{
		return distanceUnit;
	}

	/**
	 * Get the size of the results.
	 *
	 * @return the size
	 */
	public int size()
	{
		return results.size();
	}
	
	/**
	 * Sets the distance unit.
	 *
	 * @param distanceUnit
	 *            the new distance unit
	 */
	public void setDistanceUnit(DistanceUnit distanceUnit)
	{
		if (distanceUnit == null)
			throw new IllegalArgumentException("unit must not be null");
		this.distanceUnit = distanceUnit;
	}

	/**
	 * Gets the intensity unit.
	 *
	 * @return the intensity unit
	 */
	public IntensityUnit getIntensityUnit()
	{
		return intensityUnit;
	}

	/**
	 * Sets the intensity unit.
	 *
	 * @param intensityUnit
	 *            the new intensity unit
	 */
	public void setIntensityUnit(IntensityUnit intensityUnit)
	{
		if (intensityUnit == null)
			throw new IllegalArgumentException("unit must not be null");
		this.intensityUnit = intensityUnit;
	}

	/**
	 * Allocate the array to store the data based on the current results size.
	 *
	 * @param data
	 *            the data
	 * @return the float[]
	 */
	protected float[] allocate(float[] data)
	{
		return (data == null || data.length < size()) ? new float[size()] : data;
	}

	/**
	 * Allocate the array to store the data based on the current results size.
	 *
	 * @param data
	 *            the data
	 * @return the float[]
	 */
	protected double[] allocate(double[] data)
	{
		return (data == null || data.length < size()) ? new double[size()] : data;
	}
}