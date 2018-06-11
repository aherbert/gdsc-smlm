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
package gdsc.smlm.results.procedures;

import gdsc.core.data.DataException;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;

/**
 * Contains functionality to obtain the standard calibrated data for results.
 */
//@formatter:off
public class StandardResultProcedure extends UnitResultProcedure implements 
		BResultProcedure, 
        BIXYResultProcedure, 
        BIXYZResultProcedure, 
        IResultProcedure, 
        IXYResultProcedure, 
        IXYRResultProcedure, 
        IXYZResultProcedure,
		TResultProcedure, 
		TXYResultProcedure, 
		XYResultProcedure, 
		XYRResultProcedure, 
		XYZResultProcedure,
		ZResultProcedure
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

	/** The peak results. */
	public PeakResult[] peakResults;

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
	 * @param intensityUnit
	 *            the intensity unit
	 */
	public StandardResultProcedure(MemoryPeakResults results, IntensityUnit intensityUnit, DistanceUnit distanceUnit)
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
	 * Gets the B data in the configured units.
	 * 
	 * @throws DataException
	 *             if conversion to the required units is not possible
	 */
	public void getB() throws DataException
	{
		i = 0;
		allocateB();
		results.forEach(getIntensityUnit(), (BResultProcedure) this);
	}

	@Override
	public void executeB(float background)
	{
		this.background[i++] = background;
	}

	/**
	 * Gets the BIXY data in the configured units.
	 * 
	 * @throws DataException
	 *             if conversion to the required units is not possible
	 */
	public void getBIXY() throws DataException
	{
		i = 0;
		allocateB();
		allocateI();
		allocateX();
		allocateY();
		results.forEach(getIntensityUnit(), getDistanceUnit(), (BIXYResultProcedure) this);
	}

	@Override
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
	 * @throws DataException
	 *             if conversion to the required units is not possible
	 */
	public void getBIXYZ() throws DataException
	{
		i = 0;
		allocateB();
		allocateI();
		allocateX();
		allocateY();
		allocateZ();
		results.forEach(getIntensityUnit(), getDistanceUnit(), (BIXYZResultProcedure) this);
	}

	@Override
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
	 * @throws DataException
	 *             if conversion to the required units is not possible
	 */
	public void getI() throws DataException
	{
		i = 0;
		allocateI();
		results.forEach(getIntensityUnit(), (IResultProcedure) this);
	}

	@Override
	public void executeI(float intensity)
	{
		this.intensity[i] = intensity;
		i++;
	}

	/**
	 * Gets the IXY data in the configured units.
	 * 
	 * @throws DataException
	 *             if conversion to the required units is not possible
	 */
	public void getIXY() throws DataException
	{
		i = 0;
		allocateI();
		allocateX();
		allocateY();
		results.forEach(getIntensityUnit(), getDistanceUnit(), (IXYResultProcedure) this);
	}

	@Override
	public void executeIXY(float intensity, float x, float y)
	{
		this.intensity[i] = intensity;
		this.x[i] = x;
		this.y[i] = y;
		i++;
	}

	/**
	 * Gets the IXYR data in the configured units.
	 * 
	 * @throws DataException
	 *             if conversion to the required units is not possible
	 */
	public void getIXYR() throws DataException
	{
		i = 0;
		allocateI();
		allocateX();
		allocateY();
		allocateR();
		results.forEach(getIntensityUnit(), getDistanceUnit(), (IXYRResultProcedure) this);
	}

	@Override
	public void executeIXYR(float intensity, float x, float y, PeakResult result)
	{
		this.intensity[i] = intensity;
		this.x[i] = x;
		this.y[i] = y;
		peakResults[i] = result;
		i++;
	}

	/**
	 * Gets the IXYZ data in the configured units.
	 * 
	 * @throws DataException
	 *             if conversion to the required units is not possible
	 */
	public void getIXYZ() throws DataException
	{
		i = 0;
		allocateI();
		allocateX();
		allocateY();
		allocateZ();
		results.forEach(getIntensityUnit(), getDistanceUnit(), (IXYZResultProcedure) this);
	}

	@Override
	public void executeIXYZ(float intensity, float x, float y, float z)
	{
		this.intensity[i] = intensity;
		this.x[i] = x;
		this.y[i] = y;
		this.z[i] = z;
		i++;
	}

	/**
	 * Gets the T data in the configured units.
	 */
	public void getT() throws DataException
	{
		i = 0;
		allocateT();
		results.forEach(this);
	}

	@Override
	public void executeT(int frame)
	{
		this.frame[i++] = frame;
	}

	/**
	 * Gets the TXY data in the configured units.
	 * 
	 * @throws DataException
	 *             if conversion to the required units is not possible
	 */
	public void getTXY() throws DataException
	{
		i = 0;
		allocateT();
		allocateX();
		allocateY();
		results.forEach(getDistanceUnit(), (TXYResultProcedure) this);
	}

	@Override
	public void executeTXY(int frame, float x, float y)
	{
		this.frame[i] = frame;
		this.x[i] = x;
		this.y[i] = y;
		i++;
	}

	/**
	 * Gets the XY data in the configured units.
	 * 
	 * @throws DataException
	 *             if conversion to the required units is not possible
	 */
	public void getXY() throws DataException
	{
		i = 0;
		allocateX();
		allocateY();
		results.forEach(getDistanceUnit(), (XYResultProcedure) this);
	}

	@Override
	public void executeXY(float x, float y)
	{
		this.x[i] = x;
		this.y[i] = y;
		i++;
	}

	/**
	 * Gets the XYR data in the configured units.
	 * 
	 * @throws DataException
	 *             if conversion to the required units is not possible
	 */
	public void getXYR() throws DataException
	{
		i = 0;
		allocateX();
		allocateY();
		allocateR();
		results.forEach(getDistanceUnit(), (XYRResultProcedure) this);
	}

	@Override
	public void executeXYR(float x, float y, PeakResult result)
	{
		this.x[i] = x;
		this.y[i] = y;
		peakResults[i] = result;
		i++;
	}

	/**
	 * Gets the XYZ data in the configured units.
	 * 
	 * @throws DataException
	 *             if conversion to the required units is not possible
	 */
	public void getXYZ() throws DataException
	{
		i = 0;
		allocateX();
		allocateY();
		allocateZ();
		results.forEach(getDistanceUnit(), (XYZResultProcedure) this);
	}

	@Override
	public void executeXYZ(float x, float y, float z)
	{
		this.x[i] = x;
		this.y[i] = y;
		this.z[i] = z;
		i++;
	}

	/**
	 * Gets the Z data in the configured units.
	 * 
	 * @throws DataException
	 *             if conversion to the required units is not possible
	 */
	public void getZ() throws DataException
	{
		i = 0;
		allocateZ();
		results.forEach(getDistanceUnit(), (ZResultProcedure) this);
	}

	@Override
	public void executeZ(float z)
	{
		this.z[i++] = z;
	}

	private void allocateT()
	{
		this.frame = allocate(this.frame);
	}

	private void allocateB()
	{
		this.background = allocate(this.background);
	}

	private void allocateI()
	{
		this.intensity = allocate(this.intensity);
	}

	private void allocateX()
	{
		this.x = allocate(this.x);
	}

	private void allocateY()
	{
		this.y = allocate(this.y);
	}

	private void allocateZ()
	{
		this.z = allocate(this.z);
	}

	private void allocateR()
	{
		this.peakResults = allocate(this.peakResults);
	}
}
