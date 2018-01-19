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
 * Contains functionality to obtain the localisation precision for results.
 */
//@formatter:off
public class PrecisionResultProcedure extends AbstractResultProcedure implements 
	StoredPrecisionProcedure, 
	LSEPrecisionProcedure, 
	LSEPrecisionBProcedure, 
	MLEPrecisionProcedure, 
	MLEPrecisionBProcedure
//@formatter:on
{
	/** The precision. */
	public double[] precision;

	/**
	 * Instantiates a new precision result procedure.
	 *
	 * @param results
	 *            the results
	 */
	public PrecisionResultProcedure(MemoryPeakResults results)
	{
		super(results);
	}

	/**
	 * Gets the precision for the results.
	 * <p>
	 * If the results contain stored precision values then these are used. Otherwise an attempt is made to compute the
	 * precision using {@link #getLSEPrecision()}. If no exception is thrown then the precision has been computed.
	 * 
	 * @throws DataException
	 *             if conversion to the required units for precision is not possible
	 */
	public void getPrecision()
	{
		if (results.hasPrecision())
			getStoredPrecision();
		else
			getLSEPrecision();
	}

	/**
	 * Gets the precision stored in the results
	 */
	public void getStoredPrecision()
	{
		i = 0;
		precision = allocate(precision);
		results.forEach((StoredPrecisionProcedure) this);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.procedures.StoredPrecisionProcedure#executeStoredPrecision(double)
	 */
	public void executeStoredPrecision(double precision)
	{
		this.precision[i++] = precision;
	}

	/**
	 * Gets the precision assuming a Gaussian 2D PSF and a Least Squares Estimator and a local noise estimate.
	 * 
	 * @throws DataException
	 *             if conversion to the required units for precision is not possible
	 */
	public void getLSEPrecision() throws DataException
	{
		i = 0;
		precision = allocate(precision);
		results.forEach((LSEPrecisionProcedure) this);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.procedures.LSEPrecisionProcedure#executeLSEPrecision(double)
	 */
	public void executeLSEPrecision(double precision)
	{
		this.precision[i++] = precision;
	}

	/**
	 * Gets the precision assuming a Gaussian 2D PSF and a Least Squares Estimator and a local background estimate.
	 * 
	 * @throws DataException
	 *             if conversion to the required units for precision is not possible
	 */
	public void getLSEPrecisionB() throws DataException
	{
		i = 0;
		precision = allocate(precision);
		results.forEach((LSEPrecisionBProcedure) this);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.procedures.LSEPrecisionBProcedure#executeLSEPrecisionB(double)
	 */
	public void executeLSEPrecisionB(double precision)
	{
		this.precision[i++] = precision;
	}

	/**
	 * Gets the precision assuming a Gaussian 2D PSF and a Maximum Likelihood Estimator and a local noise estimate.
	 * 
	 * @throws DataException
	 *             if conversion to the required units for precision is not possible
	 */
	public void getMLEPrecision() throws DataException
	{
		i = 0;
		precision = allocate(precision);
		results.forEach((MLEPrecisionProcedure) this);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.procedures.MLEPrecisionProcedure#executeMLEPrecision(double)
	 */
	public void executeMLEPrecision(double precision)
	{
		this.precision[i++] = precision;
	}

	/**
	 * Gets the precision assuming a Gaussian 2D PSF and a Maximum Likelihood Estimator and a local background estimate.
	 * 
	 * @throws DataException
	 *             if conversion to the required units for precision is not possible
	 */
	public void getMLEPrecisionB() throws DataException
	{
		i = 0;
		precision = allocate(precision);
		results.forEach((MLEPrecisionBProcedure) this);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.procedures.MLEPrecisionBProcedure#executeMLEPrecisionB(double)
	 */
	public void executeMLEPrecisionB(double precision)
	{
		this.precision[i++] = precision;
	}
}