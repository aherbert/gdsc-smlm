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
	 * Gets the precision assuming a Least Squares Estimator and a local noise estimate.
	 *
	 * @return the precision
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
	 * Gets the precision assuming a Least Squares Estimator and a local background estimate.
	 *
	 * @return the precision
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
	 * Gets the precision assuming a Maximum Likelihood Estimator and a local noise estimate.
	 *
	 * @return the precision
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
	 * Gets the precision assuming a Maximum Likelihood Estimator and a local background estimate.
	 *
	 * @return the precision
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