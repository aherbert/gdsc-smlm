package gdsc.smlm.results.procedures;

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
public class PrecisionResultProcedure extends AbstractResultProcedure implements LSEPrecisionProcedure
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
	 * Gets the precision.
	 *
	 * @return the precision
	 * @throws ConversionException
	 *             if conversion to the required units for precision is not possible
	 */
	public void getPrecision()
	{
		i = 0;
		precision = allocate(precision);
		results.forEach(this);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.procedures.LSEPrecisionProcedure#getLSEPrecision(double)
	 */
	public void executeLSEPrecision(double precision)
	{
		this.precision[i++] = precision;
	}
}