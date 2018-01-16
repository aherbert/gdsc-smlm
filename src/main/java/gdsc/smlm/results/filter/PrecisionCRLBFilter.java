package gdsc.smlm.results.filter;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

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

import gdsc.smlm.results.PeakResult;

/**
 * Filter results using a precision threshold. Calculates the precision using the Cramér-Rao lower bound (CRLB) of the
 * variance of estimators of the fit parameter. The variance for the fitted X and Y position is averaged to produce a
 * localisation precision.
 */
public class PrecisionCRLBFilter extends DirectFilter implements IMultiFilter
{
	@XStreamAsAttribute
	final double precision;
	@XStreamOmitField
	float variance;

	public PrecisionCRLBFilter(double precision)
	{
		this.precision = Math.max(0, precision);
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		// Add the 2-fold scale factor here:
		// (varX + varY)/2 < precision^2
		// (varX + varY) < precision^2 * 2
		variance = Filter.getUpperSquaredLimit(precision) * 2f;
	}

	@Override
	public boolean accept(PeakResult peak)
	{
		// Use the estimated parameter deviations for the peak
		if (peak.hasParameterDeviations())
		{
			float vx = peak.getParameterDeviation(PeakResult.X);
			float vy = peak.getParameterDeviation(PeakResult.Y);
			return (vx * vx + vy * vy) <= variance;
		}
		return true;
	}

	public int getValidationFlags()
	{
		return V_LOCATION_VARIANCE_CRLB;
	}

	@Override
	public int validate(final PreprocessedPeakResult peak)
	{
		if (peak.getLocationVarianceCRLB() > variance)
			return V_LOCATION_VARIANCE_CRLB;
		return 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getDescription()
	 */
	@Override
	public String getDescription()
	{
		return "Filter results using an upper precision threshold (uses fitted parameter variance).";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#requiresParameterDeviations()
	 */
	@Override
	public boolean requiresParameterDeviations()
	{
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getNumberOfParameters()
	 */
	@Override
	public int getNumberOfParameters()
	{
		return 1;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getParameterValueInternal(int)
	 */
	@Override
	protected double getParameterValueInternal(int index)
	{
		return precision;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getParameterIncrement(int)
	 */
	@Override
	public double getParameterIncrement(int index)
	{
		checkIndex(index);
		return PrecisionFilter.DEFAULT_INCREMENT;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getParameterType(int)
	 */
	@Override
	public ParameterType getParameterType(int index)
	{
		checkIndex(index);
		return ParameterType.PRECISION_CRLB;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#adjustParameter(int, double)
	 */
	@Override
	public Filter adjustParameter(int index, double delta)
	{
		checkIndex(index);
		return new PrecisionCRLBFilter(updateParameter(precision, delta, PrecisionFilter.DEFAULT_RANGE));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#create(double[])
	 */
	@Override
	public Filter create(double... parameters)
	{
		return new PrecisionCRLBFilter(parameters[0]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#weakestParameters(double[])
	 */
	@Override
	public void weakestParameters(double[] parameters)
	{
		setMax(parameters, 0, precision);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.DirectFilter#lowerBoundOrientation(int)
	 */
	@Override
	public int lowerBoundOrientation(int index)
	{
		return 1;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#upperLimit()
	 */
	@Override
	public double[] upperLimit()
	{
		return new double[] { PrecisionFilter.UPPER_LIMIT };
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#mutationStepRange()
	 */
	public double[] mutationStepRange()
	{
		return new double[] { PrecisionFilter.DEFAULT_RANGE };
	}

	public double getSignal()
	{
		return 0;
	}

	public double getSNR()
	{
		return 0;
	}

	public double getMinWidth()
	{
		return 0;
	}

	public double getMaxWidth()
	{
		return 0;
	}

	public double getShift()
	{
		return 0;
	}

	public double getEShift()
	{
		return 0;
	}

	public double getPrecision()
	{
		return precision;
	}

	public PrecisionType getPrecisionType()
	{
		return PrecisionType.CRLB;
	}
}