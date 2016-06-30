package gdsc.smlm.results.filter;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

/**
 * Filter results using a precision threshold. Any results below the lower precision limit are included. Any
 * results above the upper precision limit are excluded. Any results between the limits are included only if they can be
 * traced through time, optionally via other candidates, to a valid result.
 */
public class PrecisionHysteresisFilter extends HysteresisFilter
{
	@XStreamAsAttribute
	final double strictPrecision;
	@XStreamAsAttribute
	final double range;
	@XStreamOmitField
	double lowerVariance;
	@XStreamOmitField
	double upperVariance;
	@XStreamOmitField
	double nmPerPixel;
	@XStreamOmitField
	double gain;
	@XStreamOmitField
	boolean emCCD = true;

	/**
	 * @param searchDistance
	 * @param searchDistanceMode
	 *            0 = relative to the precision of the candidates; 1 = Absolute (in nm)
	 * @param timeThreshold
	 * @param timeThresholdMode
	 *            0 = frames; 1 = seconds
	 * @param strictPrecision
	 * @param range
	 */
	public PrecisionHysteresisFilter(double searchDistance, int searchDistanceMode, double timeThreshold,
			int timeThresholdMode, double strictPrecision, double range)
	{
		super(searchDistance, searchDistanceMode, timeThreshold, timeThresholdMode);
		this.strictPrecision = Math.max(0, strictPrecision);
		this.range = Math.max(0, range);
	}

	@Override
	protected String generateName()
	{
		return String.format("Precision Hysteresis %.2f +%.2f (%s)", strictPrecision, range, getTraceParameters());
	}

	@Override
	protected String generateType()
	{
		return "Precision Hysteresis";
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		lowerVariance = Filter.getDUpperSquaredLimit(strictPrecision);
		upperVariance = Filter.getDUpperSquaredLimit(strictPrecision + range);
		nmPerPixel = peakResults.getNmPerPixel();
		gain = peakResults.getGain();
		emCCD = peakResults.isEMCCD();
		super.setup(peakResults);
	}

	@Override
	protected PeakStatus getStatus(PeakResult result)
	{
		// Use the background noise to estimate precision 
		final double variance = result.getVariance(nmPerPixel, gain, emCCD);
		if (variance <= lowerVariance)
			return PeakStatus.OK;
		else if (variance <= upperVariance)
			return PeakStatus.CANDIDATE;
		return PeakStatus.REJECT;
	}

	@Override
	public double getNumericalValue()
	{
		return strictPrecision;
	}

	@Override
	public String getNumericalValueName()
	{
		return "Precision +" + range;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getDescription()
	 */
	@Override
	public String getDescription()
	{
		return "Filter results using a precision threshold. Any results below the lower precision " +
				"limit are included. Any results above the upper precision limit are excluded. " +
				super.getDescription();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getNumberOfParameters()
	 */
	@Override
	public int getNumberOfParameters()
	{
		return 2 + super.getNumberOfParameters();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getParameterValue(int)
	 */
	@Override
	public double getParameterValue(int index)
	{
		checkIndex(index);
		if (index < super.getNumberOfParameters())
			return super.getParameterValue(index);
		index -= super.getNumberOfParameters();
		switch (index)
		{
			case 0:
				return strictPrecision;
			default:
				return range;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getParameterName(int)
	 */
	@Override
	public String getParameterName(int index)
	{
		checkIndex(index);
		if (index < super.getNumberOfParameters())
			return super.getParameterName(index);
		index -= super.getNumberOfParameters();
		switch (index)
		{
			case 0:
				return "Strict Precision";
			default:
				return "Range";
		}
	}

	static double[] defaultRange = new double[] { 0, 0, 0, 0, PrecisionFilter.DEFAULT_RANGE,
			PrecisionFilter.DEFAULT_RANGE };

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#adjustParameter(int, double)
	 */
	@Override
	public Filter adjustParameter(int index, double delta)
	{
		checkIndex(index);
		// No adjustment of the mode parameters
		if (index == 1 || index == 3)
			return this;
		double[] parameters = new double[] { searchDistance, searchDistanceMode, timeThreshold, timeThresholdMode,
				strictPrecision, range };
		if (index == 0)
			parameters[0] = updateParameter(parameters[0], delta, getDefaultSearchRange());
		else if (index == 2)
			parameters[2] = updateParameter(parameters[2], delta, getDefaultTimeRange());
		else
			parameters[index] = updateParameter(parameters[index], delta, defaultRange[index]);
		return create(parameters);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#create(double[])
	 */
	@Override
	public Filter create(double... parameters)
	{
		return new PrecisionHysteresisFilter(parameters[0], (int) parameters[1], parameters[2], (int) parameters[3],
				parameters[4], parameters[5]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#weakestParameters(double[])
	 */
	@Override
	public void weakestParameters(double[] parameters)
	{
		super.weakestParameters(parameters);

		// Hysteresis filters require all the potential candidates, so disable hysteresis above the candidate threshold  
		setMax(parameters, 4, strictPrecision + range);
		parameters[5] = 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#length()
	 */
	public int length()
	{
		return 4;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#upperLimit()
	 */
	@Override
	public double[] upperLimit()
	{
		return new double[] { Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, PrecisionFilter.UPPER_LIMIT,
				PrecisionFilter.UPPER_LIMIT };
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#sequence()
	 */
	public double[] sequence()
	{
		// Ignore the mode parameters
		return new double[] { searchDistance, timeThreshold, strictPrecision, range };
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#mutationStepRange()
	 */
	public double[] mutationStepRange()
	{
		return new double[] { getDefaultSearchRange(), getDefaultTimeRange(), PrecisionFilter.DEFAULT_RANGE,
				PrecisionFilter.DEFAULT_RANGE };
	}
}