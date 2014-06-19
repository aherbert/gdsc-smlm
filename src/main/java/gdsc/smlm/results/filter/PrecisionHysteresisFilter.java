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
	float lowerPrecision;
	@XStreamAsAttribute
	float range;
	@XStreamOmitField
	float upperPrecision;
	@XStreamOmitField
	double nmPerPixel;
	@XStreamOmitField
	double gain;

	public PrecisionHysteresisFilter(double searchDistance, float lowerPrecision, float range)
	{
		super(searchDistance);
		this.lowerPrecision = lowerPrecision;
		this.range = Math.abs(range);
	}

	@Override
	protected String generateName()
	{
		return String.format("Precision Hysteresis %.2f +%.2f (@%.2fx)", lowerPrecision, range, searchDistance);
	}
	
	@Override
	protected String generateType()
	{
		return "Precision Hysteresis";
	}
	
	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		upperPrecision = lowerPrecision + range;
		nmPerPixel = peakResults.getNmPerPixel();
		gain = peakResults.getGain();
		super.setup(peakResults);
	}

	@Override
	protected PeakStatus getStatus(PeakResult result)
	{
		double p = result.getPrecision(nmPerPixel, gain);
		if (p <= lowerPrecision)
			return PeakStatus.OK;
		else if (p <= upperPrecision)
			return PeakStatus.Candidate;
		return PeakStatus.Reject;
	}

	@Override
	public double getNumericalValue()
	{
		return lowerPrecision;
	}

	@Override
	public String getNumericalValueName()
	{
		return "Precision + " + (upperPrecision - lowerPrecision);
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
		return 3;
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
		switch (index)
		{
			case 0:
				return searchDistance;
			case 1:
				return lowerPrecision;
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
		switch (index)
		{
			case 0:
				return "Search distance";
			case 1:
				return "Lower Precision";
			default:
				return "Range";
		}
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
		switch (index)
		{
			case 0:
				return new PrecisionHysteresisFilter(updateParameter(searchDistance, delta), lowerPrecision, range);
			case 1:
				return new PrecisionHysteresisFilter(searchDistance, updateParameter(lowerPrecision, delta), range);
			default:
				return new PrecisionHysteresisFilter(searchDistance, lowerPrecision, updateParameter(range, delta));
		}
	}
}