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

import gdsc.smlm.fitting.function.Gaussian2DFunction;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

/**
 * Filter results using a precision threshold
 */
public class PrecisionFilter extends Filter
{
	@XStreamAsAttribute
	float precision;
	@XStreamOmitField
	double nmPerPixel = 100;
	@XStreamOmitField
	boolean emCCD = true;
	@XStreamOmitField
	double bias = 0;
	@XStreamOmitField
	double gain = 1;

	public PrecisionFilter(float precision)
	{
		this.precision = precision;
	}

	@Override
	protected String generateName()
	{
		return "Precision " + precision;
	}

	@Override
	protected String generateType()
	{
		return "Precision";
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		nmPerPixel = peakResults.getNmPerPixel();
		gain = peakResults.getGain();
		emCCD = peakResults.isEMCCD();
		if (peakResults.getCalibration() != null)
		{
			bias = peakResults.getCalibration().bias;
		}
	}

	@Override
	public boolean accept(PeakResult peak)
	{
		if (bias != 0)
		{
			// Use the estimated background for the peak
			final double s = nmPerPixel * peak.getSD();
			final double N = peak.getSignal();
			return PeakResult.getPrecisionX(nmPerPixel, s, N,
					Math.max(0, peak.params[Gaussian2DFunction.BACKGROUND] - bias) / gain, emCCD) <= precision;
		}
		// Use the background noise to estimate precision 
		return peak.getPrecision(nmPerPixel, gain, emCCD) <= precision;
	}

	@Override
	public double getNumericalValue()
	{
		return precision;
	}

	@Override
	public String getNumericalValueName()
	{
		return "Precision";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getDescription()
	 */
	@Override
	public String getDescription()
	{
		return "Filter results using an upper precision threshold.";
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
	 * @see gdsc.smlm.results.filter.Filter#getParameterValue(int)
	 */
	@Override
	public double getParameterValue(int index)
	{
		checkIndex(index);
		return precision;
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
		return "Precision";
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
		return new PrecisionFilter(updateParameter(precision, delta));
	}
}