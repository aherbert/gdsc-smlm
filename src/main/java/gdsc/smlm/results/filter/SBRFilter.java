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
 * Filter results using a signal-to-background ratio (SBR) threshold.
 * <p>
 * Requires the bias to be configured at or above zero. If the background is below the configured bias or there is no
 * bias then the filter resorts to a signal-to-noise filter. If there is a background level above the bias then this is
 * assumed to be the variance of the photon shot noise and the noise is taken at the square root of the background
 * level.
 */
public class SBRFilter extends Filter
{
	@XStreamAsAttribute
	final float sbr;
	@XStreamOmitField
	double bias = -1;

	public SBRFilter(float sbr)
	{
		this.sbr = Math.max(0, sbr);
	}

	@Override
	protected String generateName()
	{
		return "SBR " + sbr;
	}

	@Override
	protected String generateType()
	{
		return "SBR";
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		bias = -1;
		if (peakResults.getCalibration() != null)
		{
			if (peakResults.getCalibration().bias >= 0)
				bias = peakResults.getCalibration().bias;
		}
	}

	@Override
	public boolean accept(PeakResult peak)
	{
		if (bias != -1)
		{
			final double background = peak.getBackground() - bias;
			if (background > 0)
				return peak.getSignal() / Math.sqrt(background) >= this.sbr;
		}
		return SNRFilter.getSNR(peak) >= this.sbr;
	}

	@Override
	public double getNumericalValue()
	{
		return sbr;
	}

	@Override
	public String getNumericalValueName()
	{
		return "SBR";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getDescription()
	 */
	@Override
	public String getDescription()
	{
		return "Filter results using a lower SBR threshold.";
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
		return sbr;
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
		return "SBR";
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
		return new SBRFilter(updateParameter(sbr, delta, SNRFilter.DEFAULT_RANGE));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#create(double[])
	 */
	@Override
	public Filter create(double... parameters)
	{
		return new SBRFilter((float) parameters[0]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#weakestParameters(double[])
	 */
	@Override
	public void weakestParameters(double[] parameters)
	{
		setMin(parameters, 0, sbr);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#length()
	 */
	@Override
	public int length()
	{
		return 1;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#sequence()
	 */
	@Override
	public double[] sequence()
	{
		return new double[] { sbr };
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#mutationStepRange()
	 */
	@Override
	public double[] mutationStepRange()
	{
		return new double[] { SNRFilter.DEFAULT_RANGE };
	}
}