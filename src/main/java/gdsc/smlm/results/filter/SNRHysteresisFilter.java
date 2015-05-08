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
 * Filter results using a signal-to-noise (SNR) threshold. Any results above the upper SNR limit are
 * included. Any results below the lower SNR limit are excluded. Any results between the limits are included only if
 * they can be traced through time, optionally via other candidates, to a valid result.
 */
public class SNRHysteresisFilter extends HysteresisFilter
{
	@XStreamAsAttribute
	final float strictSnr;
	@XStreamAsAttribute
	final float range;
	@XStreamOmitField
	float weakSnr;

	/**
	 * @param searchDistance
	 * @param searchDistanceMode
	 *            0 = relative to the precision of the candidates; 1 = Absolute (in nm)
	 * @param timeThreshold
	 * @param timeThresholdMode
	 *            0 = frames; 1 = seconds
	 * @param strictSnr
	 * @param range
	 */
	public SNRHysteresisFilter(double searchDistance, int searchDistanceMode, double timeThreshold,
			int timeThresholdMode, float strictSnr, float range)
	{
		super(searchDistance, searchDistanceMode, timeThreshold, timeThresholdMode);
		this.strictSnr = Math.max(0, strictSnr);
		this.range = Math.max(0, range);
	}

	@Override
	protected String generateName()
	{
		return String.format("SNR Hysteresis %.2f -%.2f (%s)", strictSnr, range, getTraceParameters());
	}

	@Override
	protected String generateType()
	{
		return "SNR Hysteresis";
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		weakSnr = strictSnr - range;
		super.setup(peakResults);
	}

	@Override
	protected PeakStatus getStatus(PeakResult result)
	{
		final float snr = SNRFilter.getSNR(result);
		if (snr >= strictSnr)
			return PeakStatus.OK;
		else if (snr >= weakSnr)
			return PeakStatus.CANDIDATE;
		return PeakStatus.REJECT;
	}

	@Override
	public double getNumericalValue()
	{
		return strictSnr;
	}

	@Override
	public String getNumericalValueName()
	{
		return "SNR -" + range;
	}

	@Override
	public String getDescription()
	{
		return "Filter results using a signal-to-noise (SNR) threshold. Any results above the upper SNR " +
				"limit are included. Any results below the lower SNR limit are excluded. " + super.getDescription();
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
				return strictSnr;
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
				return "Strict SNR";
			default:
				return "Range";
		}
	}

	static double[] defaultRange = new double[] { 0, 0, 0, 0, SNRFilter.DEFAULT_RANGE, SNRFilter.DEFAULT_RANGE };

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
				strictSnr, range };
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
		return new SNRHysteresisFilter(parameters[0], (int) parameters[1], parameters[2], (int) parameters[3],
				(float) parameters[4], (float) parameters[5]);
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
		setMin(parameters, 4, strictSnr);
		parameters[5] = 0;
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#length()
	 */
	@Override
	public int length()
	{
		return 4;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#sequence()
	 */
	@Override
	public double[] sequence()
	{
		// Ignore the mode parameters
		return new double[] { searchDistance, timeThreshold, strictSnr, range };
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#mutationStepRange()
	 */
	@Override
	public double[] mutationStepRange()
	{
		return new double[] { getDefaultSearchRange(), getDefaultTimeRange(), SNRFilter.DEFAULT_RANGE,
				SNRFilter.DEFAULT_RANGE };
	}
}