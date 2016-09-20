package gdsc.smlm.results.filter;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

import gdsc.smlm.results.PeakResult;

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

/**
 * Filter results using an amplitude-to-noise ratio (ANR) threshold and width range
 */
public class ANRFilter2 extends DirectFilter
{
	@XStreamAsAttribute
	final float anr;
	@XStreamAsAttribute
	final double minWidth;
	@XStreamAsAttribute
	final double maxWidth;
	@XStreamOmitField
	float lowerSigmaThreshold;
	@XStreamOmitField
	float upperSigmaThreshold;
	@XStreamOmitField
	boolean widthEnabled;

	public ANRFilter2(float anr, double minWidth, double maxWidth)
	{
		this.anr = Math.max(0, anr);
		if (maxWidth < minWidth)
		{
			double f = maxWidth;
			maxWidth = minWidth;
			minWidth = f;
		}
		this.minWidth = Math.max(0, minWidth);
		this.maxWidth = Math.max(0, maxWidth);
	}

	@Override
	protected String generateName()
	{
		return "ANR " + anr + ", width " + minWidth + "-" + maxWidth;
	}

	@Override
	protected String generateType()
	{
		return "ANR2";
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		// Set the width limit
		lowerSigmaThreshold = 0;
		upperSigmaThreshold = Float.POSITIVE_INFINITY;
		Pattern pattern = Pattern.compile("initialSD0>([\\d\\.]+)");
		Matcher match = pattern.matcher(peakResults.getConfiguration());
		if (match.find())
		{
			double s = Double.parseDouble(match.group(1));
			lowerSigmaThreshold = (float) (s * minWidth);
			upperSigmaThreshold = (float) (s * maxWidth);
		}
	}

	@Override
	public void setup()
	{
		setup(true);
	}

	@Override
	public void setup(int flags)
	{
		setup(!areSet(flags, DirectFilter.NO_WIDTH));
	}

	private void setup(final boolean widthEnabled)
	{
		this.widthEnabled = widthEnabled;
		if (widthEnabled)
		{
			lowerSigmaThreshold = (float) minWidth;
			upperSigmaThreshold = Filter.getUpperLimit(maxWidth);
		}
	}

	@Override
	public boolean accept(PeakResult peak)
	{
		return ANRFilter.getANR(peak) >= this.anr && peak.getSD() <= upperSigmaThreshold &&
				peak.getSD() >= lowerSigmaThreshold;
	}

	@Override
	public int validate(final PreprocessedPeakResult peak)
	{
		if (ANRFilter.getANR(peak) < this.anr)
			return V_AMPLITUDE | V_NOISE;
		if (widthEnabled)
		{
			if (peak.getXSDFactor() > upperSigmaThreshold && peak.getXSDFactor() < lowerSigmaThreshold)
				return V_X_SD_FACTOR;
		}
		return 0;
	}

	@Override
	public double getNumericalValue()
	{
		return anr;
	}

	@Override
	public String getNumericalValueName()
	{
		return "ANR";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getDescription()
	 */
	@Override
	public String getDescription()
	{
		return "Filter results using a lower ANR threshold and width range. (Width is relative to initial peak width.)";
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
				return anr;
			case 1:
				return minWidth;
			default:
				return maxWidth;
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
				return "ANR";
			case 1:
				return "Min width";
			default:
				return "Max width";
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
				return new ANRFilter2(updateParameter(anr, delta, SNRFilter.DEFAULT_RANGE), minWidth, maxWidth);
			case 1:
				return new ANRFilter2(anr, updateParameter(minWidth, delta, WidthFilter2.DEFAULT_MIN_RANGE), maxWidth);
			default:
				return new ANRFilter2(anr, minWidth, updateParameter(maxWidth, delta, WidthFilter.DEFAULT_RANGE));
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#create(double[])
	 */
	@Override
	public Filter create(double... parameters)
	{
		return new ANRFilter2((float) parameters[0], parameters[1], parameters[2]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#weakestParameters(double[])
	 */
	@Override
	public void weakestParameters(double[] parameters)
	{
		setMin(parameters, 0, anr);
		setMin(parameters, 1, minWidth);
		setMax(parameters, 2, maxWidth);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#length()
	 */
	public int length()
	{
		return 3;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#sequence()
	 */
	public double[] sequence()
	{
		return new double[] { anr, minWidth, maxWidth };
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#mutationStepRange()
	 */
	public double[] mutationStepRange()
	{
		return new double[] { SNRFilter.DEFAULT_RANGE, WidthFilter2.DEFAULT_MIN_RANGE, WidthFilter.DEFAULT_RANGE };
	}
}