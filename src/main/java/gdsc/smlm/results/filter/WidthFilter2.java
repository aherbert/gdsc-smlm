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

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

/**
 * Filter results using a width range
 */
public class WidthFilter2 extends Filter
{
	@XStreamAsAttribute
	final double minWidth;
	@XStreamAsAttribute
	final double maxWidth;
	@XStreamOmitField
	float lowerSigmaThreshold;
	@XStreamOmitField
	float upperSigmaThreshold;

	public WidthFilter2(double minWidth, double maxWidth)
	{
		if (maxWidth < minWidth)
		{
			double f = maxWidth;
			maxWidth = minWidth;
			minWidth = f;
		}
		this.minWidth = minWidth;
		this.maxWidth = maxWidth;
	}

	@Override
	protected String generateName()
	{
		return "Width " + minWidth + "-" + maxWidth;
	}

	@Override
	protected String generateType()
	{
		return "Width2";
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
	public boolean accept(PeakResult peak)
	{
		return peak.getSD() >= lowerSigmaThreshold && peak.getSD() <= upperSigmaThreshold;
	}

	@Override
	public double getNumericalValue()
	{
		return minWidth;
	}

	@Override
	public String getNumericalValueName()
	{
		return "Width";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getDescription()
	 */
	@Override
	public String getDescription()
	{
		return "Filter results using a width range. (Width is relative to initial peak width.)";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getNumberOfParameters()
	 */
	@Override
	public int getNumberOfParameters()
	{
		return 2;
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
				return new WidthFilter2(updateParameter(minWidth, delta), maxWidth);
			default:
				return new WidthFilter2(minWidth, updateParameter(maxWidth, delta));
		}
	}
}