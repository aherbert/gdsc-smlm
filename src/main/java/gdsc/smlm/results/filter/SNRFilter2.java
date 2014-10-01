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

import gdsc.smlm.fitting.Gaussian2DFitter;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

/**
 * Filter results using a signal-to-noise (SNR) threshold and width range
 */
public class SNRFilter2 extends Filter
{
	@XStreamAsAttribute
	float snr;
	@XStreamAsAttribute
	float minWidth;
	@XStreamAsAttribute
	float maxWidth;
	@XStreamOmitField
	float lowerSigmaThreshold;
	@XStreamOmitField
	float upperSigmaThreshold;

	public SNRFilter2(float snr, float minWidth, float maxWidth)
	{
		this.snr = snr;
		if (maxWidth < minWidth)
		{
			float f = maxWidth;
			maxWidth = minWidth;
			minWidth = f;
		}
		this.minWidth = minWidth;
		this.maxWidth = maxWidth;
	}

	@Override
	protected String generateName()
	{
		return "SNR " + snr + ", width " + minWidth + "-" + maxWidth;
	}

	@Override
	protected String generateType()
	{
		return "SNR2";
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		// Set the width limit
		lowerSigmaThreshold = 0;
		upperSigmaThreshold = Float.POSITIVE_INFINITY;
		Pattern pattern = Pattern.compile("initialPeakWidth0>([\\d\\.]+)");
		Matcher match = pattern.matcher(peakResults.getConfiguration());
		if (match.find())
		{
			float s = Float.parseFloat(match.group(1));
			lowerSigmaThreshold = Gaussian2DFitter.fwhm2sd(s * minWidth);
			upperSigmaThreshold = Gaussian2DFitter.fwhm2sd(s * maxWidth);
		}
	}

	@Override
	public boolean accept(PeakResult peak)
	{
		return getSNR(peak) >= this.snr && peak.getSD() >= lowerSigmaThreshold &&
				peak.getSD() <= upperSigmaThreshold;
	}

	static float getSNR(PeakResult peak)
	{
		return (peak.noise > 0) ? peak.getSignal() / peak.noise : Float.POSITIVE_INFINITY;
	}

	@Override
	public double getNumericalValue()
	{
		return snr;
	}

	@Override
	public String getNumericalValueName()
	{
		return "SNR";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getDescription()
	 */
	@Override
	public String getDescription()
	{
		return "Filter results using a lower SNR threshold and width range. (Width is relative to initial peak width.)";
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
				return snr;
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
				return "SNR";
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
				return new SNRFilter2(updateParameter(snr, delta), minWidth, maxWidth);
			case 1:
				return new SNRFilter2(snr, updateParameter(minWidth, delta), maxWidth);
			default:
				return new SNRFilter2(snr, minWidth, updateParameter(maxWidth, delta));
		}
	}
}