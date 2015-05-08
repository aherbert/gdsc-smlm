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
 * Filter results using multiple thresholds: Signal, SNR, width, coordinate shift and precision
 */
public class MultiFilter extends Filter
{
	@XStreamAsAttribute
	final double signal;
	@XStreamAsAttribute
	final float snr;
	@XStreamAsAttribute
	final double minWidth;
	@XStreamAsAttribute
	final double maxWidth;
	@XStreamAsAttribute
	final double shift;
	@XStreamAsAttribute
	final double precision;

	@XStreamOmitField
	float signalThreshold;
	@XStreamOmitField
	float lowerSigmaThreshold;
	@XStreamOmitField
	float upperSigmaThreshold;
	@XStreamOmitField
	float offset;
	@XStreamOmitField
	double variance;
	@XStreamOmitField
	double nmPerPixel = 100;
	@XStreamOmitField
	boolean emCCD = true;
	@XStreamOmitField
	double gain = 1;

	public MultiFilter(double signal, float snr, double minWidth, double maxWidth, double shift, double precision)
	{
		this.signal = Math.max(0, signal);
		this.snr = Math.max(0, snr);
		if (maxWidth < minWidth)
		{
			double f = maxWidth;
			maxWidth = minWidth;
			minWidth = f;
		}
		this.minWidth = Math.max(0, minWidth);
		this.maxWidth = Math.max(0, maxWidth);
		this.shift = Math.max(0, shift);
		this.precision = Math.max(0, precision);
	}

	@Override
	protected String generateName()
	{
		return String.format("Multi: Signal=%.1f, SNR=%.1f, Width=%.2f-%.2f, Shift=%.2f, Precision=%.1f", signal, snr,
				minWidth, maxWidth, shift, precision);
	}

	@Override
	protected String generateType()
	{
		return "Multi";
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		// Set the signal limit using the gain
		signalThreshold = (float) (signal * peakResults.getCalibration().gain);

		// Set the width limit
		lowerSigmaThreshold = 0;
		upperSigmaThreshold = Float.POSITIVE_INFINITY;
		// Set the shift limit
		offset = Float.POSITIVE_INFINITY;
		Pattern pattern = Pattern.compile("initialSD0>([\\d\\.]+)");
		Matcher match = pattern.matcher(peakResults.getConfiguration());
		if (match.find())
		{
			double s = Double.parseDouble(match.group(1));
			lowerSigmaThreshold = (float) (s * minWidth);
			upperSigmaThreshold = (float) (s * maxWidth);
			offset = (float) (s * shift);
		}

		// Configure the precision limit
		variance = PrecisionFilter.getVarianceLimit(precision);
		nmPerPixel = peakResults.getNmPerPixel();
		gain = peakResults.getGain();
		emCCD = peakResults.isEMCCD();
	}

	@Override
	public boolean accept(PeakResult peak)
	{
		if (peak.getSignal() < signalThreshold)
			return false;
		if (SNRFilter.getSNR(peak) < this.snr)
			return false;
		final float sd = peak.getSD();
		if (sd < lowerSigmaThreshold || sd > upperSigmaThreshold)
			return false;
		if (Math.abs(peak.getXPosition()) > offset || Math.abs(peak.getYPosition()) > offset)
			return false;
		// Use the background noise to estimate precision
		return peak.getVariance(nmPerPixel, gain, emCCD) <= variance;
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
		return "Filter results using an multiple thresholds: Signal, SNR, width, shift and precision";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getNumberOfParameters()
	 */
	@Override
	public int getNumberOfParameters()
	{
		return 6;
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
				return signal;
			case 1:
				return snr;
			case 2:
				return minWidth;
			case 3:
				return maxWidth;
			case 4:
				return shift;
			default:
				return precision;
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
				return "Signal";
			case 1:
				return "SNR";
			case 2:
				return "Min width";
			case 3:
				return "Max width";
			case 4:
				return "Shift";
			default:
				return "Precision";
		}
	}

	static double[] defaultRange = new double[]{
		SignalFilter.DEFAULT_RANGE,
		SNRFilter.DEFAULT_RANGE,
		WidthFilter2.DEFAULT_MIN_RANGE,
		WidthFilter.DEFAULT_RANGE,
		ShiftFilter.DEFAULT_RANGE,
		PrecisionFilter.DEFAULT_RANGE		
	};
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#adjustParameter(int, double)
	 */
	@Override
	public Filter adjustParameter(int index, double delta)
	{
		checkIndex(index);
		double[] params = new double[] { signal, snr, minWidth, maxWidth, shift, precision };
		params[index] = updateParameter(params[index], delta, defaultRange[index]);
		return new MultiFilter(params[0], (float) params[1], params[2], params[3], params[4], params[5]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#create(double[])
	 */
	@Override
	public Filter create(double... parameters)
	{
		return new MultiFilter(parameters[0], (float) parameters[1], parameters[2], parameters[3], parameters[4],
				parameters[5]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#weakestParameters(double[])
	 */
	@Override
	public void weakestParameters(double[] parameters)
	{
		setMin(parameters, 0, signal);
		setMin(parameters, 1, snr);
		setMin(parameters, 2, minWidth);
		setMax(parameters, 3, maxWidth);
		setMax(parameters, 4, shift);
		setMax(parameters, 5, precision);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#length()
	 */
	@Override
	public int length()
	{
		return 6;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#upperLimit()
	 */
	@Override
	public double[] upperLimit()
	{
		return new double[] { Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, WidthFilter.UPPER_LIMIT,
				WidthFilter.UPPER_LIMIT, ShiftFilter.UPPER_LIMIT, PrecisionFilter.UPPER_LIMIT };
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#sequence()
	 */
	@Override
	public double[] sequence()
	{
		return new double[] { signal, snr, minWidth, maxWidth, shift, precision };
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#mutationStepRange()
	 */
	@Override
	public double[] mutationStepRange()
	{
		return defaultRange;
	}
}