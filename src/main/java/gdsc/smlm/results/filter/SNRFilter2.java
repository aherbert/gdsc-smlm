package gdsc.smlm.results.filter;

import gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;

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
 * Filter results using a signal-to-noise ratio (SNR) threshold and width range
 */
public class SNRFilter2 extends DirectFilter implements IMultiFilter
{
	@XStreamAsAttribute
	final float snr;
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

	@XStreamOmitField
	private Gaussian2DPeakResultCalculator calculator;

	public SNRFilter2(float snr, double minWidth, double maxWidth)
	{
		this.snr = Math.max(0, snr);
		// Only swap if max width is enabled
		if (maxWidth != 0 && maxWidth < minWidth)
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
		return "SNR " + snr + ", width " + minWidth + "-" + maxWidth;
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		calculator = Gaussian2DPeakResultHelper.create(peakResults.getPSF(), peakResults.getCalibration(), 0);

		// Set the width limit
		lowerSigmaThreshold = 0;
		upperSigmaThreshold = Float.POSITIVE_INFINITY;
		Pattern pattern = Pattern.compile("initialSD0>([\\d\\.]+)");
		Matcher match = pattern.matcher(peakResults.getConfiguration());
		if (match.find())
		{
			double s = Double.parseDouble(match.group(1));
			lowerSigmaThreshold = (float) (s * minWidth);
			upperSigmaThreshold = Filter.getUpperLimit(s * maxWidth);
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
		float sd = calculator.getStandardDeviation(peak.getParameters());
		return SNRFilter.getSNR(peak) >= this.snr && sd <= upperSigmaThreshold && sd >= lowerSigmaThreshold;
	}

	@Override
	public int validate(final PreprocessedPeakResult peak)
	{
		if (peak.getSNR() < this.snr)
			return V_SNR;
		if (widthEnabled)
		{
			if (peak.getXSDFactor() > upperSigmaThreshold || peak.getXSDFactor() < lowerSigmaThreshold)
				return V_X_SD_FACTOR;
		}
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
	 * @see gdsc.smlm.results.filter.Filter#getParameterValueInternal(int)
	 */
	@Override
	protected double getParameterValueInternal(int index)
	{
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
	 * @see gdsc.smlm.results.filter.Filter#getParameterIncrement(int)
	 */
	@Override
	public double getParameterIncrement(int index)
	{
		checkIndex(index);
		switch (index)
		{
			case 0:
				return SNRFilter.DEFAULT_INCREMENT;
			case 1:
				return WidthFilter2.DEFAULT_MIN_INCREMENT;
			default:
				return WidthFilter.DEFAULT_INCREMENT;
		}
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
		switch (index)
		{
			case 0:
				return ParameterType.SNR;
			case 1:
				return ParameterType.MIN_WIDTH;
			default:
				return ParameterType.MAX_WIDTH;
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
				return new SNRFilter2(updateParameter(snr, delta, SNRFilter.DEFAULT_RANGE), minWidth, maxWidth);
			case 1:
				return new SNRFilter2(snr, updateParameter(minWidth, delta, WidthFilter2.DEFAULT_MIN_RANGE), maxWidth);
			default:
				return new SNRFilter2(snr, minWidth, updateParameter(maxWidth, delta, WidthFilter.DEFAULT_RANGE));
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
		return new SNRFilter2((float) parameters[0], parameters[1], parameters[2]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#weakestParameters(double[])
	 */
	@Override
	public void weakestParameters(double[] parameters)
	{
		setMin(parameters, 0, snr);
		setMin(parameters, 1, minWidth);
		setMax(parameters, 2, maxWidth);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.DirectFilter#lowerBoundOrientation(int)
	 */
	@Override
	public int lowerBoundOrientation(int index)
	{
		return (index == 2) ? 1 : -1;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#upperLimit()
	 */
	@Override
	public double[] upperLimit()
	{
		return new double[] { Double.POSITIVE_INFINITY, WidthFilter.UPPER_LIMIT, WidthFilter.UPPER_LIMIT };
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

	public double getSignal()
	{
		return 0;
	}

	public double getSNR()
	{
		return snr;
	}

	public double getMinWidth()
	{
		return minWidth;
	}

	public double getMaxWidth()
	{
		return maxWidth;
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
		return 0;
	}

	public boolean isPrecisionUsesLocalBackground()
	{
		return false;
	}
}