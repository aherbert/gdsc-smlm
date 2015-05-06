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

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

/**
 * Filter results using multiple thresholds: Signal, SNR, width, coordinate shift and precision. Calculates the
 * precision using the true fitted background if a bias is provided.
 * <p>
 * Any results with the strict limits are included. Any results outside the weak limits are excluded. Any results
 * between the strict and weak limits are included only if they can be traced through time, optionally via other
 * candidates, to a valid result.
 */
public class MultiHysteresisFilter2 extends HysteresisFilter
{
	@XStreamAsAttribute
	final double strictSignal;
	@XStreamAsAttribute
	final float strictSnr;
	@XStreamAsAttribute
	final double strictMinWidth;
	@XStreamAsAttribute
	final double strictMaxWidth;
	@XStreamAsAttribute
	final double strictShift;
	@XStreamAsAttribute
	final double strictPrecision;
	@XStreamAsAttribute
	final double rangeSignal;
	@XStreamAsAttribute
	final float rangeSnr;
	@XStreamAsAttribute
	final double rangeMinWidth;
	@XStreamAsAttribute
	final double rangeMaxWidth;
	@XStreamAsAttribute
	final double rangeShift;
	@XStreamAsAttribute
	final double rangePrecision;

	@XStreamOmitField
	float strictSignalThreshold;
	@XStreamOmitField
	float weakSignalThreshold;
	@XStreamOmitField
	float weakSnr;
	@XStreamOmitField
	float strictMinSigmaThreshold;
	@XStreamOmitField
	float weakMinSigmaThreshold;
	@XStreamOmitField
	float strictMaxSigmaThreshold;
	@XStreamOmitField
	float weakMaxSigmaThreshold;
	@XStreamOmitField
	float strictOffset;
	@XStreamOmitField
	float weakOffset;
	@XStreamOmitField
	double strictVariance;
	@XStreamOmitField
	double weakVariance;

	@XStreamOmitField
	double nmPerPixel = 100;
	@XStreamOmitField
	boolean emCCD = true;
	@XStreamOmitField
	double gain = 1;
	@XStreamOmitField
	double bias = -1;

	public MultiHysteresisFilter2(double searchDistance, int searchDistanceMode, double strictSignal,
			double rangeSignal, float strictSnr, float rangeSnr, double strictMinWidth, double rangeMinWidth,
			double strictMaxWidth, double rangeMaxWidth, double strictShift, double rangeShift, double strictPrecision,
			double rangePrecision)
	{
		super(searchDistance, searchDistanceMode);
		this.strictSignal = Math.max(0, strictSignal);
		this.rangeSignal = Math.max(0, rangeSignal);
		this.strictSnr = Math.max(0, strictSnr);
		this.rangeSnr = Math.max(0, rangeSnr);
		this.strictMinWidth = Math.max(0, strictMinWidth);
		this.rangeMinWidth = Math.max(0, rangeMinWidth);
		this.strictMaxWidth = Math.max(0, strictMaxWidth);
		this.rangeMaxWidth = Math.max(0, rangeMaxWidth);
		this.strictShift = Math.max(0, strictShift);
		this.rangeShift = Math.max(0, rangeShift);
		this.strictPrecision = Math.max(0, strictPrecision);
		this.rangePrecision = Math.max(0, rangePrecision);
	}

	@Override
	protected String generateName()
	{
		return String
				.format("Multi2 Hysteresis: Signal=%.1f-%.1f, SNR=%.1f-%.1f, MinWidth=%.2f-%.2f, MaxWidth=%.2f+%.2f, Shift=%.2f+%.2f, Precision=%.1f+%.1f (@%.2f %s)",
						strictSignal, rangeSignal, strictSnr, rangeSnr, strictMinWidth, rangeMinWidth, strictMaxWidth,
						rangeMaxWidth, strictShift, rangeShift, strictPrecision, rangePrecision, searchDistance,
						getSearchName());
	}

	@Override
	protected String generateType()
	{
		return "Multi2 Hysteresis";
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		// Calibration
		nmPerPixel = peakResults.getNmPerPixel();
		gain = peakResults.getGain();
		emCCD = peakResults.isEMCCD();
		bias = -1;
		if (peakResults.getCalibration() != null)
		{
			bias = peakResults.getCalibration().bias;
		}

		// Set the signal limit using the gain
		strictSignalThreshold = (float) (strictSignal * gain);
		weakSignalThreshold = (float) ((strictSignal - rangeSignal) * gain);

		weakSnr = strictSnr - rangeSnr;

		// Set the width limit
		strictMinSigmaThreshold = weakMinSigmaThreshold = 0;
		strictMaxSigmaThreshold = weakMaxSigmaThreshold = Float.POSITIVE_INFINITY;
		// Set the shift limit
		strictOffset = weakOffset = Float.POSITIVE_INFINITY;

		Pattern pattern = Pattern.compile("initialSD0>([\\d\\.]+)");
		Matcher match = pattern.matcher(peakResults.getConfiguration());
		if (match.find())
		{
			double s = Double.parseDouble(match.group(1));
			strictMinSigmaThreshold = (float) (s * strictMinWidth);
			strictMaxSigmaThreshold = (float) (s * strictMaxWidth);
			weakMinSigmaThreshold = (float) (s * (strictMinWidth - rangeMinWidth));
			weakMaxSigmaThreshold = (float) (s * (strictMaxWidth + rangeMaxWidth));

			strictOffset = (float) (s * strictShift);
			weakOffset = (float) (s * (strictShift + rangeShift));
		}

		// Configure the precision limit
		strictVariance = PrecisionFilter.getVarianceLimit(strictPrecision);
		weakVariance = PrecisionFilter.getVarianceLimit(strictPrecision + rangePrecision);

		super.setup(peakResults);
	}

	@Override
	protected PeakStatus getStatus(PeakResult result)
	{
		// Check weak thresholds
		if (result.getSignal() < weakSignalThreshold)
			return PeakStatus.REJECT;
		final float snr = SNRFilter.getSNR(result);
		if (snr < weakSnr)
			return PeakStatus.REJECT;
		final float sd = result.getSD();
		if (sd < weakMinSigmaThreshold || sd > weakMaxSigmaThreshold)
			return PeakStatus.REJECT;
		if (Math.abs(result.getXPosition()) > weakOffset || Math.abs(result.getYPosition()) > weakOffset)
			return PeakStatus.REJECT;
		// Use the background directly
		final double variance;
		if (bias != -1)
		{
			// Use the estimated background for the peak
			final double s = nmPerPixel * result.getSD();
			final double N = result.getSignal();
			variance = PeakResult.getVarianceX(nmPerPixel, s, N / gain,
					Math.max(0, result.params[Gaussian2DFunction.BACKGROUND] - bias) / gain, emCCD);
		}
		else
		{
			// Use the background noise to estimate precision
			variance = result.getVariance(nmPerPixel, gain, emCCD);
		}
		if (variance > weakVariance)
			return PeakStatus.REJECT;

		// Check the strict thresholds
		if (result.getSignal() < strictSignalThreshold)
			return PeakStatus.CANDIDATE;
		if (snr < strictSnr)
			return PeakStatus.CANDIDATE;
		if (sd < strictMinSigmaThreshold || sd > strictMaxSigmaThreshold)
			return PeakStatus.CANDIDATE;
		if (Math.abs(result.getXPosition()) > strictOffset || Math.abs(result.getYPosition()) > strictOffset)
			return PeakStatus.CANDIDATE;
		if (variance > strictVariance)
			return PeakStatus.CANDIDATE;

		return PeakStatus.OK;
	}

	@Override
	public double getNumericalValue()
	{
		return strictSnr;
	}

	@Override
	public String getNumericalValueName()
	{
		return "SNR +" + rangeSnr;
	}

	@Override
	public String getDescription()
	{
		return "Filter results using a multiple thresholds: Signal, SNR, width, shift, precision (uses fitted background to set noise). Any results within the " +
				"strict limits are included. Any results outside the weak limits are excluded. " +
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
		return 14;
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
				return searchDistanceMode;
			case 2:
				return strictSignal;
			case 3:
				return rangeSignal;
			case 4:
				return strictSnr;
			case 5:
				return rangeSnr;
			case 6:
				return strictMinWidth;
			case 7:
				return rangeMinWidth;
			case 8:
				return strictMaxWidth;
			case 9:
				return rangeMaxWidth;
			case 10:
				return strictShift;
			case 11:
				return rangeShift;
			case 12:
				return strictPrecision;
			default:
				return rangePrecision;
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
				return "Search mode";
			case 2:
				return "Strict Signal";
			case 3:
				return "Range Signal";
			case 4:
				return "Strict SNR";
			case 5:
				return "Range SNR";
			case 6:
				return "Strict Min Width";
			case 7:
				return "Range Min Width";
			case 8:
				return "Strict Max Width";
			case 9:
				return "Range Max Width";
			case 10:
				return "Strict Shift";
			case 11:
				return "Range Shift";
			case 12:
				return "Strict Precision";
			default:
				return "Range Precision";
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
		if (index == 1)
			return this;
		double[] parameters = new double[] { searchDistance, searchDistanceMode, strictSignal, rangeSignal, strictSnr,
				rangeSnr, strictMinWidth, rangeMinWidth, strictMaxWidth, rangeMaxWidth, strictShift, rangeShift,
				strictPrecision, rangePrecision };
		if (index == 0)
			parameters[0] = updateParameter(parameters[0], delta, getDefaultSearchRange());
		else
			parameters[index] = updateParameter(parameters[index], delta, MultiHysteresisFilter.defaultRange[index]);
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
		return new MultiHysteresisFilter2(parameters[0], (int) parameters[1], parameters[2], parameters[3],
				(float) parameters[4], (float) parameters[5], parameters[6], parameters[7], parameters[8],
				parameters[9], parameters[10], parameters[11], parameters[12], parameters[13]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#weakestParameters(double[])
	 */
	@Override
	public void weakestParameters(double[] parameters)
	{
		setMax(parameters, 0, searchDistance);
		// Hysteresis filters require all the potential candidates, so disable hysteresis above the candidate threshold  
		setMin(parameters, 2, strictSignal - rangeSignal);
		parameters[3] = 0;
		setMin(parameters, 4, strictSnr - rangeSnr);
		parameters[5] = 0;
		setMin(parameters, 6, strictMinWidth - rangeMinWidth);
		parameters[7] = 0;
		setMax(parameters, 8, strictMaxWidth + rangeMaxWidth);
		parameters[9] = 0;
		setMax(parameters, 10, strictShift + rangeShift);
		parameters[11] = 0;
		setMax(parameters, 12, strictPrecision + rangePrecision);
		parameters[13] = 0;
	}
}