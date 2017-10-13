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

import gdsc.smlm.ga.Chromosome;
import gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

/**
 * Filter results using multiple thresholds: Signal, SNR, width, coordinate shift and precision.
 * <p>
 * Any results with the strict limits are included. Any results outside the weak limits are excluded. Any results
 * between the strict and weak limits are included only if they can be traced through time, optionally via other
 * candidates, to a valid result.
 */
public class MultiHysteresisFilter extends HysteresisFilter
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
	private Gaussian2DPeakResultCalculator calculator;

	/**
	 * @param searchDistance
	 * @param searchDistanceMode
	 *            0 = relative to the precision of the candidates; 1 = Absolute (in nm)
	 * @param timeThreshold
	 * @param timeThresholdMode
	 *            0 = frames; 1 = seconds
	 * @param strictSignal
	 * @param rangeSignal
	 * @param strictSnr
	 * @param rangeSnr
	 * @param strictMinWidth
	 * @param rangeMinWidth
	 * @param strictMaxWidth
	 * @param rangeMaxWidth
	 * @param strictShift
	 * @param rangeShift
	 * @param strictPrecision
	 * @param rangePrecision
	 */
	public MultiHysteresisFilter(double searchDistance, int searchDistanceMode, double timeThreshold,
			int timeThresholdMode, double strictSignal, double rangeSignal, float strictSnr, float rangeSnr,
			double strictMinWidth, double rangeMinWidth, double strictMaxWidth, double rangeMaxWidth,
			double strictShift, double rangeShift, double strictPrecision, double rangePrecision)
	{
		super(searchDistance, searchDistanceMode, timeThreshold, timeThresholdMode);
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
		return String.format(
				"Multi Hysteresis: Signal=%.1f-%.1f, SNR=%.1f-%.1f, MinWidth=%.2f-%.2f, MaxWidth=%.2f+%.2f, Shift=%.2f+%.2f, Precision=%.1f+%.1f (%s)",
				strictSignal, rangeSignal, strictSnr, rangeSnr, strictMinWidth, rangeMinWidth, strictMaxWidth,
				rangeMaxWidth, strictShift, rangeShift, strictPrecision, rangePrecision, getTraceParameters());
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		calculator = Gaussian2DPeakResultHelper.create(peakResults.getPSF(), peakResults.getCalibration(),
				Gaussian2DPeakResultHelper.LSE_PRECISION);

		// Set the signal limit using the gain
		strictSignalThreshold = (float) (strictSignal);
		weakSignalThreshold = (float) (strictSignal - rangeSignal);

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
			strictMaxSigmaThreshold = Filter.getUpperLimit(s * strictMaxWidth);
			weakMinSigmaThreshold = (float) (s * (strictMinWidth - rangeMinWidth));
			weakMaxSigmaThreshold = Filter.getUpperLimit(s * (strictMaxWidth + rangeMaxWidth));

			strictOffset = Filter.getUpperLimit(s * strictShift);
			weakOffset = Filter.getUpperLimit(s * (strictShift + rangeShift));
		}

		// Configure the precision limit
		strictVariance = Filter.getDUpperSquaredLimit(strictPrecision);
		weakVariance = Filter.getDUpperSquaredLimit(strictPrecision + rangePrecision);

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
		final float sd = calculator.getStandardDeviation(result.getParameters());
		if (sd < weakMinSigmaThreshold || sd > weakMaxSigmaThreshold)
			return PeakStatus.REJECT;
		if (Math.abs(result.getXPosition()) > weakOffset || Math.abs(result.getYPosition()) > weakOffset)
			return PeakStatus.REJECT;
		final double variance = calculator.getLSEVariance(result.getParameters(), result.getNoise());
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
		return ParameterType.SNR.toString() + " +" + rangeSnr;
	}

	@Override
	public String getDescription()
	{
		return "Filter results using a multiple thresholds: Signal, SNR, width, shift, precision. Any results within the " +
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
		return 12 + super.getNumberOfParameters();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getParameterValueInternal(int)
	 */
	@Override
	protected double getParameterValueInternal(int index)
	{
		if (index < super.getNumberOfParameters())
		{
			return super.getParameterValueInternal(index);
		}
		index -= super.getNumberOfParameters();
		switch (index)
		{
			case 0:
				return strictSignal;
			case 1:
				return rangeSignal;
			case 2:
				return strictSnr;
			case 3:
				return rangeSnr;
			case 4:
				return strictMinWidth;
			case 5:
				return rangeMinWidth;
			case 6:
				return strictMaxWidth;
			case 7:
				return rangeMaxWidth;
			case 8:
				return strictShift;
			case 9:
				return rangeShift;
			case 10:
				return strictPrecision;
			default:
				return rangePrecision;
		}
	}

	@Override
	public double getParameterIncrement(int index)
	{
		if (index < super.getNumberOfParameters())
		{
			return super.getParameterValueInternal(index);
		}
		index -= super.getNumberOfParameters();
		switch (index)
		{
			case 0:
				return SignalFilter.DEFAULT_INCREMENT;
			case 1:
				return SignalFilter.DEFAULT_INCREMENT;
			case 2:
				return SNRFilter.DEFAULT_INCREMENT;
			case 3:
				return SNRFilter.DEFAULT_INCREMENT;
			case 4:
				return WidthFilter2.DEFAULT_MIN_INCREMENT;
			case 5:
				return WidthFilter2.DEFAULT_MIN_INCREMENT;
			case 6:
				return WidthFilter.DEFAULT_INCREMENT;
			case 7:
				return WidthFilter.DEFAULT_INCREMENT;
			case 8:
				return ShiftFilter.DEFAULT_INCREMENT;
			case 9:
				return ShiftFilter.DEFAULT_INCREMENT;
			case 10:
				return PrecisionFilter.DEFAULT_INCREMENT;
			default:
				return PrecisionFilter.DEFAULT_INCREMENT;
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
		if (index < super.getNumberOfParameters())
		{
			return super.getParameterType(index);
		}
		index -= super.getNumberOfParameters();
		switch (index)
		{
			case 0:
				return ParameterType.SIGNAL;
			case 1:
				return ParameterType.SIGNAL_RANGE;
			case 2:
				return ParameterType.SNR;
			case 3:
				return ParameterType.SNR_RANGE;
			case 4:
				return ParameterType.MIN_WIDTH;
			case 5:
				return ParameterType.MIN_WIDTH_RANGE;
			case 6:
				return ParameterType.MAX_WIDTH;
			case 7:
				return ParameterType.MAX_WIDTH_RANGE;
			case 8:
				return ParameterType.SHIFT;
			case 9:
				return ParameterType.SHIFT_RANGE;
			case 10:
				return ParameterType.PRECISION;
			default:
				return ParameterType.PRECISION_RANGE;
		}
	}

	static double[] defaultRange = new double[] { 0, 0, 0, 0, SignalFilter.DEFAULT_RANGE, SignalFilter.DEFAULT_RANGE,
			SNRFilter.DEFAULT_RANGE, SNRFilter.DEFAULT_RANGE, WidthFilter2.DEFAULT_MIN_RANGE,
			WidthFilter2.DEFAULT_MIN_RANGE, WidthFilter.DEFAULT_RANGE, WidthFilter.DEFAULT_RANGE,
			ShiftFilter.DEFAULT_RANGE, ShiftFilter.DEFAULT_RANGE, PrecisionFilter.DEFAULT_RANGE,
			PrecisionFilter.DEFAULT_RANGE };

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
				strictSignal, rangeSignal, strictSnr, rangeSnr, strictMinWidth, rangeMinWidth, strictMaxWidth,
				rangeMaxWidth, strictShift, rangeShift, strictPrecision, rangePrecision };
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
		return new MultiHysteresisFilter(parameters[0], (int) parameters[1], parameters[2], (int) parameters[3],
				parameters[4], parameters[5], (float) parameters[6], (float) parameters[7], parameters[8],
				parameters[9], parameters[10], parameters[11], parameters[12], parameters[13], parameters[14],
				parameters[15]);
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
		setMin(parameters, 4, strictSignal - rangeSignal);
		parameters[5] = 0;
		setMin(parameters, 6, strictSnr - rangeSnr);
		parameters[7] = 0;
		setMin(parameters, 8, strictMinWidth - rangeMinWidth);
		parameters[9] = 0;
		setMax(parameters, 10, strictMaxWidth + rangeMaxWidth);
		parameters[11] = 0;
		setMax(parameters, 12, strictShift + rangeShift);
		parameters[13] = 0;
		setMax(parameters, 14, strictPrecision + rangePrecision);
		parameters[15] = 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#newChromosome(double[])
	 */
	@Override
	public Chromosome<FilterScore> newChromosome(double[] sequence)
	{
		// Override the default Hysteresis filter implementation for speed since this is the filter we
		// will most likely optimise using the genetic algorithm
		return new MultiHysteresisFilter(sequence[0], searchDistanceMode, sequence[1], timeThresholdMode, sequence[2],
				sequence[3], (float) sequence[4], (float) sequence[5], sequence[6], sequence[7], sequence[8],
				sequence[9], sequence[10], sequence[11], sequence[12], sequence[13]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#upperLimit()
	 */
	@Override
	public double[] upperLimit()
	{
		return new double[] { Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY,
				Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, WidthFilter.UPPER_LIMIT,
				WidthFilter.UPPER_LIMIT, WidthFilter.UPPER_LIMIT, WidthFilter.UPPER_LIMIT, ShiftFilter.UPPER_LIMIT,
				ShiftFilter.UPPER_LIMIT, PrecisionFilter.UPPER_LIMIT, PrecisionFilter.UPPER_LIMIT };
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ga.Chromosome#mutationStepRange()
	 */
	public double[] mutationStepRange()
	{
		return new double[] { getDefaultSearchRange(), getDefaultTimeRange(), SignalFilter.DEFAULT_RANGE,
				SignalFilter.DEFAULT_RANGE, SNRFilter.DEFAULT_RANGE, SNRFilter.DEFAULT_RANGE,
				WidthFilter2.DEFAULT_MIN_RANGE, WidthFilter2.DEFAULT_MIN_RANGE, WidthFilter.DEFAULT_RANGE,
				WidthFilter.DEFAULT_RANGE, ShiftFilter.DEFAULT_RANGE, ShiftFilter.DEFAULT_RANGE,
				PrecisionFilter.DEFAULT_RANGE, PrecisionFilter.DEFAULT_RANGE };
	}
}