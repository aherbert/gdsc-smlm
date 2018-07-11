/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm.results.filter;

import com.thoughtworks.xstream.annotations.XStreamOmitField;

import gdsc.smlm.data.config.ConfigurationException;
import gdsc.smlm.ga.Chromosome;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.results.PeakResult;

/**
 * Filter results using multiple thresholds: Signal, SNR, width, coordinate shift and precision. Calculates the
 * precision using the true fitted background if a bias is provided.
 * <p>
 * Any results with the strict limits are included. Any results outside the weak limits are excluded. Any results
 * between the strict and weak limits are included only if they can be traced through time, optionally via other
 * candidates, to a valid result.
 */
public class MultiHysteresisFilter2 extends MultiHysteresisFilter
{
	@XStreamOmitField
	boolean useBackground = false;

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
	public MultiHysteresisFilter2(double searchDistance, int searchDistanceMode, double timeThreshold,
			int timeThresholdMode, double strictSignal, double rangeSignal, float strictSnr, float rangeSnr,
			double strictMinWidth, double rangeMinWidth, double strictMaxWidth, double rangeMaxWidth,
			double strictShift, double rangeShift, double strictPrecision, double rangePrecision)
	{
		super(searchDistance, searchDistanceMode, timeThreshold, timeThresholdMode, strictSignal, rangeSignal,
				strictSnr, rangeSnr, strictMinWidth, rangeMinWidth, strictMaxWidth, rangeMaxWidth, strictShift,
				rangeShift, strictPrecision, rangePrecision);
	}

	@Override
	protected String generateName()
	{
		return String.format(
				"Multi Hysteresis2: Signal=%.1f-%.1f, SNR=%.1f-%.1f, MinWidth=%.2f-%.2f, MaxWidth=%.2f+%.2f, Shift=%.2f+%.2f, Precision2=%.1f+%.1f (%s)",
				strictSignal, rangeSignal, strictSnr, rangeSnr, strictMinWidth, rangeMinWidth, strictMaxWidth,
				rangeMaxWidth, strictShift, rangeShift, strictPrecision, rangePrecision, getTraceParameters());
	}

	@Override
	protected void setupCalculator(MemoryPeakResults peakResults)
	{
		try
		{
			calculator = Gaussian2DPeakResultHelper.create(peakResults.getPSF(), peakResults.getCalibration(),
					Gaussian2DPeakResultHelper.LSE_PRECISION_X);
			useBackground = true;
		}
		catch (ConfigurationException e)
		{
			calculator = Gaussian2DPeakResultHelper.create(peakResults.getPSF(), peakResults.getCalibration(),
					Gaussian2DPeakResultHelper.LSE_PRECISION);
			useBackground = false;
		}
	}

	@Override
	protected double getVariance(PeakResult result)
	{
		if (useBackground)
		{
			return calculator.getLSEVariance(result.getParameters());
		}
		else
		{
			return calculator.getLSEVariance(result.getParameters(), result.getNoise());
		}
	}

	@Override
	public String getDescription()
	{
		return "Filter results using a multiple thresholds: Signal, SNR, width, shift, precision (uses fitted background to set noise). Any results within the " +
				"strict limits are included. Any results outside the weak limits are excluded. " +
				super.getDescription();
	}

	@Override
	protected ParameterType getPrecisionParamaterType()
	{
		return ParameterType.PRECISION2;
	}

	@Override
	protected ParameterType getPrecisionRangeParamaterType()
	{
		return ParameterType.PRECISION2_RANGE;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.filter.Filter#create(double[])
	 */
	@Override
	public Filter create(double... parameters)
	{
		return new MultiHysteresisFilter2(parameters[0], (int) parameters[1], parameters[2], (int) parameters[3],
				parameters[4], parameters[5], (float) parameters[6], (float) parameters[7], parameters[8],
				parameters[9], parameters[10], parameters[11], parameters[12], parameters[13], parameters[14],
				parameters[15]);
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
		return new MultiHysteresisFilter2(sequence[0], searchDistanceMode, sequence[1], timeThresholdMode, sequence[2],
				sequence[3], (float) sequence[4], (float) sequence[5], sequence[6], sequence[7], sequence[8],
				sequence[9], sequence[10], sequence[11], sequence[12], sequence[13]);
	}
}
