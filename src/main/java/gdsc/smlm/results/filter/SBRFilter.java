package gdsc.smlm.results.filter;

import com.thoughtworks.xstream.annotations.XStreamAsAttribute;
import com.thoughtworks.xstream.annotations.XStreamOmitField;

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
import gdsc.smlm.results.PeakResultData;
import gdsc.smlm.results.data.Gaussian2DPeakResultDataMeanSignal;

/**
 * Filter results using a signal-to-background ratio (SBR) threshold.
 * <p>
 * Requires the bias to be configured at or above zero. If the background is below the configured bias or there is no
 * bias then the filter resorts to a signal-to-noise filter. If there is a background level above the bias then this is
 * assumed to be the variance of the photon shot noise and the noise is taken at the square root of the background
 * level.
 * <p>
 * This filter assumes the input results are Gaussian2D peak results.
 */
public class SBRFilter extends DirectFilter
{
	@XStreamAsAttribute
	final float sbr;

	@XStreamOmitField
	final PeakResultData<Float> converter = new Gaussian2DPeakResultDataMeanSignal();

	public SBRFilter(float sbr)
	{
		this.sbr = Math.max(0, sbr);
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
	}

	@Override
	public boolean accept(PeakResult peak)
	{
		final double background = peak.getBackground();
		if (background > 0)
			return converter.getValue(peak) / Math.sqrt(background) >= this.sbr;
		return SNRFilter.getSNR(peak) >= this.sbr;
	}

	public int getValidationFlags()
	{
		return V_PHOTONS | V_BACKGROUND | V_SNR;
	}

	@Override
	public int validate(final PreprocessedPeakResult peak)
	{
		final double background = peak.getBackground();
		if (background > 0)
		{
			// Get the mean signal assuming the integral / area of 1 SD of the Gaussian
			if (Gaussian2DPeakResultHelper.getMeanSignalUsingR1(peak.getSignal(), peak.getXSD(), peak.getYSD()) /
					Math.sqrt(background) < this.sbr)
				return V_PHOTONS | V_BACKGROUND;
			return 0;
		}
		if (peak.getSNR() < this.sbr)
			return V_SNR;
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
	 * @see gdsc.smlm.results.filter.Filter#getParameterValueInternal(int)
	 */
	@Override
	protected double getParameterValueInternal(int index)
	{
		return sbr;
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
		return SNRFilter.DEFAULT_INCREMENT;
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
		return ParameterType.SBR;
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
	 * @see gdsc.smlm.ga.Chromosome#mutationStepRange()
	 */
	public double[] mutationStepRange()
	{
		return new double[] { SNRFilter.DEFAULT_RANGE };
	}
}