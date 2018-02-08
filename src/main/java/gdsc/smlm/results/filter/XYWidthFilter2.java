package gdsc.smlm.results.filter;

import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
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

/**
 * Filter results using a width range. Assumes width is different on the X and Y axis and they are combined
 * using s = sqrt(s0*s1).
 */
public class XYWidthFilter2 extends WidthFilter2 implements IMultiFilter
{
	public XYWidthFilter2(double minWidth, double maxWidth)
	{
		super(minWidth, maxWidth);
	}

	@Override
	protected String generateName()
	{
		return "XYWidth " + minWidth + "-" + maxWidth;
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		calculator = Gaussian2DPeakResultHelper.create(peakResults.getPSF(), peakResults.getCalibration(), 0);

		// Set the width limit
		double[] s = PSFHelper.getGaussian2DWxWy(peakResults.getPSF());
		double s2 = s[0] * s[1];
		lowerSigmaThreshold = (float) (s2 * minWidth * minWidth);
		upperSigmaThreshold = Filter.getUpperLimit(s2 * maxWidth * maxWidth);
	}

	@Override
	void setup(final boolean widthEnabled)
	{
		this.widthEnabled = widthEnabled;
		if (widthEnabled)
		{
			lowerSigmaThreshold = (float) (minWidth * minWidth);
			upperSigmaThreshold = Filter.getUpperLimit(maxWidth * maxWidth);
		}
	}

	@Override
	public boolean accept(PeakResult peak)
	{
		final float sd2 = calculator.getStandardDeviation2(peak.getParameters());
		return sd2 <= upperSigmaThreshold && sd2 >= lowerSigmaThreshold;
	}

	public int getValidationFlags()
	{
		return V_X_SD_FACTOR | V_Y_SD_FACTOR;
	}

	@Override
	public int validate(final PreprocessedPeakResult peak)
	{
		if (widthEnabled)
		{
			final float s2 = peak.getXSDFactor() * peak.getYSDFactor();
			if (s2 > upperSigmaThreshold || s2 < lowerSigmaThreshold)
				return V_X_SD_FACTOR | V_Y_SD_FACTOR;
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
		return "Filter results using an XY width range. (Width is relative to initial peak width.)";
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
				return new XYWidthFilter2(updateParameter(minWidth, delta, DEFAULT_MIN_RANGE), maxWidth);
			default:
				return new XYWidthFilter2(minWidth, updateParameter(maxWidth, delta, WidthFilter.DEFAULT_RANGE));
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
		return new XYWidthFilter2(parameters[0], parameters[1]);
	}
}