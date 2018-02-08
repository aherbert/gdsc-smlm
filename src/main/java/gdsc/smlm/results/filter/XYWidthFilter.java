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
 * Filter results using an upper width factor. Assumes width is different on the X and Y axis and they are combined
 * using s = sqrt(s0*s1).
 */
public class XYWidthFilter extends WidthFilter implements IMultiFilter
{
	public XYWidthFilter(double width)
	{
		super(width);
	}

	@Override
	public void setup(MemoryPeakResults peakResults)
	{
		calculator = Gaussian2DPeakResultHelper.create(peakResults.getPSF(), peakResults.getCalibration(), 0);

		// Set the width limit
		double[] s = PSFHelper.getGaussian2DWxWy(peakResults.getPSF());
		upperSigmaThreshold = Filter.getUpperLimit(s[0] * s[1] * width * width);
	}

	@Override
	void setup(final boolean widthEnabled)
	{
		this.widthEnabled = widthEnabled;
		if (widthEnabled)
		{
			upperSigmaThreshold = Filter.getUpperLimit(width * width);
		}
	}

	@Override
	public boolean accept(PeakResult peak)
	{
		return calculator.getStandardDeviation2(peak.getParameters()) <= upperSigmaThreshold;
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
			if (peak.getXSDFactor() * peak.getYSDFactor() > upperSigmaThreshold)
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
		return "Filter results using an upper XY width factor. (Width is relative to initial peak width.)";
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
		return new XYWidthFilter(updateParameter(width, delta, DEFAULT_RANGE));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#create(double[])
	 */
	@Override
	public Filter create(double... parameters)
	{
		return new XYWidthFilter(parameters[0]);
	}
}