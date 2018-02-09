package gdsc.smlm.results.filter;

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

/**
 * Filter results using Width. Assume XY width are different.
 */
public class MultiFilterXYWidthComponent extends MultiFilterComponent
{
	final float lowerSigmaThreshold, upperSigmaThreshold;

	public MultiFilterXYWidthComponent(double minWidth, double maxWidth)
	{
		if (minWidth > 0 && minWidth < 1)
			lowerSigmaThreshold = (float) (minWidth * minWidth);
		else
			lowerSigmaThreshold = 0;
		upperSigmaThreshold = Filter.getUpperLimit(maxWidth * maxWidth);
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.MultiFilterComponent#fail(gdsc.smlm.results.filter.PreprocessedPeakResult)
	 */
	public boolean fail(final PreprocessedPeakResult peak)
	{
		final float s2 = peak.getXSDFactor() * peak.getYSDFactor();
		return (s2 > upperSigmaThreshold || s2 < lowerSigmaThreshold);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.MultiFilterComponent#getType()
	 */
	public int getType()
	{
		return IDirectFilter.V_X_SD_FACTOR | IDirectFilter.V_Y_SD_FACTOR;
	}
}