package gdsc.smlm.results.filter;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Filter results using Shift
 */
public class MultiFilterWidthComponent extends MultiFilterComponent
{
	final float lowerSigmaThreshold, upperSigmaThreshold;

	public MultiFilterWidthComponent(double minWidth, double maxWidth)
	{
		this.lowerSigmaThreshold = (float) minWidth;
		this.upperSigmaThreshold = Filter.getUpperLimit(maxWidth);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.MultiFilterComponent#getType()
	 */
	public boolean fail(final PreprocessedPeakResult peak)
	{
		final float xsdf = peak.getXSDFactor();
		return (xsdf > upperSigmaThreshold || xsdf < lowerSigmaThreshold);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.MultiFilterComponent#getType()
	 */
	public int getType()
	{
		return IDirectFilter.V_X_SD_FACTOR;
	}
}