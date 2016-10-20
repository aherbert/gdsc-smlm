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
public class MultiFilterShiftComponent extends MultiFilterComponent
{
	final static int type = IDirectFilter.V_X_RELATIVE_SHIFT | IDirectFilter.V_Y_RELATIVE_SHIFT;
	
	final float offset;

	public MultiFilterShiftComponent(double shift)
	{
		this.offset = Filter.getUpperSquaredLimit(shift);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.MultiFilterComponent#getType()
	 */
	public boolean fail(final PreprocessedPeakResult peak)
	{
		if (peak.getXRelativeShift2() > offset)
			return true;
		return (peak.getYRelativeShift2() > offset);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.MultiFilterComponent#getType()
	 */
	public int getType()
	{
		return type;
	}
}