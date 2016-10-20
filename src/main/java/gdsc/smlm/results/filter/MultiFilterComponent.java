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
 * Filter results using a single component of the multi filter.
 */
public abstract class MultiFilterComponent
{
	/**
	 * Validate the peak
	 *
	 * @param peak
	 *            the peak
	 * @return true, if it fails the filter
	 */
	public abstract boolean fail(final PreprocessedPeakResult peak);

	/**
	 * Gets the type of the component. The return value will match the constants defined in IDirectFilter.
	 *
	 * @return the type
	 */
	public abstract int getType();
}