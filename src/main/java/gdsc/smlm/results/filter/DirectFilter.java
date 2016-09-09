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
 * Support direct filtering of PreprocessedPeakResult objects.
 * <p>
 * The decision to support for filtering as both a DirectFilter and Filter at the same time is left to the implementing
 * class. It is not a requirement.
 */
public abstract class DirectFilter extends Filter
{
	/**
	 * Disable filtering using the width of the result
	 */
	public static final int NO_WIDTH = 1;

	/**
	 * Called before the accept method is called for PreprocessedPeakResult
	 * <p>
	 * This should be called once to initialise the filter before processing a batch of results.
	 * 
	 * @see #accept(PreprocessedPeakResult)
	 */
	public void setup()
	{
	}

	/**
	 * Called before the accept method is called for PreprocessedPeakResult. the flags can control the type of filtering
	 * requested. Filters are asked to respect the flags defined in this class.
	 * <p>
	 * This should be called once to initialise the filter before processing a batch of results.
	 * 
	 * @param flags
	 *            Flags used to control the filter
	 * @see #accept(PreprocessedPeakResult)
	 */
	public void setup(final int flags)
	{
	}

	/**
	 * Check if the given bits are set in the flags
	 * 
	 * @param flags
	 * @param bits
	 * @return True if all are set
	 */
	public static boolean areSet(final int flags, final int bits)
	{
		return (flags & bits) == bits;
	}

	/**
	 * Filter the peak result.
	 * 
	 * @param peak
	 *            The peak result
	 * @return true if the peak should be accepted, otherwise false to reject.
	 */
	public abstract boolean accept(final PreprocessedPeakResult peak);
}