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
public abstract class DirectFilter extends Filter implements IDirectFilter
{
	private int result = 0;

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#setup()
	 */
	public void setup()
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#setup(int)
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#accept(gdsc.smlm.results.filter.PreprocessedPeakResult)
	 */
	final public boolean accept(final PreprocessedPeakResult peak)
	{
		return (result = validate(peak)) == 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#validate(gdsc.smlm.results.filter.PreprocessedPeakResult)
	 */
	public abstract int validate(final PreprocessedPeakResult peak);

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.Filter#getFilterType()
	 */
	public FilterType getFilterType()
	{
		return FilterType.DIRECT;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.filter.IDirectFilter#getResult()
	 */
	public int getResult()
	{
		return result;
	}
}