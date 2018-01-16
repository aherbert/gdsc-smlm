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
 * Contains a set of components of the multi filter.
 */
public abstract class MultiFilterComponentSet
{
	/**
	 * Gets the validation flags. These are possible return flags from the {@link #validate(PreprocessedPeakResult)}
	 * method.
	 *
	 * @return the validation flags
	 */
	abstract public int getValidationFlags();

	/**
	 * Validate the peak
	 *
	 * @param peak
	 *            the peak
	 * @return the result
	 */
	public abstract int validate(final PreprocessedPeakResult peak);

	/**
	 * Replace the first component.
	 *
	 * @param c
	 *            the replacement component
	 */
	abstract void replace0(MultiFilterComponent c);
}