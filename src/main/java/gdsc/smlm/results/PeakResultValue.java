package gdsc.smlm.results;

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
 * Gets a value from a peak result.
 */
public interface PeakResultValue
{
	/**
	 * Gets the value of the result.
	 *
	 * @param result
	 *            the result
	 * @return the value
	 */
	public float getValue(PeakResult result);
}