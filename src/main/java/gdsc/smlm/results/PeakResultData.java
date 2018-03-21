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
 * Gets a data value from a peak result.
 */
public interface PeakResultData<E>
{
	/**
	 * Gets the value of the result.
	 *
	 * @param result
	 *            the result
	 * @return the value
	 */
	public E getValue(PeakResult result);
	
	/**
	 * Gets the name of the value.
	 *
	 * @return the name
	 */
	public String getValueName();

	/**
	 * Gets the class type of the value.
	 *
	 * @return the name
	 */
	public Class<?> getValueClass();
}