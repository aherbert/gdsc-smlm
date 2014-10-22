package gdsc.smlm.fitting;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Define the status of a fit result
 */
public enum FitStatus
{
	/**
	 * 
	 */
	OK,
	/**
	 * 
	 */
	SINGULAR_NON_LINEAR_MODEL,
	/**
	 * 
	 */
	SINGULAR_NON_LINEAR_SOLUTION,
	/**
	 * 
	 */
	INVALID_GRADIENTS_IN_NON_LINEAR_MODEL,
	/**
	 * 
	 */
	FAILED_TO_CONVERGE,
	/**
	 * 
	 */
	BAD_PARAMETERS,
	/**
	 * 
	 */
	COORDINATES_MOVED,
	/**
	 * 
	 */
	INSUFFICIENT_SIGNAL,
	/**
	 * 
	 */
	WIDTH_DIVERGED,
	/**
	 * 
	 */
	INSUFFICIENT_PRECISION,
	/**
	 * 
	 */
	UNKNOWN
}