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
 * 
 * This is an adaption of the C-code contained in the CcpNmr Analysis Program:
 * CCPN website (http://www.ccpn.ac.uk/)
 *---------------------------------------------------------------------------*/

/**
 * Define the criteria used to stop the fitting process
 */
public enum FitCriteria
{
	/**
	 * Stop fitting when the least-squared-error does not change significantly
	 */
	LEAST_SQUARED_ERROR,
	/**
	 * Stop fitting when the least-squared-error does not change significantly.
	 * Add an additional check for slowly improving or plateaued fitting is performed.
	 */
	LEAST_SQUARED_PLUS,
	/**
	 * Stop fitting when the XY coordinates do not change significantly
	 */
	COORDINATES,
	/**
	 * Stop fitting when the parameters do not change significantly
	 */
	PARAMETERS
}