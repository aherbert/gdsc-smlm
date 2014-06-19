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
 * Define the 2D Gaussian fitting function
 */
public enum FitFunction
{
	/**
	 * Fixed width fitting
	 */
	FIXED,
	/**
	 * Fit XY width simultaneously
	 */
	CIRCULAR,
	/**
	 * Fit XY width independently
	 */
	FREE_CIRCULAR,
	/**
	 * Fit elliptical Gaussian
	 */
	FREE
}