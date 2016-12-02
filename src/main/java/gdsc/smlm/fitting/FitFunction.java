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
 * Define the fitting function
 */
public enum FitFunction
{
	//@formatter:off
	/**
	 * Fixed width 2D Gaussian
	 */
	FIXED{ public String getName() { return "Fixed"; }},
	/**
	 * Fit 2D Gaussian with XY width simultaneously 
	 */
	CIRCULAR{ public String getName() { return "Circular"; }},
	/**
	 * Fit 2D Gaussian with XY width independently
	 */
	FREE_CIRCULAR{ public String getName() { return "Free circular"; }},
	/**
	 * Fit elliptical 2D Gaussian
	 */
	FREE{ public String getName() { return "Free"; }};
	//@formatter:on

	@Override
	public String toString()
	{
		return getName();
	}

	/**
	 * Gets the name.
	 *
	 * @return the name
	 */
	abstract public String getName();
}