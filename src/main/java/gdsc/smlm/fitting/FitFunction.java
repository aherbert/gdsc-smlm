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
	/**
	 * Fixed width 2D Gaussian
	 */
	FIXED("Fixed"),
	/**
	 * Fit 2D Gaussian with XY width simultaneously 
	 */
	CIRCULAR("Circular"),
	/**
	 * Fit 2D Gaussian with XY width independently
	 */
	FREE_CIRCULAR("Free circular"),
	/**
	 * Fit elliptical 2D Gaussian
	 */
	FREE("Free");

	private String name;

	private FitFunction(String name)
	{
		this.name = name;
	}

	@Override
	public String toString()
	{
		return name;
	}
}