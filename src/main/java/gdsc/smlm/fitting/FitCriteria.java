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
 * Define the criteria used to stop the fitting process
 */
public enum FitCriteria
{
	//@formatter:off
	/**
	 * Stop fitting when the least-squared-error does not change significantly
	 */
	LEAST_SQUARED_ERROR{ public String getName() { return "Least-squared error"; }},
	/**
	 * Stop fitting when the least-squared-error does not change significantly.
	 * Add an additional check for slowly improving or plateaued fitting is performed.
	 */
	LEAST_SQUARED_PLUS{ public String getName() { return "Least-squared error plateau"; }},
	/**
	 * Stop fitting when the XY coordinates do not change significantly
	 */
	COORDINATES{ public String getName() { return "Coordinates"; }},
	/**
	 * Stop fitting when the parameters do not change significantly
	 */
	PARAMETERS{ public String getName() { return "Parameters"; }};
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