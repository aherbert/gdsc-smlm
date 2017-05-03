package gdsc.smlm.function;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Interface for functions to produce a value, first and second partial derivatives
 */
public interface Gradient2Procedure
{
	/**
	 * Executes this procedure.
	 *
	 * @param value
	 *            the value of the function
	 * @param dy_da
	 *            Partial first derivative of function with respect to each coefficient identified by
	 *            {@link #gradientIndices()}
	 * @param dy_da
	 *            Partial second derivative of function with respect to each coefficient identified by
	 *            {@link #gradientIndices()}
	 */
	void execute(double value, double[] dy_da, double[] d2y_da2);
}