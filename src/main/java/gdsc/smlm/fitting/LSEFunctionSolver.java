package gdsc.smlm.fitting;

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
 * Defines methods to fit a function with coefficients (a) using least-squares estimation.
 */
public interface LSEFunctionSolver extends FunctionSolver
{
	/**
	 * Gets the total sum of squares.
	 *
	 * @return the total sum of squares
	 */
	public double getTotalSumOfSquares();

	/**
	 * Gets the residual sum of squares.
	 *
	 * @return the residual sum of squares
	 */
	public double getResidualSumOfSquares();

	/**
	 * Gets the coefficient of determination (R^2 = 1 - SSresiduals / SStotal).
	 *
	 * @return the coefficient of determination
	 */
	public double getCoefficientOfDetermination();

	/**
	 * Gets the adjusted coefficient of determination (Adjusted R^2 = 1 - [SSresiduals / SStotal] * [[n - 1] / [n - p - 1]])
	 *
	 * @return the adjusted coefficient of determination
	 */
	public double getAdjustedCoefficientOfDetermination();

	/**
	 * Gets the mean squared error. This is the residual sum of squares divided by the degrees of freedom.
	 *
	 * @return the mean squared error
	 */
	public double getMeanSquaredError();
}
