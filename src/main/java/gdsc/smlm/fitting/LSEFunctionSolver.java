/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 * 
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm.fitting;


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
	 * Gets the adjusted coefficient of determination (Adjusted R^2 = 1 - [SSresiduals / SStotal] * [[n - 1] / [n - p -
	 * 1]])
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
