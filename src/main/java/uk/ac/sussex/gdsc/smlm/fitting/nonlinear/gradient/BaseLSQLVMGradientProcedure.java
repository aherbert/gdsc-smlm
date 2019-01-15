/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient;

import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;

/**
 * Calculates the scaled Hessian matrix (the square matrix of second-order partial derivatives of a
 * function) and the scaled gradient vector of the function's partial first derivatives with respect
 * to the parameters. This is used within the Levenberg-Marquardt method to fit a nonlinear model
 * with coefficients (a) for a set of data points (x, y).
 *
 * <p>Note that the Hessian matrix is scaled by 1/2 and the gradient vector is scaled by -1/2 for
 * convenience in solving the non-linear model. See Numerical Recipes in C++, 2nd Ed. Equation
 * 15.5.8 for Nonlinear Models.
 */
public abstract class BaseLSQLVMGradientProcedure extends LVMGradientProcedure {

  /**
   * Instantiates a new procedure.
   *
   * @param y Data to fit
   * @param func Gradient function
   */
  public BaseLSQLVMGradientProcedure(final double[] y, final Gradient1Function func) {
    super(y, func);
  }

  /**
   * Instantiates a new procedure.
   *
   * @param y Data to fit
   * @param b Baseline pre-computed y-values
   * @param func Gradient function
   */
  public BaseLSQLVMGradientProcedure(final double[] y, final double[] b,
      final Gradient1Function func) {
    super(y, b, func);
  }

  @Override
  public void execute(double value) {
    // Produce a sum-of-squares
    final double dy = y[++yi] - value;
    this.value += dy * dy;
  }

  @Override
  protected void initialiseValue() {
    // Do nothing
  }

  @Override
  protected void finishValue() {
    // Do nothing
  }
}
