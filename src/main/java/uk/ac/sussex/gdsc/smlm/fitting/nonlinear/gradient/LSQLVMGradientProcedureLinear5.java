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

package uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient;

import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;

/**
 * Calculates the Hessian matrix (the square matrix of second-order partial derivatives of a
 * function) and the scaled gradient vector of the function's partial first derivatives with respect
 * to the parameters. This is used within the Levenberg-Marquardt method to fit a nonlinear model
 * with coefficients (a) for a set of data points (x, y).
 *
 * <p>Note that the Hessian matrix is scaled by 1/2 and the gradient vector is scaled by -1/2 for
 * convenience in solving the non-linear model. See Numerical Recipes in C++, 2nd Ed. Equation
 * 15.5.8 for Nonlinear Models.
 */
public class LSQLVMGradientProcedureLinear5 extends LSQLVMGradientProcedureLinear {
  /**
   * @param y Data to fit
   * @param b Baseline pre-computed y-values
   * @param func Gradient function
   */
  public LSQLVMGradientProcedureLinear5(final double[] y, final double[] b,
      final Gradient1Function func) {
    super(y, b, func);
    if (n != 5) {
      throw new IllegalArgumentException("Function must compute 5 gradients");
    }
  }

  /** {@inheritDoc} */
  @Override
  public void execute(double value, double[] dy_da) {
    final double dy = y[++yi] - value;

    alpha[0] += dy_da[0] * dy_da[0];
    alpha[1] += dy_da[0] * dy_da[1];
    alpha[2] += dy_da[0] * dy_da[2];
    alpha[3] += dy_da[0] * dy_da[3];
    alpha[4] += dy_da[0] * dy_da[4];
    alpha[6] += dy_da[1] * dy_da[1];
    alpha[7] += dy_da[1] * dy_da[2];
    alpha[8] += dy_da[1] * dy_da[3];
    alpha[9] += dy_da[1] * dy_da[4];
    alpha[12] += dy_da[2] * dy_da[2];
    alpha[13] += dy_da[2] * dy_da[3];
    alpha[14] += dy_da[2] * dy_da[4];
    alpha[18] += dy_da[3] * dy_da[3];
    alpha[19] += dy_da[3] * dy_da[4];
    alpha[24] += dy_da[4] * dy_da[4];

    beta[0] += dy_da[0] * dy;
    beta[1] += dy_da[1] * dy;
    beta[2] += dy_da[2] * dy;
    beta[3] += dy_da[3] * dy;
    beta[4] += dy_da[4] * dy;

    this.value += dy * dy;
  }

  @Override
  protected void initialiseGradient() {
    alpha[0] = 0;
    alpha[1] = 0;
    alpha[2] = 0;
    alpha[3] = 0;
    alpha[4] = 0;
    alpha[6] = 0;
    alpha[7] = 0;
    alpha[8] = 0;
    alpha[9] = 0;
    alpha[12] = 0;
    alpha[13] = 0;
    alpha[14] = 0;
    alpha[18] = 0;
    alpha[19] = 0;
    alpha[24] = 0;
    beta[0] = 0;
    beta[1] = 0;
    beta[2] = 0;
    beta[3] = 0;
    beta[4] = 0;
  }

  @Override
  protected void finishGradient() {
    alpha[5] = alpha[1];
    alpha[10] = alpha[2];
    alpha[15] = alpha[3];
    alpha[20] = alpha[4];
    alpha[11] = alpha[7];
    alpha[16] = alpha[8];
    alpha[21] = alpha[9];
    alpha[17] = alpha[13];
    alpha[22] = alpha[14];
    alpha[23] = alpha[19];
  }
}
