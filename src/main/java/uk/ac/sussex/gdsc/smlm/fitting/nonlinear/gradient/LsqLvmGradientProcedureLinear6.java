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
 * Calculates the Hessian matrix (the square matrix of second-order partial derivatives of a
 * function) and the scaled gradient vector of the function's partial first derivatives with respect
 * to the parameters. This is used within the Levenberg-Marquardt method to fit a nonlinear model
 * with coefficients (a) for a set of data points (x, y).
 *
 * <p>Note that the Hessian matrix is scaled by 1/2 and the gradient vector is scaled by -1/2 for
 * convenience in solving the non-linear model. See Numerical Recipes in C++, 2nd Ed. Equation
 * 15.5.8 for Nonlinear Models.
 */
public class LsqLvmGradientProcedureLinear6 extends LsqLvmGradientProcedureLinear {

  /**
   * Instantiates a new procedure.
   *
   * @param y Data to fit
   * @param baseline Baseline pre-computed y-values
   * @param func Gradient function
   */
  public LsqLvmGradientProcedureLinear6(final double[] y, final double[] baseline,
      final Gradient1Function func) {
    super(y, baseline, func);
    if (numberOfGradients != 6) {
      throw new IllegalArgumentException("Function must compute 6 gradients");
    }
  }

  @Override
  public void execute(double value, double[] dyDa) {
    final double dy = y[++yi] - value;

    alpha[0] += dyDa[0] * dyDa[0];
    alpha[1] += dyDa[0] * dyDa[1];
    alpha[2] += dyDa[0] * dyDa[2];
    alpha[3] += dyDa[0] * dyDa[3];
    alpha[4] += dyDa[0] * dyDa[4];
    alpha[5] += dyDa[0] * dyDa[5];
    alpha[7] += dyDa[1] * dyDa[1];
    alpha[8] += dyDa[1] * dyDa[2];
    alpha[9] += dyDa[1] * dyDa[3];
    alpha[10] += dyDa[1] * dyDa[4];
    alpha[11] += dyDa[1] * dyDa[5];
    alpha[14] += dyDa[2] * dyDa[2];
    alpha[15] += dyDa[2] * dyDa[3];
    alpha[16] += dyDa[2] * dyDa[4];
    alpha[17] += dyDa[2] * dyDa[5];
    alpha[21] += dyDa[3] * dyDa[3];
    alpha[22] += dyDa[3] * dyDa[4];
    alpha[23] += dyDa[3] * dyDa[5];
    alpha[28] += dyDa[4] * dyDa[4];
    alpha[29] += dyDa[4] * dyDa[5];
    alpha[35] += dyDa[5] * dyDa[5];

    beta[0] += dyDa[0] * dy;
    beta[1] += dyDa[1] * dy;
    beta[2] += dyDa[2] * dy;
    beta[3] += dyDa[3] * dy;
    beta[4] += dyDa[4] * dy;
    beta[5] += dyDa[5] * dy;

    this.value += dy * dy;
  }

  @Override
  protected void initialiseGradient() {
    alpha[0] = 0;
    alpha[1] = 0;
    alpha[2] = 0;
    alpha[3] = 0;
    alpha[4] = 0;
    alpha[5] = 0;
    alpha[7] = 0;
    alpha[8] = 0;
    alpha[9] = 0;
    alpha[10] = 0;
    alpha[11] = 0;
    alpha[14] = 0;
    alpha[15] = 0;
    alpha[16] = 0;
    alpha[17] = 0;
    alpha[21] = 0;
    alpha[22] = 0;
    alpha[23] = 0;
    alpha[28] = 0;
    alpha[29] = 0;
    alpha[35] = 0;
    beta[0] = 0;
    beta[1] = 0;
    beta[2] = 0;
    beta[3] = 0;
    beta[4] = 0;
    beta[5] = 0;
  }

  @Override
  protected void finishGradient() {
    alpha[6] = alpha[1];
    alpha[12] = alpha[2];
    alpha[18] = alpha[3];
    alpha[24] = alpha[4];
    alpha[30] = alpha[5];
    alpha[13] = alpha[8];
    alpha[19] = alpha[9];
    alpha[25] = alpha[10];
    alpha[31] = alpha[11];
    alpha[20] = alpha[15];
    alpha[26] = alpha[16];
    alpha[32] = alpha[17];
    alpha[27] = alpha[22];
    alpha[33] = alpha[23];
    alpha[34] = alpha[29];
  }
}
