/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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

import uk.ac.sussex.gdsc.core.utils.ValidationUtils;
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
public class LsqLvmGradientProcedureLinear4 extends LsqLvmGradientProcedureLinear {

  /**
   * Instantiates a new procedure.
   *
   * @param y Data to fit
   * @param baseline Baseline pre-computed y-values
   * @param func Gradient function
   */
  public LsqLvmGradientProcedureLinear4(final double[] y, final double[] baseline,
      final Gradient1Function func) {
    super(y, baseline, func);
    ValidationUtils.checkArgument(numberOfGradients == 4, "Function must compute 4 gradients");
  }

  @Override
  public void execute(double value, double[] dyDa) {
    final double dy = y[++yi] - value;

    alpha[0] += dyDa[0] * dyDa[0];
    alpha[1] += dyDa[0] * dyDa[1];
    alpha[2] += dyDa[0] * dyDa[2];
    alpha[3] += dyDa[0] * dyDa[3];
    alpha[5] += dyDa[1] * dyDa[1];
    alpha[6] += dyDa[1] * dyDa[2];
    alpha[7] += dyDa[1] * dyDa[3];
    alpha[10] += dyDa[2] * dyDa[2];
    alpha[11] += dyDa[2] * dyDa[3];
    alpha[15] += dyDa[3] * dyDa[3];

    beta[0] += dyDa[0] * dy;
    beta[1] += dyDa[1] * dy;
    beta[2] += dyDa[2] * dy;
    beta[3] += dyDa[3] * dy;

    this.value += dy * dy;
  }

  @Override
  protected void initialiseGradient() {
    alpha[0] = 0;
    alpha[1] = 0;
    alpha[2] = 0;
    alpha[3] = 0;
    alpha[5] = 0;
    alpha[6] = 0;
    alpha[7] = 0;
    alpha[10] = 0;
    alpha[11] = 0;
    alpha[15] = 0;

    beta[0] = 0;
    beta[1] = 0;
    beta[2] = 0;
    beta[3] = 0;
  }

  @Override
  protected void finishGradient() {
    alpha[4] = alpha[1];
    alpha[8] = alpha[2];
    alpha[12] = alpha[3];
    alpha[9] = alpha[6];
    alpha[13] = alpha[7];
    alpha[14] = alpha[11];
  }
}
