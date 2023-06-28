/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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
public class LsqLvmGradientProcedureMatrix extends BaseLsqLvmGradientProcedure {

  /** The scaled Hessian curvature matrix (size n*n). */
  public final double[][] alpha;

  /**
   * Instantiates a new procedure.
   *
   * @param y Data to fit
   * @param baseline Baseline pre-computed y-values
   * @param func Gradient function
   */
  public LsqLvmGradientProcedureMatrix(final double[] y, final double[] baseline,
      final Gradient1Function func) {
    super(y, baseline, func);
    alpha = new double[numberOfGradients][numberOfGradients];
  }

  @Override
  public void execute(double value, double[] dyDa) {
    final double dy = y[++yi] - value;

    // Compute:
    // - the scaled Hessian matrix (the square matrix of second-order partial derivatives of a
    // function;
    // that is, it describes the local curvature of a function of many variables.)
    // - the scaled gradient vector of the function's partial first derivatives with respect to the
    // parameters

    for (int j = 0; j < numberOfGradients; j++) {
      final double wgt = dyDa[j];

      for (int k = 0; k <= j; k++) {
        alpha[j][k] += wgt * dyDa[k];
      }

      beta[j] += wgt * dy;
    }

    this.value += dy * dy;
  }

  @Override
  protected void initialiseGradient() {
    for (int i = 0; i < numberOfGradients; i++) {
      beta[i] = 0;
      for (int j = 0; j <= i; j++) {
        alpha[i][j] = 0;
      }
    }
  }

  @Override
  protected void finishGradient() {
    // Generate symmetric matrix
    for (int i = 0; i < numberOfGradients - 1; i++) {
      for (int j = i + 1; j < numberOfGradients; j++) {
        alpha[i][j] = alpha[j][i];
      }
    }
  }

  @Override
  protected boolean checkGradients() {
    for (int i = 0; i < numberOfGradients; i++) {
      if (Double.isNaN(beta[i])) {
        return true;
      }
      for (int j = 0; j <= i; j++) {
        if (Double.isNaN(alpha[i][j])) {
          return true;
        }
      }
    }
    return false;
  }

  @Override
  public double[][] getAlphaMatrix() {
    return alpha;
  }

  @Override
  public void getAlphaMatrix(double[][] alpha) {
    for (int i = 0; i < numberOfGradients; i++) {
      System.arraycopy(this.alpha, 0, alpha, 0, numberOfGradients);
    }
  }

  @Override
  public void getAlphaLinear(double[] alpha) {
    toLinear(this.alpha, alpha);
  }
}
