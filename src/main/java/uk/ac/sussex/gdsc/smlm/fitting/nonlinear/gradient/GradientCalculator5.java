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

import uk.ac.sussex.gdsc.smlm.function.NonLinearFunction;

/**
 * Calculates the Hessian matrix (the square matrix of second-order partial derivatives of a
 * function) and the gradient vector of the function's partial first derivatives with respect to the
 * parameters. This is used within the Levenberg-Marquardt method to fit a nonlinear model with
 * coefficients (a) for a set of data points (x, y).
 */
public class GradientCalculator5 extends GradientCalculator {
  /**
   * Instantiates a new gradient calculator for a 5x5 matrix.
   */
  public GradientCalculator5() {
    super(5);
  }

  @Override
  public double findLinearised(int[] x, double[] y, double[] a, double[][] alpha, double[] beta,
      NonLinearFunction func) {
    double ssx = 0;
    final double[] dyDa = new double[5];

    alpha[0][0] = 0;
    alpha[1][0] = 0;
    alpha[1][1] = 0;
    alpha[2][0] = 0;
    alpha[2][1] = 0;
    alpha[2][2] = 0;
    alpha[3][0] = 0;
    alpha[3][1] = 0;
    alpha[3][2] = 0;
    alpha[3][3] = 0;
    alpha[4][0] = 0;
    alpha[4][1] = 0;
    alpha[4][2] = 0;
    alpha[4][3] = 0;
    alpha[4][4] = 0;

    beta[0] = 0;
    beta[1] = 0;
    beta[2] = 0;
    beta[3] = 0;
    beta[4] = 0;

    func.initialise(a);

    if (func.canComputeWeights()) {
      final double[] wgt = new double[1];
      for (int i = 0; i < x.length; i++) {
        final double dy = y[i] - func.evalw(x[i], dyDa, wgt);
        final double weight = getWeight(wgt[0]);

        alpha[0][0] += dyDa[0] * weight * dyDa[0];
        alpha[1][0] += dyDa[1] * weight * dyDa[0];
        alpha[1][1] += dyDa[1] * weight * dyDa[1];
        alpha[2][0] += dyDa[2] * weight * dyDa[0];
        alpha[2][1] += dyDa[2] * weight * dyDa[1];
        alpha[2][2] += dyDa[2] * weight * dyDa[2];
        alpha[3][0] += dyDa[3] * weight * dyDa[0];
        alpha[3][1] += dyDa[3] * weight * dyDa[1];
        alpha[3][2] += dyDa[3] * weight * dyDa[2];
        alpha[3][3] += dyDa[3] * weight * dyDa[3];
        alpha[4][0] += dyDa[4] * weight * dyDa[0];
        alpha[4][1] += dyDa[4] * weight * dyDa[1];
        alpha[4][2] += dyDa[4] * weight * dyDa[2];
        alpha[4][3] += dyDa[4] * weight * dyDa[3];
        alpha[4][4] += dyDa[4] * weight * dyDa[4];

        beta[0] += dyDa[0] * weight * dy;
        beta[1] += dyDa[1] * weight * dy;
        beta[2] += dyDa[2] * weight * dy;
        beta[3] += dyDa[3] * weight * dy;
        beta[4] += dyDa[4] * weight * dy;

        ssx += dy * dy * weight;
      }
    } else {
      for (int i = 0; i < x.length; i++) {
        final double dy = y[i] - func.eval(x[i], dyDa);

        alpha[0][0] += dyDa[0] * dyDa[0];
        alpha[1][0] += dyDa[1] * dyDa[0];
        alpha[1][1] += dyDa[1] * dyDa[1];
        alpha[2][0] += dyDa[2] * dyDa[0];
        alpha[2][1] += dyDa[2] * dyDa[1];
        alpha[2][2] += dyDa[2] * dyDa[2];
        alpha[3][0] += dyDa[3] * dyDa[0];
        alpha[3][1] += dyDa[3] * dyDa[1];
        alpha[3][2] += dyDa[3] * dyDa[2];
        alpha[3][3] += dyDa[3] * dyDa[3];
        alpha[4][0] += dyDa[4] * dyDa[0];
        alpha[4][1] += dyDa[4] * dyDa[1];
        alpha[4][2] += dyDa[4] * dyDa[2];
        alpha[4][3] += dyDa[4] * dyDa[3];
        alpha[4][4] += dyDa[4] * dyDa[4];

        beta[0] += dyDa[0] * dy;
        beta[1] += dyDa[1] * dy;
        beta[2] += dyDa[2] * dy;
        beta[3] += dyDa[3] * dy;
        beta[4] += dyDa[4] * dy;

        ssx += dy * dy;
      }
    }

    // Generate symmetric matrix
    alpha[0][1] = alpha[1][0];
    alpha[0][2] = alpha[2][0];
    alpha[0][3] = alpha[3][0];
    alpha[0][4] = alpha[4][0];
    alpha[1][2] = alpha[2][1];
    alpha[1][3] = alpha[3][1];
    alpha[1][4] = alpha[4][1];
    alpha[2][3] = alpha[3][2];
    alpha[2][4] = alpha[4][2];
    alpha[3][4] = alpha[4][3];

    return checkGradients(alpha, beta, nparams, ssx);
  }

  @Override
  public double findLinearised(int n, double[] y, double[] a, double[][] alpha, double[] beta,
      NonLinearFunction func) {
    double ssx = 0;
    final double[] dyDa = new double[5];

    alpha[0][0] = 0;
    alpha[1][0] = 0;
    alpha[1][1] = 0;
    alpha[2][0] = 0;
    alpha[2][1] = 0;
    alpha[2][2] = 0;
    alpha[3][0] = 0;
    alpha[3][1] = 0;
    alpha[3][2] = 0;
    alpha[3][3] = 0;
    alpha[4][0] = 0;
    alpha[4][1] = 0;
    alpha[4][2] = 0;
    alpha[4][3] = 0;
    alpha[4][4] = 0;

    beta[0] = 0;
    beta[1] = 0;
    beta[2] = 0;
    beta[3] = 0;
    beta[4] = 0;

    func.initialise(a);

    if (func.canComputeWeights()) {
      final double[] wgt = new double[1];
      for (int i = 0; i < n; i++) {
        final double dy = y[i] - func.evalw(i, dyDa, wgt);
        final double weight = getWeight(wgt[0]);

        alpha[0][0] += dyDa[0] * weight * dyDa[0];
        alpha[1][0] += dyDa[1] * weight * dyDa[0];
        alpha[1][1] += dyDa[1] * weight * dyDa[1];
        alpha[2][0] += dyDa[2] * weight * dyDa[0];
        alpha[2][1] += dyDa[2] * weight * dyDa[1];
        alpha[2][2] += dyDa[2] * weight * dyDa[2];
        alpha[3][0] += dyDa[3] * weight * dyDa[0];
        alpha[3][1] += dyDa[3] * weight * dyDa[1];
        alpha[3][2] += dyDa[3] * weight * dyDa[2];
        alpha[3][3] += dyDa[3] * weight * dyDa[3];
        alpha[4][0] += dyDa[4] * weight * dyDa[0];
        alpha[4][1] += dyDa[4] * weight * dyDa[1];
        alpha[4][2] += dyDa[4] * weight * dyDa[2];
        alpha[4][3] += dyDa[4] * weight * dyDa[3];
        alpha[4][4] += dyDa[4] * weight * dyDa[4];

        beta[0] += dyDa[0] * weight * dy;
        beta[1] += dyDa[1] * weight * dy;
        beta[2] += dyDa[2] * weight * dy;
        beta[3] += dyDa[3] * weight * dy;
        beta[4] += dyDa[4] * weight * dy;

        ssx += dy * dy * weight;
      }
    } else {
      for (int i = 0; i < n; i++) {
        final double dy = y[i] - func.eval(i, dyDa);

        alpha[0][0] += dyDa[0] * dyDa[0];
        alpha[1][0] += dyDa[1] * dyDa[0];
        alpha[1][1] += dyDa[1] * dyDa[1];
        alpha[2][0] += dyDa[2] * dyDa[0];
        alpha[2][1] += dyDa[2] * dyDa[1];
        alpha[2][2] += dyDa[2] * dyDa[2];
        alpha[3][0] += dyDa[3] * dyDa[0];
        alpha[3][1] += dyDa[3] * dyDa[1];
        alpha[3][2] += dyDa[3] * dyDa[2];
        alpha[3][3] += dyDa[3] * dyDa[3];
        alpha[4][0] += dyDa[4] * dyDa[0];
        alpha[4][1] += dyDa[4] * dyDa[1];
        alpha[4][2] += dyDa[4] * dyDa[2];
        alpha[4][3] += dyDa[4] * dyDa[3];
        alpha[4][4] += dyDa[4] * dyDa[4];

        beta[0] += dyDa[0] * dy;
        beta[1] += dyDa[1] * dy;
        beta[2] += dyDa[2] * dy;
        beta[3] += dyDa[3] * dy;
        beta[4] += dyDa[4] * dy;

        ssx += dy * dy;
      }
    }

    // Generate symmetric matrix
    alpha[0][1] = alpha[1][0];
    alpha[0][2] = alpha[2][0];
    alpha[0][3] = alpha[3][0];
    alpha[0][4] = alpha[4][0];
    alpha[1][2] = alpha[2][1];
    alpha[1][3] = alpha[3][1];
    alpha[1][4] = alpha[4][1];
    alpha[2][3] = alpha[3][2];
    alpha[2][4] = alpha[4][2];
    alpha[3][4] = alpha[4][3];

    return checkGradients(alpha, beta, nparams, ssx);
  }
}
