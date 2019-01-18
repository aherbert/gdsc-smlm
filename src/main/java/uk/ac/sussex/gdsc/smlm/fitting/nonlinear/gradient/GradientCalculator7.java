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
public class GradientCalculator7 extends GradientCalculator {
  /**
   * Instantiates a new gradient calculator for a 7x7 matrix.
   */
  public GradientCalculator7() {
    super(7);
  }

  @Override
  public double findLinearised(int[] x, double[] y, double[] a, double[][] alpha, double[] beta,
      NonLinearFunction func) {
    double ssx = 0;
    final double[] dyDa = new double[7];

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
    alpha[5][0] = 0;
    alpha[5][1] = 0;
    alpha[5][2] = 0;
    alpha[5][3] = 0;
    alpha[5][4] = 0;
    alpha[5][5] = 0;
    alpha[6][0] = 0;
    alpha[6][1] = 0;
    alpha[6][2] = 0;
    alpha[6][3] = 0;
    alpha[6][4] = 0;
    alpha[6][5] = 0;
    alpha[6][6] = 0;

    beta[0] = 0;
    beta[1] = 0;
    beta[2] = 0;
    beta[3] = 0;
    beta[4] = 0;
    beta[5] = 0;
    beta[6] = 0;

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
        alpha[5][0] += dyDa[5] * weight * dyDa[0];
        alpha[5][1] += dyDa[5] * weight * dyDa[1];
        alpha[5][2] += dyDa[5] * weight * dyDa[2];
        alpha[5][3] += dyDa[5] * weight * dyDa[3];
        alpha[5][4] += dyDa[5] * weight * dyDa[4];
        alpha[5][5] += dyDa[5] * weight * dyDa[5];
        alpha[6][0] += dyDa[6] * weight * dyDa[0];
        alpha[6][1] += dyDa[6] * weight * dyDa[1];
        alpha[6][2] += dyDa[6] * weight * dyDa[2];
        alpha[6][3] += dyDa[6] * weight * dyDa[3];
        alpha[6][4] += dyDa[6] * weight * dyDa[4];
        alpha[6][5] += dyDa[6] * weight * dyDa[5];
        alpha[6][6] += dyDa[6] * weight * dyDa[6];

        beta[0] += dyDa[0] * weight * dy;
        beta[1] += dyDa[1] * weight * dy;
        beta[2] += dyDa[2] * weight * dy;
        beta[3] += dyDa[3] * weight * dy;
        beta[4] += dyDa[4] * weight * dy;
        beta[5] += dyDa[5] * weight * dy;
        beta[6] += dyDa[6] * weight * dy;

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
        alpha[5][0] += dyDa[5] * dyDa[0];
        alpha[5][1] += dyDa[5] * dyDa[1];
        alpha[5][2] += dyDa[5] * dyDa[2];
        alpha[5][3] += dyDa[5] * dyDa[3];
        alpha[5][4] += dyDa[5] * dyDa[4];
        alpha[5][5] += dyDa[5] * dyDa[5];
        alpha[6][0] += dyDa[6] * dyDa[0];
        alpha[6][1] += dyDa[6] * dyDa[1];
        alpha[6][2] += dyDa[6] * dyDa[2];
        alpha[6][3] += dyDa[6] * dyDa[3];
        alpha[6][4] += dyDa[6] * dyDa[4];
        alpha[6][5] += dyDa[6] * dyDa[5];
        alpha[6][6] += dyDa[6] * dyDa[6];

        beta[0] += dyDa[0] * dy;
        beta[1] += dyDa[1] * dy;
        beta[2] += dyDa[2] * dy;
        beta[3] += dyDa[3] * dy;
        beta[4] += dyDa[4] * dy;
        beta[5] += dyDa[5] * dy;
        beta[6] += dyDa[6] * dy;

        ssx += dy * dy;
      }
    }

    // Generate symmetric matrix
    // for (int i = 0; i < m - 1; i++)
    // for (int j = i + 1; j < m; j++)
    // alpha[i][j] = alpha[j][i];

    alpha[0][1] = alpha[1][0];
    alpha[0][2] = alpha[2][0];
    alpha[0][3] = alpha[3][0];
    alpha[0][4] = alpha[4][0];
    alpha[0][5] = alpha[5][0];
    alpha[0][6] = alpha[6][0];
    alpha[1][2] = alpha[2][1];
    alpha[1][3] = alpha[3][1];
    alpha[1][4] = alpha[4][1];
    alpha[1][5] = alpha[5][1];
    alpha[1][6] = alpha[6][1];
    alpha[2][3] = alpha[3][2];
    alpha[2][4] = alpha[4][2];
    alpha[2][5] = alpha[5][2];
    alpha[2][6] = alpha[6][2];
    alpha[3][4] = alpha[4][3];
    alpha[3][5] = alpha[5][3];
    alpha[3][6] = alpha[6][3];
    alpha[4][5] = alpha[5][4];
    alpha[4][6] = alpha[6][4];
    alpha[5][6] = alpha[6][5];

    return checkGradients(alpha, beta, nparams, ssx);
  }

  @Override
  public double findLinearised(int n, double[] y, double[] a, double[][] alpha, double[] beta,
      NonLinearFunction func) {
    double ssx = 0;
    final double[] dyDa = new double[7];

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
    alpha[5][0] = 0;
    alpha[5][1] = 0;
    alpha[5][2] = 0;
    alpha[5][3] = 0;
    alpha[5][4] = 0;
    alpha[5][5] = 0;
    alpha[6][0] = 0;
    alpha[6][1] = 0;
    alpha[6][2] = 0;
    alpha[6][3] = 0;
    alpha[6][4] = 0;
    alpha[6][5] = 0;
    alpha[6][6] = 0;

    beta[0] = 0;
    beta[1] = 0;
    beta[2] = 0;
    beta[3] = 0;
    beta[4] = 0;
    beta[5] = 0;
    beta[6] = 0;

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
        alpha[5][0] += dyDa[5] * weight * dyDa[0];
        alpha[5][1] += dyDa[5] * weight * dyDa[1];
        alpha[5][2] += dyDa[5] * weight * dyDa[2];
        alpha[5][3] += dyDa[5] * weight * dyDa[3];
        alpha[5][4] += dyDa[5] * weight * dyDa[4];
        alpha[5][5] += dyDa[5] * weight * dyDa[5];
        alpha[6][0] += dyDa[6] * weight * dyDa[0];
        alpha[6][1] += dyDa[6] * weight * dyDa[1];
        alpha[6][2] += dyDa[6] * weight * dyDa[2];
        alpha[6][3] += dyDa[6] * weight * dyDa[3];
        alpha[6][4] += dyDa[6] * weight * dyDa[4];
        alpha[6][5] += dyDa[6] * weight * dyDa[5];
        alpha[6][6] += dyDa[6] * weight * dyDa[6];

        beta[0] += dyDa[0] * weight * dy;
        beta[1] += dyDa[1] * weight * dy;
        beta[2] += dyDa[2] * weight * dy;
        beta[3] += dyDa[3] * weight * dy;
        beta[4] += dyDa[4] * weight * dy;
        beta[5] += dyDa[5] * weight * dy;
        beta[6] += dyDa[6] * weight * dy;

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
        alpha[5][0] += dyDa[5] * dyDa[0];
        alpha[5][1] += dyDa[5] * dyDa[1];
        alpha[5][2] += dyDa[5] * dyDa[2];
        alpha[5][3] += dyDa[5] * dyDa[3];
        alpha[5][4] += dyDa[5] * dyDa[4];
        alpha[5][5] += dyDa[5] * dyDa[5];
        alpha[6][0] += dyDa[6] * dyDa[0];
        alpha[6][1] += dyDa[6] * dyDa[1];
        alpha[6][2] += dyDa[6] * dyDa[2];
        alpha[6][3] += dyDa[6] * dyDa[3];
        alpha[6][4] += dyDa[6] * dyDa[4];
        alpha[6][5] += dyDa[6] * dyDa[5];
        alpha[6][6] += dyDa[6] * dyDa[6];

        beta[0] += dyDa[0] * dy;
        beta[1] += dyDa[1] * dy;
        beta[2] += dyDa[2] * dy;
        beta[3] += dyDa[3] * dy;
        beta[4] += dyDa[4] * dy;
        beta[5] += dyDa[5] * dy;
        beta[6] += dyDa[6] * dy;

        ssx += dy * dy;
      }
    }

    // Generate symmetric matrix
    // for (int i = 0; i < m - 1; i++)
    // for (int j = i + 1; j < m; j++)
    // alpha[i][j] = alpha[j][i];

    alpha[0][1] = alpha[1][0];
    alpha[0][2] = alpha[2][0];
    alpha[0][3] = alpha[3][0];
    alpha[0][4] = alpha[4][0];
    alpha[0][5] = alpha[5][0];
    alpha[0][6] = alpha[6][0];
    alpha[1][2] = alpha[2][1];
    alpha[1][3] = alpha[3][1];
    alpha[1][4] = alpha[4][1];
    alpha[1][5] = alpha[5][1];
    alpha[1][6] = alpha[6][1];
    alpha[2][3] = alpha[3][2];
    alpha[2][4] = alpha[4][2];
    alpha[2][5] = alpha[5][2];
    alpha[2][6] = alpha[6][2];
    alpha[3][4] = alpha[4][3];
    alpha[3][5] = alpha[5][3];
    alpha[3][6] = alpha[6][3];
    alpha[4][5] = alpha[5][4];
    alpha[4][6] = alpha[6][4];
    alpha[5][6] = alpha[6][5];

    return checkGradients(alpha, beta, nparams, ssx);
  }

  /**
   * Zero the working region of the input matrix alpha and vector beta.
   *
   * @param alpha the alpha
   * @param beta the beta
   */
  protected void zero(final double[][] alpha, final double[] beta) {
    for (int i = 0; i < nparams; i++) {
      beta[i] = 0;
      for (int j = 0; j <= i; j++) {
        alpha[i][j] = 0;
      }
    }
  }

  /**
   * Compute the matrix alpha and vector beta.
   *
   * @param alpha the alpha
   * @param beta the beta
   * @param dfiDa the gradient of the function with respect to each parameter a
   * @param fi the function value at index i
   * @param xi the data value at index i
   */
  protected void compute(final double[][] alpha, final double[] beta, final double[] dfiDa,
      final double fi, final double xi) {
    final double xi_fi = xi / fi;
    final double xi_fi2 = xi_fi / fi;
    final double e = 1 - (xi_fi);

    // Compute:
    // Laurence &amp; Chromy (2010) Nature Methods 7, 338-339, SI
    // alpha - the Hessian matrix (the square matrix of second-order partial derivatives of a
    // function;
    // that is, it describes the local curvature of a function of many variables.)
    // beta - the gradient vector of the function's partial first derivatives with respect to the
    // parameters

    for (int k = 0; k < nparams; k++) {
      final double wgt = dfiDa[k] * xi_fi2;

      for (int l = 0; l <= k; l++) {
        // This is the non-optimised version:
        // alpha[j][k] += dyDa[j] * dyDa[k] * y[i] / (ymod * ymod);
        alpha[k][l] += wgt * dfiDa[l];
      }

      // This is the non-optimised version:
      // beta[j] -= (1 - y[i] / ymod) * dyDa[j];
      beta[k] -= e * dfiDa[k];
    }
  }

  /**
   * Generate a symmetric matrix alpha.
   *
   * @param alpha the alpha
   */
  protected void symmetric(final double[][] alpha) {
    alpha[0][1] = alpha[1][0];
    alpha[0][2] = alpha[2][0];
    alpha[0][3] = alpha[3][0];
    alpha[0][4] = alpha[4][0];
    alpha[0][5] = alpha[5][0];
    alpha[0][6] = alpha[6][0];
    alpha[1][2] = alpha[2][1];
    alpha[1][3] = alpha[3][1];
    alpha[1][4] = alpha[4][1];
    alpha[1][5] = alpha[5][1];
    alpha[1][6] = alpha[6][1];
    alpha[2][3] = alpha[3][2];
    alpha[2][4] = alpha[4][2];
    alpha[2][5] = alpha[5][2];
    alpha[2][6] = alpha[6][2];
    alpha[3][4] = alpha[4][3];
    alpha[3][5] = alpha[5][3];
    alpha[3][6] = alpha[6][3];
    alpha[4][5] = alpha[5][4];
    alpha[4][6] = alpha[6][4];
    alpha[5][6] = alpha[6][5];
  }
}
