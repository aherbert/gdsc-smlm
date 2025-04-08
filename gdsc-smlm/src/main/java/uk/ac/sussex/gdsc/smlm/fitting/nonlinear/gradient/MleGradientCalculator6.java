/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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

/**
 * Calculates the Hessian matrix (the square matrix of second-order partial derivatives of a
 * function) and the gradient vector of the function's partial first derivatives with respect to the
 * parameters. This is used within the Levenberg-Marquardt method to fit a nonlinear model with
 * coefficients (a) for a set of data points (x, y).
 *
 * <p>This calculator computes a modified Chi-squared expression to perform Maximum Likelihood
 * Estimation assuming Poisson model. See Laurence &amp; Chromy (2010) Efficient maximum likelihood
 * estimator. Nature Methods 7, 338-339. The input data must be Poisson distributed for this to be
 * relevant.
 */
public class MleGradientCalculator6 extends MleGradientCalculator {
  /**
   * Instantiates a new MLE gradient calculator.
   */
  public MleGradientCalculator6() {
    super(6);
  }

  @Override
  protected void zero(final double[][] alpha, final double[] beta) {
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

    beta[0] = 0;
    beta[1] = 0;
    beta[2] = 0;
    beta[3] = 0;
    beta[4] = 0;
    beta[5] = 0;
  }

  @Override
  protected void compute(final double[][] alpha, final double[] beta, final double[] dfiDa,
      final double fi, final double xi) {
    final double xi_fi = xi / fi;
    final double xi_fi2 = xi_fi / fi;
    final double e = 1 - (xi_fi);

    alpha[0][0] += dfiDa[0] * xi_fi2 * dfiDa[0];
    double wgt;
    wgt = dfiDa[1] * xi_fi2;
    alpha[1][0] += wgt * dfiDa[0];
    alpha[1][1] += wgt * dfiDa[1];
    wgt = dfiDa[2] * xi_fi2;
    alpha[2][0] += wgt * dfiDa[0];
    alpha[2][1] += wgt * dfiDa[1];
    alpha[2][2] += wgt * dfiDa[2];
    wgt = dfiDa[3] * xi_fi2;
    alpha[3][0] += wgt * dfiDa[0];
    alpha[3][1] += wgt * dfiDa[1];
    alpha[3][2] += wgt * dfiDa[2];
    alpha[3][3] += wgt * dfiDa[3];
    wgt = dfiDa[4] * xi_fi2;
    alpha[4][0] += wgt * dfiDa[0];
    alpha[4][1] += wgt * dfiDa[1];
    alpha[4][2] += wgt * dfiDa[2];
    alpha[4][3] += wgt * dfiDa[3];
    alpha[4][4] += wgt * dfiDa[4];
    wgt = dfiDa[5] * xi_fi2;
    alpha[5][0] += wgt * dfiDa[0];
    alpha[5][1] += wgt * dfiDa[1];
    alpha[5][2] += wgt * dfiDa[2];
    alpha[5][3] += wgt * dfiDa[3];
    alpha[5][4] += wgt * dfiDa[4];
    alpha[5][5] += wgt * dfiDa[5];

    beta[0] -= e * dfiDa[0];
    beta[1] -= e * dfiDa[1];
    beta[2] -= e * dfiDa[2];
    beta[3] -= e * dfiDa[3];
    beta[4] -= e * dfiDa[4];
    beta[5] -= e * dfiDa[5];
  }

  @Override
  protected void symmetric(final double[][] alpha) {
    alpha[0][1] = alpha[1][0];
    alpha[0][2] = alpha[2][0];
    alpha[0][3] = alpha[3][0];
    alpha[0][4] = alpha[4][0];
    alpha[0][5] = alpha[5][0];
    alpha[1][2] = alpha[2][1];
    alpha[1][3] = alpha[3][1];
    alpha[1][4] = alpha[4][1];
    alpha[1][5] = alpha[5][1];
    alpha[2][3] = alpha[3][2];
    alpha[2][4] = alpha[4][2];
    alpha[2][5] = alpha[5][2];
    alpha[3][4] = alpha[4][3];
    alpha[3][5] = alpha[5][3];
    alpha[4][5] = alpha[5][4];
  }
}
