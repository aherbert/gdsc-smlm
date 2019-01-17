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
 * <p>This procedure computes a modified Chi-squared expression to perform Maximum Likelihood
 * Estimation assuming Poisson model. See Laurence &amp; Chromy (2010) Efficient maximum likelihood
 * estimator. Nature Methods 7, 338-339. The input data must be Poisson distributed for this to be
 * relevant.
 */
public class WLSQLVMGradientProcedure6 extends WLSQLVMGradientProcedure {

  /**
   * Instantiates a new procedure.
   *
   * @param y Data to fit (must be positive)
   * @param var the base variance of each observation (must be positive)
   * @param func Gradient function
   */
  public WLSQLVMGradientProcedure6(final double[] y, final double[] var,
      final Gradient1Function func) {
    super(y, var, func);
    if (n != 6) {
      throw new IllegalArgumentException("Function must compute 6 gradients");
    }
  }

  @Override
  public void execute(double value, double[] dyDa) {
    final double dy = y[++yi] - value;
    final double w = this.w[yi];
    this.value += dy * dy * w;

    double wgt;
    wgt = dyDa[0] * w;
    alpha[0] += wgt * dyDa[0];
    beta[0] += wgt * dy;
    wgt = dyDa[1] * w;
    alpha[1] += wgt * dyDa[0];
    alpha[2] += wgt * dyDa[1];
    beta[1] += wgt * dy;
    wgt = dyDa[2] * w;
    alpha[3] += wgt * dyDa[0];
    alpha[4] += wgt * dyDa[1];
    alpha[5] += wgt * dyDa[2];
    beta[2] += wgt * dy;
    wgt = dyDa[3] * w;
    alpha[6] += wgt * dyDa[0];
    alpha[7] += wgt * dyDa[1];
    alpha[8] += wgt * dyDa[2];
    alpha[9] += wgt * dyDa[3];
    beta[3] += wgt * dy;
    wgt = dyDa[4] * w;
    alpha[10] += wgt * dyDa[0];
    alpha[11] += wgt * dyDa[1];
    alpha[12] += wgt * dyDa[2];
    alpha[13] += wgt * dyDa[3];
    alpha[14] += wgt * dyDa[4];
    beta[4] += wgt * dy;
    wgt = dyDa[5] * w;
    alpha[15] += wgt * dyDa[0];
    alpha[16] += wgt * dyDa[1];
    alpha[17] += wgt * dyDa[2];
    alpha[18] += wgt * dyDa[3];
    alpha[19] += wgt * dyDa[4];
    alpha[20] += wgt * dyDa[5];
    beta[5] += wgt * dy;
  }

  @Override
  protected void initialiseGradient() {
    GradientProcedureHelper.initialiseWorkingMatrix6(alpha);
    beta[0] = 0;
    beta[1] = 0;
    beta[2] = 0;
    beta[3] = 0;
    beta[4] = 0;
    beta[5] = 0;
  }

  @Override
  public void getAlphaMatrix(double[][] alpha) {
    GradientProcedureHelper.getMatrix6(this.alpha, alpha);
  }

  @Override
  public void getAlphaLinear(double[] alpha) {
    GradientProcedureHelper.getMatrix6(this.alpha, alpha);
  }
}
