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
 * <p>This procedure computes a modified Chi-squared expression to perform Weighted Least Squares
 * Estimation assuming a Poisson model with a Gaussian noise component. The weight per observation
 * is equal to 1/[variance + max(y, 0) + 1].
 *
 *
 * <p>See Ruisheng, et al (2017) Algorithmic corrections for localization microscopy with sCMOS
 * cameras - characterisation of a computationally efficient localization approach. Optical Express
 * 25, Issue 10, pp 11701-11716.
 */
public class WLsqLvmGradientProcedure extends LsqLvmGradientProcedure {
  /** The weights. */
  protected final double[] wgt;

  /**
   * Instantiates a new procedure.
   *
   * @param y Data to fit
   * @param var the base variance of each observation (must be positive)
   * @param func Gradient function
   */
  public WLsqLvmGradientProcedure(final double[] y, final double[] var,
      final Gradient1Function func) {
    super(y, func);
    final int n = y.length;
    wgt = new double[n];

    // From Ruisheng, et al (2017):
    // Total noise = variance + max(di, 0) + 1

    if (var != null && var.length == n) {
      // Include the variance in the weight. Assume variance is positive.
      for (int i = 0; i < n; i++) {
        wgt[i] = (y[i] > 0) ? 1.0 / (var[i] + y[i] + 1.0) : 1.0 / (var[i] + 1.0);
      }
    } else {
      for (int i = 0; i < n; i++) {
        wgt[i] = (y[i] > 0) ? 1.0 / (y[i] + 1.0) : 1.0;
      }
    }
  }

  @Override
  public void execute(double value, double[] dyDa) {
    final double dy = y[++yi] - value;
    final double weight = this.wgt[yi];
    this.value += dy * dy * weight;

    // Compute:
    // - the scaled Hessian matrix (the square matrix of second-order partial derivatives of a
    // function;
    // that is, it describes the local curvature of a function of many variables.)
    // - the scaled gradient vector of the function's partial first derivatives with respect to the
    // parameters

    for (int j = 0, i = 0; j < n; j++) {
      final double dw = dyDa[j] * weight;

      for (int k = 0; k <= j; k++) {
        alpha[i++] += dw * dyDa[k];
      }
      beta[j] += dw * dy;
    }
  }

  @Override
  public void execute(double value) {
    // Produce a sum-of-squares
    final double dy = y[++yi] - value;
    this.value += dy * dy * wgt[yi];
  }
}
