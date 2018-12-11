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
import uk.ac.sussex.gdsc.smlm.function.Gradient1Procedure;

import java.util.Arrays;

/**
 * Calculates the Fisher information matrix for a Poisson process.
 *
 * <p>This procedure is based on computation of a modified Chi-squared expression to perform
 * Weighted Least Squares Estimation assuming a Poisson model with a Gaussian noise component. The
 * weight per observation is equal to 1/[variance + max(y, 0) + 1].
 *
 * <p>See Ruisheng, et al (2017) Algorithmic corrections for localization microscopy with sCMOS
 * cameras - characterisation of a computationally efficient localization approach. Optical Express
 * 25, Issue 10, pp 11701-11716.
 */
public class WPoissonGradientProcedure implements Gradient1Procedure {
  /** The weights. */
  protected final double[] w;

  /** The function. */
  protected final Gradient1Function func;

  /**
   * The number of gradients.
   */
  public final int n;

  /**
   * Working space for the Fisher information matrix (size n * (n + 1) / 2)
   */
  protected double[] data;

  /** The y index counter. */
  protected int yi;

  /**
   * @param y Data to fit
   * @param var the base variance of each observation (must be positive)
   * @param func Gradient function
   */
  public WPoissonGradientProcedure(final double[] y, final double[] var,
      final Gradient1Function func) {
    this.func = func;
    this.n = func.getNumberOfGradients();

    final int n = y.length;
    w = new double[n];

    // From Ruisheng, et al (2017):
    // Total noise = variance + max(di, 0) + 1

    if (var != null && var.length == n) {
      // Include the variance in the weight. Assume variance is positive.
      for (int i = 0; i < n; i++) {
        w[i] = (y[i] > 0) ? 1.0 / (var[i] + y[i] + 1.0) : 1.0 / (var[i] + 1.0);
      }
    } else {
      for (int i = 0; i < n; i++) {
        w[i] = (y[i] > 0) ? 1.0 / (y[i] + 1.0) : 1.0;
      }
    }
  }

  /**
   * Compute the Fisher information matrix.
   *
   * <pre>
   * Iab = E [ ( d ln(L(x|p)) / da ) * ( d ln(L(x|p)) / db ) ]
   * p = parameters
   * x = observed values
   * L(x|p) = likelihood of X given p
   * E = expected value
   * </pre>
   *
   * Note that this is only a true Fisher information matrix if the function returns the expected
   * value for a Poisson process (see Smith, et al (2010)). In this case the equation reduces to:
   *
   * <pre>
   * Iab = sum(i) (dYi da) * (dYi db) / Yi
   * </pre>
   *
   * In this case Yi refers to the expected value at observation i. This expression was updated
   * (Ruisheng, et al (2017)) to use Yi as the observed value at observation i (Oi). To increase
   * stability for zero or small Oi a Baysian prior is added using max(0, Oi) + 1. To account for
   * Gaussian noise per observation using the variance (vari) the weights can be combined resulting
   * in:
   *
   * <pre>
   * Iab = sum(i) (dYi da) * (dYi db) / (max(0, Oi) + 1 + vari)
   * </pre>
   *
   * <p>See Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically
   * minimum uncertainty. Nature Methods 7, 373-375 (supplementary note), Eq. 9.
   *
   * <p>See Ruisheng, et al (2017) Algorithmic corrections for localization microscopy with sCMOS
   * cameras - characterisation of a computationally efficient localization approach. Optical
   * Express 25, Issue 10, pp 11701-11716.
   *
   * A call to {@link #isNaNGradients()} will indicate if the gradients were invalid.
   *
   * @param a Set of coefficients for the function (if null then the function must be
   *        pre-initialised)
   */
  public void computeFisherInformation(final double[] a) {
    yi = 0;
    if (data == null) {
      data = new double[n * (n + 1) / 2];
    } else {
      initialiseWorkingMatrix();
    }
    if (a != null) {
      func.initialise1(a);
    }
    func.forEach(this);
  }

  /**
   * Initialise for the computation of the Fisher information matrix. Zero the working matrix.
   */
  protected void initialiseWorkingMatrix() {
    Arrays.fill(data, 0.0);
  }

  /** {@inheritDoc} */
  @Override
  public void execute(double value, double[] dy_da) {
    // Note: Ignore the value
    final double w = this.w[yi++];
    for (int j = 0, i = 0; j < n; j++) {
      final double wgt = dy_da[j] * w;
      for (int k = 0; k <= j; k++) {
        data[i++] += wgt * dy_da[k];
      }
    }
  }

  /**
   * Check the gradients are NaN.
   *
   * @return true, if NaN
   */
  protected boolean checkGradients() {
    for (int i = 0, len = data.length; i < len; i++) {
      if (Double.isNaN(data[i])) {
        return true;
      }
    }
    return false;
  }

  /**
   * Get the scaled Fisher Information matrix (size n*n).
   *
   * @return the alpha
   */
  public double[][] getMatrix() {
    final double[][] a = new double[n][n];
    getMatrix(a);
    return a;
  }

  /**
   * Get the scaled Fisher Information matrix (size n*n) into the provided storage.
   *
   * @param matrix the matrix
   */
  public void getMatrix(double[][] matrix) {
    GradientProcedureHelper.getMatrix(data, matrix, n);
  }

  /**
   * Get the scaled Fisher Information matrix (size n*n).
   *
   * @return the alpha
   */
  public double[] getLinear() {
    final double[] a = new double[n * n];
    getLinear(a);
    return a;
  }

  /**
   * Get the scaled Fisher Information matrix (size n*n) into the provided storage.
   *
   * @param matrix the matrix
   */
  public void getLinear(double[] matrix) {
    GradientProcedureHelper.getMatrix(data, matrix, n);
  }

  /**
   * @return True if the last calculation produced gradients with NaN values.
   */
  public boolean isNaNGradients() {
    return checkGradients();
  }
}
