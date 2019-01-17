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
import uk.ac.sussex.gdsc.smlm.function.Gradient1Procedure;

import java.util.Arrays;

/**
 * Calculates the Fisher information matrix for a Poisson process.
 *
 * <p>Ref: Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically
 * minimum uncertainty. Nature Methods 7, 373-375 (supplementary note), Eq. 9.
 */
public class PoissonGradientProcedure implements Gradient1Procedure {
  /** The function. */
  protected final Gradient1Function func;

  /**
   * The number of gradients.
   */
  public final int n;

  /** Working space for the Fisher information matrix (size n * (n + 1) / 2). */
  protected double[] data;

  /**
   * Instantiates a new procedure.
   *
   * @param func Gradient function
   */
  public PoissonGradientProcedure(final Gradient1Function func) {
    this.func = func;
    this.n = func.getNumberOfGradients();
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
   * <p>Note that this is only a true Fisher information matrix if the function returns the expected
   * value for a Poisson process (see Smith, et al (2010)). In this case the equation reduces to:
   *
   * <pre>
   * Iab = sum(i) (dYi da) * (dYi db) / Yi
   * </pre>
   *
   * <p>This expression was extended (Huang et al, (2015)) to account for Gaussian noise per
   * observation using the variance (vari):
   *
   * <pre>
   * Iab = sum(i) (dYi da) * (dYi db) / (Yi + vari)
   * </pre>
   *
   * <p>Thus per-observation noise can be handled by wrapping the input function with a pre-computed
   * gradient function and pre-computed noise values.
   *
   * <p>See Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically
   * minimum uncertainty. Nature Methods 7, 373-375 (supplementary note), Eq. 9.
   *
   * <p>See: Huang et al, (2015). Video-rate nanoscopy using sCMOS camera–specific single-molecule
   * localization algorithms. Nature Methods 10, 653–658.
   *
   * <p>A call to {@link #isNaNGradients()} will indicate if the gradients were invalid.
   *
   * @param a Set of coefficients for the function (if null then the function must be
   *        pre-initialised)
   */
  public void computeFisherInformation(final double[] a) {
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

  @Override
  public void execute(double value, double[] dyDa) {
    if (value > 0.0) {
      final double function = 1.0 / value;
      for (int j = 0, i = 0; j < n; j++) {
        final double wgt = function * dyDa[j];
        for (int k = 0; k <= j; k++) {
          data[i++] += wgt * dyDa[k];
        }
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
   * Checks if is na N gradients.
   *
   * @return True if the last calculation produced gradients with NaN values.
   */
  public boolean isNaNGradients() {
    return checkGradients();
  }
}
