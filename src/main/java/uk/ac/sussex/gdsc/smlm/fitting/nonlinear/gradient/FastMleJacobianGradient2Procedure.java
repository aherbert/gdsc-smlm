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

import java.util.Arrays;
import uk.ac.sussex.gdsc.smlm.function.ExtendedGradient2Function;
import uk.ac.sussex.gdsc.smlm.function.ExtendedGradient2Procedure;

/**
 * Calculates the Newton-Raphson update vector for a Poisson process using the first and second
 * partial derivatives.
 *
 * <p>Computes the Jacobian matrix of the partial derivatives, dFi/dxj, for all n parameters. dFi is
 * the first partial derivative of the log likelihood function with respect to parameter i. dFi/dxj
 * is the first partial derivative of dFi with respect to parameter j.
 *
 * <p>Ref: Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically
 * minimum uncertainty. Nature Methods 7, 373-375 (supplementary note), Eq. 12.
 */
public class FastMleJacobianGradient2Procedure extends FastMleGradient2Procedure
    implements ExtendedGradient2Procedure {
  /** The function. */
  protected final ExtendedGradient2Function exfunc;

  /**
   * The Jacobian matrix of the partial derivatives, dFi/dxj, for all n parameters. dFi is the first
   * partial derivative of the log likelihood function with respect to parameter i. dFi/dxj is the
   * first partial derivative of dFi with respect to parameter j.
   */
  private final double[] jacobian;

  /**
   * Instantiates a new procedure.
   *
   * @param x Data to fit (must be positive, i.e. the value of a Poisson process)
   * @param func Gradient function (must produce a strictly positive value, i.e. the mean of a
   *        Poisson process)
   */
  public FastMleJacobianGradient2Procedure(final double[] x, final ExtendedGradient2Function func) {
    super(x, func);
    jacobian = new double[numberOfGradients * (numberOfGradients + 1) / 2];
    this.exfunc = func;
  }

  /**
   * Calculates the first and second derivative of the Poisson log likelihood with respect to each
   * parameter.
   *
   * @param a Set of coefficients for the function
   */
  public void computeJacobian(final double[] a) {
    k = 0;
    resetExtended2();
    exfunc.initialiseExtended2(a);
    exfunc.forEach((ExtendedGradient2Procedure) this);
    for (int i = 0, index = 0; i < numberOfGradients; i++, index += i + 1) {
      d2[i] = jacobian[index];
    }
  }

  /**
   * Reset the computation vectors to zero.
   */
  protected void resetExtended2() {
    Arrays.fill(d1, 0);
    Arrays.fill(jacobian, 0);
  }

  @Override
  public void executeExtended(double uk, double[] dukDt, double[] d2ukDtDs) {
    u[k] = uk;
    final double xk = x[k++];

    final double xk_uk_minus1 = xk / uk - 1.0;
    final double xk_uk2 = xk / (uk * uk);
    for (int i = 0, index = 0; i < numberOfGradients; i++) {
      d1[i] += dukDt[i] * xk_uk_minus1;

      for (int j = 0, k = i * numberOfGradients; j <= i; j++, index++, k++) {
        // This requires the partial second derivative with respect to i and j
        jacobian[index] += d2ukDtDs[k] * xk_uk_minus1 - dukDt[i] * dukDt[j] * xk_uk2;
      }
    }
  }

  @Override
  public boolean isNaNGradients() {
    for (int i = numberOfGradients; i-- > 0;) {
      if (Double.isNaN(d1[i])) {
        return true;
      }
    }
    for (int i = jacobian.length; i-- > 0;) {
      if (Double.isNaN(jacobian[i])) {
        return true;
      }
    }
    return false;
  }

  /**
   * Gets the Jacobian (size n*n).
   *
   * @return the Jacobian
   */
  public double[] getJacobianLinear() {
    final double[] jacobianMatrix = new double[numberOfGradients * numberOfGradients];
    GradientProcedureHelper.getMatrix(this.jacobian, jacobianMatrix, numberOfGradients);
    return jacobianMatrix;
  }

  /**
   * Gets the Jacobian (size n*n) into the provided storage.
   *
   * @param jacobian the Jacobian
   */
  public void getJacobianLinear(double[] jacobian) {
    GradientProcedureHelper.getMatrix(this.jacobian, jacobian, numberOfGradients);
  }
}
