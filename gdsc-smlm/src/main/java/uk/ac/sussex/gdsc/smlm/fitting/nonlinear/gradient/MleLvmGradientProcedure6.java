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
public class MleLvmGradientProcedure6 extends MleLvmGradientProcedure {

  /**
   * Instantiates a new procedure.
   *
   * @param y Data to fit (must be positive)
   * @param func Gradient function
   */
  public MleLvmGradientProcedure6(final double[] y, final Gradient1Function func) {
    super(y, func);
    ValidationUtils.checkArgument(numberOfGradients == 6, "Function must compute 6 gradients");
  }

  @Override
  public void execute(double fi, double[] dfiDa) {
    ++yi;
    if (fi > 0.0) {
      final double xi = y[yi];

      // We assume y[i] is positive but must handle zero
      if (xi > 0.0) {
        value += (fi - xi - xi * Math.log(fi / xi));

        final double xi_fi2 = xi / fi / fi;
        final double e = 1 - (xi / fi);

        beta[0] -= e * dfiDa[0];
        beta[1] -= e * dfiDa[1];
        beta[2] -= e * dfiDa[2];
        beta[3] -= e * dfiDa[3];
        beta[4] -= e * dfiDa[4];
        beta[5] -= e * dfiDa[5];

        alpha[0] += dfiDa[0] * xi_fi2 * dfiDa[0];
        double wgt;
        wgt = dfiDa[1] * xi_fi2;
        alpha[1] += wgt * dfiDa[0];
        alpha[2] += wgt * dfiDa[1];
        wgt = dfiDa[2] * xi_fi2;
        alpha[3] += wgt * dfiDa[0];
        alpha[4] += wgt * dfiDa[1];
        alpha[5] += wgt * dfiDa[2];
        wgt = dfiDa[3] * xi_fi2;
        alpha[6] += wgt * dfiDa[0];
        alpha[7] += wgt * dfiDa[1];
        alpha[8] += wgt * dfiDa[2];
        alpha[9] += wgt * dfiDa[3];
        wgt = dfiDa[4] * xi_fi2;
        alpha[10] += wgt * dfiDa[0];
        alpha[11] += wgt * dfiDa[1];
        alpha[12] += wgt * dfiDa[2];
        alpha[13] += wgt * dfiDa[3];
        alpha[14] += wgt * dfiDa[4];
        wgt = dfiDa[5] * xi_fi2;
        alpha[15] += wgt * dfiDa[0];
        alpha[16] += wgt * dfiDa[1];
        alpha[17] += wgt * dfiDa[2];
        alpha[18] += wgt * dfiDa[3];
        alpha[19] += wgt * dfiDa[4];
        alpha[20] += wgt * dfiDa[5];
      } else {
        value += fi;
        beta[0] -= dfiDa[0];
        beta[1] -= dfiDa[1];
        beta[2] -= dfiDa[2];
        beta[3] -= dfiDa[3];
        beta[4] -= dfiDa[4];
        beta[5] -= dfiDa[5];
      }
    }
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
