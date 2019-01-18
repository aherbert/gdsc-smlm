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

import uk.ac.sussex.gdsc.smlm.function.Gradient2Function;

/**
 * Calculates the Newton-Raphson update vector for a Poisson process using the first and second
 * partial derivatives.
 *
 * <p>Ref: Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically
 * minimum uncertainty. Nature Methods 7, 373-375 (supplementary note), Eq. 12.
 */
public class FastMleGradient2Procedure4 extends FastMleGradient2Procedure {

  /**
   * Instantiates a new procedure.
   *
   * @param x Data to fit (must be positive, i.e. the value of a Poisson process)
   * @param func Gradient function (must produce a strictly positive value, i.e. the mean of a
   *        Poisson process)
   */
  public FastMleGradient2Procedure4(final double[] x, final Gradient2Function func) {
    super(x, func);
    if (n != 4) {
      throw new IllegalArgumentException("Function must compute 4 gradients");
    }
  }

  @Override
  protected void reset2() {
    d1[0] = 0;
    d1[1] = 0;
    d1[2] = 0;
    d1[3] = 0;
    d2[0] = 0;
    d2[1] = 0;
    d2[2] = 0;
    d2[3] = 0;
  }

  /**
   * Reset the first derivative vector.
   */
  @Override
  protected void reset1() {
    d1[0] = 0;
    d1[1] = 0;
    d1[2] = 0;
    d1[3] = 0;
  }

  @Override
  public void execute(double uk, double[] dukDt, double[] d2ukDt2) {
    u[k] = uk;
    final double xk = x[k++];
    final double xk_uk_minus1 = xk / uk - 1.0;
    final double xk_uk2 = xk / (uk * uk);
    d1[0] += dukDt[0] * xk_uk_minus1;
    d1[1] += dukDt[1] * xk_uk_minus1;
    d1[2] += dukDt[2] * xk_uk_minus1;
    d1[3] += dukDt[3] * xk_uk_minus1;
    d2[0] += d2ukDt2[0] * xk_uk_minus1 - dukDt[0] * dukDt[0] * xk_uk2;
    d2[1] += d2ukDt2[1] * xk_uk_minus1 - dukDt[1] * dukDt[1] * xk_uk2;
    d2[2] += d2ukDt2[2] * xk_uk_minus1 - dukDt[2] * dukDt[2] * xk_uk2;
    d2[3] += d2ukDt2[3] * xk_uk_minus1 - dukDt[3] * dukDt[3] * xk_uk2;
  }

  @Override
  public void execute(double uk, double[] dukDt) {
    u[k] = uk;
    final double xk = x[k++];
    final double xk_uk_minus1 = xk / uk - 1.0;
    d1[0] += dukDt[0] * xk_uk_minus1;
    d1[1] += dukDt[1] * xk_uk_minus1;
    d1[2] += dukDt[2] * xk_uk_minus1;
    d1[3] += dukDt[3] * xk_uk_minus1;
  }
}
