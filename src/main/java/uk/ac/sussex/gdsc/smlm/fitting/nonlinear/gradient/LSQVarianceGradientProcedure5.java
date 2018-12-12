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

import uk.ac.sussex.gdsc.smlm.fitting.linear.EJMLLinearSolver;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;

import java.util.Arrays;

/**
 * Compute the variance of the parameters of the function assuming a least squares fit of a Poisson
 * process.
 *
 * <p>Uses the Mortensen formula (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 25.
 */
public class LSQVarianceGradientProcedure5 extends LSQVarianceGradientProcedure {
  /**
   * Instantiates a new procedure.
   *
   * @param func Gradient function
   */
  public LSQVarianceGradientProcedure5(final Gradient1Function func) {
    super(func);
    if (n != 5) {
      throw new IllegalArgumentException("Function must compute 5 gradients");
    }
  }

  /**
   * Instantiates a new procedure.
   *
   * @param func Gradient function
   * @param solver The solver used to invert the Fisher information matrix to find the Cramér–Rao
   *        lower bound (CRLB).
   * @throws IllegalArgumentException if the solver is null
   */
  public LSQVarianceGradientProcedure5(final Gradient1Function func, EJMLLinearSolver solver) {
    super(func, solver);
    if (n != 5) {
      throw new IllegalArgumentException("Function must compute 5 gradients");
    }
  }

  @Override
  protected void initialise() {
    I[0] = 0;
    E[0] = 0;
    I[5] = 0;
    E[5] = 0;
    I[6] = 0;
    E[6] = 0;
    I[10] = 0;
    E[10] = 0;
    I[11] = 0;
    E[11] = 0;
    I[12] = 0;
    E[12] = 0;
    I[15] = 0;
    E[15] = 0;
    I[16] = 0;
    E[16] = 0;
    I[17] = 0;
    E[17] = 0;
    I[18] = 0;
    E[18] = 0;
    I[20] = 0;
    E[20] = 0;
    I[21] = 0;
    E[21] = 0;
    I[22] = 0;
    E[22] = 0;
    I[23] = 0;
    E[23] = 0;
    I[24] = 0;
    E[24] = 0;
    Arrays.fill(variance, 0);
  }

  @Override
  protected boolean finish() {
    if (I[0] != I[0] || I[5] != I[5] || I[6] != I[6] || I[10] != I[10] || I[11] != I[11]
        || I[12] != I[12] || I[15] != I[15] || I[16] != I[16] || I[17] != I[17] || I[18] != I[18]
        || I[20] != I[20] || I[21] != I[21] || I[22] != I[22] || I[23] != I[23] || I[24] != I[24]) {
      return true;
    }

    I[1] = I[5];
    E[1] = E[5];
    I[2] = I[10];
    E[2] = E[10];
    I[7] = I[11];
    E[7] = E[11];
    I[3] = I[15];
    E[3] = E[15];
    I[8] = I[16];
    E[8] = E[16];
    I[13] = I[17];
    E[13] = E[17];
    I[4] = I[20];
    E[4] = E[20];
    I[9] = I[21];
    E[9] = E[21];
    I[14] = I[22];
    E[14] = E[22];
    I[19] = I[23];
    E[19] = E[23];
    return false;
  }

  @Override
  protected void computeVariance() {
    variance[0] = I[0] * E[0] * I[0] + I[0] * E[1] * I[5] + I[0] * E[2] * I[10]
        + I[0] * E[3] * I[15] + I[0] * E[4] * I[20] + I[1] * E[5] * I[0] + I[1] * E[6] * I[5]
        + I[1] * E[7] * I[10] + I[1] * E[8] * I[15] + I[1] * E[9] * I[20] + I[2] * E[10] * I[0]
        + I[2] * E[11] * I[5] + I[2] * E[12] * I[10] + I[2] * E[13] * I[15] + I[2] * E[14] * I[20]
        + I[3] * E[15] * I[0] + I[3] * E[16] * I[5] + I[3] * E[17] * I[10] + I[3] * E[18] * I[15]
        + I[3] * E[19] * I[20] + I[4] * E[20] * I[0] + I[4] * E[21] * I[5] + I[4] * E[22] * I[10]
        + I[4] * E[23] * I[15] + I[4] * E[24] * I[20];
    variance[1] = I[5] * E[0] * I[1] + I[5] * E[1] * I[6] + I[5] * E[2] * I[11]
        + I[5] * E[3] * I[16] + I[5] * E[4] * I[21] + I[6] * E[5] * I[1] + I[6] * E[6] * I[6]
        + I[6] * E[7] * I[11] + I[6] * E[8] * I[16] + I[6] * E[9] * I[21] + I[7] * E[10] * I[1]
        + I[7] * E[11] * I[6] + I[7] * E[12] * I[11] + I[7] * E[13] * I[16] + I[7] * E[14] * I[21]
        + I[8] * E[15] * I[1] + I[8] * E[16] * I[6] + I[8] * E[17] * I[11] + I[8] * E[18] * I[16]
        + I[8] * E[19] * I[21] + I[9] * E[20] * I[1] + I[9] * E[21] * I[6] + I[9] * E[22] * I[11]
        + I[9] * E[23] * I[16] + I[9] * E[24] * I[21];
    variance[2] = I[10] * E[0] * I[2] + I[10] * E[1] * I[7] + I[10] * E[2] * I[12]
        + I[10] * E[3] * I[17] + I[10] * E[4] * I[22] + I[11] * E[5] * I[2] + I[11] * E[6] * I[7]
        + I[11] * E[7] * I[12] + I[11] * E[8] * I[17] + I[11] * E[9] * I[22] + I[12] * E[10] * I[2]
        + I[12] * E[11] * I[7] + I[12] * E[12] * I[12] + I[12] * E[13] * I[17]
        + I[12] * E[14] * I[22] + I[13] * E[15] * I[2] + I[13] * E[16] * I[7]
        + I[13] * E[17] * I[12] + I[13] * E[18] * I[17] + I[13] * E[19] * I[22]
        + I[14] * E[20] * I[2] + I[14] * E[21] * I[7] + I[14] * E[22] * I[12]
        + I[14] * E[23] * I[17] + I[14] * E[24] * I[22];
    variance[3] = I[15] * E[0] * I[3] + I[15] * E[1] * I[8] + I[15] * E[2] * I[13]
        + I[15] * E[3] * I[18] + I[15] * E[4] * I[23] + I[16] * E[5] * I[3] + I[16] * E[6] * I[8]
        + I[16] * E[7] * I[13] + I[16] * E[8] * I[18] + I[16] * E[9] * I[23] + I[17] * E[10] * I[3]
        + I[17] * E[11] * I[8] + I[17] * E[12] * I[13] + I[17] * E[13] * I[18]
        + I[17] * E[14] * I[23] + I[18] * E[15] * I[3] + I[18] * E[16] * I[8]
        + I[18] * E[17] * I[13] + I[18] * E[18] * I[18] + I[18] * E[19] * I[23]
        + I[19] * E[20] * I[3] + I[19] * E[21] * I[8] + I[19] * E[22] * I[13]
        + I[19] * E[23] * I[18] + I[19] * E[24] * I[23];
    variance[4] = I[20] * E[0] * I[4] + I[20] * E[1] * I[9] + I[20] * E[2] * I[14]
        + I[20] * E[3] * I[19] + I[20] * E[4] * I[24] + I[21] * E[5] * I[4] + I[21] * E[6] * I[9]
        + I[21] * E[7] * I[14] + I[21] * E[8] * I[19] + I[21] * E[9] * I[24] + I[22] * E[10] * I[4]
        + I[22] * E[11] * I[9] + I[22] * E[12] * I[14] + I[22] * E[13] * I[19]
        + I[22] * E[14] * I[24] + I[23] * E[15] * I[4] + I[23] * E[16] * I[9]
        + I[23] * E[17] * I[14] + I[23] * E[18] * I[19] + I[23] * E[19] * I[24]
        + I[24] * E[20] * I[4] + I[24] * E[21] * I[9] + I[24] * E[22] * I[14]
        + I[24] * E[23] * I[19] + I[24] * E[24] * I[24];
  }
}
