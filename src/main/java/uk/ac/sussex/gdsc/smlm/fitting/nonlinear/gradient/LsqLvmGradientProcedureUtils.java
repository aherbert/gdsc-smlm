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

import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;

/**
 * Create a gradient procedure.
 */
public final class LsqLvmGradientProcedureUtils {
  /** No public constructor. */
  private LsqLvmGradientProcedureUtils() {}

  /**
   * Create a new gradient procedure.
   *
   * @param y Data to fit
   * @param func Gradient function
   * @return the gradient procedure
   */
  public static LsqLvmGradientProcedure create(final double[] y, final Gradient1Function func) {
    switch (func.getNumberOfGradients()) {
      case 5:
        return new LsqLvmGradientProcedure5(y, func);
      case 4:
        return new LsqLvmGradientProcedure4(y, func);
      case 6:
        return new LsqLvmGradientProcedure6(y, func);

      default:
        return new LsqLvmGradientProcedure(y, func);
    }
  }

  /**
   * Create a new gradient procedure.
   *
   * @param y Data to fit
   * @param baseline Baseline pre-computed y-values
   * @param func Gradient function
   * @return the gradient procedure
   */
  public static LsqLvmGradientProcedure create(final double[] y, final double[] baseline,
      final Gradient1Function func) {
    switch (func.getNumberOfGradients()) {
      case 5:
        return new LsqLvmGradientProcedure5(y, baseline, func);
      case 4:
        return new LsqLvmGradientProcedure4(y, baseline, func);
      case 6:
        return new LsqLvmGradientProcedure6(y, baseline, func);

      default:
        return new LsqLvmGradientProcedure(y, baseline, func);
    }
  }
}
