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

/**
 * Create a gradient procedure.
 */
public class LSQVarianceGradientProcedureFactory {
  /**
   * Create a new gradient procedure.
   *
   * @param func Gradient function
   * @return the gradient procedure
   */
  public static LSQVarianceGradientProcedure create(final Gradient1Function func) {
    switch (func.getNumberOfGradients()) {
      case 5:
        return new LSQVarianceGradientProcedure5(func);
      case 4:
        return new LSQVarianceGradientProcedure4(func);
      case 6:
        return new LSQVarianceGradientProcedure6(func);
      default:
        return new LSQVarianceGradientProcedure(func);
    }
  }

  /**
   * Create a new gradient procedure.
   *
   * @param func Gradient function
   * @param solver The solver used to invert the Fisher information matrix to find the Cramér–Rao
   *        lower bound (CRLB).
   * @return the gradient procedure
   */
  public static LSQVarianceGradientProcedure create(final Gradient1Function func,
      EJMLLinearSolver solver) {
    switch (func.getNumberOfGradients()) {
      case 5:
        return new LSQVarianceGradientProcedure5(func, solver);
      case 4:
        return new LSQVarianceGradientProcedure4(func, solver);
      case 6:
        return new LSQVarianceGradientProcedure6(func, solver);
      default:
        return new LSQVarianceGradientProcedure(func, solver);
    }
  }
}
