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

import uk.ac.sussex.gdsc.smlm.function.Gradient2Function;

/**
 * Create a Fast MLE gradient procedure.
 */
public class FastMLEGradient2ProcedureFactory {
  /**
   * Create a new gradient procedure.
   *
   * @param x Data to fit (must be positive, i.e. the value of a Poisson process)
   * @param func Gradient function
   * @return the gradient procedure
   */
  public static FastMLEGradient2Procedure create(final double[] x, final Gradient2Function func) {
    return new FastMLEGradient2Procedure(x, func);
    // Note:
    // JUnit speed tests show the unrolled version are slower, i.e. the JVM is able to
    // efficiently optimise the single for loops in the procedure. So just return the
    // default implementation.
    // return createUnrolled(x, func);
  }

  /**
   * Create a new gradient procedure that has the loops unrolled.
   *
   * @param x Data to fit (must be positive, i.e. the value of a Poisson process)
   * @param func Gradient function
   * @return the gradient procedure
   */
  static FastMLEGradient2Procedure createUnrolled(final double[] x, final Gradient2Function func) {
    switch (func.getNumberOfGradients()) {
      case 5:
        return new FastMLEGradient2Procedure5(x, func);
      case 4:
        return new FastMLEGradient2Procedure4(x, func);
      case 6:
        return new FastMLEGradient2Procedure6(x, func);
      default:
        return new FastMLEGradient2Procedure(x, func);
    }
  }
}
