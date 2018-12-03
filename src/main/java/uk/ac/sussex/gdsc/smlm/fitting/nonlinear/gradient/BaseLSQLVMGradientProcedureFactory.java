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

/**
 * Create a gradient procedure.
 */
abstract class BaseLSQLVMGradientProcedureFactory {
  /**
   * Create a LSQ LVM procedure. <p> Instance methods for testing.
   *
   * @param y the y
   * @param func the function
   * @return the LSQ LVM gradient procedure
   */
  BaseLSQLVMGradientProcedure createProcedure(final double[] y, final Gradient1Function func) {
    return createProcedure(y, null, func);
  }

  /**
   * Create a LSQ LVM procedure. <p> Instance methods for testing.
   *
   * @param y the y
   * @param b the pre-computed background (can be null)
   * @param func the function
   * @return the LSQ LVM gradient procedure
   */
  abstract BaseLSQLVMGradientProcedure createProcedure(final double[] y, final double[] b,
      final Gradient1Function func);
}
