/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.function;

/**
 * Wraps a value function to add a pre-computed offset to the value during the forEach procedure.
 */
public final class OffsetFunctionFactory {

  /** No public construction. */
  private OffsetFunctionFactory() {}

  /**
   * Wrap a function with a pre-computed offset.
   *
   * @param func the function
   * @param baseline Baseline pre-computed offset for the values
   * @return the wrapped function (or the original if pre-computed values are null or wrong length)
   */
  public static ValueFunction wrapFunction(final ValueFunction func, final double[] baseline) {
    if (baseline != null && baseline.length == func.size()) {
      // Wrap appropriately
      if (func instanceof ExtendedGradient2Function) {
        return OffsetExtendedGradient2Function
            .wrapExtendedGradient2Function((ExtendedGradient2Function) func, baseline);
      }
      if (func instanceof Gradient2Function) {
        return OffsetGradient2Function.wrapGradient2Function((Gradient2Function) func, baseline);
      }
      if (func instanceof Gradient1Function) {
        return OffsetGradient1Function.wrapGradient1Function((Gradient1Function) func, baseline);
      }
      return OffsetValueFunction.wrapValueFunction(func, baseline);
    }
    return func;
  }
}
