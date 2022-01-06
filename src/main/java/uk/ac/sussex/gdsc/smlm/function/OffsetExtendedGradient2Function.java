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

package uk.ac.sussex.gdsc.smlm.function;

/**
 * Wraps a value function to add a pre-computed offset to the value during the forEach procedure.
 */
public class OffsetExtendedGradient2Function extends OffsetGradient2Function
    implements ExtendedGradient2Function, ExtendedGradient2Procedure {
  /** The extended gradient2 function. */
  protected final ExtendedGradient2Function ef2;

  /** The procedure. */
  protected ExtendedGradient2Procedure procedureEx2;

  /**
   * Instantiates a new offset extended gradient2 function.
   *
   * @param function the function
   * @param values the pre-computed offset values
   * @throws IllegalArgumentException if the values length does not match the function size
   */
  protected OffsetExtendedGradient2Function(ExtendedGradient2Function function, double[] values) {
    super(function, values);
    ef2 = function;
  }

  /**
   * Instantiates a new offset extended gradient2 function.
   *
   * @param pre the pre
   * @param values the pre-computed offset values
   */
  protected OffsetExtendedGradient2Function(OffsetExtendedGradient2Function pre, double[] values) {
    super(pre, values);
    ef2 = (ExtendedGradient2Function) vf;
  }

  /**
   * Gets the extended gradient 2 function.
   *
   * @return the extended gradient 2 function
   */
  public ExtendedGradient2Function getExtendedGradient2Function() {
    return ef2;
  }

  @Override
  public void initialiseExtended2(double[] a) {
    ef2.initialiseExtended2(a);
  }

  @Override
  public void forEach(ExtendedGradient2Procedure procedure) {
    this.procedureEx2 = procedure;
    index = 0;
    ef2.forEach((ExtendedGradient2Procedure) this);
  }

  @Override
  public void executeExtended(double value, double[] dyda, double[] d2ydadb) {
    procedureEx2.executeExtended(value + values[index++], dyda, d2ydadb);
  }

  /**
   * Wrap a function with pre-computed values.
   *
   * @param func the function
   * @param baseline Baseline pre-computed y-values
   * @return the wrapped function (or the original if pre-computed values are null or wrong length)
   */
  public static ExtendedGradient2Function
      wrapExtendedGradient2Function(final ExtendedGradient2Function func, final double[] baseline) {
    if (baseline != null && baseline.length == func.size()) {
      // Avoid multiple wrapping
      if (func instanceof OffsetExtendedGradient2Function) {
        return new OffsetExtendedGradient2Function((OffsetExtendedGradient2Function) func,
            baseline);
      }
      return new OffsetExtendedGradient2Function(func, baseline);
    }
    return func;
  }
}
