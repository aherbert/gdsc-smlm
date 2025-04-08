/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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
public class OffsetGradient2Function extends OffsetGradient1Function
    implements Gradient2Function, Gradient2Procedure {
  /** The gradient2 function. */
  protected final Gradient2Function f2;

  /** The procedure. */
  protected Gradient2Procedure procedure2;

  /**
   * Instantiates a new offset gradient2 function.
   *
   * @param function the function
   * @param values the precomputed values
   * @throws IllegalArgumentException if the values length does not match the function size
   */
  protected OffsetGradient2Function(Gradient2Function function, double[] values) {
    super(function, values);
    f2 = function;
  }

  /**
   * Instantiates a new offset gradient2 function.
   *
   * @param pre the function
   * @param values the precomputed values
   * @throws IllegalArgumentException if the values length does not match the function size
   */
  protected OffsetGradient2Function(OffsetGradient2Function pre, double[] values) {
    super(pre, values);
    f2 = (Gradient2Function) vf;
  }

  /**
   * Gets the gradient 2 function.
   *
   * @return the gradient 2 function
   */
  public Gradient2Function getGradient2Function() {
    return f2;
  }

  @Override
  public void initialise2(double[] a) {
    f2.initialise2(a);
  }

  @Override
  public void forEach(Gradient2Procedure procedure) {
    this.procedure2 = procedure;
    index = 0;
    f2.forEach((Gradient2Procedure) this);
  }

  @Override
  public void execute(double value, double[] dyda, double[] d2yda2) {
    procedure2.execute(value + values[index++], dyda, d2yda2);
  }

  /**
   * Wrap a function with pre-computed values.
   *
   * @param func the function
   * @param baseline Baseline pre-computed y-values
   * @return the wrapped function (or the original if pre-computed values are null or wrong length)
   */
  public static Gradient2Function wrapGradient2Function(final Gradient2Function func,
      final double[] baseline) {
    if (baseline != null && baseline.length == func.size()) {
      // Avoid multiple wrapping
      if (func instanceof OffsetGradient2Function) {
        return new OffsetGradient2Function((OffsetGradient2Function) func, baseline);
      }
      return new OffsetGradient2Function(func, baseline);
    }
    return func;
  }
}
