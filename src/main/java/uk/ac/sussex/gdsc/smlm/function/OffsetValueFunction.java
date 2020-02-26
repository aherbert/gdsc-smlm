/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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
public class OffsetValueFunction extends PrecomputedValueFunction
    implements ValueProcedure, NamedFunction {
  /** The value function. */
  protected final ValueFunction vf;
  /** The counter i. */
  protected int index;
  /** The procedure. */
  protected ValueProcedure procedure;

  /**
   * Instantiates a new offset value function.
   *
   * @param function the function
   * @param values the precomputed values
   * @throws IllegalArgumentException if the values length does not match the function size
   */
  protected OffsetValueFunction(ValueFunction function, double[] values) {
    super(values);
    if (function.size() != values.length) {
      throw new IllegalArgumentException("Length of precomputed values must match function size");
    }
    this.vf = function;
  }

  /**
   * Instantiates a new offset value function by combining the current precomputed values with more
   * precomputed values. This is used internally and so no checks are made on the size of values
   * arrays (which must match).
   *
   * @param pre the pre-computed function
   * @param values the second set of precomputed values
   */
  protected OffsetValueFunction(OffsetValueFunction pre, double[] values) {
    // Clone the values as they will be modified
    super(values.clone());
    this.vf = pre.vf;
    final int n = vf.size();
    final double[] values1 = pre.values;
    for (int i = 0; i < n; i++) {
      this.values[i] += values1[i];
    }
  }

  /**
   * Gets the value function.
   *
   * @return the value function
   */
  public ValueFunction getValueFunction() {
    return vf;
  }

  @Override
  public void initialise0(double[] a) {
    vf.initialise0(a);
  }

  @Override
  public void forEach(ValueProcedure procedure) {
    this.procedure = procedure;
    index = 0;
    vf.forEach(this);
  }

  @Override
  public void execute(double value) {
    procedure.execute(value + values[index++]);
  }

  /**
   * Wrap a function with pre-computed values.
   *
   * @param func the function
   * @param baseline Baseline pre-computed y-values
   * @return the wrapped function (or the original if pre-computed values are null or wrong length)
   */
  public static ValueFunction wrapValueFunction(final ValueFunction func, final double[] baseline) {
    if (baseline != null && baseline.length == func.size()) {
      // Avoid multiple wrapping
      if (func instanceof OffsetValueFunction) {
        return new OffsetValueFunction((OffsetValueFunction) func, baseline);
      }
      return new OffsetValueFunction(func, baseline);
    }
    return func;
  }

  @Override
  public String getParameterName(int index) {
    if (vf instanceof NamedFunction) {
      return ((NamedFunction) vf).getParameterName(index);
    }
    return "Unknown";
  }
}
