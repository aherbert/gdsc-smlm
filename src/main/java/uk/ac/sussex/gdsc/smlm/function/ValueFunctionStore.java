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

package uk.ac.sussex.gdsc.smlm.function;

/**
 * Wrap a function and store the values from the procedure.
 */
public class ValueFunctionStore implements ValueFunction, ValueProcedure {
  private final ValueFunction function;
  private ValueProcedure procedure;

  /** The counter to use as a result index during the procedure. */
  protected int index;
  /**
   * The values from the last call to {@link #forEach(ValueProcedure)}.
   */
  public double[] values;

  /**
   * Instantiates a new value function store.
   *
   * @param function the function
   */
  public ValueFunctionStore(ValueFunction function) {
    this(function, null);
  }

  /**
   * Instantiates a new value function store with storage.
   *
   * @param function the function
   * @param values the values
   */
  public ValueFunctionStore(ValueFunction function, double[] values) {
    this.function = function;
    this.values = values;
  }

  /** {@inheritDoc} */
  @Override
  public int size() {
    return function.size();
  }

  /** {@inheritDoc} */
  @Override
  public void initialise0(double[] a) {
    function.initialise0(a);
  }

  /** {@inheritDoc} */
  @Override
  public void forEach(ValueProcedure procedure) {
    index = 0;
    createValues();
    this.procedure = procedure;
    function.forEach(this);
  }

  /**
   * Creates the {@link #values} array.
   */
  protected void createValues() {
    if (values == null || values.length != function.size()) {
      values = new double[function.size()];
    }
  }

  /** {@inheritDoc} */
  @Override
  public void execute(double value) {
    values[index++] = value;
    procedure.execute(value);
  }
}
