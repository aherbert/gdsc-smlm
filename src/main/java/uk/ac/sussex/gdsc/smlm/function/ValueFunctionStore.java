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

package uk.ac.sussex.gdsc.smlm.function;

/**
 * Wrap a function and store the values from the procedure.
 */
public class ValueFunctionStore implements ValueFunction, ValueProcedure {
  private final ValueFunction f;
  private ValueProcedure procedure;

  /** The counter i. */
  protected int i;
  /** The values from the last call to {@link #forEach(ValueProcedure)}. */
  public double[] values;

  /**
   * Instantiates a new value function store.
   *
   * @param f the f
   */
  public ValueFunctionStore(ValueFunction f) {
    this(f, null);
  }

  /**
   * Instantiates a new value function store with storage.
   *
   * @param f the f
   * @param values the values
   */
  public ValueFunctionStore(ValueFunction f, double[] values) {
    this.f = f;
    this.values = values;
  }

  /** {@inheritDoc} */
  @Override
  public int size() {
    return f.size();
  }

  /** {@inheritDoc} */
  @Override
  public void initialise0(double[] a) {
    f.initialise0(a);
  }

  /** {@inheritDoc} */
  @Override
  public void forEach(ValueProcedure procedure) {
    i = 0;
    createValues();
    this.procedure = procedure;
    f.forEach(this);
  }

  /**
   * Creates the {@link #values} array.
   */
  protected void createValues() {
    if (values == null || values.length != f.size()) {
      values = new double[f.size()];
    }
  }

  /** {@inheritDoc} */
  @Override
  public void execute(double value) {
    values[i++] = value;
    procedure.execute(value);
  }
}
