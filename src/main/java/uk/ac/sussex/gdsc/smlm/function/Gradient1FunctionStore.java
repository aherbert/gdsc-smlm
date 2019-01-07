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
public class Gradient1FunctionStore extends ValueFunctionStore
    implements Gradient1Function, Gradient1Procedure {
  private final Gradient1Function function1;
  private Gradient1Procedure procedure;

  /** The number of gradients. */
  protected final int length;

  /**
   * The gradients from the last call to {@link #forEach(Gradient1Procedure)}.
   */
  public double[][] dyda;

  /**
   * Instantiates a new gradient 1 function store.
   *
   * @param function the function
   */
  public Gradient1FunctionStore(Gradient1Function function) {
    this(function, null, null);
  }

  /**
   * Instantiates a new gradient 1 function store with storage.
   *
   * @param function the function
   * @param values the values
   * @param dyda the dyda
   */
  public Gradient1FunctionStore(Gradient1Function function, double[] values, double[][] dyda) {
    super(function, values);
    this.function1 = function;
    this.dyda = dyda;
    length = function.getNumberOfGradients();
  }

  /** {@inheritDoc} */
  @Override
  public void initialise(double[] a) {
    function1.initialise(a);
  }

  /** {@inheritDoc} */
  @Override
  public void initialise1(double[] a) {
    function1.initialise(a);
  }

  /** {@inheritDoc} */
  @Override
  public int[] gradientIndices() {
    return function1.gradientIndices();
  }

  /** {@inheritDoc} */
  @Override
  public int getNumberOfGradients() {
    return function1.getNumberOfGradients();
  }

  /** {@inheritDoc} */
  @Override
  public void forEach(Gradient1Procedure procedure) {
    index = 0;
    createValues();
    createDYDA();
    this.procedure = procedure;
    function1.forEach((Gradient1Procedure) this);
  }

  /**
   * Creates the {@link #dyda} matrix.
   */
  protected void createDYDA() {
    if (dyda == null || dyda.length != function1.size()) {
      dyda = new double[values.length][length];
    }
  }

  /** {@inheritDoc} */
  @Override
  public void execute(double value, double[] dy_da) {
    values[index] = value;
    System.arraycopy(dy_da[index], 0, dyda[index], 0, length);
    index++;
    procedure.execute(value, dy_da);
  }
}
