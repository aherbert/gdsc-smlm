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
public class Gradient2FunctionStore extends Gradient1FunctionStore
    implements Gradient2Function, Gradient2Procedure {
  private final Gradient2Function function2;
  private Gradient2Procedure procedure;

  /**
   * The gradients from the last call to {@link #forEach(Gradient2Procedure)}.
   */
  public double[][] d2yda2;

  /**
   * Instantiates a new gradient 2 function store.
   *
   * @param function the function
   */
  public Gradient2FunctionStore(Gradient2Function function) {
    this(function, null, null, null);
  }

  /**
   * Instantiates a new gradient 2 function store with storage.
   *
   * @param function the function
   * @param values the values
   * @param dyda the dyda
   * @param d2yda2 the d2yda2
   */
  public Gradient2FunctionStore(Gradient2Function function, double[] values, double[][] dyda,
      double[][] d2yda2) {
    super(function, values, dyda);
    this.function2 = function;
    this.d2yda2 = d2yda2;
  }

  /** {@inheritDoc} */
  @Override
  public void initialise2(double[] a) {
    function2.initialise2(a);
  }

  /** {@inheritDoc} */
  @Override
  public void forEach(Gradient2Procedure procedure) {
    index = 0;
    createValues();
    createDYDA();
    createD2YDA2();
    this.procedure = procedure;
    function2.forEach((Gradient2Procedure) this);
  }

  /**
   * Creates the {@link #d2yda2} matrix.
   */
  protected void createD2YDA2() {
    if (d2yda2 == null || d2yda2.length != function2.size()) {
      d2yda2 = new double[values.length][length];
    }
  }

  /** {@inheritDoc} */
  @Override
  public void execute(double value, double[] dy_da, double[] d2y_da2) {
    values[index] = value;
    System.arraycopy(dy_da[index], 0, dyda[index], 0, length);
    System.arraycopy(d2y_da2[index], 0, d2yda2[index], 0, length);
    index++;
    procedure.execute(value, dy_da, d2y_da2);
  }
}
