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
 * Wraps a value function to add a pre-computed offset to the value during the forEach procedure.
 */
public class OffsetGradient1Function extends OffsetValueFunction
    implements Gradient1Function, Gradient1Procedure, NonLinearFunction {
  /** The gradient1 function. */
  protected final Gradient1Function f1;

  /** The procedure. */
  protected Gradient1Procedure procedure1;

  /**
   * Class for evaluating a function and storing the values and gradients.
   */
  protected class FunctionStore implements Gradient1Procedure {
    private int index;

    /** The values. */
    public final double[] values;

    /** The dyda gradients. */
    public final double[][] dyda;

    /** The number of gradients. */
    public final int length;

    /**
     * Instantiates a new function store.
     *
     * @param values the values
     * @param dyda the dyda
     */
    public FunctionStore(double[] values, double[][] dyda) {
      length = f1.getNumberOfGradients();
      if (values == null) {
        values = new double[f1.size()];
        dyda = new double[values.length][length];
      }
      this.values = values;
      this.dyda = dyda;
    }

    /**
     * Gets the values.
     */
    public void getValues() {
      index = 0;
      f1.forEach(this);
    }

    /** {@inheritDoc} */
    @Override
    public void execute(double value, double[] dyda) {
      values[index] = value;
      System.arraycopy(dyda, 0, this.dyda[index], 0, length);
      index++;
    }
  }

  /** Used to store all the values and gradients for the NonLinearFunction interface. */
  protected FunctionStore store;

  /** All the values for the NonLinearFunction interface. */
  protected double[] allValues;

  /** All the gradients for the NonLinearFunction interface. */
  protected double[][] allDyda;

  /**
   * Instantiates a new offset gradient1 function.
   *
   * @param function the function
   * @param values the precomputed values
   * @throws IllegalArgumentException if the values length does not match the function size
   */
  protected OffsetGradient1Function(Gradient1Function function, double[] values) {
    super(function, values);
    f1 = function;
  }

  /**
   * Instantiates a new offset gradient1 function.
   *
   * @param pre the function
   * @param values the precomputed values
   * @throws IllegalArgumentException if the values length does not match the function size
   */
  protected OffsetGradient1Function(OffsetGradient1Function pre, double[] values) {
    super(pre, values);
    f1 = (Gradient1Function) vf;
  }

  /**
   * Gets the gradient 1 function.
   *
   * @return the gradient 1 function
   */
  public Gradient1Function getGradient1Function() {
    return f1;
  }

  @Override
  public void initialise(double[] a) {
    store = null;
    f1.initialise(a);
    index = 0;
  }

  @Override
  public void initialise1(double[] a) {
    f1.initialise1(a);
  }

  @Override
  public int[] gradientIndices() {
    return f1.gradientIndices();
  }

  @Override
  public int getNumberOfGradients() {
    return f1.getNumberOfGradients();
  }

  @Override
  public void forEach(Gradient1Procedure procedure) {
    this.procedure1 = procedure;
    index = 0;
    f1.forEach((Gradient1Procedure) this);
  }

  @Override
  public void execute(double value, double[] dyda) {
    procedure1.execute(value + values[index++], dyda);
  }

  /**
   * Wrap a function with pre-computed values.
   *
   * @param func the function
   * @param baseline Baseline pre-computed y-values
   * @return the wrapped function (or the original if pre-computed values are null or wrong length)
   */
  public static Gradient1Function wrapGradient1Function(final Gradient1Function func,
      final double[] baseline) {
    if (baseline != null && baseline.length == func.size()) {
      // Avoid multiple wrapping
      if (func instanceof OffsetGradient1Function) {
        return new OffsetGradient1Function((OffsetGradient1Function) func, baseline);
      }
      return new OffsetGradient1Function(func, baseline);
    }
    return func;
  }

  @Override
  public double eval(int x, double[] dyda) {
    createStore();
    System.arraycopy(allDyda[index], 0, dyda, 0, store.length);
    return allValues[x];
  }

  @Override
  public double eval(int x) {
    createStore();
    return store.values[x];
  }

  @Override
  public double eval(int x, double[] dyda, double[] weight) {
    weight[0] = 1;
    return eval(x, dyda);
  }

  @Override
  public double evalw(int x, double[] weight) {
    weight[0] = 1;
    return eval(x);
  }

  @Override
  public boolean canComputeWeights() {
    return false;
  }

  private void createStore() {
    if (store == null) {
      store = new FunctionStore(allValues, allDyda);
      store.getValues();
      // Re-use space
      allValues = store.values;
      allDyda = store.dyda;
    }
  }
}
