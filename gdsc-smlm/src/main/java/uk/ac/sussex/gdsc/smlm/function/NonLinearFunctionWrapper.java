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

import org.apache.commons.lang3.tuple.Pair;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;

/**
 * Wrap the NonLinearFunction to remove the parameters that are fixed from the evaluation methods.
 */
public class NonLinearFunctionWrapper implements ExtendedNonLinearFunction {
  private final NonLinearFunction fun;
  private final double[] params;
  private final int dataSize;
  private final int[] gradientIndices;

  /**
   * Create a new instance using the full set of parameters for the function and the number of
   * points the function evaluates. The parameters that are not within the function gradient indices
   * array will be fixed.
   *
   * @param fun The function
   * @param params The parameters
   * @param dataSize The number of data points to evaluate
   */
  public NonLinearFunctionWrapper(NonLinearFunction fun, double[] params, int dataSize) {
    this.fun = fun;
    this.params = params.clone();
    this.dataSize = dataSize;
    // This wrapper will evaluate all the indices that are not fixed
    gradientIndices = SimpleArrayUtils.natural(fun.getNumberOfGradients());
  }

  /**
   * Set the predictor coefficients for the function that are not fixed (i.e. those corresponding to
   * the gradient indices in the wrapped function). The fixed coefficients are set in the
   * constructor.
   */
  @Override
  public void initialise(double[] variables) {
    final int[] gi = fun.gradientIndices();
    for (int i = 0; i < gi.length; i++) {
      params[gi[i]] = variables[i];
    }
    fun.initialise(params);
  }

  @Override
  public int[] gradientIndices() {
    return gradientIndices;
  }

  @Override
  public int getNumberOfGradients() {
    return gradientIndices.length;
  }

  @Override
  public double eval(int x, double[] dyda) {
    return fun.eval(x, dyda);
  }

  @Override
  public double eval(int x) {
    return fun.eval(x);
  }

  @Override
  public double evalw(int x, double[] dyda, double[] weight) {
    return fun.evalw(x, dyda, weight);
  }

  @Override
  public double evalw(int x, double[] weight) {
    return fun.eval(x, weight);
  }

  @Override
  public boolean canComputeWeights() {
    return fun.canComputeWeights();
  }

  @Override
  public double[] computeValues(double[] variables) {
    initialise(variables);
    final double[] values = new double[dataSize];
    for (int i = 0; i < values.length; i++) {
      // Assume linear X from 0..N
      values[i] = fun.eval(i);
    }
    return values;
  }

  @Override
  public double[][] computeJacobian(double[] variables) {
    initialise(variables);

    final double[][] jacobian = new double[dataSize][];

    for (int i = 0; i < dataSize; ++i) {
      // Assume linear X from 0..N
      final double[] dyda = new double[variables.length];
      fun.eval(i, dyda);
      jacobian[i] = dyda;
    }

    return jacobian;
  }

  @Override
  public boolean canComputeValuesAndJacobian() {
    return true;
  }

  @Override
  public Pair<double[], double[][]> computeValuesAndJacobian(double[] variables) {
    initialise(variables);

    final double[][] jacobian = new double[dataSize][];
    final double[] values = new double[dataSize];

    for (int i = 0; i < dataSize; ++i) {
      // Assume linear X from 0..N
      final double[] dyda = new double[variables.length];
      values[i] = fun.eval(i, dyda);
      jacobian[i] = dyda;
    }

    return Pair.of(values, jacobian);
  }
}
