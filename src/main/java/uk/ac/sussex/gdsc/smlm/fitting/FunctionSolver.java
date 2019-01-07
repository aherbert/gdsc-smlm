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

package uk.ac.sussex.gdsc.smlm.fitting;

/**
 * Defines methods to fit a function with coefficients (a).
 */
public interface FunctionSolver {
  /**
   * Gets the type of function solver.
   *
   * @return the type
   */
  FunctionSolverType getType();

  /**
   * Fit a parameterised function using a set of data points, {@code f(x|parameters) = data}.
   *
   * <p>Note: No input x data is specified. It can be assumed that the data points x[i]
   * corresponding to observed data[i] are consecutive integers from zero.
   *
   * @param data Set of data points to fit (input)
   * @param fx The evaluated function data points (output)
   * @param parameters The function parameters (input/output)
   * @param parameterVariances Variance of the fitted parameters (output)
   * @return The fit status
   */
  FitStatus fit(final double[] data, final double[] fx, final double[] parameters,
      final double[] parameterVariances);

  /**
   * Gets the number of fitted parameters.
   *
   * @return the number of fitted parameters
   */
  int getNumberOfFittedParameters();

  /**
   * Gets the number of fitted points.
   *
   * @return the number of fitted points
   */
  int getNumberOfFittedPoints();

  /**
   * Gets the number of iterations used to solve the function.
   *
   * @return The number of iterations used to solve the function.
   */
  int getIterations();

  /**
   * Gets the number of function evaluations used to solve the function.
   *
   * @return The number of function evaluations used to solve the function.
   */
  int getEvaluations();

  /**
   * Specifies if the function solver supports a bounded search, that is a search of parameter space
   * within the total allowed space of valid parameters. If true then the bounds can be set before a
   * call to the fit(...) method.
   *
   * @return True if the function solver supports a bounded search
   */
  boolean isBounded();

  /**
   * Specifies if the function solver supports constraints on the parameters. If true then the
   * constraints can be set before a call to the fit(...) method.
   *
   * <p>Note that constraints are to be used to specify the values that are absolutely not allowed.
   * They are not meant to be as restrictive as the bounds for a solver that supports a bounded
   * search. For example the constraints on a parameter may be 0 - Infinity but the bounds may be 5
   * - 15. A bounded solver can be used to search within the expected range for a parameter.
   *
   * @return True if the function solver supports a constrained search
   */
  boolean isConstrained();

  /**
   * Specifies if the function solver supports per observation weights.
   *
   * @return True if the function solver supports per observation weights
   */
  boolean isWeighted();

  /**
   * Checks if the function solver requires a strictly positive function.
   *
   * @return true, if the function solver requires a strictly positive function
   */
  boolean isStrictlyPositiveFunction();

  /**
   * Set the bounds for each of the parameters. If a subset of the parameters are fitted then the
   * bounds can be ignored for the fixed parameters.
   *
   * <p>The bounds can be used to set the expected range for a parameter.
   *
   * @param lower the lower bounds
   * @param upper the upper bounds
   */
  void setBounds(double[] lower, double[] upper);

  /**
   * Set the constraints for each of the parameters. If a subset of the parameters are fitted then
   * the bounds can be ignored for the fixed parameters.
   *
   * @param lower the lower bounds
   * @param upper the upper bounds
   */
  void setConstraints(double[] lower, double[] upper);

  /**
   * Sets the weights for each of the observations. The weights must match the length of the
   * observations passed to {@link #fit(double[], double[], double[], double[])}.
   *
   * <p>The weights are the variances of each observation.
   *
   * @param weights the new weights
   */
  void setWeights(double[] weights);

  /**
   * The optimised function value for the solution.
   *
   * @return the value
   */
  double getValue();

  /**
   * Evaluate a parameterised function using a set of data points, {@code f(x|parameters) = data}.
   *
   * <p>Note: No input x data is specified. It can be assumed that the data points x[i]
   * corresponding to observed data[i] are consecutive integers from zero.
   *
   * <p>The evaluated data points should be the same values as the result from the
   * {@link #fit(double[], double[], double[], double[])} method if this is called with the fitted
   * parameters.
   *
   * @param data Set of data points (input)
   * @param fx The evaluated function data points (output)
   * @param parameters The function parameters (input)
   * @return True if evaluation was performed
   */
  boolean evaluate(final double[] data, final double[] fx, final double[] parameters);

  /**
   * Compute the deviations for a parameterised function using a set of data points,
   * {@code f(x|parameters) = data}.
   *
   * <p>Note: No input x data is specified. It can be assumed that the data points x[i]
   * corresponding to observed data[i] are consecutive integers from zero.
   *
   * <p>The deviations should be the same values as the result from the
   * {@link #fit(double[], double[], double[], double[])} method if this is called with the fitted
   * parameters.
   *
   * @param data Set of data points (input)
   * @param parameters The function parameters (input/output)
   * @param parameterVariances Variance of the fitted parameters (output)
   * @return True if computation was performed
   */
  boolean computeDeviations(final double[] data, final double[] parameters,
      final double[] parameterVariances);

  /**
   * Gets the name of the parameter for the specified index.
   *
   * @param index the parameter index
   * @return the name
   */
  String getName(int index);
}
