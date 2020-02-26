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

package uk.ac.sussex.gdsc.smlm.fitting.nonlinear;

import java.util.Arrays;
import uk.ac.sussex.gdsc.core.utils.ValidationUtils;
import uk.ac.sussex.gdsc.smlm.fitting.FisherInformationMatrix;
import uk.ac.sussex.gdsc.smlm.fitting.FitStatus;
import uk.ac.sussex.gdsc.smlm.fitting.FunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.FunctionSolverType;
import uk.ac.sussex.gdsc.smlm.function.GradientFunction;
import uk.ac.sussex.gdsc.smlm.function.NamedFunction;

/**
 * Abstract class with utility methods for the FunctionSolver interface.
 */
public abstract class BaseFunctionSolver implements FunctionSolver {
  private FunctionSolverType type;

  /** The gradient function. */
  protected GradientFunction function;

  private int maxEvaluations = 20;

  /** The number of fitted points. */
  protected int numberOfFittedPoints;

  /** The iterations. */
  protected int iterations;

  /** The evaluations. */
  protected int evaluations;

  /** The value. */
  protected double value;

  // Cache the data to fit on success

  /** The data Y from the last successful fit. */
  protected double[] lastY;

  /** The parameters A from the last successful fit. */
  protected double[] lastA;

  /**
   * Default constructor.
   *
   * @param type the type
   * @param function the function
   * @throws NullPointerException if the function or type is null
   */
  public BaseFunctionSolver(FunctionSolverType type, GradientFunction function) {
    this.function = ValidationUtils.checkNotNull(function, "Function must not be null");
    this.type = ValidationUtils.checkNotNull(type, "Type must not be null");
  }

  /**
   * Sets the type.
   *
   * @param type the new type
   */
  protected void setType(FunctionSolverType type) {
    this.type = ValidationUtils.checkNotNull(type, "Type must not be null");
  }

  @Override
  public FunctionSolverType getType() {
    return type;
  }

  @Override
  public FitStatus fit(double[] y, double[] fx, double[] a, double[] parametersVariance) {
    // Reset the results
    numberOfFittedPoints = y.length;
    iterations = 0;
    evaluations = 0;
    value = 0;
    lastY = null;
    lastA = null;
    preProcess();
    final FitStatus status = computeFit(y, fx, a, parametersVariance);
    if (status == FitStatus.OK) {
      if (lastY == null) {
        lastY = y;
      }
      if (lastA == null) {
        lastA = a;
      }
      postProcess();
    }
    return status;
  }

  /**
   * Run before fit/evaluate.
   */
  protected void preProcess() {
    // To be optionally over-ridden
  }

  /**
   * Run if the fit/evaluate was successful.
   */
  protected void postProcess() {
    // To be optionally over-ridden
  }

  @Override
  public boolean evaluate(double[] y, double[] fx, double[] a) {
    // Reset the results
    numberOfFittedPoints = y.length;
    iterations = 0;
    evaluations = 0;
    value = 0;
    lastY = null;
    lastA = null;
    preProcess();
    final boolean status = computeValue(y, fx, a);
    if (status) {
      if (lastY == null) {
        lastY = y;
      }
      if (lastA == null) {
        lastA = a;
      }
      postProcess();
    }
    return status;
  }

  @Override
  public boolean computeDeviations(double[] y, double[] a, double[] parametersVariance) {
    // Use a dedicated solver optimised for inverting the matrix diagonal.
    final FisherInformationMatrix m = computeFisherInformationMatrix(y, a);

    setDeviations(parametersVariance, m);

    return true;
  }

  /**
   * Compute fit.
   *
   * @param y the y values to fit
   * @param fx the final fitted y values
   * @param a the parameters a
   * @param parametersVariance the deviations for the parameters a
   * @return the fit status
   */
  protected abstract FitStatus computeFit(double[] y, double[] fx, double[] a,
      double[] parametersVariance);

  /**
   * Evaluate the function.
   *
   * @param y the y values to fit
   * @param fx the final fitted y values
   * @param a the parameters a
   * @return true if evaluated
   */
  protected abstract boolean computeValue(double[] y, double[] fx, double[] a);

  /**
   * Compute the Fisher Information matrix. This can be used to set the deviations for each of the
   * fitted parameters.
   *
   * <p>Alternatively a sub-class can override
   * {@link #computeDeviations(double[], double[], double[])} directly and provide a dummy
   * implementation of this function as it will not be used, e.g. throw an exception.
   *
   * @param y the y values
   * @param a the parameters
   * @return the Fisher Information matrix
   */
  protected abstract FisherInformationMatrix computeFisherInformationMatrix(double[] y, double[] a);

  /**
   * Copy the parameter values into an initial solution at positions defined by the
   * {@link GradientFunction#gradientIndices()}.
   *
   * @param params the parameters
   * @return the initial solution
   */
  public double[] getInitialSolution(double[] params) {
    final int[] indices = function.gradientIndices();
    final double[] initialSolution = new double[indices.length];
    for (int i = 0; i < indices.length; i++) {
      initialSolution[i] = params[indices[i]];
    }
    return initialSolution;
  }

  /**
   * Copy the solution values into the parameters at positions defined by the
   * {@link GradientFunction#gradientIndices()}.
   *
   * @param params the parameters
   * @param solution the solution
   */
  public void setSolution(double[] params, double[] solution) {
    final int[] indices = function.gradientIndices();
    for (int i = 0; i < indices.length; i++) {
      params[indices[i]] = solution[i];
    }
  }

  /**
   * Copy the covariance matrix diagonal values into the deviations at positions defined by the
   * {@link GradientFunction#gradientIndices()}.
   *
   * @param deviations the deviations
   * @param covar the covariance matrix (assumed to be NxN with N = gradientIndices().length)
   */
  public void setDeviationsFromMatrix(double[] deviations, double[][] covar) {
    Arrays.fill(deviations, 0);
    final int[] indices = function.gradientIndices();
    for (int i = 0; i < indices.length; i++) {
      deviations[indices[i]] = checkVariance(covar[i][i]);
    }
  }

  /**
   * Copy the covariance matrix diagonal values into the deviations at positions defined by the
   * {@link GradientFunction#gradientIndices()}.
   *
   * @param deviations the deviations
   * @param covar the covariance matrix (assumed to be NxN with N = gradientIndices().length)
   */
  public void setDeviationsFromLinearMatrix(double[] deviations, double[] covar) {
    Arrays.fill(deviations, 0);
    final int[] indices = function.gradientIndices();
    final int n = indices.length;
    for (int i = 0, j = 0; i < n; i++, j += (n + 1)) {
      deviations[indices[i]] = checkVariance(covar[j]);
    }
  }

  /**
   * Copy the covariance values into the deviations at positions defined by the
   * {@link GradientFunction#gradientIndices()}.
   *
   * @param deviations the deviations
   * @param covar the covariance values (assumed to be length N with N = gradientIndices().length)
   */
  public void setDeviations(double[] deviations, double[] covar) {
    Arrays.fill(deviations, 0);
    final int[] indices = function.gradientIndices();
    final int n = indices.length;
    for (int i = 0; i < n; i++) {
      deviations[indices[i]] = checkVariance(covar[i]);
    }
  }

  /**
   * Invert the Fisher information matrix to determine the Cramér-Rao Lower Bound and copy the CRLB
   * values into the deviations at positions defined by the
   * {@link GradientFunction#gradientIndices()}.
   *
   * @param deviations the deviations
   * @param fim the Fisher information matrix
   */
  public void setDeviations(double[] deviations, FisherInformationMatrix fim) {
    // Use this method for robustness, i.e. it will not fail
    setDeviations(deviations, fim.crlb(true));
  }

  private static double checkVariance(double value) {
    return (value > 0) ? value : 0;
  }

  @Override
  public int getNumberOfFittedParameters() {
    return function.getNumberOfGradients();
  }

  @Override
  public int getNumberOfFittedPoints() {
    return numberOfFittedPoints;
  }

  @Override
  public int getIterations() {
    return iterations;
  }

  @Override
  public int getEvaluations() {
    return evaluations;
  }

  /**
   * Gets the max evaluations.
   *
   * @return the max evaluations
   */
  public int getMaxEvaluations() {
    return maxEvaluations;
  }

  /**
   * Sets the max evaluations.
   *
   * @param maxEvaluations the new max evaluations
   */
  public void setMaxEvaluations(int maxEvaluations) {
    this.maxEvaluations = maxEvaluations;
  }

  @Override
  public boolean isBounded() {
    return false;
  }

  @Override
  public boolean isConstrained() {
    return false;
  }

  @Override
  public boolean isWeighted() {
    return false;
  }

  @Override
  public boolean isStrictlyPositiveFunction() {
    // Provide a default implementation based on the type.
    // This can be overridden if the solver can handle negative function data.
    switch (type) {
      case LSE:
        // Assume least-square estimation can handle any function data
        return false;
      case MLE:
        // Assume maximum likelihood estimation requires a positive function value
        // for computing the likelihood
        return true;
      case WLSE:
        // Assume least-square estimation can handle any function data
        return false;
      default:
        // Leave this as disabled by default
        return false;
    }
  }

  @Override
  public void setBounds(double[] lower, double[] upper) {
    // To be over-ridden
  }

  @Override
  public void setConstraints(double[] lower, double[] upper) {
    // To be over-ridden
  }

  @Override
  public void setWeights(double[] weights) {
    // To be over-ridden
  }

  @Override
  public double getValue() {
    return value;
  }

  /**
   * Update the function.
   *
   * @param function the new function
   * @throws NullPointerException if the function is null
   */
  public void setGradientFunction(GradientFunction function) {
    if (function == null) {
      throw new NullPointerException("Function must not be null");
    }
    this.function = function;
  }

  /**
   * Gets the gradient function.
   *
   * @return The function.
   */
  public GradientFunction getGradientFunction() {
    return function;
  }

  @Override
  public String getName(int index) {
    if (function instanceof NamedFunction) {
      return ((NamedFunction) function).getParameterName(index);
    }
    return "Unknown";
  }

  /**
   * Ensure positive values. If values are negative a copy is made with those values set to zero.
   *
   * @param y the y
   * @return the positive y values
   */
  public static double[] ensurePositive(double[] y) {
    return ensurePositive(y.length, y);
  }

  /**
   * Ensure positive values. If values are negative a copy is made with those values set to zero.
   *
   * @param n the number of values to check
   * @param y the y
   * @return the positive y values
   */
  public static double[] ensurePositive(final int n, double[] y) {
    for (int i = 0; i < n; i++) {
      if (y[i] < 0) {
        // Not positive so create a clone

        final double[] y2 = new double[n];
        if (i != 0) {
          // Copy the values that were positive
          System.arraycopy(y, 0, y2, 0, i);
        }

        // Note that java initialises the array to zero so only copy the positives.
        // Skip current index i as it was not positive.
        while (++i < n) {
          if (y[i] > 0) {
            // Copy positive values
            y2[i] = y[i];
          }
        }
        return y2;
      }
    }
    return y;
  }
}
