//@formatter:off
/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package uk.ac.sussex.gdsc.smlm.math3.analysis.integration;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.exception.NumberIsTooLargeException;
import org.apache.commons.math3.exception.NumberIsTooSmallException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;

/**
 * Implements <a href="http://mathworld.wolfram.com/SimpsonsRule.html"> Simpson's Rule</a> for
 * integration of real univariate functions. For reference, see <b>Introduction to Numerical
 * Analysis</b>, ISBN 038795452X, chapter 3.
 *
 * <p>This implementation employs the basic trapezoid rule to calculate Simpson's rule.
 *
 * <p>Extends the default SimpsonIntegrator to allow the last computed sum to be returned even upon
 * failure. Exceptions can be suppressed and the reason for failure can be obtained.
 *
 * <p>Note: Integration allows a full iteration to complete before the max function evaluations are
 * checked. This allows the iterations to be limited using an approximate maximum number of function
 * evaluations, and the result of the last iteration can be returned. The number of evaluations to
 * complete iteration i is equal to 1 + 2<sup>i+1</sup>.
 *
 * <p>If the maximum evaluations is reached or exceeded at the <strong>end</strong> of a
 * non-convereged iteration then a {@link TooManyEvaluationsException} is raised. This can be
 * suppressed using a property.
 *
 * <p>If the maximum iterations is exceeded then a {@link MaxCountExceededException} is raised.
 * This can be suppressed using a property.
 *
 * <p>The status of integration when using suppression of exceptions can be checked using
 * {@link #getStatus()}. The last computed result is available using {@link #getLastSum()}.
 */
public class CustomSimpsonIntegrator extends SimpsonIntegrator {
  // CHECKSTYLE.OFF: MemberName
  private long n;
  // CHECKSTYLE.ON: MemberName
  private double lastSum;

  /** Function to integrate. */
  private UnivariateFunction function;
  /** Maximum allowed evaluations. */
  private int maxEval;
  /** Function evaluations. */
  private int eval;
  /** Itegration iterations. */
  private int iter;
  /** Flag to suppress max evaluations exception. */
  private boolean noThrowOnMax;
  /** Integration status. */
  private Status status;

  /** Maximal number of iterations for Simpson. This is set lower than
   * {@link SimpsonIntegrator#SIMPSON_MAX_ITERATIONS_COUNT} as the iterations are limited
   * by the maximum integer value for the max evaluations limit. */
  public static final int SIMPSON_MAX_ITERATIONS_COUNT = 30;

  /**
   * Status of the integrator.
   */
  public enum Status {
    /** Initialising. Set during initialisation of integration. */
    INIT,
    /** Integration achieved the convergence criteria. */
    OK,
    /** Integration reached or exceeded the maximum evaluations without convergence. */
    MAX_EVAL,
    /** Integration exceeded the maximum iterations without convergence. */
    MAX_ITER;
  }

  /**
   * Build a Simpson integrator with given accuracies and iterations counts.
   *
   * @param relativeAccuracy relative accuracy of the result
   * @param absoluteAccuracy absolute accuracy of the result
   * @param minimalIterationCount minimum number of iterations
   * @param maximalIterationCount maximum number of iterations (must be less than or equal to
   *        {@link #SIMPSON_MAX_ITERATIONS_COUNT})
   * @throws NotStrictlyPositiveException if minimal number of iterations is not strictly positive
   * @throws NumberIsTooSmallException if maximal number of iterations is lesser than or equal to
   *         the minimal number of iterations
   * @throws NumberIsTooLargeException if maximal number of iterations is greater than
   *         {@link #SIMPSON_MAX_ITERATIONS_COUNT}
   */
  public CustomSimpsonIntegrator(final double relativeAccuracy, final double absoluteAccuracy,
      final int minimalIterationCount, final int maximalIterationCount) {
    super(relativeAccuracy, absoluteAccuracy, minimalIterationCount, maximalIterationCount);
    if (maximalIterationCount > SIMPSON_MAX_ITERATIONS_COUNT) {
      throw new NumberIsTooLargeException(maximalIterationCount, SIMPSON_MAX_ITERATIONS_COUNT,
          false);
    }
  }

  /**
   * Build a Simpson integrator with given iteration counts.
   *
   * @param minimalIterationCount minimum number of iterations
   * @param maximalIterationCount maximum number of iterations (must be less than or equal to
   *        {@link #SIMPSON_MAX_ITERATIONS_COUNT})
   * @throws NotStrictlyPositiveException if minimal number of iterations is not strictly positive
   * @throws NumberIsTooSmallException if maximal number of iterations is lesser than or equal to
   *         the minimal number of iterations
   * @throws NumberIsTooLargeException if maximal number of iterations is greater than
   *         {@link #SIMPSON_MAX_ITERATIONS_COUNT}
   */
  public CustomSimpsonIntegrator(final int minimalIterationCount, final int maximalIterationCount) {
    super(minimalIterationCount, maximalIterationCount);
    if (maximalIterationCount > SIMPSON_MAX_ITERATIONS_COUNT) {
      throw new NumberIsTooLargeException(maximalIterationCount, SIMPSON_MAX_ITERATIONS_COUNT,
          false);
    }
  }

  /**
   * Construct an integrator with default settings. (max iteration count set to
   * {@link #SIMPSON_MAX_ITERATIONS_COUNT})
   */
  public CustomSimpsonIntegrator() {
    // This can now be 1 since at least 1 refinement is performed.
    super(1, SIMPSON_MAX_ITERATIONS_COUNT);
  }

  @Override
  protected double doIntegrate() {
    // This is a modification from the base SimpsonIntegrator.
    // That only computed a single iteration if getMinimalIterationCount() == 1.
    //
    // The code has been modified to remove the use of a TrapezoidIntegrator
    // to compute the trapezoid integral. This is done using a iterative
    // function call.

    // Simpson's rule requires at least two trapezoid stages.
    // So we set the first sum using two trapezoid stages.

    // First stage of trapezoid rule
    final double min = getMin();
    double spacing = getMax() - min;
    double s0 = 0.5 * spacing * (computeObjectiveValue(getMax()) + computeObjectiveValue(min));

    // Update trapezoid stage
    int np = 1;
    double oldt = updateTrapezoid(s0, min, spacing, np);
    double olds = (4 * oldt - s0) / 3.0;
    n = 2L;
    for (;;) {
      // The first iteration is now the first refinement of the sum.
      // This matches how the TrapezoidIntegrator works.
      if (++iter > getMaximalIterationCount()) {
        status = Status.MAX_ITER;
        // Optional suppression of exception
        if (noThrowOnMax) {
          return lastSum;
        }
        throw new MaxCountExceededException(getMaximalIterationCount());
      }

      // Update trapezoid stage
      spacing *= 0.5;
      np <<= 1;
      final double t = updateTrapezoid(oldt, min, spacing, np);
      final double s = (4 * t - oldt) / 3.0;
      n <<= 1;
      lastSum = s;
      if (iter >= getMinimalIterationCount()) {
        final double delta = Math.abs(s - olds);
        final double rLimit = getRelativeAccuracy() * (Math.abs(olds) + Math.abs(s)) * 0.5;
        if ((delta <= rLimit) || (delta <= getAbsoluteAccuracy())) {
          status = Status.OK;
          return s;
        }
      }
      // Only compare evaluations after a full iteration has not converged
      if (Integer.compareUnsigned(eval, maxEval) >= 0) {
        status = Status.MAX_EVAL;
        // Optional suppression of exception
        if (noThrowOnMax) {
          return lastSum;
        }
        throw new TooManyEvaluationsException(maxEval);
      }
      olds = s;
      oldt = t;
    }
  }

  /**
   * Compute the next stage integral of trapezoid rule. The first point is evaluated at the min plus
   * half of the spacing. The remaining points are evaluated at increments thereafter of the
   * spacing.
   *
   * <pre>
   * sum_i f(min + (i+0.5) * spacing)
   * </pre>
   *
   * <p>This method can be iteratively called with a doubling of the number points and a halving of
   * the spacing. This incrementally evaluates the function using new points between the min and max
   * of the range spaced at midpoints between the last set of evaluated points.
   *
   * @param s previous trapezoid sum
   * @param min the minimum of the range
   * @param spacing the spacing between points
   * @param np the number of points to evaluate
   * @return the new trapezoid integral
   */
  double updateTrapezoid(double s, double min, double spacing, final int np) {
    double sum = 0;
    for (long i = 0; i < np; i++) {
      sum += computeObjectiveValue(min + (i + 0.5) * spacing);
    }
    // Update the result.
    // The sum should be the sum of *all* the points in the [min, max] range
    // computed at half-spacing intervals multiplied by 0.5 of the half-spacing.
    // We skipped computing points that were previously computed.
    // The previous result is multiplied by 0.5 to effectively compute
    // the sum of points we missed rescaled to the smaller spacing.
    return 0.5 * (s + sum * spacing);
  }

  @Override
  protected void setup(int maxEval, UnivariateFunction f, double lower, double upper)
      throws NullArgumentException, MathIllegalArgumentException {
    // Call super simply for the argument validation
    super.setup(maxEval, f, lower, upper);
    // Capture the variables to allow a full iteration to complete before the
    // max evaluations is triggered
    this.function = f;
    this.maxEval = maxEval;
    eval = 0;
    iter = 0;
    status = Status.INIT;
  }

  @Override
  public double getMax() {
    // This is a legacy exposure of the property and is not required for integration function.
    // It remains for binary compatibility.
    return super.getMax();
  }

  @Override
  public double getMin() {
    // This is a legacy exposure of the property and is not required for integration function.
    // It remains for binary compatibility.
    return super.getMin();
  }

  @Override
  public double computeObjectiveValue(double point) {
    // This does not call super so does not increment the evaluations counter
    // and trigger TooManyEvaluationsException
    eval++;
    return function.value(point);
  }

  @Override
  public int getEvaluations() {
    return eval;
  }

  @Override
  public int getIterations() {
    return iter;
  }

  /**
   * Gets the last sum successfully computed.
   *
   * @return the last sum
   */
  public double getLastSum() {
    return lastSum;
  }

  /**
   * Gets the number of increments between the max and min of the range from the last successful
   * computation.
   *
   * @return the n
   */
  public long getN() {
    return n;
  }

  /**
   * Sets to {@code true} to suppress throwing an exception if the maximum function evaluations or
   * maximum iterations will be exceeded before convergence.
   *
   * @param value the value
   */
  public void setNoThrowOnMax(boolean value) {
    this.noThrowOnMax = value;
  }

  /**
   * Gets the status from the last call to
   * {@link #integrate(int, UnivariateFunction, double, double)}.
   *
   * @return the status
   */
  public Status getStatus() {
    return status;
  }
}
