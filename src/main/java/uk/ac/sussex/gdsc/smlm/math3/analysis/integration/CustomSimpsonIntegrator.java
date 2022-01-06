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

import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.NumberIsTooLargeException;
import org.apache.commons.math3.exception.NumberIsTooSmallException;

/**
 * Implements <a href="http://mathworld.wolfram.com/SimpsonsRule.html"> Simpson's Rule</a> for
 * integration of real univariate functions. For reference, see <b>Introduction to Numerical
 * Analysis</b>, ISBN 038795452X, chapter 3.
 *
 * <p>This implementation employs the basic trapezoid rule to calculate Simpson's rule.
 *
 * <p>Extends the default SimpsonIntegrator to allow the last computed sum to be returned even upon
 * failure.
 */
public class CustomSimpsonIntegrator extends SimpsonIntegrator {
  // CHECKSTYLE.OFF: MemberName
  private long n;
  // CHECKSTYLE.ON: MemberName
  private double lastSum;

  /** Maximal number of iterations for Simpson. */
  public static final int SIMPSON_MAX_ITERATIONS_COUNT = 63;

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

    // Simpson's rule requires at least two trapezoid stages.
    // So we set the first sum using two trapezoid stages.
    final TrapezoidIntegratorCopy qtrap = new TrapezoidIntegratorCopy();

    final double s0 = qtrap.stage(this, 0);
    double oldt = qtrap.stage(this, 1);
    double olds = (4 * oldt - s0) / 3.0;
    n = 2L;
    for (;;) {
      // The first iteration is now the first refinement of the sum.
      // This matches how the TrapezoidIntegrator works.
      incrementCount();
      final int i = getIterations();
      final double t = qtrap.stage(this, i + 1); // 1-stage ahead of the iteration
      final double s = (4 * t - oldt) / 3.0;
      n *= 2L;
      lastSum = s;
      if (i >= getMinimalIterationCount()) {
        final double delta = Math.abs(s - olds);
        final double rLimit = getRelativeAccuracy() * (Math.abs(olds) + Math.abs(s)) * 0.5;
        if ((delta <= rLimit) || (delta <= getAbsoluteAccuracy())) {
          return s;
        }
      }
      olds = s;
      oldt = t;
    }
  }

  // The following overrides are changed from protected to public to
  // allow access by the TrapezoidIntegratorCopy.

  @Override
  public double getMax() {
    return super.getMax();
  }

  @Override
  public double getMin() {
    return super.getMin();
  }

  @Override
  public double computeObjectiveValue(double point) {
    return super.computeObjectiveValue(point);
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
}
