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

import org.apache.commons.math3.analysis.integration.BaseAbstractUnivariateIntegrator;
import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.NumberIsTooLargeException;
import org.apache.commons.math3.exception.NumberIsTooSmallException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;

/**
 * Implements the <a href="http://mathworld.wolfram.com/TrapezoidalRule.html">
 * Trapezoid Rule</a> for integration of real univariate functions. For
 * reference, see <b>Introduction to Numerical Analysis</b>, ISBN 038795452X,
 * chapter 3.
 * <p>
 * The function should be integrable.</p>
 *
 * <p><strong>Note:</strong> This class has been copied from
 * {@code org.apache.commons.math3.analysis.integration} for use in the
 * {@link CustomSimpsonIntegrator} which requires access to package-private functions.
 * The source version was 3.6.1.
 */
class TrapezoidIntegratorCopy extends BaseAbstractUnivariateIntegrator {

    /** Maximum number of iterations for trapezoid. */
    public static final int TRAPEZOID_MAX_ITERATIONS_COUNT = 64;

    /** Intermediate result. */
    private double s;

    /**
     * Build a trapezoid integrator with given accuracies and iterations counts.
     * @param relativeAccuracy relative accuracy of the result
     * @param absoluteAccuracy absolute accuracy of the result
     * @param minimalIterationCount minimum number of iterations
     * @param maximalIterationCount maximum number of iterations
     * (must be less than or equal to {@link #TRAPEZOID_MAX_ITERATIONS_COUNT}
     * @throws NotStrictlyPositiveException if minimal number of iterations
     * is not strictly positive
     * @throws NumberIsTooSmallException if maximal number of iterations
     * is lesser than or equal to the minimal number of iterations
     * @throws NumberIsTooLargeException if maximal number of iterations
     * is greater than {@link #TRAPEZOID_MAX_ITERATIONS_COUNT}
     */
    public TrapezoidIntegratorCopy(final double relativeAccuracy,
                               final double absoluteAccuracy,
                               final int minimalIterationCount,
                               final int maximalIterationCount)
        throws NotStrictlyPositiveException, NumberIsTooSmallException, NumberIsTooLargeException {
        super(relativeAccuracy, absoluteAccuracy, minimalIterationCount, maximalIterationCount);
        if (maximalIterationCount > TRAPEZOID_MAX_ITERATIONS_COUNT) {
            throw new NumberIsTooLargeException(maximalIterationCount,
                                                TRAPEZOID_MAX_ITERATIONS_COUNT, false);
        }
    }

    /**
     * Build a trapezoid integrator with given iteration counts.
     * @param minimalIterationCount minimum number of iterations
     * @param maximalIterationCount maximum number of iterations
     * (must be less than or equal to {@link #TRAPEZOID_MAX_ITERATIONS_COUNT}
     * @throws NotStrictlyPositiveException if minimal number of iterations
     * is not strictly positive
     * @throws NumberIsTooSmallException if maximal number of iterations
     * is lesser than or equal to the minimal number of iterations
     * @throws NumberIsTooLargeException if maximal number of iterations
     * is greater than {@link #TRAPEZOID_MAX_ITERATIONS_COUNT}
     */
    public TrapezoidIntegratorCopy(final int minimalIterationCount,
                               final int maximalIterationCount)
        throws NotStrictlyPositiveException, NumberIsTooSmallException, NumberIsTooLargeException {
        super(minimalIterationCount, maximalIterationCount);
        if (maximalIterationCount > TRAPEZOID_MAX_ITERATIONS_COUNT) {
            throw new NumberIsTooLargeException(maximalIterationCount,
                                                TRAPEZOID_MAX_ITERATIONS_COUNT, false);
        }
    }

    /**
     * Construct a trapezoid integrator with default settings.
     * (max iteration count set to {@link #TRAPEZOID_MAX_ITERATIONS_COUNT})
     */
    public TrapezoidIntegratorCopy() {
        super(DEFAULT_MIN_ITERATIONS_COUNT, TRAPEZOID_MAX_ITERATIONS_COUNT);
    }

    /**
     * Compute the n-th stage integral of trapezoid rule. This function
     * should only be called by API {@code integrate()} in the package.
     * To save time it does not verify arguments - caller does.
     * <p>
     * The interval is divided equally into 2^n sections rather than an
     * arbitrary m sections because this configuration can best utilize the
     * already computed values.</p>
     *
     * @param integrator integrator holding integration parameters
     * @param n the stage of 1/2 refinement, n = 0 is no refinement
     * @return the value of n-th stage integral
     * @throws TooManyEvaluationsException if the maximal number of evaluations
     * is exceeded.
     */
    double stage(final CustomSimpsonIntegrator integrator, final int n)
        throws TooManyEvaluationsException {

        if (n == 0) {
            final double max = integrator.getMax();
            final double min = integrator.getMin();
            s = 0.5 * (max - min) *
                      (integrator.computeObjectiveValue(min) +
                       integrator.computeObjectiveValue(max));
            return s;
        }

        final long np = 1L << (n-1);           // number of new points in this stage
        double sum = 0;
        final double max = integrator.getMax();
        final double min = integrator.getMin();
        // spacing between adjacent new points
        final double spacing = (max - min) / np;
        double x = min + 0.5 * spacing;    // the first new point
        for (long i = 0; i < np; i++) {
            sum += integrator.computeObjectiveValue(x);
            x += spacing;
        }
        // add the new sum to previously calculated result
        s = 0.5 * (s + sum * spacing);
        return s;
    }

    /**
     * Compute the n-th stage integral of trapezoid rule. This function
     * should only be called by API {@code integrate()} in the package.
     * To save time it does not verify arguments - caller does.
     * <p>
     * The interval is divided equally into 2^n sections rather than an
     * arbitrary m sections because this configuration can best utilize the
     * already computed values.</p>
     *
     * @param integrator integrator holding integration parameters
     * @param n the stage of 1/2 refinement, n = 0 is no refinement
     * @return the value of n-th stage integral
     * @throws TooManyEvaluationsException if the maximal number of evaluations
     * is exceeded.
     */
    double stage(final TrapezoidIntegratorCopy integrator, final int n)
        throws TooManyEvaluationsException {

        if (n == 0) {
            final double max = integrator.getMax();
            final double min = integrator.getMin();
            s = 0.5 * (max - min) *
                      (integrator.computeObjectiveValue(min) +
                       integrator.computeObjectiveValue(max));
            return s;
        }

        final long np = 1L << (n-1);           // number of new points in this stage
        double sum = 0;
        final double max = integrator.getMax();
        final double min = integrator.getMin();
        // spacing between adjacent new points
        final double spacing = (max - min) / np;
        double x = min + 0.5 * spacing;    // the first new point
        for (long i = 0; i < np; i++) {
            sum += integrator.computeObjectiveValue(x);
            x += spacing;
        }
        // add the new sum to previously calculated result
        s = 0.5 * (s + sum * spacing);
        return s;
    }

    /** {@inheritDoc} */
    @Override
    protected double doIntegrate()
        throws MathIllegalArgumentException, TooManyEvaluationsException,
               MaxCountExceededException {

        double oldt = stage(this, 0);
        incrementCount();
        for (;;) {
            final int i = getIterations();
            final double t = stage(this, i);
            if (i >= getMinimalIterationCount()) {
                final double delta = Math.abs(t - oldt);
                final double rLimit =
                    getRelativeAccuracy() * (Math.abs(oldt) + Math.abs(t)) * 0.5;
                if ((delta <= rLimit) || (delta <= getAbsoluteAccuracy())) {
                    return t;
                }
            }
            oldt = t;
            incrementCount();
        }
    }
}
