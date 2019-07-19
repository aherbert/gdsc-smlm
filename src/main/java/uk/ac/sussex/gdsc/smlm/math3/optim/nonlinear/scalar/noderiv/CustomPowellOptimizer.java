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
/*
 * Licensed to the Apache Software Foundation (ASF) under one or more contributor license
 * agreements. See the NOTICE file distributed with this work for additional information regarding
 * copyright ownership. The ASF licenses this file to You under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance with the License. You may obtain a
 * copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software distributed under the License
 * is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
 * or implied. See the License for the specific language governing permissions and limitations under
 * the License.
 */

package uk.ac.sussex.gdsc.smlm.math3.optim.nonlinear.scalar.noderiv;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.smlm.math3.optim.PositionChecker;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.MathUnsupportedOperationException;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.NumberIsTooSmallException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.univariate.BracketFinder;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.SimpleUnivariateValueChecker;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.util.FastMath;

import java.util.Arrays;

/**
 * Powell's algorithm.
 *
 * <p>The class is based on the
 * org.apache.commons.math3.optim.nonlinear.scalar.noderiv.PowellOptimizer but updated to support:
 * (a) convergence on the position; (b) convergence when using the original basis vectors; (c)
 * support bounds checking on the current point within the optimisation space.
 */
public class CustomPowellOptimizer extends MultivariateOptimizer {
  /**
   * Minimum relative tolerance.
   */
  private static final double MIN_RELATIVE_TOLERANCE = 2 * FastMath.ulp(1d);
  /**
   * Relative threshold.
   */
  private final double relativeThreshold;
  /**
   * Absolute threshold.
   */
  private final double absoluteThreshold;
  /**
   * Line search.
   */
  private final LineSearch line;
  /** Convergence tolerance on position. */
  private PositionChecker positionChecker;
  /** Only allow convergence when using initial basis vectors. */
  private final boolean basisConvergence;
  /** Allow custom basis search direction. */
  private double[] basis;

  /** Flags to indicate if bounds are present. */
  private boolean isLower;
  private boolean isUpper;
  private double[] lower;
  private double[] upper;

  /**
   * The initial step is used to construct the basis vectors for the search direction. By default
   * the identity matrix is used for the search. The magnitude of each diagonal position can be set
   * using this data object.
   */
  public static class BasisStep implements OptimizationData {
    /** Initial step in each direction. */
    private final double[] step;

    /**
     * Instantiates a new basis step.
     *
     * @param step Initial step for the bracket search.
     */
    public BasisStep(double[] step) {
      this.step = step;
    }

    /**
     * Gets the initial step.
     *
     * @return the initial step.
     */
    public double[] getStep() {
      return step;
    }
  }

  /**
   * This constructor allows to specify a user-defined convergence checker, in addition to the
   * parameters that control the default convergence checking procedure. <br> The internal line
   * search tolerances are set to the square-root of their corresponding value in the multivariate
   * optimizer.
   *
   * @param rel Relative threshold.
   * @param abs Absolute threshold.
   * @param checker Convergence checker.
   * @param basisConvergence Only allow convergence when using initial basis vectors
   * @throws NotStrictlyPositiveException if {@code abs <= 0}.
   * @throws NumberIsTooSmallException if {@code rel < 2 * Math.ulp(1d)}.
   */
  public CustomPowellOptimizer(double rel, double abs, ConvergenceChecker<PointValuePair> checker,
      boolean basisConvergence) {
    this(rel, abs, FastMath.sqrt(rel), FastMath.sqrt(abs), checker, basisConvergence);
  }

  /**
   * This constructor allows to specify a user-defined convergence checker, in addition to the
   * parameters that control the default convergence checking procedure and the line search
   * tolerances.
   *
   * @param rel Relative threshold for this optimizer.
   * @param abs Absolute threshold for this optimizer.
   * @param lineRel Relative threshold for the internal line search optimizer.
   * @param lineAbs Absolute threshold for the internal line search optimizer.
   * @param checker Convergence checker.
   * @param basisConvergence Only allow convergence when using initial basis vectors. If true then
   *        the vectors are reset if they have been modified and the search continues.
   * @throws NotStrictlyPositiveException if {@code abs <= 0}.
   * @throws NumberIsTooSmallException if {@code rel < 2 * Math.ulp(1d)}.
   */
  public CustomPowellOptimizer(double rel, double abs, double lineRel, double lineAbs,
      ConvergenceChecker<PointValuePair> checker, boolean basisConvergence) {
    super(checker);

    if (rel < MIN_RELATIVE_TOLERANCE) {
      throw new NumberIsTooSmallException(rel, MIN_RELATIVE_TOLERANCE, true);
    }
    if (abs <= 0) {
      throw new NotStrictlyPositiveException(abs);
    }
    relativeThreshold = rel;
    absoluteThreshold = abs;
    this.basisConvergence = basisConvergence;

    // Create the line search optimizer.
    line = new LineSearch(lineRel, lineAbs);
  }

  /**
   * The parameters control the default convergence checking procedure. <br> The internal line
   * search tolerances are set to the square-root of their corresponding value in the multivariate
   * optimizer.
   *
   * @param rel Relative threshold.
   * @param abs Absolute threshold.
   * @throws NotStrictlyPositiveException if {@code abs <= 0}.
   * @throws NumberIsTooSmallException if {@code rel < 2 * Math.ulp(1d)}.
   */
  public CustomPowellOptimizer(double rel, double abs) {
    this(rel, abs, null, false);
  }

  /**
   * Builds an instance with the default convergence checking procedure.
   *
   * @param rel Relative threshold.
   * @param abs Absolute threshold.
   * @param lineRel Relative threshold for the internal line search optimizer.
   * @param lineAbs Absolute threshold for the internal line search optimizer.
   * @throws NotStrictlyPositiveException if {@code abs <= 0}.
   * @throws NumberIsTooSmallException if {@code rel < 2 * Math.ulp(1d)}.
   */
  public CustomPowellOptimizer(double rel, double abs, double lineRel, double lineAbs) {
    this(rel, abs, lineRel, lineAbs, null, false);
  }

  // @CHECKSTYLE.OFF: LocalVariableName
  // @CHECKSTYLE.OFF: ParameterName

  @Override
  protected PointValuePair doOptimize() {
    final GoalType goal = getGoalType();
    final double[] guess = getStartPoint();
    final int n = guess.length;

    // Mark when we have modified the basis vectors
    boolean nonBasis = false;
    double[][] direc = createBasisVectors(n);

    final ConvergenceChecker<PointValuePair> checker = getConvergenceChecker();

    double[] x = guess;
    // Ensure the point is within bounds
    applyBounds(x);
    double functionValue = computeObjectiveValue(x);
    double[] x1 = x.clone();
    for (;;) {
      incrementIterationCount();

      final double fX = functionValue;
      double fX2 = 0;
      double delta = 0;
      int bigInd = 0;

      for (int i = 0; i < n; i++) {
        fX2 = functionValue;

        final UnivariatePointValuePair optimum = line.search(x, direc[i]);
        functionValue = optimum.getValue();
        x = newPoint(x, direc[i], optimum.getPoint());

        if ((fX2 - functionValue) > delta) {
          delta = fX2 - functionValue;
          bigInd = i;
        }
      }

      final PointValuePair previous = new PointValuePair(x1, fX, false);
      final PointValuePair current = new PointValuePair(x, functionValue, false);
      boolean stop = false;
      if (positionChecker != null) {
        // Check for convergence on the position
        stop = positionChecker.converged(getIterations(), previous, current);
      }
      if (!stop) {
        // Check if we have improved from an impossible position
        if (Double.isInfinite(fX) || Double.isNaN(fX)) {
          if (Double.isInfinite(functionValue) || Double.isNaN(functionValue)) {
            // Nowhere to go
            stop = true;
          }
        } else {
          stop = DoubleEquality.almostEqualRelativeOrAbsolute(fX, functionValue, relativeThreshold,
              absoluteThreshold);
        }
      }

      if (!stop && checker != null) {
        stop = checker.converged(getIterations(), previous, current);
      }

      boolean reset = false;
      if (stop) {
        // Only allow convergence using the basis vectors, i.e. we cannot move along any dimension
        if (basisConvergence && nonBasis) {
          // Reset to the basis vectors and continue
          reset = true;
        } else {
          final PointValuePair answer;
          if (goal == GoalType.MINIMIZE) {
            answer = (functionValue < fX) ? current : previous;
          } else {
            answer = (functionValue > fX) ? current : previous;
          }
          return answer;
        }
      }

      if (reset) {
        direc = createBasisVectors(n);
        nonBasis = false;
      }

      final double[] d = new double[n];
      final double[] x2 = new double[n];
      for (int i = 0; i < n; i++) {
        d[i] = x[i] - x1[i];
        x2[i] = x[i] + d[i];
      }
      applyBounds(x2);

      x1 = x.clone();
      fX2 = computeObjectiveValue(x2);

      // See if we can continue along the overall search direction to find a better value
      if (fX > fX2) {
        // Check if:
        // 1. The decrease along the average direction was not due to any single direction's
        // decrease
        // 2. There is a substantial second derivative along the average direction and we are close
        // to it minimum
        double t = 2 * (fX + fX2 - 2 * functionValue);
        double temp = fX - functionValue - delta;
        t *= temp * temp;
        temp = fX - fX2;
        t -= delta * temp * temp;

        if (t < 0.0) {
          final UnivariatePointValuePair optimum = line.search(x, d);
          functionValue = optimum.getValue();
          if (reset) {
            x = newPoint(x, d, optimum.getPoint());
            continue;
          }

          final double[][] result = newPointAndDirection(x, d, optimum.getPoint());
          x = result[0];

          final int lastInd = n - 1;
          direc[bigInd] = direc[lastInd];
          direc[lastInd] = result[1];
          nonBasis = true;
        }
      }
    }
  }

  private double[][] createBasisVectors(final int n) {
    final double[][] direc = new double[n][n];
    double[] step;
    if (basis != null && basis.length == n) {
      step = basis;
    } else {
      step = new double[n];
      Arrays.fill(step, 1);
    }
    for (int i = 0; i < n; i++) {
      direc[i][i] = step[i];
    }
    return direc;
  }

  /**
   * Compute a new point (in the original space) and a new direction vector, resulting from the line
   * search.
   *
   * @param p Point used in the line search.
   * @param d Direction used in the line search.
   * @param optimum Optimum found by the line search.
   * @return a 2-element array containing the new point (at index 0) and the new direction (at index
   *         1).
   */
  private double[][] newPointAndDirection(final double[] p, final double[] d,
      final double optimum) {
    final int n = p.length;
    final double[] nP = new double[n];
    final double[] nD = new double[n];
    for (int i = 0; i < n; i++) {
      nD[i] = d[i] * optimum;
      nP[i] = p[i] + nD[i];
    }
    applyBounds(nP);

    final double[][] result = new double[2][];
    result[0] = nP;
    result[1] = nD;

    return result;
  }

  /**
   * Compute a new point (in the original space) resulting from the line search.
   *
   * @param p Point used in the line search.
   * @param d Direction used in the line search.
   * @param optimum Optimum found by the line search.
   * @return array containing the new point.
   */
  private double[] newPoint(final double[] p, final double[] d, final double optimum) {
    final int n = p.length;
    final double[] nP = new double[n];
    for (int i = 0; i < n; i++) {
      nP[i] = p[i] + d[i] * optimum;
    }
    applyBounds(nP);
    return nP;
  }

  /**
   * Value that will pass the precondition check for {@link BrentOptimizer} but will not pass the
   * convergence check, so that the custom checker will always decide when to stop the line search.
   */
  private static final double REL_TOL_UNUSED;

  static {
    REL_TOL_UNUSED = 2 * FastMath.ulp(1d);
  }

  /**
   * Class for finding the minimum of the objective function along a given direction.
   */
  private class LineSearch extends BrentOptimizer {
    /**
     * Value that will pass the precondition check for {@link BrentOptimizer} but will not pass the
     * convergence check, so that the custom checker will always decide when to stop the line
     * search.
     */
    private static final double ABS_TOL_UNUSED = Double.MIN_VALUE;
    /**
     * Automatic bracketing.
     */
    private final BracketFinder bracket = new BracketFinder();

    /**
     * The "BrentOptimizer" default stopping criterion uses the tolerances to check the domain
     * (point) values, not the function values. We thus create a custom checker to use function
     * values.
     *
     * @param rel Relative threshold.
     * @param abs Absolute threshold.
     */
    LineSearch(double rel, double abs) {
      super(REL_TOL_UNUSED, ABS_TOL_UNUSED, new SimpleUnivariateValueChecker(rel, abs));
    }

    /**
     * Find the minimum of the function {@code f(p + alpha * d)}.
     *
     * @param p Starting point.
     * @param d Search direction.
     * @return the optimum.
     * @throws org.apache.commons.math3.exception.TooManyEvaluationsException if the number of
     *         evaluations is exceeded.
     */
    public UnivariatePointValuePair search(final double[] p, final double[] d) {
      final int n = p.length;
      final UnivariateFunction f = new UnivariateFunction() {
        final double[] x = new double[n];

        @Override
        public double value(double alpha) {
          for (int i = 0; i < n; i++) {
            x[i] = p[i] + alpha * d[i];
          }
          // Ensure the point is within bounds
          applyBounds(x);
          return CustomPowellOptimizer.this.computeObjectiveValue(x);
        }
      };

      final GoalType goal = CustomPowellOptimizer.this.getGoalType();
      bracket.search(f, goal, 0, 1);
      // Passing "MAX_VALUE" as a dummy value because it is the enclosing
      // class that counts the number of evaluations (and will eventually
      // generate the exception).
      return optimize(new MaxEval(Integer.MAX_VALUE), new UnivariateObjectiveFunction(f), goal,
          new SearchInterval(bracket.getLo(), bracket.getHi(), bracket.getMid()));
    }
  }

  /**
   * Scans the list of (required and optional) optimization data that characterize the problem.
   *
   * @param optData Optimization data. The following data will be looked for: <ul>
   *        <li>{@link PositionChecker}</li> <li>{@link BasisStep}</li> </ul>
   */
  @Override
  protected void parseOptimizationData(OptimizationData... optData) {
    // Allow base class to register its own data.
    super.parseOptimizationData(optData);

    // The existing values (as set by the previous call) are reused if
    // not provided in the argument list.
    for (final OptimizationData data : optData) {
      if (data instanceof PositionChecker) {
        positionChecker = (PositionChecker) data;
      } else if (data instanceof BasisStep) {
        basis = ((BasisStep) data).getStep();
      }
    }

    checkParameters();
  }

  /**
   * Check the parameters.
   *
   * @throws MathUnsupportedOperationException if the basis step passed to the
   *         {@link #optimize(OptimizationData[]) optimize} method is zero for any dimension
   */
  private void checkParameters() {
    lower = getLowerBound();
    upper = getUpperBound();
    isLower = checkArray(lower, Double.NEGATIVE_INFINITY);
    isUpper = checkArray(upper, Double.POSITIVE_INFINITY);
    // Check that the upper bound is above the lower bound
    if (isUpper && isLower) {
      for (int i = 0; i < lower.length; i++) {
        if (lower[i] > upper[i]) {
          throw new MathUnsupportedOperationException(LocalizedFormats.CONSTRAINT);
        }
      }
    }
    if (basis != null) {
      for (final double d : basis) {
        if (d == 0) {
          throw new MathUnsupportedOperationException(LocalizedFormats.CONSTRAINT);
        }
      }
    }
  }

  /**
   * Check if the array contains anything other than value.
   *
   * @param array the array
   * @param value the value
   * @return True if the array has another value
   */
  private static boolean checkArray(double[] array, double value) {
    if (array == null) {
      return false;
    }
    for (int i = 0; i < array.length; i++) {
      if (value != array[i]) {
        return true;
      }
    }
    return false;
  }

  /**
   * Check the point falls within the configured bounds truncating if necessary.
   *
   * @param point the point
   */
  private void applyBounds(double[] point) {
    if (isUpper) {
      for (int i = 0; i < point.length; i++) {
        if (point[i] > upper[i]) {
          point[i] = upper[i];
        }
      }
    }
    if (isLower) {
      for (int i = 0; i < point.length; i++) {
        if (point[i] < lower[i]) {
          point[i] = lower[i];
        }
      }
    }
  }
}
