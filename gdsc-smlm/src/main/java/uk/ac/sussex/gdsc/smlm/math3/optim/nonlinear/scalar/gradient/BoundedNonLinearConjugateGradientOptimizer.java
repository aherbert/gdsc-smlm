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

package uk.ac.sussex.gdsc.smlm.math3.optim.nonlinear.scalar.gradient;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BrentSolver;
import org.apache.commons.math3.analysis.solvers.UnivariateSolver;
import org.apache.commons.math3.exception.MathIllegalStateException;
import org.apache.commons.math3.exception.MathInternalError;
import org.apache.commons.math3.exception.MathUnsupportedOperationException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.SimpleValueChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.GradientMultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunctionGradient;
import org.apache.commons.math3.optim.nonlinear.scalar.gradient.Preconditioner;
import org.apache.commons.math3.optim.univariate.BracketFinder;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.SimpleUnivariateValueChecker;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;

/**
 * Non-linear conjugate gradient optimizer.
 *
 * <p>This class supports both the Fletcher-Reeves and the Polak-Ribière update formulas for the
 * conjugate search directions. It also supports optional preconditioning.
 *
 * <p>The class is based on the Apache Commons Math 3 class
 * {@link org.apache.commons.math3.optim.nonlinear.scalar.gradient.NonLinearConjugateGradientOptimizer
 * NonLinearConjugateGradientOptimizer}. It has been updated to support bounds checking on the
 * current point within the optimisation space.
 */
public class BoundedNonLinearConjugateGradientOptimizer extends GradientMultivariateOptimizer {
  /** Update formula for the beta parameter. */
  private final Formula updateFormula;
  /** Preconditioner (may be null). */
  private final Preconditioner preconditioner;
  /** solver to use in the line search (may be null). */
  private final UnivariateSolver solver;
  /** Initial step used to bracket the optimum in line search. */
  private double initialStep = 1;

  /** Flags to indicate if bounds are present. */
  private boolean isLower;
  private boolean isUpper;
  private double[] lower;
  private double[] upper;

  private boolean useGradientLineSearch = true;
  private boolean noBracket;
  private double sign;

  /**
   * Constructor with default {@link BrentSolver line search solver} and
   * {@link IdentityPreconditioner preconditioner}.
   *
   * @param updateFormula formula to use for updating the &beta; parameter, must be one of
   *        {@link Formula#FLETCHER_REEVES} or {@link Formula#POLAK_RIBIERE}.
   * @param checker Convergence checker.
   */
  public BoundedNonLinearConjugateGradientOptimizer(final Formula updateFormula,
      ConvergenceChecker<PointValuePair> checker) {
    this(updateFormula, checker, new BrentSolver(), new IdentityPreconditioner());
  }

  /**
   * Available choices of update formulas for the updating the parameter that is used to compute the
   * successive conjugate search directions. For non-linear conjugate gradients, there are two
   * formulas: <ul> <li>Fletcher-Reeves formula</li> <li>Polak-Ribière formula</li> </ul>
   *
   * <p>On the one hand, the Fletcher-Reeves formula is guaranteed to converge if the start point is
   * close enough of the optimum whether the Polak-Ribière formula may not converge in rare cases.
   * On the other hand, the Polak-Ribière formula is often faster when it does converge.
   * Polak-Ribière is often used.
   */
  public enum Formula {
    /** Fletcher-Reeves formula. */
    FLETCHER_REEVES,
    /** Polak-Ribière formula. */
    POLAK_RIBIERE
  }

  /**
   * The initial step is a factor with respect to the search direction (which itself is roughly
   * related to the gradient of the function). <br> It is used to find an interval that brackets the
   * optimum in line search.
   */
  public static class BracketingStep implements OptimizationData {
    /** Initial step. */
    private final double initialStep;

    /**
     * Instantiates a new bracketing step.
     *
     * @param step Initial step for the bracket search.
     */
    public BracketingStep(double step) {
      initialStep = step;
    }

    /**
     * Gets the initial step.
     *
     * @return the initial step.
     */
    public double getBracketingStep() {
      return initialStep;
    }
  }

  /**
   * Constructor with default {@link IdentityPreconditioner preconditioner}.
   *
   * @param updateFormula formula to use for updating the &beta; parameter, must be one of
   *        {@link Formula#FLETCHER_REEVES} or {@link Formula#POLAK_RIBIERE}.
   * @param checker Convergence checker.
   * @param lineSearchSolver Solver to use during line search.
   */
  public BoundedNonLinearConjugateGradientOptimizer(final Formula updateFormula,
      ConvergenceChecker<PointValuePair> checker, final UnivariateSolver lineSearchSolver) {
    this(updateFormula, checker, lineSearchSolver, new IdentityPreconditioner());
  }

  /**
   * Instantiates a new bounded non linear conjugate gradient optimizer.
   *
   * @param updateFormula formula to use for updating the &beta; parameter, must be one of
   *        {@link Formula#FLETCHER_REEVES} or {@link Formula#POLAK_RIBIERE}.
   * @param checker Convergence checker.
   * @param lineSearchSolver Solver to use during line search.
   * @param preconditioner Preconditioner.
   */
  public BoundedNonLinearConjugateGradientOptimizer(final Formula updateFormula,
      ConvergenceChecker<PointValuePair> checker, final UnivariateSolver lineSearchSolver,
      final Preconditioner preconditioner) {
    super(checker);

    this.updateFormula = updateFormula;
    solver = lineSearchSolver;
    this.preconditioner = preconditioner;
    initialStep = 1;
  }

  /**
   * {@inheritDoc}
   *
   * @param optData Optimization data. This method will register the following data: <ul>
   *        <li>{@link MaxEval}</li> <li>{@link MaxIter}</li> <li>{@link InitialGuess}</li>
   *        <li>{@link SimpleBounds}</li> <li>{@link ObjectiveFunction}</li>
   *        <li>{@link GoalType}</li> <li>{@link ObjectiveFunctionGradient}</li>
   *        <li>{@link BracketingStep}</li> </ul>
   * @return {@inheritDoc}
   * @throws TooManyEvaluationsException if the maximal number of evaluations (of the objective
   *         function) is exceeded.
   */
  @Override
  public PointValuePair optimize(OptimizationData... optData) throws TooManyEvaluationsException {
    // Set up base class and perform computation.
    return super.optimize(optData);
  }

  @Override
  protected PointValuePair doOptimize() {
    final ConvergenceChecker<PointValuePair> checker = getConvergenceChecker();
    final double[] point = getStartPoint();
    final GoalType goal = getGoalType();
    final int n = point.length;

    sign = (goal == GoalType.MINIMIZE) ? -1 : 1;
    double[] unbounded = point.clone();
    applyBounds(point);
    double[] gradient = computeObjectiveGradient(point);
    checkGradients(gradient, unbounded);

    if (goal == GoalType.MINIMIZE) {
      for (int i = 0; i < n; i++) {
        gradient[i] = -gradient[i];
      }
    }

    // Initial search direction.
    double[] steepestDescent = preconditioner.precondition(point, gradient);
    double[] searchDirection = steepestDescent.clone();

    double delta = 0;
    for (int i = 0; i < n; ++i) {
      delta += gradient[i] * searchDirection[i];
    }

    // Used for non-gradient based line search
    LineSearch line = null;
    double rel = 1e-6;
    double abs = 1e-10;
    if (getConvergenceChecker() instanceof SimpleValueChecker) {
      rel = ((SimpleValueChecker) getConvergenceChecker()).getRelativeThreshold();
      abs = ((SimpleValueChecker) getConvergenceChecker()).getRelativeThreshold();
    }
    line = new LineSearch(Math.sqrt(rel), Math.sqrt(abs));

    PointValuePair current = null;
    int maxEval = getMaxEvaluations();
    for (;;) {
      incrementIterationCount();

      final double objective = computeObjectiveValue(point);
      final PointValuePair previous = current;
      current = new PointValuePair(point, objective);
      if (previous != null && checker.converged(getIterations(), previous, current)) {
        // We have found an optimum.
        return current;
      }

      double step;
      if (useGradientLineSearch) {
        // Classic code using the gradient function for the line search:

        // Find the optimal step in the search direction.
        final UnivariateFunction lsf = new LineSearchFunction(point, searchDirection);

        final double uB;
        try {
          uB = findUpperBound(lsf, 0, initialStep);

          // Use this method to do some more extensive searching with checking for invalid gradients
          // uB = findUpperBoundWithCheck(lsf, 0, initialStep);

          // Check if the bracket found a minimum. Otherwise just move to the new point.
          if (noBracket) {
            step = uB;
          } else {
            // XXX Last parameters is set to a value close to zero in order to
            // work around the divergence problem in the "testCircleFitting"
            // unit test (see MATH-439).
            // System.out.printf("Bracket %f - %f - %f\n", 0., 1e-15, uB);
            step = solver.solve(maxEval, lsf, 0, uB, 1e-15);
            maxEval -= solver.getEvaluations(); // Subtract used up evaluations.
          }
        } catch (final MathIllegalStateException ex) {
          // System.out.printf("Failed to bracket %s @ %s\n", Arrays.toString(point),
          // Arrays.toString(searchDirection));

          // Line search without gradient (as per Powell optimiser)
          final UnivariatePointValuePair optimum = line.search(point, searchDirection);
          step = optimum.getPoint();

          // throw ex;
        }
      } else {
        // Line search without gradient (as per Powell optimiser)
        final UnivariatePointValuePair optimum = line.search(point, searchDirection);
        step = optimum.getPoint();
      }

      // Validate new point.
      // System.out.printf("Step = %f x %s\n", step, Arrays.toString(searchDirection));
      for (int i = 0; i < point.length; ++i) {
        point[i] += step * searchDirection[i];
      }
      unbounded = point.clone();
      applyBounds(point);
      gradient = computeObjectiveGradient(point);
      checkGradients(gradient, unbounded);

      if (goal == GoalType.MINIMIZE) {
        for (int i = 0; i < n; ++i) {
          gradient[i] = -gradient[i];
        }
      }

      // Compute beta.
      final double deltaOld = delta;
      final double[] newSteepestDescent = preconditioner.precondition(point, gradient);
      delta = 0;
      for (int i = 0; i < n; ++i) {
        delta += gradient[i] * newSteepestDescent[i];
      }

      if (delta == 0) {
        return new PointValuePair(point, computeObjectiveValue(point));
      }

      final double beta;
      switch (updateFormula) {
        case FLETCHER_REEVES:
          beta = delta / deltaOld;
          break;
        case POLAK_RIBIERE:
          double deltaMid = 0;
          for (int i = 0; i < gradient.length; ++i) {
            deltaMid += gradient[i] * steepestDescent[i];
          }
          beta = (delta - deltaMid) / deltaOld;
          break;
        default:
          // Should never happen.
          throw new MathInternalError();
      }
      steepestDescent = newSteepestDescent;

      // Compute conjugate search direction.
      if (getIterations() % n == 0 || beta < 0) {
        // Break conjugation: reset search direction.
        searchDirection = steepestDescent.clone();
      } else {
        // Compute new conjugate search direction.
        for (int i = 0; i < n; ++i) {
          searchDirection[i] = steepestDescent[i] + beta * searchDirection[i];
        }
      }

      // The gradient has already been adjusted for the search direction
      checkGradients(searchDirection, unbounded, -sign);
    }
  }

  /**
   * Scans the list of (required and optional) optimization data that characterize the problem.
   *
   * @param optData Optimization data. The following data will be looked for: <ul>
   *        <li>{@link BracketingStep}</li> </ul>
   */
  @Override
  protected void parseOptimizationData(OptimizationData... optData) {
    // Allow base class to register its own data.
    super.parseOptimizationData(optData);

    // The existing values (as set by the previous call) are reused if
    // not provided in the argument list.
    for (final OptimizationData data : optData) {
      if (data instanceof BracketingStep) {
        initialStep = ((BracketingStep) data).getBracketingStep();
        // If more data must be parsed, this statement _must_ be
        // changed to "continue".
        break;
      }
    }

    checkParameters();
  }

  /**
   * Finds the upper bound b ensuring bracketing of a root between a and b.
   *
   * @param func function whose root must be bracketed.
   * @param a lower bound of the interval.
   * @param initialStep initial step to try.
   * @return b such that f(a) and f(b) have opposite signs.
   * @throws MathIllegalStateException if no bracket can be found.
   */
  private static double findUpperBound(final UnivariateFunction func, final double a,
      final double initialStep) {
    final double yA = func.value(a);
    double yB = yA;
    for (double step = initialStep; step < Double.MAX_VALUE; step *= Math.max(2, yA / yB)) {
      final double b = a + step;
      yB = func.value(b);
      if (yA * yB <= 0) {
        return b;
      }
      if (Double.isNaN(yB)) {
        throw new MathIllegalStateException(
            LocalizedFormats.UNABLE_TO_BRACKET_OPTIMUM_IN_LINE_SEARCH);
      }
    }
    throw new MathIllegalStateException(LocalizedFormats.UNABLE_TO_BRACKET_OPTIMUM_IN_LINE_SEARCH);
  }

  /**
   * Finds the upper bound b ensuring bracketing of a root between a and b.
   *
   * @param func function whose root must be bracketed.
   * @param lower lower bound of the interval.
   * @param initialStep initial step to try.
   * @return b such that f(a) and f(b) have opposite signs.
   * @throws MathIllegalStateException if no bracket can be found.
   */
  @SuppressWarnings("unused")
  private double findUpperBoundWithChecks(final UnivariateFunction func, final double lower,
      final double initialStep) {
    noBracket = false;
    final double yA = func.value(lower);

    // Check we have a gradient. This should be true unless something slipped by.
    if (Double.isNaN(yA)) {
      throw new MathIllegalStateException(
          LocalizedFormats.UNABLE_TO_BRACKET_OPTIMUM_IN_LINE_SEARCH);
    }

    double yB = yA;
    double lastB = Double.NaN;
    for (double step = initialStep; step < Double.MAX_VALUE; step *= Math.max(2, yA / yB)) {
      double upper = lower + step;
      yB = func.value(upper);
      if (yA * yB <= 0) {
        return upper;
      }

      if (Double.isNaN(yB)) {
        // We have moved along the vector to a point where we have no gradient.
        // Bracketing is impossible.
        noBracket = true;

        // Check we made at least one step to a place with a new gradient
        if (!Double.isNaN(lastB)) {
          // Return the point we reached as the minimum
          return lastB;
        }

        // We made no valid steps. Do an inner loop reducing the step size until we find a point
        // with a valid gradient
        for (step *= 0.1; step > Double.MIN_VALUE; step *= 0.1) {
          upper = lower + step;
          yB = func.value(upper);
          if (yA * yB <= 0) {
            return upper;
          }
          if (!Double.isNaN(yB)) {
            lastB = upper;
          }
        }
        if (!Double.isNaN(lastB)) {
          // Return the point we reached as the minimum
          return lastB;
        }
        throw new MathIllegalStateException(
            LocalizedFormats.UNABLE_TO_BRACKET_OPTIMUM_IN_LINE_SEARCH);
      }
      lastB = upper;
    }
    throw new MathIllegalStateException(LocalizedFormats.UNABLE_TO_BRACKET_OPTIMUM_IN_LINE_SEARCH);
  }

  /** Default identity preconditioner. */
  public static class IdentityPreconditioner implements Preconditioner {
    @Override
    public double[] precondition(double[] variables, double[] direction) {
      return direction.clone();
    }
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
    private static final double REL_TOL_UNUSED = 1e-15;
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
     * @param point Starting point.
     * @param direction Search direction.
     * @return the optimum.
     * @throws org.apache.commons.math3.exception.TooManyEvaluationsException if the number of
     *         evaluations is exceeded.
     */
    public UnivariatePointValuePair search(final double[] point, final double[] direction) {
      final int n = point.length;
      final UnivariateFunction f = alpha -> {
        final double[] x = new double[n];
        for (int i = 0; i < n; i++) {
          x[i] = point[i] + alpha * direction[i];
        }
        applyBounds(x);
        return BoundedNonLinearConjugateGradientOptimizer.this.computeObjectiveValue(x);
      };

      final GoalType goal = BoundedNonLinearConjugateGradientOptimizer.this.getGoalType();
      bracket.search(f, goal, 0, initialStep);
      // System.out.printf("Bracket %f (%f) - %f (%f) - %f (%f)\n", bracket.getHi(),
      // bracket.getFHi(),
      // bracket.getMid(), bracket.getFMid(), bracket.getLo(), bracket.getFLo());
      // Passing "MAX_VALUE" as a dummy value because it is the enclosing
      // class that counts the number of evaluations (and will eventually
      // generate the exception).
      return optimize(new MaxEval(Integer.MAX_VALUE), new UnivariateObjectiveFunction(f), goal,
          new SearchInterval(bracket.getLo(), bracket.getHi(), bracket.getMid()));
    }
  }

  /**
   * Internal class for line search.
   *
   * <p>The function represented by this class is the dot product of the objective function gradient
   * and the search direction. Its value is zero when the gradient is orthogonal to the search
   * direction, i.e. when the objective function value is a local extremum along the search
   * direction.
   */
  private class LineSearchFunction implements UnivariateFunction {
    /** Current point. */
    private final double[] currentPoint;
    /** Search direction. */
    private final double[] searchDirection;

    /**
     * Instantiates a new line search function.
     *
     * @param point Current point.
     * @param direction Search direction.
     */
    public LineSearchFunction(double[] point, double[] direction) {
      currentPoint = point.clone();
      searchDirection = direction.clone();
    }

    @Override
    public double value(double x) {
      // current point in the search direction
      final double[] shiftedPoint = currentPoint.clone();
      for (int i = 0; i < shiftedPoint.length; ++i) {
        shiftedPoint[i] += x * searchDirection[i];
      }

      final double[] unbounded = shiftedPoint.clone();

      // Ensure the point is within bounds
      applyBounds(shiftedPoint);

      // gradient of the objective function
      final double[] gradient = computeObjectiveGradient(shiftedPoint);

      // Check for bad gradients here, e.g. the gradients can be NaN or infinity.
      if (checkGradients(gradient, unbounded)) {
        return Double.NaN;
      }

      // dot product with the search direction
      double dotProduct = 0;
      for (int i = 0; i < gradient.length; ++i) {
        dotProduct += gradient[i] * searchDirection[i];
      }

      return dotProduct;
    }
  }

  /**
   * Checks if there are lower or upper bounds that are not -Infinity or +Infinity.
   *
   * @throws MathUnsupportedOperationException if invalid bounds were passed to the
   *         {@link #optimize(OptimizationData[]) optimize} method.
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
    for (final double v : array) {
      if (v != value) {
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

  /**
   * Check if the point falls outside configured bounds truncating the gradient to zero if it is
   * moving further outside the bounds.
   *
   * @param gradient the gradient
   * @param point the point
   * @return true if NaN gradients
   */
  private boolean checkGradients(double[] gradient, double[] point) {
    return checkGradients(gradient, point, sign);
  }

  /**
   * Check if the point falls outside configured bounds truncating the gradient to zero if it is
   * moving further outside the bounds (defined by the sign of the search direction).
   *
   * @param gradient the gradient
   * @param point the point
   * @param sign the sign
   * @return true if NaN gradients
   */
  private boolean checkGradients(double[] gradient, double[] point, final double sign) {
    if (isUpper) {
      for (int i = 0; i < point.length; i++) {
        if (point[i] >= upper[i] && Math.signum(gradient[i]) == sign) {
          gradient[i] = 0;
        }
      }
    }
    if (isLower) {
      for (int i = 0; i < point.length; i++) {
        if (point[i] <= lower[i] && Math.signum(gradient[i]) == -sign) {
          gradient[i] = 0;
        }
      }
    }
    boolean isNaN = false;
    for (int i = 0; i < point.length; i++) {
      if (Double.isNaN(gradient[i])) {
        isNaN = true;
        gradient[i] = 0;
      }
    }
    return isNaN;
  }

  /**
   * Checks if using gradient line search.
   *
   * @return true, if using gradient line search
   */
  public boolean isUseGradientLineSearch() {
    return useGradientLineSearch;
  }

  /**
   * Sets whether to use gradient line search.
   *
   * @param useGradientLineSearch the new use gradient line search
   */
  public void setUseGradientLineSearch(boolean useGradientLineSearch) {
    this.useGradientLineSearch = useGradientLineSearch;
  }
}
