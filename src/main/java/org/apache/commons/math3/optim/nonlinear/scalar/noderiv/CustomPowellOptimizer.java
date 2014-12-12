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
package org.apache.commons.math3.optim.nonlinear.scalar.noderiv;

import gdsc.smlm.utils.DoubleEquality;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.MathUnsupportedOperationException;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.NumberIsTooSmallException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.PositionChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.univariate.BracketFinder;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.SimpleUnivariateValueChecker;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.util.FastMath;

/**
 * Powell's algorithm.
 * <p>
 * The class is based on the org.apache.commons.math3.optim.nonlinear.scalar.noderiv.PowellOptimizer but updated to
 * support: (a) convergence on the position; (b) convergence when using the original basis vectors.
 */
public class CustomPowellOptimizer extends MultivariateOptimizer
{
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
	/** Convergence tolerance on position */
	private PositionChecker positionChecker = null;
	/** Only allow convergence when using initial basis vectors */
	private final boolean basisConvergence;

	/**
	 * This constructor allows to specify a user-defined convergence checker,
	 * in addition to the parameters that control the default convergence
	 * checking procedure. <br/>
	 * The internal line search tolerances are set to the square-root of their
	 * corresponding value in the multivariate optimizer.
	 *
	 * @param rel
	 *            Relative threshold.
	 * @param abs
	 *            Absolute threshold.
	 * @param checker
	 *            Convergence checker.
	 * @param basisConvergence
	 *            Only allow convergence when using initial basis vectors
	 * @throws NotStrictlyPositiveException
	 *             if {@code abs <= 0}.
	 * @throws NumberIsTooSmallException
	 *             if {@code rel < 2 * Math.ulp(1d)}.
	 */
	public CustomPowellOptimizer(double rel, double abs, ConvergenceChecker<PointValuePair> checker,
			boolean basisConvergence)
	{
		this(rel, abs, FastMath.sqrt(rel), FastMath.sqrt(abs), checker, basisConvergence);
	}

	/**
	 * This constructor allows to specify a user-defined convergence checker,
	 * in addition to the parameters that control the default convergence
	 * checking procedure and the line search tolerances.
	 *
	 * @param rel
	 *            Relative threshold for this optimizer.
	 * @param abs
	 *            Absolute threshold for this optimizer.
	 * @param lineRel
	 *            Relative threshold for the internal line search optimizer.
	 * @param lineAbs
	 *            Absolute threshold for the internal line search optimizer.
	 * @param checker
	 *            Convergence checker.
	 * @param basisConvergence
	 *            Only allow convergence when using initial basis vectors. If true then the vectors are reset if they
	 *            have been modified and the search continues.
	 * @throws NotStrictlyPositiveException
	 *             if {@code abs <= 0}.
	 * @throws NumberIsTooSmallException
	 *             if {@code rel < 2 * Math.ulp(1d)}.
	 */
	public CustomPowellOptimizer(double rel, double abs, double lineRel, double lineAbs,
			ConvergenceChecker<PointValuePair> checker, boolean basisConvergence)
	{
		super(checker);

		if (rel < MIN_RELATIVE_TOLERANCE)
		{
			throw new NumberIsTooSmallException(rel, MIN_RELATIVE_TOLERANCE, true);
		}
		if (abs <= 0)
		{
			throw new NotStrictlyPositiveException(abs);
		}
		relativeThreshold = rel;
		absoluteThreshold = abs;
		this.basisConvergence = basisConvergence;

		// Create the line search optimizer.
		line = new LineSearch(lineRel, lineAbs);
	}

	/**
	 * The parameters control the default convergence checking procedure. <br/>
	 * The internal line search tolerances are set to the square-root of their
	 * corresponding value in the multivariate optimizer.
	 *
	 * @param rel
	 *            Relative threshold.
	 * @param abs
	 *            Absolute threshold.
	 * @throws NotStrictlyPositiveException
	 *             if {@code abs <= 0}.
	 * @throws NumberIsTooSmallException
	 *             if {@code rel < 2 * Math.ulp(1d)}.
	 */
	public CustomPowellOptimizer(double rel, double abs)
	{
		this(rel, abs, null, false);
	}

	/**
	 * Builds an instance with the default convergence checking procedure.
	 *
	 * @param rel
	 *            Relative threshold.
	 * @param abs
	 *            Absolute threshold.
	 * @param lineRel
	 *            Relative threshold for the internal line search optimizer.
	 * @param lineAbs
	 *            Absolute threshold for the internal line search optimizer.
	 * @throws NotStrictlyPositiveException
	 *             if {@code abs <= 0}.
	 * @throws NumberIsTooSmallException
	 *             if {@code rel < 2 * Math.ulp(1d)}.
	 */
	public CustomPowellOptimizer(double rel, double abs, double lineRel, double lineAbs)
	{
		this(rel, abs, lineRel, lineAbs, null, false);
	}

	/** {@inheritDoc} */
	@Override
	protected PointValuePair doOptimize()
	{
		final GoalType goal = getGoalType();
		final double[] guess = getStartPoint();
		final int n = guess.length;

		// Mark when we have modified the basis vectors
		boolean nonBasis = false;
		double[][] direc = createBasisVectors(n);

		final ConvergenceChecker<PointValuePair> checker = getConvergenceChecker();

		//int resets = 0;

		double[] x = guess;
		double fVal = computeObjectiveValue(x);
		double[] x1 = x.clone();
		while (true)
		{
			incrementIterationCount();

			final double fX = fVal;
			double fX2 = 0;
			double delta = 0;
			int bigInd = 0;

			for (int i = 0; i < n; i++)
			{
				fX2 = fVal;

				final UnivariatePointValuePair optimum = line.search(x, direc[i]);
				fVal = optimum.getValue();
				x = newPoint(x, direc[i], optimum.getPoint());

				if ((fX2 - fVal) > delta)
				{
					delta = fX2 - fVal;
					bigInd = i;
				}
			}

			boolean stop;
			if (positionChecker != null)
			{
				// Check for convergence on the position
				stop = positionChecker.converged(x1, x);
			}
			else
			{
				// Default convergence check on value
				//stop = 2 * (fX - fVal) <= (relativeThreshold * (FastMath.abs(fX) + FastMath.abs(fVal)) + absoluteThreshold);
				stop = DoubleEquality.almostEqualRelativeOrAbsolute(fX, fVal, relativeThreshold, absoluteThreshold);
			}

			final PointValuePair previous = new PointValuePair(x1, fX);
			final PointValuePair current = new PointValuePair(x, fVal);
			if (!stop && checker != null)
			{ // User-defined stopping criteria.
				stop = checker.converged(getIterations(), previous, current);
			}

			boolean reset = false;
			if (stop)
			{
				// Only allow convergence using the basis vectors, i.e. we cannot move along any dimension
				if (basisConvergence && nonBasis)
				{
					// Reset to the basis vectors and continue
					reset = true;
					//resets++;
				}
				else
				{
					//System.out.printf("Resets = %d\n", resets);
					if (goal == GoalType.MINIMIZE)
					{
						return (fVal < fX) ? current : previous;
					}
					else
					{
						return (fVal > fX) ? current : previous;
					}
				}
			}

			if (reset)
			{
				direc = createBasisVectors(n);
				nonBasis = false;
			}

			final double[] d = new double[n];
			final double[] x2 = new double[n];
			for (int i = 0; i < n; i++)
			{
				d[i] = x[i] - x1[i];
				x2[i] = x[i] + d[i];
			}

			x1 = x.clone();
			fX2 = computeObjectiveValue(x2);

			// See if we can continue along the overall search direction to find a better value
			if (fX > fX2)
			{
				// Check if:
				// 1. The decrease along the average direction was not due to any single direction's decrease
				// 2. There is a substantial second derivative along the average direction and we are close to
				// it minimum
				double t = 2 * (fX + fX2 - 2 * fVal);
				double temp = fX - fVal - delta;
				t *= temp * temp;
				temp = fX - fX2;
				t -= delta * temp * temp;

				if (t < 0.0)
				{
					final UnivariatePointValuePair optimum = line.search(x, d);
					fVal = optimum.getValue();
					if (reset)
					{
						x = newPoint(x, d, optimum.getPoint());
						continue;
					}
					else
					{
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
	}

	private double[][] createBasisVectors(final int n)
	{
		double[][] direc = new double[n][n];
		for (int i = 0; i < n; i++)
		{
			direc[i][i] = 1;
		}
		return direc;
	}

	/**
	 * Compute a new point (in the original space) and a new direction
	 * vector, resulting from the line search.
	 *
	 * @param p
	 *            Point used in the line search.
	 * @param d
	 *            Direction used in the line search.
	 * @param optimum
	 *            Optimum found by the line search.
	 * @return a 2-element array containing the new point (at index 0) and
	 *         the new direction (at index 1).
	 */
	private double[][] newPointAndDirection(final double[] p, final double[] d, final double optimum)
	{
		final int n = p.length;
		final double[] nP = new double[n];
		final double[] nD = new double[n];
		for (int i = 0; i < n; i++)
		{
			nD[i] = d[i] * optimum;
			nP[i] = p[i] + nD[i];
		}

		final double[][] result = new double[2][];
		result[0] = nP;
		result[1] = nD;

		return result;
	}

	/**
	 * Compute a new point (in the original space) resulting from the line search.
	 *
	 * @param p
	 *            Point used in the line search.
	 * @param d
	 *            Direction used in the line search.
	 * @param optimum
	 *            Optimum found by the line search.
	 * @return array containing the new point.
	 */
	private double[] newPoint(final double[] p, final double[] d, final double optimum)
	{
		final int n = p.length;
		final double[] nP = new double[n];
		for (int i = 0; i < n; i++)
		{
			nP[i] = p[i] + d[i] * optimum;
		}
		return nP;
	}

	/**
	 * Value that will pass the precondition check for {@link BrentOptimizer} but will not pass the convergence
	 * check, so that the custom checker
	 * will always decide when to stop the line search.
	 */
	private static final double REL_TOL_UNUSED;
	static {
		REL_TOL_UNUSED = 2 * FastMath.ulp(1d);
	}
	
	/**
	 * Class for finding the minimum of the objective function along a given
	 * direction.
	 */
	private class LineSearch extends BrentOptimizer
	{
		/**
		 * Value that will pass the precondition check for {@link BrentOptimizer} but will not pass the convergence
		 * check, so that the custom checker
		 * will always decide when to stop the line search.
		 */
		private static final double ABS_TOL_UNUSED = Double.MIN_VALUE;
		/**
		 * Automatic bracketing.
		 */
		private final BracketFinder bracket = new BracketFinder();

		/**
		 * The "BrentOptimizer" default stopping criterion uses the tolerances
		 * to check the domain (point) values, not the function values.
		 * We thus create a custom checker to use function values.
		 *
		 * @param rel
		 *            Relative threshold.
		 * @param abs
		 *            Absolute threshold.
		 */
		LineSearch(double rel, double abs)
		{
			super(REL_TOL_UNUSED, ABS_TOL_UNUSED, new SimpleUnivariateValueChecker(rel, abs));
		}

		/**
		 * Find the minimum of the function {@code f(p + alpha * d)}.
		 *
		 * @param p
		 *            Starting point.
		 * @param d
		 *            Search direction.
		 * @return the optimum.
		 * @throws org.apache.commons.math3.exception.TooManyEvaluationsException
		 *             if the number of evaluations is exceeded.
		 */
		public UnivariatePointValuePair search(final double[] p, final double[] d)
		{
			final int n = p.length;
			final UnivariateFunction f = new UnivariateFunction()
			{
				final double[] x = new double[n];

				public double value(double alpha)
				{
					for (int i = 0; i < n; i++)
					{
						x[i] = p[i] + alpha * d[i];
					}
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
	 * Scans the list of (required and optional) optimization data that
	 * characterize the problem.
	 *
	 * @param optData
	 *            Optimization data.
	 *            The following data will be looked for:
	 *            <ul>
	 *            <li>{@link PositionChecker}</li>
	 *            </ul>
	 */
	@Override
	protected void parseOptimizationData(OptimizationData... optData)
	{
		// Allow base class to register its own data.
		super.parseOptimizationData(optData);

		// The existing values (as set by the previous call) are reused if
		// not provided in the argument list.
		for (OptimizationData data : optData)
		{
			if (data instanceof PositionChecker)
			{
				positionChecker = (PositionChecker) data;
				break;
			}
		}

		checkParameters();
	}

	/**
	 * @throws MathUnsupportedOperationException
	 *             if bounds were passed to the {@link #optimize(OptimizationData[]) optimize} method.
	 */
	private void checkParameters()
	{
		if (getLowerBound() != null || getUpperBound() != null)
		{
			throw new MathUnsupportedOperationException(LocalizedFormats.CONSTRAINT);
		}
	}
}
