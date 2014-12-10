/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2014 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This code is based on the ideas expressed in Numerical Recipes in C++, 
 * The Art of Scientific Computing, Second Edition, W.H. Press, 
 * S.A. Teukolsky, W.T. Vetterling, B.P. Flannery (Cambridge University Press, 
 * Cambridge, 2002).
 *---------------------------------------------------------------------------*/

package org.apache.commons.math3.optim.nonlinear.scalar.gradient;

import java.util.Locale;

import org.apache.commons.math3.exception.MathUnsupportedOperationException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.exception.util.Localizable;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.PositionChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.GradientMultivariateOptimizer;
import org.apache.commons.math3.util.FastMath;

/**
 * Implementation of the Broyden-Fletcher-Goldfarb-Shanno (BFGS) variant of the Davidson-Fletcher-Powell (DFP)
 * minimisation.
 * <p>
 * This is not part of the Apache Commons Math library but extends the same base classes to allow an easy swap with
 * existing code based on the Apache library.
 * <p>
 * Note that although rare, it may happen that the algorithm converges since the search direction no longer leads
 * downhill. In case of doubt restarting the algorithm should overcome this issue.
 * <p>
 * The implementation is based upon that presented in: Numerical Recipes in C++, The Art of Scientific Computing, Second
 * Edition, W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery (Cambridge University Press, Cambridge, 2002).
 * The algorithm has been updated to support a bounded search and convergence checking on position and gradient.
 */
public class BFGSOptimizer extends GradientMultivariateOptimizer
{
	/** Maximum step length used in line search. */
	private double[] maximumStepLength = null;

	/** Convergence tolerance on gradient */
	private double gradientTolerance;

	/** Maximum number of restarts on the convergence point */
	private int restarts = 0;

	/** Maximum number of restarts in the event of roundoff error */
	private int roundoffRestarts = 3;

	/** Convergence tolerance on position */
	private PositionChecker positionChecker = null;

	/** Flags to indicate if bounds are present */
	private boolean isLower, isUpper;
	private double[] lower, upper;

	private double sign = 0;

	/**
	 * Specify the maximum step length in each dimension
	 */
	public static class StepLength implements OptimizationData
	{
		private double[] step;

		/**
		 * Build an instance
		 * 
		 * @param step
		 *            The maximum step size in each dimension
		 */
		public StepLength(double[] step)
		{
			this.step = step;
		}

		public double[] getStep()
		{
			return step;
		}
	}

	/**
	 * Specify the tolerance on the gradient convergence with zero
	 */
	public static class GradientTolerance implements OptimizationData
	{
		private double tolerance;

		/**
		 * Build an instance
		 * 
		 * @param tolerance
		 *            The tolerance on the gradient
		 */
		public GradientTolerance(double tolerance)
		{
			this.tolerance = tolerance;
		}

		public double getTolerance()
		{
			return tolerance;
		}
	}

	/**
	 * Specify the maximum number of restarts on the converged point in the event that the gradient has not yet
	 * converged on zero.
	 */
	public static class MaximumRestarts implements OptimizationData
	{
		private int restarts;

		/**
		 * Build an instance
		 * 
		 * @param restarts
		 *            The restarts on the gradient
		 */
		public MaximumRestarts(int restarts)
		{
			this.restarts = restarts;
		}

		public int getRestarts()
		{
			return restarts;
		}
	}


	/**
	 * Specify the maximum number of restarts in the event of roundoff error.
	 */
	public static class MaximumRoundoffRestarts extends MaximumRestarts
	{
		/**
		 * Build an instance
		 * 
		 * @param restarts
		 *            The restarts on the gradient
		 */
		public MaximumRoundoffRestarts(int restarts)
		{
			super(restarts);
		}
	}
	
	/**
	 * Constructor
	 */
	public BFGSOptimizer()
	{
		super(null);
	}

	/**
	 * @param checker
	 *            Convergence checker.
	 */
	public BFGSOptimizer(ConvergenceChecker<PointValuePair> checker)
	{
		super(checker);
	}

	/**
	 * {@inheritDoc}
	 *
	 * @param optData
	 *            Optimization data. In addition to those documented in
	 *            {@link GradientMultivariateOptimizer#parseOptimizationData(OptimizationData[])
	 *            GradientMultivariateOptimizer}, this method will register the following data:
	 *            <ul>
	 *            <li>{@link MaximumStepLength}</li>
	 *            </ul>
	 * @return {@inheritDoc}
	 * @throws TooManyEvaluationsException
	 *             if the maximal number of
	 *             evaluations (of the objective function) is exceeded.
	 */
	@Override
	public PointValuePair optimize(OptimizationData... optData) throws TooManyEvaluationsException
	{
		// Set up base class and perform computation.
		return super.optimize(optData);
	}

	private int converged;
	private static final int CHECKER = 0;
	private static final int POSITION = 1;
	private static final int GRADIENT = 2;
	private static final int ROUNDOFF_ERROR = 3;

	/** {@inheritDoc} */
	@Override
	protected PointValuePair doOptimize()
	{
		final ConvergenceChecker<PointValuePair> checker = getConvergenceChecker();
		double[] p = getStartPoint();

		// Assume minimisation
		sign = -1;

		LineStepSearch lineSearch = new LineStepSearch();

		// In case there are no restarts
		if (restarts <= 0)
			return bfgsWithRoundoffCheck(checker, p, lineSearch);

		PointValuePair lastResult = null;
		PointValuePair result = null;
		//int lastConverge = 0;
		int iteration = 0;
		//int initialConvergenceIteration = 0;
		//int[] count = new int[3];
		while (iteration <= restarts)
		{
			iteration++;
			result = bfgsWithRoundoffCheck(checker, p, lineSearch);
			//count[converged]++;

			//if (lastResult == null)
			//	initialConvergenceIteration = getIterations();

			if (converged == GRADIENT)
			{
				// If no gradient remains then we cannot move anywhere so return
				break;
			}

			if (lastResult != null)
			{
				//// Check if the optimum was improved using the last convergence criteria
				//if (lastConverge == CHECKER)
				//{
				//	if (checker.converged(getIterations(), lastResult, result))
				//	{
				//		break;
				//	}
				//}
				//else
				//{
				//	if (positionChecker.converged(lastResult.getPointRef(), result.getPointRef()))
				//	{
				//		break;
				//	}
				//}

				// Check if the optimum was improved using the convergence criteria
				if (checker != null && checker.converged(getIterations(), lastResult, result))
				{
					break;
				}
				if (positionChecker.converged(lastResult.getPointRef(), result.getPointRef()))
				{
					break;
				}
			}

			// Store the new optimum and repeat
			lastResult = result;
			//lastConverge = converged;
			p = lastResult.getPointRef();
		}

		//System.out.printf("Iter=%d (%d > %d): %s\n", iteration, initialConvergenceIteration, getIterations(),
		//		java.util.Arrays.toString(count));

		return result;
	}

	/**
	 * Repeat the BFGS algorithm until it converges without roundoff error on the search direction
	 * 
	 * @param checker
	 * @param p
	 * @param lineSearch
	 * @return
	 */
	protected PointValuePair bfgsWithRoundoffCheck(ConvergenceChecker<PointValuePair> checker, double[] p,
			LineStepSearch lineSearch)
	{
		// Note: Position might converge if the hessian becomes singular or non-positive-definite
		// In this case the simple check is to restart the algorithm.
		int iteration = 0;

		PointValuePair result = bfgs(checker, p, lineSearch);
		
		// Allow restarts in the case of roundoff convergence
		while (converged == ROUNDOFF_ERROR && iteration < roundoffRestarts)
		{
			iteration++;
			p = result.getPointRef();
			result = bfgs(checker, p, lineSearch);
		}
		
		// If restarts did not work then this is a failure
		if (converged == ROUNDOFF_ERROR)
			throw new LineSearchRoundoffException();
		
		//if (iteration > 0)
		//	System.out.printf("Restarts for roundoff error = %d\n", iteration);
		
		return result;
	}

	protected PointValuePair bfgs(ConvergenceChecker<PointValuePair> checker, double[] p, LineStepSearch lineSearch)
	{
		final int n = p.length;

		final double EPS = epsilon;

		double[] hdg = new double[n];
		double[] xi = new double[n];
		double[][] hessian = new double[n][n];

		// Get the gradient for the the bounded point
		applyBounds(p);
		double[] g = computeObjectiveGradient(p);
		checkGradients(g, p);

		// Initialise the hessian and search direction
		for (int i = 0; i < n; i++)
		{
			hessian[i][i] = 1.0;
			xi[i] = -g[i];
		}

		PointValuePair current = null;

		while (true)
		{
			incrementIterationCount();

			// Get the value of the point
			double fp = computeObjectiveValue(p);

			if (checker != null)
			{
				PointValuePair previous = current;
				current = new PointValuePair(p, fp);
				if (previous != null && checker.converged(getIterations(), previous, current))
				{
					// We have found an optimum.
					converged = CHECKER;
					return current;
				}
			}

			// Move along the search direction.
			final double[] pnew;
			try
			{
				pnew = lineSearch.lineSearch(p, fp, g, xi);
			}
			catch (LineSearchRoundoffException e)
			{
				// This can happen if the Hessian is nearly singular or non-positive-definite.
				// In this case the algorithm should be restarted.
				converged = ROUNDOFF_ERROR;
				//System.out.printf("Roundoff error, iter=%d\n", getIterations());
				return new PointValuePair(p, fp);
			}

			// We assume the new point is on/within the bounds since the line search is constrained
			double fret = lineSearch.f;

			// Test for convergence on change in position
			if (positionChecker.converged(p, pnew))
			{
				converged = POSITION;
				return new PointValuePair(pnew, fret);
			}

			// Update the line direction
			for (int i = 0; i < n; i++)
			{
				xi[i] = pnew[i] - p[i];
			}
			p = pnew;

			// Save the old gradient
			double[] dg = g;

			// Get the gradient for the new point
			g = computeObjectiveGradient(p);
			checkGradients(g, p);

			// If necessary recompute the function value. 
			// Doing this after the gradient evaluation allows the value to be cached when 
			// computing the objective gradient
			fp = fret;

			// Test for convergence on zero gradient.
			double test = 0;
			for (int i = 0; i < n; i++)
			{
				final double temp = Math.abs(g[i]) * FastMath.max(Math.abs(p[i]), 1);
				//final double temp = Math.abs(g[i]);
				if (test < temp)
					test = temp;
			}
			// Compute the biggest gradient relative to the objective function
			test /= FastMath.max(Math.abs(fp), 1);
			if (test < gradientTolerance)
			{
				converged = GRADIENT;
				return new PointValuePair(p, fp);
			}

			for (int i = 0; i < n; i++)
				dg[i] = g[i] - dg[i];
			for (int i = 0; i < n; i++)
			{
				hdg[i] = 0.0;
				for (int j = 0; j < n; j++)
					hdg[i] += hessian[i][j] * dg[j];
			}
			double fac = 0, fae = 0, sumdg = 0, sumxi = 0;
			for (int i = 0; i < n; i++)
			{
				fac += dg[i] * xi[i];
				fae += dg[i] * hdg[i];
				sumdg += dg[i] * dg[i];
				sumxi += xi[i] * xi[i];
			}
			if (fac > Math.sqrt(EPS * sumdg * sumxi))
			{
				fac = 1.0 / fac;
				final double fad = 1.0 / fae;
				for (int i = 0; i < n; i++)
					dg[i] = fac * xi[i] - fad * hdg[i];
				for (int i = 0; i < n; i++)
				{
					for (int j = i; j < n; j++)
					{
						hessian[i][j] += fac * xi[i] * xi[j] - fad * hdg[i] * hdg[j] + fae * dg[i] * dg[j];
						hessian[j][i] = hessian[i][j];
					}
				}
			}
			for (int i = 0; i < n; i++)
			{
				xi[i] = 0.0;
				for (int j = 0; j < n; j++)
					xi[i] -= hessian[i][j] * g[j];
			}
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
	 *            <li>{@link GradientChecker}</li>
	 *            <li>{@link PositionChecker}</li>
	 *            <li>{@link MaximumStepLength}</li>
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
			}
			else if (data instanceof StepLength)
			{
				maximumStepLength = ((StepLength) data).getStep();
			}
			else if (data instanceof GradientTolerance)
			{
				gradientTolerance = ((GradientTolerance) data).getTolerance();
			}
			else if (data instanceof MaximumRestarts)
			{
				restarts = ((MaximumRestarts) data).getRestarts();
			}
			else if (data instanceof MaximumRoundoffRestarts)
			{
				roundoffRestarts = ((MaximumRoundoffRestarts) data).getRestarts();
			}
		}

		checkParameters();
	}

	/**
	 * The minimum value between two doubles
	 */
	private static double epsilon = calculateMachineEpsilonDouble();

	/**
	 * @return The minimum value between two doubles
	 * @see http://en.wikipedia.org/wiki/Machine_epsilon#Approximation_using_C.2B.2B
	 */
	private static double calculateMachineEpsilonDouble()
	{
		double machEps = 1.0;

		do
			machEps /= 2.0;
		while ((1.0 + (machEps / 2.0)) != 1.0);

		// ISO standard is 2^-52 = 2.220446049e-16

		return machEps;
	}

	public static class LineSearchRoundoffException extends RuntimeException
	{
		private static final long serialVersionUID = -8974644703023090107L;
		private final double slope;

		public LineSearchRoundoffException(double slope)
		{
			super();
			this.slope = slope;
		}

		public LineSearchRoundoffException()
		{
			super();
			this.slope = 0;
		}

		@Override
		public String getMessage()
		{
			return (slope != 0) ? "Round-off problem. Slope = " + slope : "Round-off problem";
		}
	}

	/**
	 * Internal class for a line search with backtracking
	 * <p>
	 * Adapted from NR::lnsrch, as discussed in Numerical Recipes section 9.7. The algorithm has been changed to support
	 * bounds on the point, limits on the search direction in all dimensions and checking for bad function evaluations
	 * when backtracking.
	 */
	private class LineStepSearch
	{
		/**
		 * Set to true when the the new point is too close to the old point. In a minimisation algorithm this signifies
		 * convergence.
		 */
		@SuppressWarnings("unused")
		boolean check;
		/**
		 * The function value at the new point
		 */
		double f;

		/**
		 * Given an n-dimension point, the function value and gradient at that point find a new point
		 * along the given search direction so that the function value has decreased sufficiently.
		 * 
		 * @param xOld
		 *            The old point
		 * @param fOld
		 *            The old point function value
		 * @param gradient
		 *            The old point function gradient
		 * @param searchDirection
		 *            The search direction
		 * @return The new point
		 * @throws LineSearchRoundoffException
		 *             if the slope of the line search is positive
		 */
		double[] lineSearch(double[] xOld, final double fOld, double[] gradient, double[] searchDirection)
				throws LineSearchRoundoffException
		{
			final double ALF = 1.0e-4, TOLX = epsilon;
			double alam2 = 0.0, f2 = 0.0;

			// New point
			double[] x = new double[xOld.length];

			final int n = xOld.length;
			check = false;

			// Limit the search step size for each dimension
			if (maximumStepLength != null)
			{
				double scale = 1;
				for (int i = 0; i < n; i++)
				{
					if (Math.abs(searchDirection[i]) * scale > maximumStepLength[i])
						scale = maximumStepLength[i] / Math.abs(searchDirection[i]);
				}
				if (scale < 1)
				{
					// Scale the entire search direction
					for (int i = 0; i < n; i++)
						searchDirection[i] *= scale;
				}
			}

			double slope = 0.0;
			for (int i = 0; i < n; i++)
				slope += gradient[i] * searchDirection[i];
			if (slope >= 0.0)
			{
				throw new LineSearchRoundoffException(slope);
			}

			// Compute lambda min
			double test = 0.0;
			for (int i = 0; i < n; i++)
			{
				final double temp = Math.abs(searchDirection[i]) / FastMath.max(Math.abs(xOld[i]), 1.0);
				if (temp > test)
					test = temp;
			}
			double alamin = TOLX / test;

			// Always try the full step first
			double alam = 1.0;
			// Count the number of backtracking steps
			int backtracking = 0;
			for (;;)
			{
				if (alam < alamin)
				{
					// Convergence (insignificant step).
					// Since we use the old f and x then we do not need to compute the objective value
					check = true;
					f = fOld;
					//System.out.printf("alam %f < alamin %f\n", alam, alamin);
					return xOld;
				}

				for (int i = 0; i < n; i++)
					x[i] = xOld[i] + alam * searchDirection[i];
				applyBounds(x);
				f = BFGSOptimizer.this.computeObjectiveValue(x);
				//System.out.printf("f=%f @ %f : %s\n", f, alam, java.util.Arrays.toString(x));
				if (f <= fOld + ALF * alam * slope)
				{
					// Sufficient function decrease
					//System.out.printf("f=%f < %f\n", f, fOld + ALF * alam * slope);
					return x;
				}
				else
				{
					// Check for bad function evaluation
					if (f == Double.POSITIVE_INFINITY)
					{
						// Reset backtracking
						backtracking = 0;

						alam *= 0.1;
						continue;
					}

					// Backtrack
					double tmplam;
					if (backtracking++ == 0)
					{
						// First backtrack iteration
						tmplam = -slope / (2.0 * (f - fOld - slope));
						// Ensure the lambda is reduced, i.e. we take a step smaller than last time
						if (tmplam > 0.9 * alam)
							tmplam = 0.9 * alam;
					}
					else
					{
						// Subsequent backtracks
						final double rhs1 = f - fOld - alam * slope;
						final double rhs2 = f2 - fOld - alam2 * slope;
						final double a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
						final double b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) /
								(alam - alam2);
						if (a == 0.0)
							tmplam = -slope / (2.0 * b);
						else
						{
							final double disc = b * b - 3.0 * a * slope;
							if (disc < 0.0)
								tmplam = 0.5 * alam;
							else if (b <= 0.0)
								tmplam = (-b + Math.sqrt(disc)) / (3.0 * a);
							else
								tmplam = -slope / (b + Math.sqrt(disc));
						}
						// Ensure the lambda is <= 0.5 lamda1, i.e. we take a step smaller than last time
						if (tmplam > 0.5 * alam)
							tmplam = 0.5 * alam;
					}

					alam2 = alam;
					f2 = f;
					// Ensure the lambda is >= 0.1 lamda1, i.e. we take reasonable step
					alam = FastMath.max(tmplam, 0.1 * alam);
				}
			}
		}
	}

	/**
	 * Checks if there are lower or upper bounds that are not -Infinity or +Infinity
	 * 
	 * @throws MathUnsupportedOperationException
	 *             if invalid bounds were passed to the {@link #optimize(OptimizationData[]) optimize} method.
	 */
	private void checkParameters()
	{
		lower = getLowerBound();
		upper = getUpperBound();
		isLower = checkArray(lower, Double.NEGATIVE_INFINITY);
		isUpper = checkArray(upper, Double.POSITIVE_INFINITY);
		// Check that the upper bound is above the lower bound
		if (isUpper && isLower)
		{
			for (int i = 0; i < lower.length; i++)
				if (lower[i] > upper[i])
					throw new MathUnsupportedOperationException(createError("Lower bound must be below upper bound"));
		}

		// Numerical Recipes set the position convergence very low
		if (positionChecker == null)
			positionChecker = new PositionChecker(4 * epsilon, 0);

		// Ensure that the step length is strictly positive
		if (maximumStepLength == null)
		{
			for (int i = 0; i < maximumStepLength.length; i++)
			{
				if (maximumStepLength[i] <= 0)
					throw new MathUnsupportedOperationException(
							createError("Maximum step length must be strictly positive"));
			}
		}

		// Set a tolerance? If not then the routine will iterate until position convergence
		//if (gradientTolerance == 0)
		//	gradientTolerance = 1e-6;
	}

	private Localizable createError(final String message)
	{
		return new Localizable()
		{
			private static final long serialVersionUID = 1L;

			@Override
			public String getSourceString()
			{
				return message;
			}

			@Override
			public String getLocalizedString(Locale locale)
			{
				return message;
			}
		};
	}

	/**
	 * Check if the array contains anything other than value
	 * 
	 * @param array
	 * @param value
	 * @return True if the array has another value
	 */
	private boolean checkArray(double[] array, double value)
	{
		if (array == null)
			return false;
		for (double v : array)
			if (v != value)
				return true;
		return false;
	}

	/**
	 * Check the point falls within the configured bounds truncating if necessary
	 * 
	 * @param point
	 * @return true if the point was truncated
	 */
	private boolean applyBounds(double[] point)
	{
		boolean truncated = false;
		if (isUpper)
		{
			for (int i = 0; i < point.length; i++)
				if (point[i] > upper[i])
				{
					point[i] = upper[i];
					truncated = true;
				}
		}
		if (isLower)
		{
			for (int i = 0; i < point.length; i++)
				if (point[i] < lower[i])
				{
					point[i] = lower[i];
					truncated = true;
				}
		}
		return truncated;
	}

	/**
	 * Check if the point falls on or outside configured bounds truncating the gradient to zero
	 * if it is moving further outside the bounds
	 * 
	 * @param r
	 * @param point
	 */
	private void checkGradients(double[] r, double[] point)
	{
		checkGradients(r, point, sign);
	}

	/**
	 * Check if the point falls on or outside configured bounds truncating the gradient to zero
	 * if it is moving further outside the bounds (defined by the sign of the search direction)
	 * 
	 * @param r
	 * @param point
	 * @param sign
	 */
	private void checkGradients(double[] r, double[] point, final double sign)
	{
		if (isUpper)
		{
			for (int i = 0; i < point.length; i++)
				if (point[i] >= upper[i] && Math.signum(r[i]) == sign)
					r[i] = 0;
		}
		if (isLower)
		{
			for (int i = 0; i < point.length; i++)
				if (point[i] <= lower[i] && Math.signum(r[i]) == -sign)
					r[i] = 0;
		}
		//boolean isNaN = false;
		//for (int i = 0; i < point.length; i++)
		//	if (Double.isNaN(r[i]))
		//	{
		//		isNaN = true;
		//		r[i] = 0;
		//	}
	}
}
