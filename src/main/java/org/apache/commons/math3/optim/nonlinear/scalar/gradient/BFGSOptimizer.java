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

import org.apache.commons.math3.exception.MathUnsupportedOperationException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.PointValuePair;
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
	private double maximumStepLength = 100;

	/** Convergence tolerance on gradient */
	private GradientChecker gradientChecker = null;

	/** Convergence tolerance on position */
	private PositionChecker positionChecker = null;

	/** Flags to indicate if bounds are present */
	private boolean isLower, isUpper;
	private double[] lower, upper;

	private double sign = 0;

	/**
	 * Check if the coordinates have converged
	 */
	public abstract static class CoordinateChecker implements OptimizationData
	{
		final double relative, absolute;

		/**
		 * Build an instance with specified thresholds.
		 *
		 * In order to perform only relative checks, the absolute tolerance
		 * must be set to a negative value. In order to perform only absolute
		 * checks, the relative tolerance must be set to a negative value.
		 *
		 * @param relativeThreshold
		 *            relative tolerance threshold
		 * @param absoluteThreshold
		 *            absolute tolerance threshold
		 */
		public CoordinateChecker(double relative, double absolute)
		{
			this.relative = relative;
			this.absolute = absolute;
		}

		/**
		 * Check if the coordinates have converged
		 * 
		 * @param p
		 *            Previous
		 * @param c
		 *            Current
		 * @return True if converged
		 */
		public boolean converged(final double[] p, final double[] c)
		{
			for (int i = 0; i < p.length; ++i)
			{
				final double pi = p[i];
				final double ci = c[i];
				final double difference = FastMath.abs(pi - ci);
				final double size = FastMath.max(FastMath.abs(pi), FastMath.abs(ci));
				if (difference > size * relative && difference > absolute)
				{
					return false;
				}
			}
			return true;
		}
	}

	/**
	 * Check if the gradient has converged
	 */
	public static class GradientChecker extends CoordinateChecker
	{
		/**
		 * Build an instance with specified thresholds.
		 *
		 * In order to perform only relative checks, the absolute tolerance
		 * must be set to a negative value. In order to perform only absolute
		 * checks, the relative tolerance must be set to a negative value.
		 *
		 * @param relativeThreshold
		 *            relative tolerance threshold
		 * @param absoluteThreshold
		 *            absolute tolerance threshold
		 */
		public GradientChecker(double relative, double absolute)
		{
			super(relative, absolute);
		}
	}

	/**
	 * Check if the position has converged
	 */
	public static class PositionChecker extends CoordinateChecker
	{
		/**
		 * Build an instance with specified thresholds.
		 *
		 * In order to perform only relative checks, the absolute tolerance
		 * must be set to a negative value. In order to perform only absolute
		 * checks, the relative tolerance must be set to a negative value.
		 *
		 * @param relativeThreshold
		 *            relative tolerance threshold
		 * @param absoluteThreshold
		 *            absolute tolerance threshold
		 */
		public PositionChecker(double relative, double absolute)
		{
			super(relative, absolute);
		}
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

	/** {@inheritDoc} */
	@Override
	protected PointValuePair doOptimize()
	{
		final ConvergenceChecker<PointValuePair> checker = getConvergenceChecker();
		double[] p = getStartPoint();
		final int n = p.length;

		// Assume minimisation
		sign = -1;

		LineStepSearch lineSearch = new LineStepSearch();

		final double EPS = epsilon;

		double[] dg = new double[n];
		double[] hdg = new double[n];
		double[] xi = new double[n];
		double[][] hessin = new double[n][n];

		// Get the gradient for the the bounded point
		double[] unbounded = p.clone();
		applyBounds(p);
		double[] g = computeObjectiveGradient(p);
		checkGradients(g, unbounded);

		double sum = 0;
		for (int i = 0; i < n; i++)
		{
			hessin[i][i] = 1.0;
			xi[i] = -g[i];
			sum += p[i] * p[i];
		}
		final double stpmax = maximumStepLength * FastMath.max(FastMath.sqrt(sum), n);

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
					return current;
				}
			}

			// Move along the search direction.
			final double[] pnew = lineSearch.lineSearch(p, fp, g, xi, stpmax);

			// TODO - Can we assume the new point is within the bounds and remove the 
			// bounds check on the point and the gradient?
			unbounded = pnew.clone();
			final double fret;
			if (applyBounds(pnew))
			{
				// Recompute the objective value within the bounds
				fret = computeObjectiveValue(pnew);
			}
			else
			{
				fret = lineSearch.f;
			}

			// Test for convergence on change in position
			if (positionChecker.converged(p, pnew))
			{
				//System.out.printf("Position converged\n");
				return new PointValuePair(pnew, fret);
			}

			fp = fret;
			// Update the line direction
			for (int i = 0; i < n; i++)
			{
				xi[i] = pnew[i] - p[i];
			}
			p = pnew;

			// Save the old gradient
			dg = g;

			// Get the gradient for the the bounded point
			g = computeObjectiveGradient(p);
			checkGradients(g, unbounded);

			// Test for convergence on zero gradient. 
			if (gradientChecker.converged(dg, g))
			{
				//System.out.printf("Gradient converged\n");
				return new PointValuePair(p, fp);
			}

			for (int i = 0; i < n; i++)
				dg[i] = g[i] - dg[i];
			for (int i = 0; i < n; i++)
			{
				hdg[i] = 0.0;
				for (int j = 0; j < n; j++)
					hdg[i] += hessin[i][j] * dg[j];
			}
			double fac = 0, fae = 0, sumdg = 0, sumxi = 0;
			for (int i = 0; i < n; i++)
			{
				fac += dg[i] * xi[i];
				fae += dg[i] * hdg[i];
				sumdg += dg[i] * dg[i];
				sumxi += xi[i] * xi[i];
			}
			if (fac > FastMath.sqrt(EPS * sumdg * sumxi))
			{
				fac = 1.0 / fac;
				final double fad = 1.0 / fae;
				for (int i = 0; i < n; i++)
					dg[i] = fac * xi[i] - fad * hdg[i];
				for (int i = 0; i < n; i++)
				{
					for (int j = i; j < n; j++)
					{
						hessin[i][j] += fac * xi[i] * xi[j] - fad * hdg[i] * hdg[j] + fae * dg[i] * dg[j];
						hessin[j][i] = hessin[i][j];
					}
				}
			}
			for (int i = 0; i < n; i++)
			{
				xi[i] = 0.0;
				for (int j = 0; j < n; j++)
					xi[i] -= hessin[i][j] * g[j];
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
			if (data instanceof GradientChecker)
			{
				gradientChecker = (GradientChecker) data;
			}
			else if (data instanceof PositionChecker)
			{
				positionChecker = (PositionChecker) data;
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

	/**
	 * Internal class for a line search with backtracking
	 * <p>
	 * Adapted from NR::lnsrch, as discussed in Numerical Recipes section 9.7. The algorithm has been changed to support
	 * bounds on the point and checking for bad function evaluations when backtracking.
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
		 * @param stepMax
		 *            The maximum step to take in the search direction
		 * @return The new point
		 */
		double[] lineSearch(double[] xOld, final double fOld, double[] gradient, double[] searchDirection,
				final double stepMax)
		{
			final double ALF = 1.0e-4, TOLX = epsilon;
			double alam2 = 0.0, f2 = 0.0;

			// New point
			double[] x = new double[xOld.length];

			final int n = xOld.length;
			check = false;

			// Limit the search step
			// TODO - add support for limiting the step size for each dimension			
			double sum = 0.0;
			for (int i = 0; i < n; i++)
				sum += searchDirection[i] * searchDirection[i];
			sum = FastMath.sqrt(sum);
			if (sum > stepMax)
				for (int i = 0; i < n; i++)
					// Scale if attempted step is too big
					searchDirection[i] *= stepMax / sum;
					
			double slope = 0.0;
			for (int i = 0; i < n; i++)
				slope += gradient[i] * searchDirection[i];
			if (slope >= 0.0)
			{
				throw new RuntimeException("Roundoff problem in line search");
			}
			// Compute lambda min
			double test = 0.0;
			for (int i = 0; i < n; i++)
			{
				final double temp = FastMath.abs(searchDirection[i]) / FastMath.max(FastMath.abs(xOld[i]), 1.0);
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
					for (int i = 0; i < n; i++)
						x[i] = xOld[i];
					check = true;
					f = fOld;
					//System.out.printf("alam %f < alamin %f\n", alam, alamin);
					return x;
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
								tmplam = (-b + FastMath.sqrt(disc)) / (3.0 * a);
							else
								tmplam = -slope / (b + FastMath.sqrt(disc));
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
					throw new MathUnsupportedOperationException(LocalizedFormats.CONSTRAINT);
		}
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
	 * Check if the point falls outside configured bounds truncating the gradient to zero
	 * if it is moving further outside the bounds
	 * 
	 * @param r
	 * @param point
	 * @return true if NaN gradients
	 */
	private boolean checkGradients(double[] r, double[] point)
	{
		return checkGradients(r, point, sign);
	}

	/**
	 * Check if the point falls outside configured bounds truncating the gradient to zero
	 * if it is moving further outside the bounds (defined by the sign of the search direction)
	 * 
	 * @param r
	 * @param point
	 * @param sign
	 * @return true if NaN gradients
	 */
	private boolean checkGradients(double[] r, double[] point, final double sign)
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
		boolean isNaN = false;
		for (int i = 0; i < point.length; i++)
			if (Double.isNaN(r[i]))
			{
				isNaN = true;
				r[i] = 0;
			}
		return isNaN;
	}

	/**
	 * @return The maximum step length used in line search
	 */
	public double getMaximumStepLength()
	{
		return maximumStepLength;
	}

	/**
	 * Set the maximum step length used in line search
	 * 
	 * @param maximumStepLength
	 *            the maximum step length used in line search
	 */
	public void setMaximumStepLength(double maximumStepLength)
	{
		if (maximumStepLength < epsilon)
			throw new IllegalArgumentException("Maximum step length must be higher than machine epsilon: " + epsilon);
		this.maximumStepLength = maximumStepLength;
	}
}
