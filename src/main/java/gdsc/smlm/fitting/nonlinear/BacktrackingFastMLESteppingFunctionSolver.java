package gdsc.smlm.fitting.nonlinear;

import gdsc.core.utils.Maths;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.function.Gradient2Function;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Uses the Fast MLE method to fit a gradient function with coefficients (a).
 * <p>
 * Calculates the Newton-Raphson update vector for a Poisson process using the first and second partial derivatives.
 * <p>
 * Ref: Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically minimum uncertainty.
 * Nature Methods 7, 373-375 (supplementary note), Eq. 12.
 * <p>
 * Ref: Huang et al, (2015). Video-rate nanoscopy using sCMOS camera–specific single-molecule localization algorithms.
 * Nature Methods 10, 653–658.
 */
public class BacktrackingFastMLESteppingFunctionSolver extends FastMLESteppingFunctionSolver
{
	private LineStepSearch lineSearch = new LineStepSearch();
	/** Maximum step length used in line search. */
	private double[] maximumStepLength = null;
	private double maximumStepSize = 0;

	private double[] aOld, searchDirection;
	private boolean firstEvaluation;

	/**
	 * The minimum value between two doubles.
	 * ISO standard is 2^-52 = 2.220446049e-16.
	 * This computes 2.220446049250313E-16.
	 */
	private static double epsilon = Math.ulp(1.0);

	/**
	 * Create a new stepping function solver.
	 *
	 * @param f
	 *            the function
	 * @param maxRelativeError
	 *            the max relative error
	 * @param maxAbsoluteError
	 *            the max absolute error
	 * @throws NullPointerException
	 *             if the function is null
	 */
	public BacktrackingFastMLESteppingFunctionSolver(Gradient2Function f, double maxRelativeError,
			double maxAbsoluteError)
	{
		super(f, maxRelativeError, maxAbsoluteError);
	}

	/**
	 * Create a new stepping function solver.
	 *
	 * @param f
	 *            the function
	 * @param tc
	 *            the tolerance checker
	 * @param bounds
	 *            the bounds
	 * @throws NullPointerException
	 *             if the function or tolerance checker is null
	 */
	public BacktrackingFastMLESteppingFunctionSolver(Gradient2Function f, ToleranceChecker tc, ParameterBounds bounds)
	{
		super(f, tc, bounds);
	}

	// Extend the Newton-Raphson method by implementing Line Search and Backtracking
	// (see Numerical Recipes in C++, 2nd Ed, page 388-389, function lnsrch).
	// This can still be done in the context of the stepping function solver.
	// Adjustments:
	// - Always ensure we compute the function value. This is needed to determine
	// if the step computed a better value.
	// - The first call to computeFitValue() computes derivatives. The value is obtained from 
	// computePseudoLogLikelihood() 
	// - All subsequent calls to computeFitValue() only call computeValue() and implement
	// backtracking if the value does not improve with the full Newton step.
	// - When a suitable value is achieved then this should be returned from computeFitValue()
	// - The next call to compute step must evaluate the derivatives.
	// - The function evaluations counter should be appropriately incremented

	@Override
	protected double[] prepareFitValue(double[] y, double[] a)
	{
		y = super.prepareFitValue(y, a);
		// We always compute the pseudolikelihood
		isPseudoLogLikelihood = true;
		firstEvaluation = true;

		// Configure maximum step length for each dimension using the bounds.
		// This is a simple check that can prevent wild Newton Raphson steps 
		// that require a lot of backtracking.
		maximumStepLength = null;
		if (maximumStepSize > 0)
		{
			double[] lower = bounds.getLower();
			double[] upper = bounds.getUpper();
			if (lower != null && upper != null)
			{
				maximumStepLength = new double[lower.length];
				for (int i = 0; i < maximumStepLength.length; i++)
				{
					maximumStepLength[i] = (upper[i] - lower[i]) * maximumStepSize;
					if (maximumStepLength[i] <= 0)
						maximumStepLength[i] = Double.POSITIVE_INFINITY;
				}
			}
		}

		return y;
	}

	@Override
	protected double computeFitValue(double[] a)
	{
		if (firstEvaluation)
		{
			firstEvaluation = false;

			// The first call to computeFitValue() computes derivatives. The value is obtained from 
			// computePseudoLogLikelihood() 
			computeGradients(a);

			evaluations++;
			ll = gradientProcedure.computePseudoLogLikelihood();

			searchDirection = new double[a.length];
			aOld = a.clone();

			return ll;
		}

		// All subsequent calls to computeFitValue() only evaluate the value and implement
		// backtracking if the value does not improve with the full Newton step.
		for (int i = 0; i < searchDirection.length; i++)
			// Configure the search direction with the full Newton step
			searchDirection[i] = a[i] - aOld[i];

		aOld = lineSearch.lineSearch(aOld, ll, gradientProcedure.d1, searchDirection);
		ll = lineSearch.f;

		// Update the parameters to reflect any backtracking
		System.arraycopy(aOld, 0, a, 0, a.length);

		return ll;
	}

	@Override
	protected void computeStep(double[] step)
	{
		if (tc.getIterations() > 0)
		{
			// After backtracking we must compute the derivatives.
			// Note we leave it to here (and not after the line search) so it 
			// can be skipped if convergence is achieved.
			computeGradients(aOld);
		}
		super.computeStep(step);
	}

	/**
	 * Internal class for a line search with backtracking
	 * <p>
	 * Adapted from NR::lnsrch, as discussed in Numerical Recipes section 9.7. The algorithm
	 * has been changed to find the maximum, support limits on the search direction in all
	 * dimensions and check for bad function evaluations when backtracking.
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
		 *            The old point function gradient (only for the gradient indices)
		 * @param searchDirection
		 *            The search direction
		 * @return The new point
		 * @throws FunctionSolverException
		 *             if the slope of the line search is positive
		 */
		double[] lineSearch(double[] xOld, final double fOld, double[] gradient, double[] searchDirection)
				throws FunctionSolverException
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

			// XXX
			// Note that the gradient is raw but the search direction may be after clamping
			// and bounds have been applied. This should be fixed.
			
			double slope = 0.0;
			final int[] gradientIndices = BacktrackingFastMLESteppingFunctionSolver.this.f.gradientIndices();
			for (int i = 0; i < gradient.length; i++)
				slope += gradient[i] * searchDirection[gradientIndices[i]];
			if (slope <= 0.0)
			{
				for (int i = 0; i < gradient.length; i++)
				{
					System.out.printf("[%d:%g %s] %g (%g/-(%g)=%g) x %g = %g\n", i, 
							xOld[gradientIndices[i]],
							Gaussian2DFunction.getName(gradientIndices[i]),
							gradient[i], 
							gradientProcedure.d1[i],
							gradientProcedure.d2[i],
							gradientProcedure.d1[i]/-gradientProcedure.d2[i],
							searchDirection[gradientIndices[i]],
							gradient[i] * searchDirection[gradientIndices[i]]);
					slope += gradient[i] * searchDirection[gradientIndices[i]];
				}
				throw new FunctionSolverException(FitStatus.LINE_SEARCH_ERROR, "Slope is negative: " + slope);
			}

			// Compute lambda min
			double test = 0.0;
			for (int i = 0; i < n; i++)
			{
				final double temp = Math.abs(searchDirection[i]) / max(Math.abs(xOld[i]), 1.0);
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
				// Compute the pseudoLikelihood
				evaluations++;
				gradientProcedure.computeValue(x);
				f = gradientProcedure.computePseudoLogLikelihood();
				//System.out.printf("f=%f @ %f : %s\n", f, alam, java.util.Arrays.toString(x));
				if (f >= fOld + ALF * alam * slope)
				{
					// Sufficient function decrease
					//System.out.printf("f=%f > %f\n", f, fOld + ALF * alam * slope);
					return x;
				}
				else
				{
					// Check for bad function evaluation
					if (!Maths.isFinite(f))
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
					alam = max(tmplam, 0.1 * alam);
				}
			}
		}
	}

	private static double max(double a, double b)
	{
		return (a > b) ? a : b;
	}

	/**
	 * Gets the maximum step size.
	 *
	 * @return the maximum step size
	 */
	public double getMaximumStepSize()
	{
		return maximumStepSize;
	}

	/**
	 * Sets the maximum step size. This is expressed as a factor of the upper-lower bound. Steps beyond this will cause
	 * the step to be truncated. Ignored if above 1. Set to below 1 to disable.
	 * <p>
	 * Note: Restriction of the steps is better controlled using the ParaneterBounds object.
	 *
	 * @param maximumStepSize
	 *            the new maximum step size
	 */
	public void setMaximumStepSize(double maximumStepSize)
	{
		this.maximumStepSize = (maximumStepSize >= 1) ? 0 : maximumStepSize;
	}

}
