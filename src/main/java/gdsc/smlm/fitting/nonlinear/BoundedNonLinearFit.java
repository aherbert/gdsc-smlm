package gdsc.smlm.fitting.nonlinear;

import java.util.Arrays;

import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.function.NonLinearFunction;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This is an adaption of the C-code contained in the CcpNmr Analysis Program:
 *   CCPN website (http://www.ccpn.ac.uk/). 
 * The CCPN code was based on Numerical Recipes. 
 *---------------------------------------------------------------------------*/

/**
 * Uses Levenberg-Marquardt method to fit a nonlinear model with coefficients (a) for a
 * set of data points (x, y).
 * <p>
 * Support bounded parameters using a hard-stop limit.
 * <p>
 * Support parameter clamping to prevent large parameter shifts. Optionally update the clamping when the search
 * direction changes.
 * <p>
 * Support ignoring an update to the LVM lambda parameter when the accepted step was not local (relative to the initial
 * parameter clamp values)
 */
public class BoundedNonLinearFit extends NonLinearFit
{
	private boolean isLower = false, isUpper = false;
	private double[] lower, upper;
	private boolean isClamped = false;
	private boolean nonLocalSearch = false;
	private double localSearch = 0;
	private double[] clampInitial, clamp;
	private int[] dir;
	private boolean dynamicClamp = false;

	/**
	 * Default constructor
	 * 
	 * @param func
	 *            The function to fit
	 */
	public BoundedNonLinearFit(NonLinearFunction func)
	{
		super(func, null);
	}

	/**
	 * Default constructor
	 * 
	 * @param func
	 *            The function to fit
	 * @param sc
	 *            The stopping criteria
	 */
	public BoundedNonLinearFit(NonLinearFunction func, StoppingCriteria sc)
	{
		super(func, sc);
	}

	/**
	 * Default constructor
	 * 
	 * @param func
	 *            The function to fit
	 * @param sc
	 *            The stopping criteria
	 * @param significantDigits
	 *            Validate the Levenberg-Marquardt fit solution to the specified number of significant digits
	 * @param maxAbsoluteError
	 *            Validate the Levenberg-Marquardt fit solution using the specified maximum absolute error
	 */
	public BoundedNonLinearFit(NonLinearFunction func, StoppingCriteria sc, int significantDigits,
			double maxAbsoluteError)
	{
		super(func, sc, significantDigits, maxAbsoluteError);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.NonLinearFit#solve(double[], int)
	 */
	protected boolean solve(double[] a, final int m)
	{
		if (super.solve(a, m))
			return true;

		// If using a bounded LVM is there a chance that the gradient against the bounds will 
		// be very large and effect the linear decomposition of the matrix? 
		// If decomposition fails try again but set the bounded params to zero (these are 
		// ignored by the solver), thus skipping these params for this iteration.

		if (atBounds(a))
		{
			//System.out.printf("Failed when point was at the bounds\n");
			createLinearProblem(m);
			ignoreAtBounds(a);

			// This handles the case when the entire set of params have been excluded
			if (solve(covar, da))
			{
				// TODO
				// See if this ever helps. It may just add overhead.
				// Add counters for:
				// - how often this occurs, 
				// - how often it allows a solution, and 
				// - how often the solution was accepted.
				System.out.printf("Ignoring parameters at the bounds allowed a solution to be found\n");
				return true;
			}
		}

		return false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.NonLinearFit#updateFitParameters(double[], int[], int, double[], double[])
	 */
	@Override
	protected void updateFitParameters(double[] a, int[] gradientIndices, int m, double[] da, double[] ap)
	{
		nonLocalSearch = false;

		if (isClamped)
		{
			for (int j = m; j-- > 0;)
			{
				if (clamp[j] == 0)
				{
					// Use the update parameter directly
					ap[gradientIndices[j]] = a[gradientIndices[j]] + da[j];
				}
				else
				{
					// This parameter is clamped
					ap[gradientIndices[j]] = a[gradientIndices[j]] + da[j] / clamp(da[j], j);
				}
			}
			applyBounds(ap, gradientIndices);

			// If using clamping should we can optionally only update lambda if we 
			// are close to the correct solution.
			if (localSearch != 0)
				nonLocalSearch = checkForNonLocalSearch(a, gradientIndices, m, ap);
		}
		else
		{
			for (int j = m; j-- > 0;)
			{
				// Use the update parameter directly
				ap[gradientIndices[j]] = a[gradientIndices[j]] + da[j];
			}
			applyBounds(ap, gradientIndices);
		}
	}

	/**
	 * Produce the clamping value.
	 * <p>
	 * See Stetson PB (1987) DAOPHOT: A compute program for crowded-field stellar photometry. Publ Astrom Soc Pac
	 * 99:191-222. pp207-208
	 *
	 * @param u
	 *            the update parameter
	 * @param k
	 *            the parameter index
	 * @return the clamping value
	 */
	private double clamp(double u, int k)
	{
		if (u == 0)
			// Nothing to clamp
			return 1;

		double ck = clamp[k];
		if (dynamicClamp)
		{
			// If the sign has changed then reduce the clamp factor
			final int newDir = (u > 0) ? 1 : -1;

			// This addition overcomes the issue when the direction vector is new (i.e. zero filled)
			if (newDir + dir[k] == 0)
			{
				// Note: By reducing the size of the clamping factor we are restricting the movement
				ck *= 0.5;
			}

			// Note: We do not update the clamp[k] array yet as the move may be rejected. 
		}

		// Denominator for clamping function
		return 1 + (Math.abs(u) / ck);
	}

	/**
	 * Check the parameter updates are within the local search parameter relative to the initial clamp values
	 * 
	 * @param a
	 *            the current fit parameters
	 * @param gradientIndices
	 *            the gradient indices (maps the fit parameter index to the parameter array)
	 * @param m
	 *            the number of fit parameters
	 * @param da
	 *            the parameter shift
	 * @param ap
	 *            the new fit parameters
	 * @return True if the search is non-local
	 */
	private boolean checkForNonLocalSearch(double[] a, int[] gradientIndices, int m, double[] ap)
	{
		// Check each update
		for (int j = m; j-- > 0;)
		{
			if (localSearch * Math.abs(ap[gradientIndices[j]] - a[gradientIndices[j]]) > clampInitial[j])
				return true;
		}
		return false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.NonLinearFit#accepted(double[], double[], int)
	 */
	@Override
	protected void accepted(double[] a, double[] ap, int m)
	{
		if (isClamped && dynamicClamp)
		{
			// Get the direction and update the clamp parameter if the direction has changed
			final int[] gradientIndices = f.gradientIndices();
			for (int k = m; k-- > 0;)
			{
				if (clamp[k] != 0)
				{
					final double u = ap[gradientIndices[k]] - a[gradientIndices[k]];
					if (u == 0)
						continue;
					final int newDir = (u > 0) ? 1 : -1;
					// This addition overcomes the issue when the direction vector is new (i.e. zero filled)
					if (newDir + dir[k] == 0)
					{
						// Note: By reducing the size of the clamping factor we are restricting the movement
						clamp[k] *= 0.5;
					}
					dir[k] = newDir;
				}
			}

		}
		if (nonLocalSearch)
		{
			// do not update the lambda parameter
			return;
		}
		super.accepted(a, ap, m);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.NonLinearFit#computeFit(int, double[], double[], double[], double[], double[],
	 * double)
	 */
	@Override
	public FitStatus computeFit(int n, double[] y, double[] y_fit, double[] a, double[] a_dev, double[] error,
			double noise)
	{
		// Initialise for clamping
		if (isClamped)
		{
			// Prevent the clamping value being destroyed by dynamic updates
			if (dynamicClamp)
			{
				final int m = f.gradientIndices().length;
				clamp = Arrays.copyOf(clampInitial, m);
				for (int i = m; i-- > 0;)
					dir[i] = 0;
			}
			else
			{
				clamp = clampInitial;
			}
		}
		return super.computeFit(n, y, y_fit, a, a_dev, error, noise);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#isBounded()
	 */
	@Override
	public boolean isBounded()
	{
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#isConstrained()
	 */
	@Override
	public boolean isConstrained()
	{
		return false;
	}

	/**
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#setBounds(double[], double[])
	 * @throws IllegalArgumentException
	 *             If the lower bound is above the upper bound
	 */
	@Override
	public void setBounds(double[] lowerB, double[] upperB)
	{
		// Extract the bounds for the parameters we are fitting
		if (lowerB == null)
		{
			lower = null;
		}
		else
		{
			final int[] indices = f.gradientIndices();
			lower = new double[indices.length];
			for (int i = 0; i < indices.length; i++)
			{
				lower[i] = lowerB[indices[i]];
			}
		}
		if (upperB == null)
		{
			upper = null;
		}
		else
		{
			final int[] indices = f.gradientIndices();
			upper = new double[indices.length];
			for (int i = 0; i < indices.length; i++)
			{
				upper[i] = upperB[indices[i]];
			}
		}
		isLower = checkArray(lower, Double.NEGATIVE_INFINITY);
		isUpper = checkArray(upper, Double.POSITIVE_INFINITY);
		// Check that the upper bound is above the lower bound
		if (isUpper && isLower)
		{
			for (int i = 0; i < lower.length; i++)
				if (lower[i] > upper[i])
					throw new IllegalArgumentException(
							"Lower bound is above upper bound: " + lower[i] + " > " + upper[i]);
		}
	}

	/**
	 * Check if the array contains anything other than value.
	 *
	 * @param array
	 *            the array
	 * @param value
	 *            the value
	 * @return True if the array has another value
	 */
	private static boolean checkArray(double[] array, double value)
	{
		if (array == null)
			return false;
		for (int i = 0; i < array.length; i++)
			if (value != array[i])
				return true;
		return false;
	}

	/**
	 * Check the point falls within the configured bounds truncating if necessary.
	 *
	 * @param point
	 *            the point
	 */
	private void applyBounds(double[] point, int[] gradientIndices)
	{
		if (isUpper)
		{
			for (int i = 0; i < gradientIndices.length; i++)
				if (point[gradientIndices[i]] > upper[i])
				{
					point[gradientIndices[i]] = upper[i];
				}
		}
		if (isLower)
		{
			for (int i = 0; i < gradientIndices.length; i++)
				if (point[gradientIndices[i]] < lower[i])
				{
					point[gradientIndices[i]] = lower[i];
				}
		}
	}

	/**
	 * Determine if the current solution (a) is at the the bounds
	 *
	 * @param a
	 *            the current parameters
	 * @return true, if the point is at the bounds
	 */
	private boolean atBounds(double[] a)
	{
		if (isUpper)
		{
			final int[] gradientIndices = f.gradientIndices();
			for (int i = 0; i < gradientIndices.length; i++)
				if (a[gradientIndices[i]] == upper[i])
				{
					return true;
				}
		}
		if (isLower)
		{
			final int[] gradientIndices = f.gradientIndices();
			for (int i = 0; i < gradientIndices.length; i++)
				if (a[gradientIndices[i]] == lower[i])
				{
					return true;
				}
		}
		return false;
	}

	/**
	 * If the current solution (a) is at the the bounds then set the gradient parameter to be solved (da) to zero. It
	 * will
	 * then be ignored.
	 *
	 * @param a
	 *            the current parameters
	 */
	private void ignoreAtBounds(double[] a)
	{
		if (isUpper)
		{
			final int[] gradientIndices = f.gradientIndices();
			for (int i = 0; i < gradientIndices.length; i++)
				if (a[gradientIndices[i]] == upper[i])
				{
					da[i] = 0;
				}
		}
		if (isLower)
		{
			final int[] gradientIndices = f.gradientIndices();
			for (int i = 0; i < gradientIndices.length; i++)
				if (a[gradientIndices[i]] == lower[i])
				{
					da[i] = 0;
				}
		}
	}

	/**
	 * Sets the parameter specific clamp values. This is the maximum permissible update to the parameter.
	 * <p>
	 * See Stetson PB (1987) DAOPHOT: A compute program for crowded-field stellar photometry. Publ Astrom Soc Pac
	 * 99:191-222.
	 * <p>
	 * Warning: If the function is changed then the clamp values may require updating. However setting a new function
	 * does not set the clamp values to null to allow caching when the clamp values are unchanged.
	 *
	 * @param clampValues
	 *            the new clamp values
	 */
	public void setClampValues(double[] clampValues)
	{
		// Extract the bounds for the parameters we are fitting
		final int[] indices = f.gradientIndices();

		if (clampValues == null)
		{
			clampInitial = null;
		}
		else
		{
			clampInitial = new double[indices.length];
			for (int i = 0; i < indices.length; i++)
			{
				final double v = clampValues[indices[i]];
				if (Double.isNaN(v) || Double.isInfinite(v))
					continue;
				clampInitial[i] = Math.abs(v);
			}
		}
		isClamped = checkArray(clampInitial, 0);
		if (isClamped && (dir == null || dir.length < clampInitial.length))
			dir = new int[clampInitial.length];
	}

	/**
	 * Checks if is dynamic clamping. The clamping factor will be reduced by a factor of 2 when the direction changes.
	 *
	 * @return true, if is dynamic clamping
	 */
	public boolean isDynamicClamp()
	{
		return dynamicClamp;
	}

	/**
	 * Set to true to reduce the clamp factor by a factor of when the direction changes.
	 *
	 * @param dynamicClamp
	 *            the new dynamic clamp
	 */
	public void setDynamicClamp(boolean dynamicClamp)
	{
		this.dynamicClamp = dynamicClamp;
	}

	/**
	 * @return the local search parameter
	 */
	public double getLocalSearch()
	{
		return localSearch;
	}

	/**
	 * When using clamping, if [update * local search parameter] > [initial clamp value] then the search is deemed to be
	 * non-local and lambda is not updated. This preserves the steepest descent search from the previous step.
	 * <p>
	 * Set to zero to disable.
	 * 
	 * @param localSearch
	 *            the local search parameter
	 */
	public void setLocalSearch(double localSearch)
	{
		this.localSearch = localSearch;
	}

	/**
	 * Warning: If the function is changed then the clamp values may require updating. However setting a new function
	 * does not set the clamp values to null to allow caching when the clamp values are unchanged, e.g. evaluation of a
	 * different function in the same parameter space.
	 * <p>
	 * Setting a new function removes the current bounds.
	 * 
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.NonLinearFit#setNonLinearFunction(gdsc.smlm.function.NonLinearFunction)
	 */
	@Override
	public void setNonLinearFunction(NonLinearFunction func)
	{
		// Do not do this to allow caching
		//setClampValues(null);

		setBounds(null, null);

		super.setNonLinearFunction(func);
	}
}
