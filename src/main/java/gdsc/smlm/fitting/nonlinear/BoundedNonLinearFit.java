package gdsc.smlm.fitting.nonlinear;

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
 */
public class BoundedNonLinearFit extends NonLinearFit
{
	private boolean isLower, isUpper;
	private double[] lower, upper;

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
	 * @see gdsc.smlm.fitting.nonlinear.NonLinearFit#updateFitParameters(double[], int[], int, double[], double[])
	 */
	@Override
	protected boolean updateFitParameters(double[] a, int[] gradientIndices, int m, double[] da, double[] ap)
	{
		for (int j = m; j-- > 0;)
			ap[gradientIndices[j]] = a[gradientIndices[j]] + da[j];
		return applyBounds(ap, gradientIndices);
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
		int[] indices = f.gradientIndices();

		if (lowerB == null)
		{
			lower = null;
		}
		else
		{
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
	 * @return true if truncated
	 */
	private boolean applyBounds(double[] point, int[] gradientIndices)
	{
		boolean truncated = false;
		if (isUpper)
		{
			for (int i = 0; i < gradientIndices.length; i++)
				if (point[gradientIndices[i]] > upper[i])
				{
					point[gradientIndices[i]] = upper[i];
					truncated = true;
				}
		}
		if (isLower)
		{
			for (int i = 0; i < gradientIndices.length; i++)
				if (point[gradientIndices[i]] < lower[i])
				{
					point[gradientIndices[i]] = lower[i];
					truncated = true;
				}
		}
		return truncated;
	}
}
