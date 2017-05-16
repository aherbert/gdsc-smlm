package gdsc.smlm.fitting.nonlinear;

import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.function.GradientFunction;
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
 */
public class BoundedNonLinearFit extends NonLinearFit
{
	private ParameterBounds bounds;

	/**
	 * Default constructor
	 * 
	 * @param func
	 *            The function to fit
	 */
	public BoundedNonLinearFit(NonLinearFunction func)
	{
		super(func, null);
		bounds = new ParameterBounds(func);
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
		bounds = new ParameterBounds(func);
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
		bounds = new ParameterBounds(func);
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

		if (bounds.atBounds(a))
		{
			//System.out.printf("Failed when point was at the bounds\n");

			// This functionality has been removed since the solver no longer
			// extracts rows/columns from the matrix if the gradient vector is zero.
			// TODO - Determine if this support is needed.
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
		bounds.applyBounds(a, da, ap);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.NonLinearFit#accepted(double[], double[], int)
	 */
	@Override
	protected void accepted(double[] a, double[] ap, int m)
	{
		bounds.accepted(a, ap);
		super.accepted(a, ap, m);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.NonLinearFit#computeFit(double[], double[], double[], double[])
	 */
	@Override
	public FitStatus computeFit(double[] y, double[] y_fit, double[] a, double[] a_dev)
	{
		bounds.initialise();
		return super.computeFit(y, y_fit, a, a_dev);
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
		bounds.setBounds(lowerB, upperB);
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
		bounds.setClampValues(clampValues);
	}

	/**
	 * Checks if is dynamic clamping. The clamping factor will be reduced by a factor of 2 when the direction changes.
	 *
	 * @return true, if is dynamic clamping
	 */
	public boolean isDynamicClamp()
	{
		return bounds.isDynamicClamp();
	}

	/**
	 * Set to true to reduce the clamp factor by a factor of when the direction changes.
	 *
	 * @param dynamicClamp
	 *            the new dynamic clamp
	 */
	public void setDynamicClamp(boolean dynamicClamp)
	{
		bounds.setDynamicClamp(dynamicClamp);
	}

	/**
	 * Warning: If the function is changed then the clamp values may require updating. However setting a new function
	 * does not set the clamp values to null to allow caching when the clamp values are unchanged, e.g. evaluation of a
	 * different function in the same parameter space.
	 * <p>
	 * Setting a new function removes the current bounds.
	 * 
	 * @param f
	 *            the new gradient function
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#setGradientFunction(gdsc.smlm.function.GradientFunction)
	 */
	@Override
	public void setGradientFunction(GradientFunction f)
	{
		super.setGradientFunction(f);
		if (bounds != null)
			bounds.setGradientFunction(f);
	}
}
