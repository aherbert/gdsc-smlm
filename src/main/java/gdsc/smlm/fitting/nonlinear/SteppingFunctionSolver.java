package gdsc.smlm.fitting.nonlinear;

import gdsc.core.utils.BitFlags;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolverType;
import gdsc.smlm.function.Gradient1Function;
import gdsc.smlm.function.GradientFunction;

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
 * Abstract class for FunctionSolvers that use update steps to the current parameters.
 */
public abstract class SteppingFunctionSolver extends BaseFunctionSolver
{
	protected int[] gradientIndices;
	protected final ToleranceChecker tc;
	protected ParameterBounds bounds;

	/**
	 * Create a new stepping function solver
	 *
	 * @param type
	 *            the type
	 * @param f
	 *            the function
	 * @throws NullPointerException
	 *             if the function is null
	 */
	public SteppingFunctionSolver(FunctionSolverType type, Gradient1Function f)
	{
		this(type, f, new ToleranceChecker(1e-3, 1e-6), null);
	}

	/**
	 * Create a new stepping function solver
	 *
	 * @param type
	 *            the type
	 * @param f
	 *            the function
	 * @param tc
	 *            the tolerance checker
	 * @param bounds
	 *            the bounds
	 * @throws NullPointerException
	 *             if the function or tolerance checker is null
	 */
	public SteppingFunctionSolver(FunctionSolverType type, Gradient1Function f, ToleranceChecker tc,
			ParameterBounds bounds)
	{
		super(type, f);
		if (tc == null)
			throw new NullPointerException("Null tolerance checker");
		this.tc = tc;
		this.bounds = bounds;
	}

	/**
	 * Compute fit.
	 *
	 * @param y
	 *            the y
	 * @param y_fit
	 *            the y_fit
	 * @param a
	 *            the a
	 * @param a_dev
	 *            the a_dev
	 * @return the fit status
	 */
	public FitStatus computeFit(double[] y, double[] y_fit, double[] a, double[] a_dev)
	{
		// Lay out a simple iteration loop for a stepping solver.
		// The sub-class must compute the next step.
		// This class handles attenuation of the step.
		// The sub-class determines if the step is accepted or rejected.

		gradientIndices = f.gradientIndices();
		final double[] step = new double[gradientIndices.length];
		final double[] newA = a.clone();

		// Initialise for fitting
		if (bounds != null)
			bounds.initialise();

		try
		{
			lastY = y = prepareFitValue(y, a);

			// First evaluation
			double currentValue = computeFitValue(y, a);

			int status = 0;
			while (true)
			{
				// Compute next step
				computeStep(step);

				// Apply bounds to the step
				if (bounds != null)
					bounds.applyBounds(a, step, newA);

				// Evaluate
				double newValue = computeFitValue(y, newA);

				// Check stopping criteria
				status = tc.converged(currentValue, a, newValue, newA);
				if (status != 0)
				{
					value = newValue;
					System.arraycopy(newA, 0, a, 0, a.length);
					break;
				}

				// Check if the step was an improvement
				if (accept(currentValue, a, newValue, newA))
				{
					currentValue = newValue;
					System.arraycopy(newA, 0, a, 0, a.length);
					if (bounds != null)
						bounds.accepted(a, newA);
				}
			}

			if (BitFlags.anySet(status, ToleranceChecker.STATUS_CONVERGED))
			{
				if (a_dev != null)
					computeDeviations(a_dev);
				if (y_fit != null)
					computeFunctionValue(y_fit);
				return FitStatus.OK;
			}

			// Check the iterations
			if (BitFlags.areSet(status, ToleranceChecker.STATUS_MAX_ITERATIONS))
				return FitStatus.TOO_MANY_ITERATIONS;

			// We should not reach here unless we missed something
			return FitStatus.FAILED_TO_CONVERGE;
		}
		catch (FunctionSolverException e)
		{
			// XXX - debugging
			System.out.printf("%s failed: %s\n", getClass().getSimpleName(), e.fitStatus.getName());
			return e.fitStatus;
		}
	}

	/**
	 * Prepare y for fitting, e.g. ensure strictly positive values.
	 *
	 * @param y
	 *            the y
	 * @param a
	 *            the parameters
	 * @return the new y
	 */
	protected abstract double[] prepareFitValue(double[] y, double[] a);

	/**
	 * Compute the fit value using the parameters. This method is followed by a call to {@link #computeStep(double[])}
	 * so the step could be pre-computed here.
	 *
	 * @param y
	 *            the y
	 * @param a
	 *            the parameters
	 * @return the fit value
	 */
	protected abstract double computeFitValue(double[] y, double[] a);

	/**
	 * Compute the update step for the current parameters.
	 *
	 * @param step
	 *            the step
	 */
	protected abstract void computeStep(double[] step);

	/**
	 * Determine if the step should be accepted. If accepted then the current parameters and function value are updated
	 * and any bounds on the step size may be updated.
	 * <p>
	 * Note that although this class handles convergence on the value/parameters it is left to the sub-class to
	 * determine if each step should be accepted.
	 *
	 * @param currentValue
	 *            the current value
	 * @param a
	 *            the current parameters
	 * @param newValue
	 *            the new value
	 * @param newA
	 *            the new parameters
	 * @return true, if successful
	 */
	protected abstract boolean accept(double currentValue, double[] a, double newValue, double[] newA);

	/**
	 * Compute the deviations for the parameters a from the last call to
	 * {@link #computeFitValue(double[], double[])}.
	 *
	 * @param a_dev
	 *            the parameter deviations
	 */
	protected abstract void computeDeviations(double[] a_dev);

	/**
	 * Compute function value using the y and parameters a from the last call to
	 * {@link #computeFitValue(double[], double[])}.
	 *
	 * @param y_fit
	 *            the y fit
	 */
	protected abstract void computeFunctionValue(double[] y_fit);

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#computeValue(double[], double[], double[])
	 */
	public boolean computeValue(double[] y, double[] y_fit, double[] a)
	{
		gradientIndices = f.gradientIndices();
		lastY = y = prepareFunctionValue(y, a);
		if (y_fit == null)
			y_fit = new double[y.length];
		value = computeFunctionValue(y, y_fit, a);
		return true;
	}

	/**
	 * Prepare y for computing the function value, e.g. ensure strictly positive values.
	 *
	 * @param y
	 *            the y
	 * @param a
	 *            the parameters
	 * @return the new y
	 */
	protected abstract double[] prepareFunctionValue(double[] y, double[] a);

	/**
	 * Compute the function value.
	 *
	 * @param y
	 *            the y
	 * @param y_fit
	 *            the y fit (this will not be null)
	 * @param a
	 *            the parameters
	 * @return the function value
	 */
	protected abstract double computeFunctionValue(double[] y, double[] y_fit, double[] a);

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#isBounded()
	 */
	public boolean isBounded()
	{
		// Bounds are tighter than constraints and we support those
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#setBounds(double[], double[])
	 */
	public void setBounds(double[] lower, double[] upper)
	{
		if (bounds == null)
		{
			if (lower == null && upper == null)
				return;
			bounds = new ParameterBounds(f);
		}
		bounds.setBounds(lower, upper);
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
