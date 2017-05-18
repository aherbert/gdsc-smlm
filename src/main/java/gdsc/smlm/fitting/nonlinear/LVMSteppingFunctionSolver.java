package gdsc.smlm.fitting.nonlinear;

import gdsc.core.utils.DoubleEquality;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolverType;
import gdsc.smlm.fitting.linear.EJMLLinearSolver;
import gdsc.smlm.fitting.nonlinear.gradient.LVMGradientProcedure;
import gdsc.smlm.function.Gradient1Function;
import gdsc.smlm.function.ValueFunction;
import gdsc.smlm.function.ValueProcedure;

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
 * Uses the Levenberg-Marquardt method to fit a gradient function with coefficients (a).
 */
public abstract class LVMSteppingFunctionSolver extends SteppingFunctionSolver
{
	protected EJMLLinearSolver solver = new EJMLLinearSolver();
	protected LVMGradientProcedure gradientProcedure;

	protected double initialLambda = 0.01;
	protected double lambda;

	// Alpha = Scaled Hessian matrix
	// beta  = Gradient vector
	// We want to solve: A x = b to find the update x

	// Current best alpha and beta 
	protected double[] alpha, beta;
	// Working alpha and beta 
	protected double[] walpha, wbeta;

	// TODO - Determine what a good solution tolerance would be. 
	// We may not need to be that strict to accept the solution.

	public static final double DEFAULT_MAX_RELATIVE_ERROR = 1e-3;
	public static final double DEFAULT_MAX_ABSOLUTE_ERROR = 1e-4;

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
	public LVMSteppingFunctionSolver(FunctionSolverType type, Gradient1Function f)
	{
		this(type, f, DEFAULT_MAX_RELATIVE_ERROR, DEFAULT_MAX_ABSOLUTE_ERROR);
	}

	/**
	 * Create a new stepping function solver.
	 *
	 * @param type
	 *            the type
	 * @param f
	 *            the function
	 * @param maxRelativeError
	 *            Validate the Levenberg-Marquardt fit solution using the specified maximum relative error
	 * @param maxAbsoluteError
	 *            Validate the Levenberg-Marquardt fit solution using the specified maximum absolute error
	 * @throws NullPointerException
	 *             if the function is null
	 */
	public LVMSteppingFunctionSolver(FunctionSolverType type, Gradient1Function f, double maxRelativeError,
			double maxAbsoluteError)
	{
		super(type, f);
		solver.setEqual(new DoubleEquality(maxRelativeError, maxAbsoluteError));
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
	public LVMSteppingFunctionSolver(FunctionSolverType type, Gradient1Function f, ToleranceChecker tc,
			ParameterBounds bounds)
	{
		this(type, f, tc, bounds, DEFAULT_MAX_RELATIVE_ERROR, DEFAULT_MAX_ABSOLUTE_ERROR);
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
	 * @param maxRelativeError
	 *            Validate the Levenberg-Marquardt fit solution using the specified maximum relative error
	 * @param maxAbsoluteError
	 *            Validate the Levenberg-Marquardt fit solution using the specified maximum absolute error
	 * @throws NullPointerException
	 *             if the function or tolerance checker is null
	 */
	public LVMSteppingFunctionSolver(FunctionSolverType type, Gradient1Function f, ToleranceChecker tc,
			ParameterBounds bounds, double maxRelativeError, double maxAbsoluteError)
	{
		super(type, f, tc, bounds);
		solver.setEqual(new DoubleEquality(maxRelativeError, maxAbsoluteError));
	}

	@Override
	protected double[] prepareFitValue(double[] y, double[] a)
	{
		// Ensure the gradient procedure is created
		y = prepareY(y);
		gradientProcedure = createGradientProcedure(y);

		// Set up the current best Hessian matrix and gradient parameter
		lambda = initialLambda;
		int n = gradientProcedure.n;
		alpha = null;
		beta = null;
		walpha = new double[n * n];
		wbeta = new double[n];

		return y;
	}

	/**
	 * Prepare Y for the gradient procedure, e.g. ensure positive values.
	 *
	 * @param y
	 *            the y
	 * @return the new y
	 */
	protected double[] prepareY(double[] y)
	{
		return y;
	}

	/**
	 * Creates the gradient procedure.
	 *
	 * @param y
	 *            the y
	 * @return the LVM gradient procedure
	 */
	protected abstract LVMGradientProcedure createGradientProcedure(double[] y);

	@Override
	protected double computeFitValue(double[] a)
	{
		gradientProcedure.gradient(a);

		if (gradientProcedure.isNaNGradients())
			throw new FunctionSolverException(FitStatus.INVALID_GRADIENTS);

		if (alpha == null)
		{
			// This is the first computation:
			// Set the current alpha and beta
			alpha = gradientProcedure.getAlphaLinear();
			beta = gradientProcedure.beta.clone();
		}
		else
		{
			// This is a subsequent computation:
			// Store the working alpha and beta which may be accepted
			gradientProcedure.getAlphaLinear(walpha);
			gradientProcedure.getBeta(wbeta);
		}

		return gradientProcedure.value;
	}

	@Override
	protected void computeStep(double[] step)
	{
		// Alpha = Scaled Hessian matrix
		// beta  = Gradient vector
		// We want to solve: A x = b to find the update x

		final int n = gradientProcedure.n;
		System.arraycopy(beta, 0, step, 0, n);
		System.arraycopy(alpha, 0, walpha, 0, alpha.length);
		final double scale = (1.0 + lambda);
		for (int i = 0; i < n; i += (n + 1))
		{
			// Scale the diagonal of the Hessian to favour direct descent
			walpha[i] *= scale;
		}
		if (!solver.solve(walpha, step))
			throw new FunctionSolverException(FitStatus.SINGULAR_NON_LINEAR_MODEL);
	}

	@Override
	protected boolean accept(double currentValue, double[] a, double newValue, double[] newA)
	{
		if (newValue <= currentValue)
		{
			// Update the current alpha and beta:
			// We can do this by swapping storage.
			double[] tmp = alpha;
			alpha = walpha;
			walpha = tmp;
			tmp = beta;
			beta = wbeta;
			wbeta = tmp;

			// Decrease Lambda
			lambda *= 0.1;
			return true;
		}
		else
		{
			// Increase Lambda
			lambda *= 10.0;
			return false;
		}
	}

	@Override
	protected double[] prepareFunctionValue(double[] y, double[] a)
	{
		// Ensure the gradient procedure is created 
		y = prepareY(y);
		gradientProcedure = createGradientProcedure(y);
		return y;
	}

	@Override
	protected double computeFunctionValue(double[] y_fit, double[] a)
	{
		gradientProcedure.value(a);
		if (y_fit != null)
			// This simple implementation causes a double evaluation of the function 
			// but only a single initialisation. The following method can be overriden
			// to use cached function values.
			computeValues(y_fit);
		return gradientProcedure.value;
	}

	private static class SimpleValueProcedure implements ValueProcedure
	{
		int i = 0;
		double[] y_fit;

		SimpleValueProcedure(double[] y_fit)
		{
			this.y_fit = y_fit;
		};

		public void execute(double value)
		{
			y_fit[i++] = value;
		}
	}
	
	/**
	 * Utility method to compute the function values using the preinitialised function.
	 * Sub-classes may override this if they have cached the function values from the 
	 * last execution of a forEach procedure.

	 * @see gdsc.smlm.fitting.nonlinear.SteppingFunctionSolver#computeValues(double[])
	 */
	@Override
	protected void computeValues(final double[] y_fit)
	{
		ValueFunction f = (ValueFunction) this.f;
		f.forEach(new SimpleValueProcedure(y_fit));
	}

	/**
	 * @param initialLambda
	 *            the initial lambda for the Levenberg-Marquardt fitting routine
	 */
	public void setInitialLambda(double initialLambda)
	{
		this.initialLambda = initialLambda;
	}

	/**
	 * @return the initialLambda
	 */
	public double getInitialLambda()
	{
		return initialLambda;
	}
}
