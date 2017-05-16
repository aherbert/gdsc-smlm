package gdsc.smlm.fitting.nonlinear;

import gdsc.core.utils.DoubleEquality;
import gdsc.smlm.fitting.FunctionSolverType;
import gdsc.smlm.fitting.linear.EJMLLinearSolver;
import gdsc.smlm.fitting.nonlinear.gradient.LVMGradientProcedure;
import gdsc.smlm.function.Gradient1Function;

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
public class LVMSteppingFunctionSolver extends SteppingFunctionSolver
{
	protected EJMLLinearSolver solver = new EJMLLinearSolver();
	protected LVMGradientProcedure calculator;

	// TODO - Determine what a good solution tolerance would be. 
	// We may not need to be that strict to accept the solution.

	public static final int DEFAULT_SIGNIFICANT_DIGITS = 3;
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
		this(type, f, DEFAULT_SIGNIFICANT_DIGITS, DEFAULT_MAX_ABSOLUTE_ERROR);
	}

	/**
	 * Create a new stepping function solver
	 *
	 * @param type
	 *            the type
	 * @param f
	 *            the function
	 * @param significantDigits
	 *            Validate the Levenberg-Marquardt fit solution to the specified number of significant digits
	 * @param maxAbsoluteError
	 *            Validate the Levenberg-Marquardt fit solution using the specified maximum absolute error
	 * @throws NullPointerException
	 *             if the function is null
	 */
	public LVMSteppingFunctionSolver(FunctionSolverType type, Gradient1Function f, int significantDigits,
			double maxAbsoluteError)
	{
		super(type, f);
		solver.setEqual(new DoubleEquality(significantDigits, maxAbsoluteError));
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
		this(type, f, tc, bounds, DEFAULT_SIGNIFICANT_DIGITS, DEFAULT_MAX_ABSOLUTE_ERROR);
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
	 * @param significantDigits
	 *            Validate the Levenberg-Marquardt fit solution to the specified number of significant digits
	 * @param maxAbsoluteError
	 *            Validate the Levenberg-Marquardt fit solution using the specified maximum absolute error
	 * @throws NullPointerException
	 *             if the function or tolerance checker is null
	 */
	public LVMSteppingFunctionSolver(FunctionSolverType type, Gradient1Function f, ToleranceChecker tc,
			ParameterBounds bounds, int significantDigits, double maxAbsoluteError)
	{
		super(type, f, tc, bounds);
		solver.setEqual(new DoubleEquality(significantDigits, maxAbsoluteError));
	}

	@Override
	protected double computeFitValue(double[] y, double[] a)
	{
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	protected void computeStep(double[] step)
	{
		// TODO Auto-generated method stub

	}

	@Override
	protected boolean accept(double currentValue, double[] a, double newValue, double[] newA)
	{
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	protected void computeDeviations(double[] a_dev)
	{
		// TODO Auto-generated method stub

	}

	@Override
	protected void computeFunctionValue(double[] y_fit)
	{
		// TODO Auto-generated method stub

	}

	@Override
	protected double computeFunctionValue(double[] y, double[] y_fit, double[] a)
	{
		// TODO Auto-generated method stub
		return 0;
	}
}
