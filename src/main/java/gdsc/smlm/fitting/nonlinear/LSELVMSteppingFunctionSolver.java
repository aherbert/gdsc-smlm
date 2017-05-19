package gdsc.smlm.fitting.nonlinear;

import gdsc.core.utils.Maths;
import gdsc.smlm.fitting.FisherInformationMatrix;
import gdsc.smlm.fitting.FunctionSolverType;
import gdsc.smlm.fitting.LSEFunctionSolver;
import gdsc.smlm.fitting.nonlinear.gradient.LSQLVMGradientProcedureFactory;
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
 * Uses the Levenberg-Marquardt method to fit a gradient function with coefficients (a) using least squares estimation.
 */
public class LSELVMSteppingFunctionSolver extends LVMSteppingFunctionSolver implements LSEFunctionSolver
{
	/** The total sum of squares. */
	protected double totalSumOfSquares = Double.NaN;

	/**
	 * Create a new stepping function solver.
	 *
	 * @param f
	 *            the function
	 * @throws NullPointerException
	 *             if the function is null
	 */
	public LSELVMSteppingFunctionSolver(Gradient1Function f)
	{
		super(FunctionSolverType.LSE, f);
	}

	/**
	 * Create a new stepping function solver.
	 *
	 * @param f
	 *            the function
	 * @param maxRelativeError
	 *            Validate the Levenberg-Marquardt fit solution using the specified maximum relative error
	 * @param maxAbsoluteError
	 *            Validate the Levenberg-Marquardt fit solution using the specified maximum absolute error
	 * @throws NullPointerException
	 *             if the function is null
	 */
	public LSELVMSteppingFunctionSolver(Gradient1Function f, double maxRelativeError, double maxAbsoluteError)
	{
		super(FunctionSolverType.LSE, f, maxRelativeError, maxAbsoluteError);
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
	public LSELVMSteppingFunctionSolver(Gradient1Function f, ToleranceChecker tc, ParameterBounds bounds)
	{
		super(FunctionSolverType.LSE, f, tc, bounds);
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
	 * @param maxRelativeError
	 *            Validate the Levenberg-Marquardt fit solution using the specified maximum relative error
	 * @param maxAbsoluteError
	 *            Validate the Levenberg-Marquardt fit solution using the specified maximum absolute error
	 * @throws NullPointerException
	 *             if the function or tolerance checker is null
	 */
	public LSELVMSteppingFunctionSolver(Gradient1Function f, ToleranceChecker tc, ParameterBounds bounds,
			double maxRelativeError, double maxAbsoluteError)
	{
		super(FunctionSolverType.LSE, f, tc, bounds, maxRelativeError, maxAbsoluteError);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#preProcess()
	 */
	@Override
	protected void preProcess()
	{
		totalSumOfSquares = Double.NaN;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.LVMSteppingFunctionSolver#createGradientProcedure(double[])
	 */
	@Override
	protected LVMGradientProcedure createGradientProcedure(double[] y)
	{
		return LSQLVMGradientProcedureFactory.create(y, (Gradient1Function) f);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.SteppingFunctionSolver#computeFisherInformationMatrix()
	 */
	@Override
	protected FisherInformationMatrix computeFisherInformationMatrix()
	{
		// The last Hessian matrix should be stored in the working alpha.
		return new FisherInformationMatrix(walpha, beta.length);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.LSEFunctionSolver#getTotalSumOfSquares()
	 */
	public double getTotalSumOfSquares()
	{
		if (Double.isNaN(totalSumOfSquares) && lastY != null)
		{
			totalSumOfSquares = LSEBaseFunctionSolver.getTotalSumOfSquares(lastY);
		}
		return totalSumOfSquares;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.LSEFunctionSolver#getResidualSumOfSquares()
	 */
	public double getResidualSumOfSquares()
	{
		return value;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.LSEFunctionSolver#getCoefficientOfDetermination()
	 */
	public double getCoefficientOfDetermination()
	{
		return 1.0 - (value / getTotalSumOfSquares());
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.LSEFunctionSolver#getAdjustedCoefficientOfDetermination()
	 */
	public double getAdjustedCoefficientOfDetermination()
	{
		return Maths.getAdjustedCoefficientOfDetermination(value, getTotalSumOfSquares(), getNumberOfFittedPoints(),
				getNumberOfFittedParameters());
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.LSEFunctionSolver#getMeanSquaredError()
	 */
	public double getMeanSquaredError()
	{
		return value / (getNumberOfFittedPoints() - getNumberOfFittedParameters());
	}
}
