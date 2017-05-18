package gdsc.smlm.fitting.nonlinear;

import gdsc.core.utils.Maths;
import gdsc.smlm.fitting.FunctionSolverType;
import gdsc.smlm.fitting.LSEFunctionSolver;
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
 * Abstract class with utility methods for the LSEFunctionSolver interface.
 */
public abstract class LSEBaseFunctionSolver extends BaseFunctionSolver implements LSEFunctionSolver
{
	protected double totalSumOfSquares = Double.NaN;

	/**
	 * Default constructor
	 * 
	 * @throws NullPointerException
	 *             if the function is null
	 */
	public LSEBaseFunctionSolver(GradientFunction f)
	{
		super(FunctionSolverType.LSE, f);
	}

	@Override
	protected void preProcess()
	{
		totalSumOfSquares = Double.NaN;
	}

	/**
	 * Gets the total sum of squares.
	 *
	 * @param y
	 *            the y
	 * @return the total sum of squares
	 */
	public static double getTotalSumOfSquares(double[] y)
	{
		double sx = 0, ssx = 0;
		for (int i = y.length; i-- > 0;)
		{
			sx += y[i];
			ssx += y[i] * y[i];
		}
		final double sumOfSquares = ssx - (sx * sx) / (y.length);
		return sumOfSquares;
	}

	/**
	 * Compute the error
	 * 
	 * @param value
	 * @param noise
	 * @param numberOfFittedPoints
	 * @param numberOfFittedParameters
	 * @return the error
	 */
	public static double getError(double value, double noise, int numberOfFittedPoints, int numberOfFittedParameters)
	{
		double error = value;

		// Divide by the uncertainty in the individual measurements yi to get the chi-squared
		if (noise > 0)
		{
			error /= numberOfFittedPoints * noise * noise;
		}

		// This updates the chi-squared value to the average error for a single fitted
		// point using the degrees of freedom (N-M)?
		// Note: This matches the mean squared error output from the MatLab fitting code.
		// If a noise estimate was provided for individual measurements then this will be the
		// reduced chi-square (see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2892436/)
		if (numberOfFittedPoints > numberOfFittedParameters)
			error /= (numberOfFittedPoints - numberOfFittedParameters);
		else
			error = 0;

		return error;
	}

	/**
	 * Update the total and residual sum-of-squares using the given data.
	 *
	 * @param y
	 *            The data
	 * @param residuals
	 *            The fit residuals
	 */
	public void updateSumOfSquares(double[] y, double[] residuals)
	{
		totalSumOfSquares = getTotalSumOfSquares(y);
		value = 0;
		for (int i = y.length; i-- > 0;)
		{
			value += residuals[i] * residuals[i];
		}
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
			totalSumOfSquares = getTotalSumOfSquares(lastY);
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
