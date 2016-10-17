package gdsc.smlm.fitting.nonlinear;

import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolver;
import gdsc.smlm.function.NonLinearFunction;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Abstract class with utility methods for the FunctionSolver interface.
 */
public abstract class BaseFunctionSolver implements FunctionSolver
{
	protected NonLinearFunction f;

	private int maxEvaluations = 20;

	protected double totalSumOfSquares;
	protected double residualSumOfSquares;
	protected int numberOfFittedPoints;
	protected int iterations;
	protected int evaluations;
	protected double value;

	/**
	 * Default constructor
	 * 
	 * @throws NullPointerException
	 *             if the function is null
	 */
	public BaseFunctionSolver(NonLinearFunction f)
	{
		setNonLinearFunction(f);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#fit(int, double[], double[], double[], double[], double[], double)
	 */
	public FitStatus fit(int n, double[] y, double[] y_fit, double[] a, double[] a_dev, double[] error, double noise)
	{
		// Reset the results
		residualSumOfSquares = 0;
		numberOfFittedPoints = n;
		iterations = 0;
		evaluations = 0;
		value = 0;
		FitStatus status = computeFit(n, y, y_fit, a, a_dev, error, noise);
		// Compute the total sum of squares if a good fit
		if (status == FitStatus.OK)
			totalSumOfSquares = getSumOfSquares(n, y);
		return status;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#evaluate(int, double[], double[], double[])
	 */
	public boolean evaluate(int n, double[] y, double[] y_fit, double[] a)
	{
		// Reset the results
		residualSumOfSquares = 0;
		numberOfFittedPoints = n;
		iterations = 0;
		evaluations = 0;
		value = 0;
		boolean status = computeValue(n, y, y_fit, a);
		// Compute the total sum of squares if a good fit
		if (status)
			totalSumOfSquares = getSumOfSquares(n, y);
		return status;
	}

	/**
	 * Compute fit.
	 *
	 * @param n
	 *            the n
	 * @param y
	 *            the y
	 * @param y_fit
	 *            the y_fit
	 * @param a
	 *            the a
	 * @param a_dev
	 *            the a_dev
	 * @param error
	 *            the error
	 * @param noise
	 *            the noise
	 * @return the fit status
	 */
	public abstract FitStatus computeFit(int n, double[] y, double[] y_fit, double[] a, double[] a_dev, double[] error,
			double noise);

	/**
	 * Evaluate the function.
	 *
	 * @param n
	 *            the n
	 * @param y
	 *            the y
	 * @param y_fit
	 *            the y_fit
	 * @param a
	 *            the a
	 * @return true if evaluated
	 */
	public abstract boolean computeValue(int n, double[] y, double[] y_fit, double[] a);

	public double[] getInitialSolution(double[] params)
	{
		final int[] indices = f.gradientIndices();
		final double[] initialSolution = new double[indices.length];
		for (int i = 0; i < indices.length; i++)
			initialSolution[i] = params[indices[i]];
		return initialSolution;
	}

	public void setSolution(double[] params, double[] solution)
	{
		final int[] indices = f.gradientIndices();
		for (int i = 0; i < indices.length; i++)
			params[indices[i]] = solution[i];
	}

	public void setDeviations(double[] deviations, double[][] covar)
	{
		final int[] indices = f.gradientIndices();
		for (int i = 0; i < indices.length; i++)
			deviations[indices[i]] = Math.sqrt(covar[i][i]);
	}

	public static double getSumOfSquares(final int n, double[] y)
	{
		double sx = 0, ssx = 0;
		for (int i = n; i-- > 0;)
		{
			sx += y[i];
			ssx += y[i] * y[i];
		}
		final double sumOfSquares = ssx - (sx * sx) / (n);
		return sumOfSquares;
	}

	/**
	 * Compute the error
	 * 
	 * @param residualSumOfSquares
	 * @param noise
	 * @param numberOfFittedPoints
	 * @param numberOfFittedParameters
	 * @return the error
	 */
	public static double getError(double residualSumOfSquares, double noise, int numberOfFittedPoints,
			int numberOfFittedParameters)
	{
		double error = residualSumOfSquares;

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

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#getTotalSumOfSquares()
	 */
	public double getTotalSumOfSquares()
	{
		return totalSumOfSquares;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#getNumberOfFittedParameters()
	 */
	public int getNumberOfFittedParameters()
	{
		return f.gradientIndices().length;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#getNumberOfFittedPoints()
	 */
	public int getNumberOfFittedPoints()
	{
		return numberOfFittedPoints;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#getResidualSumOfSquares()
	 */
	public double getResidualSumOfSquares()
	{
		return residualSumOfSquares;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#getIterations()
	 */
	public int getIterations()
	{
		return iterations;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#getEvaluations()
	 */
	public int getEvaluations()
	{
		return evaluations;
	}

	/**
	 * @return the maxEvaluations
	 */
	public int getMaxEvaluations()
	{
		return maxEvaluations;
	}

	/**
	 * @param maxEvaluations
	 *            the maxEvaluations to set
	 */
	public void setMaxEvaluations(int maxEvaluations)
	{
		this.maxEvaluations = maxEvaluations;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#isBounded()
	 */
	public boolean isBounded()
	{
		return false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#isConstrained()
	 */
	public boolean isConstrained()
	{
		return false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#setBounds(double[], double[])
	 */
	public void setBounds(double[] lower, double[] upper)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#setConstraints(double[], double[])
	 */
	public void setConstraints(double[] lower, double[] upper)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#getValue()
	 */
	public double getValue()
	{
		return value;
	}

	/**
	 * Update the function.
	 *
	 * @param f
	 *            the new function
	 * @throws NullPointerException
	 *             if the function is null
	 */
	public void setNonLinearFunction(NonLinearFunction f)
	{
		if (f == null)
			throw new NullPointerException("Function must not be null");
		this.f = f;
	}

	/**
	 * @return The function
	 */
	public NonLinearFunction getNonLinearFunction()
	{
		return f;
	}

	/**
	 * Update the total and residual sum-of-squares using the given data
	 * 
	 * @param n
	 *            The number of data points
	 * @param y
	 *            The data
	 * @param residuals
	 *            The fit residuals
	 */
	public void updateSumOfSquares(int n, double[] y, double[] residuals)
	{
		totalSumOfSquares = getSumOfSquares(n, y);
		residualSumOfSquares = 0;
		for (int i = 0; i < n; i++)
		{
			residualSumOfSquares += residuals[i] * residuals[i];
		}
	}
}
