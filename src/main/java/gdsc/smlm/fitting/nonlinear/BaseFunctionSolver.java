package gdsc.smlm.fitting.nonlinear;

import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolver;
import gdsc.smlm.fitting.function.NonLinearFunction;

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

	/**
	 * Default constructor
	 */
	public BaseFunctionSolver(NonLinearFunction f)
	{
		this.f = f;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#fit(int, double[], double[], double[], double[], double[], double)
	 */
	public abstract FitStatus fit(int n, double[] y, double[] y_fit, double[] a, double[] a_dev, double[] error,
			double noise);

	public double[] getInitialSolution(double[] params)
	{
		int[] indices = f.gradientIndices();
		double[] initialSolution = new double[indices.length];
		for (int i = 0; i < indices.length; i++)
			initialSolution[i] = params[indices[i]];
		return initialSolution;
	}

	public void setSolution(double[] params, double[] solution)
	{
		int[] indices = f.gradientIndices();
		for (int i = 0; i < indices.length; i++)
			params[indices[i]] = (double) solution[i];
	}

	public void setDeviations(double[] deviations, double[][] covar)
	{
		int[] indices = f.gradientIndices();
		for (int i = 0; i < indices.length; i++)
			deviations[indices[i]] = (double) Math.sqrt(covar[i][i]);
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
	@Override
	public double getTotalSumOfSquares()
	{
		return totalSumOfSquares;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#getNumberOfFittedParameters()
	 */
	@Override
	public int getNumberOfFittedParameters()
	{
		return f.gradientIndices().length;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#getNumberOfFittedPoints()
	 */
	@Override
	public int getNumberOfFittedPoints()
	{
		return numberOfFittedPoints;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#getResidualSumOfSquares()
	 */
	@Override
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
	@Override
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
	@Override
	public boolean isBounded()
	{
		return false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#isConstrained()
	 */
	@Override
	public boolean isConstrained()
	{
		return false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#setBounds(double[], double[])
	 */
	@Override
	public void setBounds(double[] lower, double[] upper)
	{
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#setConstraints(double[], double[])
	 */
	@Override
	public void setConstraints(double[] lower, double[] upper)
	{
	}
}
