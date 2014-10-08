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
 * 
 * This is an adaption of the C-code contained in the CcpNmr Analysis Program:
 *   CCPN website (http://www.ccpn.ac.uk/). 
 * The CCPN code was based on Numerical Recipes. 
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
	 * @see gdsc.smlm.fitting.FunctionSolver#fit(int, float[], float[], float[], float[], double[], double)
	 */
	public abstract FitStatus fit(int n, float[] y, float[] y_fit, float[] a, float[] a_dev, double[] error,
			double noise);

	public double[] getInitialSolution(float[] params)
	{
		int[] indices = f.gradientIndices();
		double[] initialSolution = new double[indices.length];
		for (int i = 0; i < indices.length; i++)
			initialSolution[i] = params[indices[i]];
		return initialSolution;
	}

	public void setSolution(float[] params, double[] solution)
	{
		int[] indices = f.gradientIndices();
		for (int i = 0; i < indices.length; i++)
			params[indices[i]] = (float) solution[i];
	}

	public void setDeviations(float[] deviations, double[][] covar)
	{
		int[] indices = f.gradientIndices();
		for (int i = 0; i < indices.length; i++)
			deviations[indices[i]] = (float) covar[i][i];
	}

	public static double getSumOfSquares(final int n, float[] y)
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
}
