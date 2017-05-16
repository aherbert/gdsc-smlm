package gdsc.smlm.fitting.nonlinear;

import java.util.Arrays;

import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolver;
import gdsc.smlm.fitting.FunctionSolverType;
import gdsc.smlm.function.GradientFunction;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;

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
	protected FunctionSolverType type;
	protected GradientFunction f;

	private int maxEvaluations = 20;

	protected int numberOfFittedPoints;
	protected int iterations;
	protected int evaluations;
	protected double value;

	// Cache the data to fit on success
	protected double[] lastY, lastA;

	/**
	 * Default constructor.
	 *
	 * @param type
	 *            the type
	 * @param f
	 *            the f
	 * @throws NullPointerException
	 *             if the function is null
	 */
	public BaseFunctionSolver(FunctionSolverType type, GradientFunction f)
	{
		setGradientFunction(f);
		this.type = type;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#getType()
	 */
	public FunctionSolverType getType()
	{
		return type;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#fit(double[], double[], double[], double[])
	 */
	public FitStatus fit(double[] y, double[] y_fit, double[] a, double[] a_dev)
	{
		// Reset the results
		numberOfFittedPoints = y.length;
		iterations = 0;
		evaluations = 0;
		value = 0;
		lastY = null;
		lastA = null;
		preProcess();
		FitStatus status = computeFit(y, y_fit, a, a_dev);
		if (status == FitStatus.OK)
		{
			if (lastY == null)
				lastY = y;
			if (lastA == null)
				lastA = a;
			postProcess();
		}
		return status;

	}

	/**
	 * Run before fit/evaluate
	 */
	protected void preProcess()
	{
		
	}
	
	/**
	 * Run if the fit/evaluate was successful
	 */
	protected void postProcess()
	{
		
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#evaluate(double[], double[], double[])
	 */
	public boolean evaluate(double[] y, double[] y_fit, double[] a)
	{
		// Reset the results
		numberOfFittedPoints = y.length;
		iterations = 0;
		evaluations = 0;
		value = 0;
		lastY = null;
		lastA = null;
		preProcess();
		boolean status = computeValue(y, y_fit, a);
		if (status)
		{
			if (lastY == null)
				lastY = y;
			if (lastA == null)
				lastA = a;
			postProcess();
		}
		return status;
	}

	/**
	 * Compute fit.
	 *
	 * @param y
	 *            the y values to fit
	 * @param y_fit
	 *            the final fitted y values
	 * @param a
	 *            the parameters a
	 * @param a_dev
	 *            the deviations for the parameters a
	 * @return the fit status
	 */
	public abstract FitStatus computeFit(double[] y, double[] y_fit, double[] a, double[] a_dev);

	/**
	 * Evaluate the function.
	 *
	 * @param y
	 *            the y values to fit
	 * @param y_fit
	 *            the final fitted y values
	 * @param a
	 *            the parameters a
	 * @return true if evaluated
	 */
	public abstract boolean computeValue(double[] y, double[] y_fit, double[] a);

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
		Arrays.fill(deviations, 0);
		final int[] indices = f.gradientIndices();
		for (int i = 0; i < indices.length; i++)
			deviations[indices[i]] = Math.sqrt(covar[i][i]);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#getNumberOfFittedParameters()
	 */
	public int getNumberOfFittedParameters()
	{
		return f.getNumberOfGradients();
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
	public void setGradientFunction(GradientFunction f)
	{
		if (f == null)
			throw new NullPointerException("Function must not be null");
		this.f = f;
	}

	/**
	 * @return The function
	 */
	public GradientFunction getGradientFunction()
	{
		return f;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#getName(int)
	 */
	public String getName(int i)
	{
		if (f instanceof Gaussian2DFunction)
		{
			((Gaussian2DFunction) f).getName(i);
		}
		return "Unknown";
	}
}
