/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm.fitting.nonlinear;

import java.util.Arrays;

import gdsc.smlm.fitting.FisherInformationMatrix;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolver;
import gdsc.smlm.fitting.FunctionSolverType;
import gdsc.smlm.function.GradientFunction;
import gdsc.smlm.function.NamedFunction;

/**
 * Abstract class with utility methods for the FunctionSolver interface.
 */
public abstract class BaseFunctionSolver implements FunctionSolver
{
	private FunctionSolverType type;
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
	 *             if the function or type is null
	 */
	public BaseFunctionSolver(FunctionSolverType type, GradientFunction f)
	{
		if (f == null)
			throw new NullPointerException("Function must not be null");
		this.f = f;
		setType(type);
	}

	/**
	 * Sets the type.
	 *
	 * @param type
	 *            the new type
	 */
	protected void setType(FunctionSolverType type)
	{
		if (type == null)
			throw new NullPointerException("Type must not be null");
		this.type = type;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.FunctionSolver#getType()
	 */
	@Override
	public FunctionSolverType getType()
	{
		return type;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.FunctionSolver#fit(double[], double[], double[], double[])
	 */
	@Override
	public FitStatus fit(double[] y, double[] yFit, double[] a, double[] aDev)
	{
		// Reset the results
		numberOfFittedPoints = y.length;
		iterations = 0;
		evaluations = 0;
		value = 0;
		lastY = null;
		lastA = null;
		preProcess();
		FitStatus status = computeFit(y, yFit, a, aDev);
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
	@Override
	public boolean evaluate(double[] y, double[] yFit, double[] a)
	{
		// Reset the results
		numberOfFittedPoints = y.length;
		iterations = 0;
		evaluations = 0;
		value = 0;
		lastY = null;
		lastA = null;
		preProcess();
		boolean status = computeValue(y, yFit, a);
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

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.FunctionSolver#computeDeviations(double[], double[], double[])
	 */
	@Override
	public boolean computeDeviations(double[] y, double[] a, double[] aDev)
	{
		// Use a dedicated solver optimised for inverting the matrix diagonal.
		final FisherInformationMatrix m = computeFisherInformationMatrix(y, a);

		setDeviations(aDev, m);

		return true;
	}

	/**
	 * Compute fit.
	 *
	 * @param y
	 *            the y values to fit
	 * @param yFit
	 *            the final fitted y values
	 * @param a
	 *            the parameters a
	 * @param aDev
	 *            the deviations for the parameters a
	 * @return the fit status
	 */
	protected abstract FitStatus computeFit(double[] y, double[] yFit, double[] a, double[] aDev);

	/**
	 * Evaluate the function.
	 *
	 * @param y
	 *            the y values to fit
	 * @param yFit
	 *            the final fitted y values
	 * @param a
	 *            the parameters a
	 * @return true if evaluated
	 */
	protected abstract boolean computeValue(double[] y, double[] yFit, double[] a);

	/**
	 * Compute the Fisher Information matrix. This can be used to set the deviations for each of the fitted parameters.
	 * <p>
	 * Alternatively a sub-class can override {@link #computeDeviations(double[], double[], double[])} directly and
	 * provide a dummy implementation of this function as it will not be used, e.g. throw an exception.
	 *
	 * @param y
	 *            the y values
	 * @param a
	 *            the parameters
	 * @return the Fisher Information matrix
	 */
	protected abstract FisherInformationMatrix computeFisherInformationMatrix(double[] y, double[] a);

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

	public void setDeviationsFromMatrix(double[] deviations, double[][] covar)
	{
		Arrays.fill(deviations, 0);
		final int[] indices = f.gradientIndices();
		for (int i = 0; i < indices.length; i++)
			deviations[indices[i]] = checkVariance(covar[i][i]);
	}

	public void setDeviationsFromLinearMatrix(double[] deviations, double[] covar)
	{
		Arrays.fill(deviations, 0);
		final int[] indices = f.gradientIndices();
		final int n = indices.length;
		for (int i = 0, j = 0; i < n; i++, j += (n + 1))
			deviations[indices[i]] = checkVariance(covar[j]);
	}

	public void setDeviations(double[] deviations, double[] covar)
	{
		Arrays.fill(deviations, 0);
		final int[] indices = f.gradientIndices();
		final int n = indices.length;
		for (int i = 0; i < n; i++)
			deviations[indices[i]] = checkVariance(covar[i]);
	}

	public void setDeviations(double[] deviations, FisherInformationMatrix m)
	{
		//// This may fail if the matrix cannot be inverted
		//final double[] crlb = m.crlb();
		//if (crlb == null)
		//	throw new FunctionSolverException(FitStatus.SINGULAR_NON_LINEAR_SOLUTION);
		//setDeviations(aDev, crlb);

		// Use this method for robustness, i.e. it will not fail
		setDeviations(deviations, m.crlb(true));
	}

	private static double checkVariance(double d)
	{
		return (d > 0) ? d : 0;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.FunctionSolver#getNumberOfFittedParameters()
	 */
	@Override
	public int getNumberOfFittedParameters()
	{
		return f.getNumberOfGradients();
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
	 * @see gdsc.smlm.fitting.FunctionSolver#getIterations()
	 */
	@Override
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
	 * @see gdsc.smlm.fitting.FunctionSolver#isWeighted()
	 */
	@Override
	public boolean isWeighted()
	{
		return false;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.FunctionSolver#isStrictlyPositiveFunction()
	 */
	@Override
	public boolean isStrictlyPositiveFunction()
	{
		// Provide a default implementation based on the type.
		// This can be overridden if the solver can handle negative function data.
		switch (type)
		{
			case LSE:
				// Assume least-square estimation can handle any function data
				return false;
			case MLE:
				// Assume maximum likelihood estimation requires a positive function value
				// for computing the likelihood
				return true;
			case WLSE:
				// Assume least-square estimation can handle any function data
				return false;
			default:
				// Leave this as disabled by default
				return false;
		}
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

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.FunctionSolver#setWeights(double[])
	 */
	@Override
	public void setWeights(double[] weights)
	{
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.fitting.FunctionSolver#getValue()
	 */
	@Override
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
	@Override
	public String getName(int i)
	{
		if (f instanceof NamedFunction)
		{
			return ((NamedFunction) f).getParameterName(i);
		}
		return "Unknown";
	}

	/**
	 * Ensure positive values. If values are negative a copy is made with those values set to zero.
	 *
	 * @param y
	 *            the y
	 * @return the positive y values
	 */
	public static double[] ensurePositive(double[] y)
	{
		return ensurePositive(y.length, y);
	}

	/**
	 * Ensure positive values. If values are negative a copy is made with those values set to zero.
	 *
	 * @param n
	 *            the number of values to check
	 * @param y
	 *            the y
	 * @return the positive y values
	 */
	public static double[] ensurePositive(final int n, double[] y)
	{
		for (int i = 0; i < n; i++)
		{
			if (y[i] < 0)
			{
				// Not positive so create a clone

				final double[] y2 = new double[n];
				if (i != 0)
					// Copy the values that were positive
					System.arraycopy(y, 0, y2, 0, i);

				// Note that java initialises the array to zero so only copy the positives

				//y2[i] = 0; // We know this was not positive so skip it
				while (++i < n)
				{
					//y2[i] = (y[i] < 0) ? 0 : y[i];
					if (y[i] > 0)
						y2[i] = y[i];
				}
				return y2;
			}
		}
		return y;
	}
}
