package gdsc.smlm.fitting.nonlinear;

import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.FunctionSolver;
import gdsc.smlm.fitting.function.DifferentiableGaussian2DFunction;
import gdsc.smlm.fitting.function.Gaussian2DFunction;

import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optimization.PointVectorValuePair;
import org.apache.commons.math3.optimization.general.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.util.Precision;

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
 * Uses Apache Commons Math Levenberg-Marquardt method to fit a nonlinear model with coefficients (a) for a
 * set of data points (x, y).
 */
public class ApacheNonLinearFit implements FunctionSolver
{
	private Gaussian2DFunction gf;

	private int maxEvaluations = 20;
	
	private double totalSumOfSquares;
	private double residualSumOfSquares;
	private int numberOfFittedPoints;
	private int iterations;

	/**
	 * Default constructor
	 */
	public ApacheNonLinearFit(Gaussian2DFunction gf)
	{
		this.gf = gf;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#fit(int, float[], float[], float[], float[], double[], double)
	 */
	public FitStatus fit(int n, float[] y, float[] y_fit, float[] a, float[] a_dev, double[] error, double noise)
	{
		numberOfFittedPoints = n;

		DifferentiableGaussian2DFunction func = new DifferentiableGaussian2DFunction(gf, a);
		func.addData(y);

		try
		{
			// Different convergence thresholds seem to have no effect on the resulting fit, only the number of
			// iterations for convergence
			double initialStepBoundFactor = 100;
			double costRelativeTolerance = 1e-10;
			double parRelativeTolerance = 1e-10;
			double orthoTolerance = 1e-10;
			double threshold = Precision.SAFE_MIN;

			// Extract the parameters to be fitted
			int[] gradientIndices = gf.gradientIndices();
			double[] initialSolution = new double[gradientIndices.length];
			for (int i = 0; i < gradientIndices.length; i++)
				initialSolution[i] = a[gradientIndices[i]];

			// TODO - Pass in stopping criteria.
			// Do this if the fitting seems to be successful
			
			LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer(initialStepBoundFactor,
					costRelativeTolerance, parRelativeTolerance, orthoTolerance, threshold);
			PointVectorValuePair optimum = optimizer.optimize(maxEvaluations, func, func.getY(), func.getWeights(),
					initialSolution);

			double[] parameters = optimum.getPoint();
			for (int i = 0; i < gradientIndices.length; i++)
				a[gradientIndices[i]] = (float) parameters[i];
			iterations = optimizer.getEvaluations();
			error[0] = optimizer.getChiSquare();
			if (a_dev != null)
			{
				double[][] covar = optimizer.getCovariances();
				for (int i = 0; i < gradientIndices.length; i++)
					a_dev[gradientIndices[i]] = (float) covar[i][i];
			}
			residualSumOfSquares = optimizer.getRMS() * optimizer.getRMS() * n;
			totalSumOfSquares = getSumOfSquares(n, y);
		}
		catch (TooManyEvaluationsException e)
		{
			return FitStatus.FAILED_TO_CONVERGE;
		}
		catch (ConvergenceException e)
		{
			// Occurs when QR decomposition fails - mark as a singular non-linear model (no solution)
			return FitStatus.SINGULAR_NON_LINEAR_MODEL;
		}
		catch (Exception e)
		{
			// TODO - Find out the other exceptions from the fitter and add return values to match. 
			return FitStatus.UNKNOWN;
		}

		return FitStatus.OK;
	}

	private double getSumOfSquares(final int n, float[] y)
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
		return gf.gradientIndices().length;
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
	 * @param maxEvaluations the maxEvaluations to set
	 */
	public void setMaxEvaluations(int maxEvaluations)
	{
		this.maxEvaluations = maxEvaluations;
	}
}
