package gdsc.smlm.fitting.nonlinear;

import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.function.Gaussian2DFunction;
import gdsc.smlm.fitting.function.NonLinearFunction;

import java.util.Arrays;

import org.apache.commons.math3.analysis.DifferentiableMultivariateFunction;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.exception.TooManyIterationsException;
import org.apache.commons.math3.optimization.BaseMultivariateOptimizer;
import org.apache.commons.math3.optimization.ConvergenceChecker;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.InitialGuess;
import org.apache.commons.math3.optimization.PointValuePair;
import org.apache.commons.math3.optimization.SimpleBounds;
import org.apache.commons.math3.optimization.SimpleValueChecker;
import org.apache.commons.math3.optimization.direct.BOBYQAOptimizer;
import org.apache.commons.math3.optimization.direct.CMAESOptimizer;
import org.apache.commons.math3.optimization.direct.PowellOptimizer;
import org.apache.commons.math3.optimization.general.ConjugateGradientFormula;
import org.apache.commons.math3.optimization.general.NonLinearConjugateGradientOptimizer;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

import com.sun.tools.javac.comp.Todo;

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
 * Uses Maximum Likelihood Estimation (MLE) to fit a nonlinear model with coefficients (a) for a
 * set of data points (x, y).
 * <p>
 * The probability mass function for observed value k is modelled as a Poisson process:<br/>
 * pmf = e^-k.(l^k / k!) <br/>
 * where: <br/>
 * k = Observed number of occurrences <br/>
 * l = Expected number of occurrences (the mean)
 * <p>
 * MLE = Max [ sum (ln(e^-k.(l^k / k!)) ] <br/>
 * = Max [ sum (k.ln(l) - l) ]
 * <p>
 * The expected number of occurrences can be modelled using any parameterised function, for example the Gaussian 2D
 * function.
 */
public class MaximumLikelihoodFitter extends BaseFunctionSolver
{
	/**
	 * Wrapper for any function to allow use of the Apache Commons optimiser for Maximum Likelihood Estimation.
	 * <p>
	 * Uses the deprecated API since the new API for version 4.0 is not a fully documented final release.
	 */
	public class PoissonLikelihoodFunction implements DifferentiableMultivariateFunction
	{
		private NonLinearFunction f;
		private float[] a, y;
		int n;

		/**
		 * @param f
		 *            The function to be used to calculated the expected values
		 * @param a
		 *            The initial parameters for the gaussian function
		 * @param y
		 *            The observed values
		 * @param n
		 *            The number of observed values
		 */
		public PoissonLikelihoodFunction(NonLinearFunction f, float[] a, float[] y, int n)
		{
			this.f = f;
			this.a = Arrays.copyOf(a, a.length);
			this.y = y;
			this.n = n;
		}

		/**
		 * Copy the variables into the appropriate parameter positions for the NonLinearFunction
		 * 
		 * @param variables
		 */
		private void initialiseFunction(double[] variables)
		{
			int[] gradientIndices = f.gradientIndices();
			for (int i = 0; i < gradientIndices.length; i++)
				a[gradientIndices[i]] = (float) variables[i];
			// Do not allow negative values to be computed (since log(-x) is undefined)
			if (a[Gaussian2DFunction.AMPLITUDE] <= 0)
				a[Gaussian2DFunction.AMPLITUDE] = 0.1f;
			if (a[Gaussian2DFunction.BACKGROUND] < 0)
				a[Gaussian2DFunction.BACKGROUND] = 0;

			f.initialise(a);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see org.apache.commons.math3.analysis.MultivariateFunction#value(double[])
		 */
		public double value(double[] variables)
		{
			initialiseFunction(variables);

			// Compute the negative log-likelihood to be minimised
			double ll = 0;
			for (int i = 0; i < n; i++)
			{
				final double l = f.eval(i);
				final double k = y[i];
				ll += l - k * Math.log(l);
			}
			return ll;
		}

		@Override
		public MultivariateFunction partialDerivative(int k)
		{
			// This is not required for the AbstractDifferentiableOptimizer classes
			return null;
		}

		@Override
		public MultivariateVectorFunction gradient()
		{
			// TODO - This is required
			return null;
		}
	}

	// Maximum iterations for the Powell optimiser
	private int maxIterations;

	public enum SearchMethod
	{
		/**
		 * Search using Powell's conjugate direction method
		 */
		POWELL,
		/**
		 * Search using Powell's Bound Optimization BY Quadratic Approximation (BOBYQA) algorithm.
		 * <p>
		 * BOBYQA could also be considered as a replacement of any derivative-based optimizer when the derivatives are
		 * approximated by finite differences.
		 */
		BOBYQA,
		/**
		 * Search using active Covariance Matrix Adaptation Evolution Strategy (CMA-ES).
		 * <p>
		 * The CMA-ES is a reliable stochastic optimization method which should be applied if derivative-based methods,
		 * e.g. conjugate gradient, fail due to a rugged search landscape
		 */
		CMAES
	}

	private SearchMethod searchMethod;
	private double[] lower, upper;

	/**
	 * Default constructor
	 */
	public MaximumLikelihoodFitter(NonLinearFunction gf)
	{
		super(gf);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#fit(int, float[], float[], float[], float[], double[], double)
	 */
	@SuppressWarnings("deprecation")
	public FitStatus fit(int n, float[] y, float[] y_fit, float[] a, float[] a_dev, double[] error, double noise)
	{
		numberOfFittedPoints = n;

		BaseMultivariateOptimizer optimizer = null;

		try
		{
			double[] startPoint = getInitialSolution(a);

			PointValuePair optimum;
			if (searchMethod == SearchMethod.POWELL)
			{
				// Non-differentiable version using Powell Optimiser
				double rel = 1e-4; // Relative threshold.
				double abs = 1e-10; // Absolute threshold.
				PowellOptimizer o = new PowellOptimizer(rel, abs, new ConvergenceChecker<PointValuePair>()
				{
					@Override
					public boolean converged(int iteration, PointValuePair previous, PointValuePair current)
					{
						if (getMaxIterations() > 0 && iteration > getMaxIterations())
							throw new TooManyIterationsException(getMaxIterations());
						return false;
					}
				});

				optimizer = o;
				optimum = o.optimize(getMaxEvaluations(), new PoissonLikelihoodFunction(f, a, y, n), GoalType.MINIMIZE,
						new InitialGuess(startPoint));
			}
			else if (searchMethod == SearchMethod.BOBYQA)
			{
				// Differentiable approximation using Powell's BOBYQA algorithm.
				// This is slower than the Powell optimiser
				int numberOfInterpolationPoints = this.getNumberOfFittedParameters() + 2;

				BOBYQAOptimizer o = new BOBYQAOptimizer(numberOfInterpolationPoints);
				optimizer = o;
				optimum = o.optimize(getMaxEvaluations(), new PoissonLikelihoodFunction(f, a, y, n), GoalType.MINIMIZE,
						new InitialGuess(startPoint), new SimpleBounds(lower, upper));
			}
			else if (searchMethod == SearchMethod.CMAES)
			{
				// CMAESOptimiser based on Matlab code:
				// https://www.lri.fr/~hansen/cmaes.m
				// Take the defaults from the Matlab documentation
				double stopFitness = 0; //Double.NEGATIVE_INFINITY;
				boolean isActiveCMA = true;
				int diagonalOnly = 0;
				int checkFeasableCount = 1;
				RandomGenerator random = new Well19937c();
				boolean generateStatistics = false;
				// The sigma determines the search range for the variables. It should be 1/3 of the initial search region.
				double[] sigma = new double[lower.length];
				for (int i = 0; i < sigma.length; i++)
					sigma[i] = (upper[i] - lower[i]) / 3;
				int popSize = (int) (4 + Math.floor(3 * Math.log(n)));

				CMAESOptimizer o = new CMAESOptimizer(popSize, sigma, getMaxIterations(), stopFitness, isActiveCMA,
						diagonalOnly, checkFeasableCount, random, generateStatistics, new SimpleValueChecker(1e-6,
								1e-10));
				optimizer = o;
				optimum = o.optimize(getMaxEvaluations(), new PoissonLikelihoodFunction(f, a, y, n), GoalType.MINIMIZE,
						new InitialGuess(startPoint), new SimpleBounds(lower, upper));
			}
			else
			{
				// TODO - differentiable version using NonLinearConjugateGradientOptimiser
				// Need to differentiate the log-likelihood function

				NonLinearConjugateGradientOptimizer o = new NonLinearConjugateGradientOptimizer(
						ConjugateGradientFormula.POLAK_RIBIERE);
				optimizer = o;
				optimum = o.optimize(getMaxEvaluations(), new PoissonLikelihoodFunction(f, a, y, n), GoalType.MINIMIZE,
						new InitialGuess(startPoint));
			}

			setSolution(a, optimum.getPointRef());
			iterations = optimizer.getEvaluations();

			//System.out.printf("Iter = %d, %g @ %s\n", iterations, optimum.getValue(),
			//		Arrays.toString(optimum.getPoint()));

			// Compute residuals for the FunctionSolver interface
			if (y_fit == null || y_fit.length < n)
				y_fit = new float[n];
			f.initialise(a);
			residualSumOfSquares = 0;
			for (int i = 0; i < n; i++)
			{
				y_fit[i] = f.eval(i);
				final double residual = y[i] - y_fit[i];
				residualSumOfSquares += residual * residual;
			}

			if (a_dev != null)
			{
				// Covariance not supported for this fitter
			}

			error[0] = NonLinearFit.getError(residualSumOfSquares, noise, n, f.gradientIndices().length);
			totalSumOfSquares = getSumOfSquares(n, y);
		}
		catch (TooManyIterationsException e)
		{
			System.out.printf("Iterations = %d\n", e.getMax());
			//e.printStackTrace();
			return FitStatus.FAILED_TO_CONVERGE;
		}
		catch (TooManyEvaluationsException e)
		{
			System.out.printf("Evaluations = %d\n", optimizer.getEvaluations());
			//e.printStackTrace();
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

	/**
	 * @return the max iterations for the Powell search method
	 */
	public int getMaxIterations()
	{
		return maxIterations;
	}

	/**
	 * @param maxIterations
	 *            the max iterations for the Powell search method
	 */
	public void setMaxIterations(int maxIterations)
	{
		this.maxIterations = maxIterations;
	}

	/**
	 * @return the search method
	 */
	public SearchMethod getSearchMethod()
	{
		return searchMethod;
	}

	/**
	 * @param searchMethod
	 *            the search method
	 */
	public void setSearchMethod(SearchMethod searchMethod)
	{
		this.searchMethod = searchMethod;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#isBounded()
	 */
	@Override
	public boolean isBounded()
	{
		return searchMethod == SearchMethod.BOBYQA || searchMethod == SearchMethod.CMAES;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#setBounds(float[], float[])
	 */
	@Override
	public void setBounds(float[] lowerB, float[] upperB)
	{
		// Extract the bounds for the parameters we are fitting
		int[] indices = f.gradientIndices();

		lower = new double[indices.length];
		upper = new double[indices.length];
		for (int i = 0; i < indices.length; i++)
		{
			lower[i] = lowerB[indices[i]];
			upper[i] = upperB[indices[i]];
		}
	}
}
