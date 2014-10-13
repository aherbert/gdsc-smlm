package gdsc.smlm.fitting.nonlinear;

import java.util.Arrays;

import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.function.NonLinearFunction;
import gdsc.smlm.fitting.function.PoissonLikelihoodFunction;

import org.apache.commons.math3.analysis.FunctionUtils;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.differentiation.MultivariateDifferentiableFunction;
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
import org.apache.commons.math3.optimization.direct.MultivariateFunctionMappingAdapter;
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
	// Maximum iterations for the Powell optimiser
	private int maxIterations;

	public enum SearchMethod
	{
		/**
		 * Search using Powell's conjugate direction method
		 */
		POWELL,
		/**
		 * Search using Powell's conjugate direction method using a mapping adapter to ensure a bounded search
		 */
		POWELL_BOUNDED,
		/**
		 * Search using Powell's Bound Optimization BY Quadratic Approximation (BOBYQA) algorithm.
		 * <p>
		 * BOBYQA could also be considered as a replacement of any derivative-based optimizer when the derivatives are
		 * approximated by finite differences. This is a bounded search.
		 */
		BOBYQA,
		/**
		 * Search using active Covariance Matrix Adaptation Evolution Strategy (CMA-ES).
		 * <p>
		 * The CMA-ES is a reliable stochastic optimization method which should be applied if derivative-based methods,
		 * e.g. conjugate gradient, fail due to a rugged search landscape. This is a bounded search.
		 */
		CMAES,
		/**
		 * Search using a non-linear conjugate gradient optimiser. Use the Fletcher-Reeves update formulas for the
		 * conjugate search directions.
		 */
		CONJUGATE_GRADIENT_FR,
		/**
		 * Search using a non-linear conjugate gradient optimiser. Use the Polak-Ribière update formulas for the
		 * conjugate search directions.
		 */
		CONJUGATE_GRADIENT_PR
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

			double ll;
			double[] solution;
			if (searchMethod == SearchMethod.POWELL || searchMethod == SearchMethod.POWELL_BOUNDED)
			{
				// Non-differentiable version using Powell Optimiser

				// This is as per the method in Numerical Recipes 10.5 (Direction Set (Powell's) method)
				// I could extend the optimiser and implement bounds on the directions moved. However the mapping
				// adapter seems to work OK.

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

				// Try using the mapping adapter
				MultivariateFunction fun = new PoissonLikelihoodFunction(f, a, y, n);
				if (searchMethod == SearchMethod.POWELL)
				{
					PointValuePair optimum = o.optimize(getMaxEvaluations(), fun, GoalType.MINIMIZE, new InitialGuess(
							startPoint));
					solution = optimum.getPointRef();
					ll = optimum.getValue();
				}
				else
				{
					MultivariateFunctionMappingAdapter adapter = new MultivariateFunctionMappingAdapter(fun, lower,
							upper);
					PointValuePair optimum = o.optimize(getMaxEvaluations(), adapter, GoalType.MINIMIZE,
							new InitialGuess(adapter.boundedToUnbounded(startPoint)));
					solution = adapter.unboundedToBounded(optimum.getPointRef());
					ll = optimum.getValue();
				}
			}
			else if (searchMethod == SearchMethod.BOBYQA)
			{
				// Differentiable approximation using Powell's BOBYQA algorithm.
				// This is slower than the Powell optimiser
				int numberOfInterpolationPoints = this.getNumberOfFittedParameters() + 2;

				BOBYQAOptimizer o = new BOBYQAOptimizer(numberOfInterpolationPoints);
				optimizer = o;
				PointValuePair optimum = o.optimize(getMaxEvaluations(), new PoissonLikelihoodFunction(f, a, y, n),
						GoalType.MINIMIZE, new InitialGuess(startPoint), new SimpleBounds(lower, upper));
				solution = optimum.getPointRef();
				ll = optimum.getValue();
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
				PointValuePair optimum = o.optimize(getMaxEvaluations(), new PoissonLikelihoodFunction(f, a, y, n),
						GoalType.MINIMIZE, new InitialGuess(startPoint), new SimpleBounds(lower, upper));
				solution = optimum.getPointRef();
				ll = optimum.getValue();
			}
			else
			{
				// TODO - Check if this works. Are the gradients correct?

				// The line search algorithm often fails. This is due to searching into a region where the 
				// function evaluates to a negative so has been clipped. This means the upper bound of the line
				// cannot be found. Maybe I can write a custom line search routine that copes with the bounds.
				// Note that running it on an easy problem (200 photons with fixed fitting (no background)) the algorithm
				// does sometimes produces results better than the Powell algorithm but it is slower.

				// Rearrange to use the non-deprecated classes in org.apache.commons.math3.optim.nonlinear.scalar.gradient
				NonLinearConjugateGradientOptimizer o = new NonLinearConjugateGradientOptimizer(
						(searchMethod == SearchMethod.CONJUGATE_GRADIENT_FR) ? ConjugateGradientFormula.FLETCHER_REEVES
								: ConjugateGradientFormula.POLAK_RIBIERE);
				optimizer = o;

				// TODO - Can I apply a bounds adapter to the differentiable function. I will have to write one.
				// However the gradients may become unstable at the edge of the bounds. 
				// Or they will not change direction if the true solution is on the bounds since the gradient will
				// always continue towards the bounds. This is key to the
				// conjugate gradient method. It searches along a vector until the direction of the gradient
				// is in the opposite direction (using dot products, i.e. cosine of angle between them)

				// It may be better to write a bounded version of the conjugate gradient optimiser that truncates
				// search directions at the bounds.

				// NR 10.7 states there is no advantage of the variable metric DFP or BFGS methods over
				// conjugate gradient methods.

				// Try this:
				// Adapt the conjugate gradient optimiser to use the gradient to pick the search direction
				// and then for the line minimisation. However if the function is out of bounds then:
				// - a) just stop at that point (no line minimisation) and switch direction
				// - b) clip the variables at the bounds and continue 
				// However the gradient is invalid when the function evaluates <0.
				// Maybe just use the gradient for the search direction then use the line minimisation/rest
				// as per the Powell optimiser. The bounds should limit the search.

				MultivariateDifferentiableFunction fun = FunctionUtils
						.toMultivariateDifferentiableFunction(new PoissonLikelihoodFunction(f, a, y, n));
				o.setInitialStep(0.1);
				PointValuePair optimum = o.optimize(getMaxEvaluations(), fun, GoalType.MINIMIZE, startPoint);
				solution = optimum.getPointRef();
				ll = optimum.getValue();
			}

			setSolution(a, solution);
			iterations = optimizer.getEvaluations();

			System.out.printf("Iter = %d, %g @ %s\n", iterations, ll, Arrays.toString(solution));

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
			//System.out.printf("Too many iterations = %d\n", e.getMax());
			//e.printStackTrace();
			return FitStatus.FAILED_TO_CONVERGE;
		}
		catch (TooManyEvaluationsException e)
		{
			//System.out.printf("Too many evaluations = %d\n", e.getMax());
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
			//e.printStackTrace();

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
		return searchMethod == SearchMethod.POWELL_BOUNDED || searchMethod == SearchMethod.BOBYQA ||
				searchMethod == SearchMethod.CMAES;
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
