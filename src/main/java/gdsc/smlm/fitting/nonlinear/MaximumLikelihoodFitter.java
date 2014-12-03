package gdsc.smlm.fitting.nonlinear;

import java.util.Arrays;

import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator;
import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculatorFactory;
import gdsc.smlm.function.NonLinearFunction;
import gdsc.smlm.function.PoissonLikelihoodWrapper;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.exception.TooManyIterationsException;
import org.apache.commons.math3.optim.GradientChecker;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.PositionChecker;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.SimpleValueChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateFunctionMappingAdapter;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunctionGradient;
import org.apache.commons.math3.optim.nonlinear.scalar.gradient.BFGSOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.gradient.BoundedNonLinearConjugateGradientOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.gradient.BoundedNonLinearConjugateGradientOptimizer.Formula;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CustomPowellOptimizer;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.util.FastMath;

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
	 * Wrap the PoissonLikelihoodFunction with classes that implement the required interfaces
	 */
	private class PoissonWrapper
	{
		PoissonLikelihoodWrapper fun;

		public PoissonWrapper(PoissonLikelihoodWrapper fun)
		{
			this.fun = fun;
		}
	}

	private class MultivariatePoisson extends PoissonWrapper implements MultivariateFunction
	{
		public MultivariatePoisson(PoissonLikelihoodWrapper fun)
		{
			super(fun);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see org.apache.commons.math3.analysis.MultivariateFunction#value(double[])
		 */
		public double value(double[] point)
		{
			return fun.value(point);
		}

		public boolean isMapped()
		{
			return false;
		}
	}

	/**
	 * Map the specified indices using the sqrt function for use with the Powell optimiser
	 */
	private class MappedMultivariatePoisson extends MultivariatePoisson
	{
		final int[] map;

		public MappedMultivariatePoisson(PoissonLikelihoodWrapper fun, int[] map)
		{
			super(fun);
			this.map = map;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see org.apache.commons.math3.analysis.MultivariateFunction#value(double[])
		 */
		public double value(double[] point)
		{
			return fun.value(unmap(point));
		}

		/**
		 * Convert the unmapped point to the mapped equivalent. The mapped point is used by the Powell optimiser.
		 * <p>
		 * This is done by square rooting the value of the mapped indices.
		 * 
		 * @param point
		 * @return The mapped point
		 */
		public double[] map(double[] point)
		{
			point = point.clone();
			for (int i : map)
			{
				point[i] = Math.sqrt(FastMath.abs(point[i])) * FastMath.signum(point[i]);
			}
			return point;
		}

		/**
		 * Convert the mapped point to the unmapped equivalent. The unmapped point is used to evaluate the function.
		 * <p>
		 * This is done by squaring the value of the mapped indices.
		 * 
		 * @param point
		 * @return The unmapped point
		 */
		public double[] unmap(double[] point)
		{
			point = point.clone();
			for (int i : map)
			{
				//point[i] = point[i] * point[i] * FastMath.signum(point[i]);
				if (point[i] > 0)
					point[i] = point[i] * point[i];
				else
					point[i] = -(point[i] * point[i]);
			}
			return point;
		}

		public boolean isMapped()
		{
			return true;
		}
	}

	private class MultivariateVectorPoisson extends PoissonWrapper implements MultivariateVectorFunction
	{
		public MultivariateVectorPoisson(PoissonLikelihoodWrapper fun)
		{
			super(fun);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see org.apache.commons.math3.analysis.MultivariateFunction#value(double[])
		 */
		public double[] value(double[] point)
		{
			double[] gradient = new double[point.length];
			fun.value(point, gradient);
			return gradient;
		}
	}

	// Maximum iterations for the Powell optimiser
	private int maxIterations;

	public enum SearchMethod
	{
		/**
		 * Search using Powell's conjugate direction method
		 */
		POWELL("Powell"),
		/**
		 * Search using Powell's conjugate direction method using a mapping adapter to ensure a bounded search
		 */
		POWELL_BOUNDED("Powell (bounded)"),
		/**
		 * Search using Powell's Bound Optimization BY Quadratic Approximation (BOBYQA) algorithm.
		 * <p>
		 * BOBYQA could also be considered as a replacement of any derivative-based optimizer when the derivatives are
		 * approximated by finite differences. This is a bounded search.
		 */
		BOBYQA("BOBYQA"),
		/**
		 * Search using active Covariance Matrix Adaptation Evolution Strategy (CMA-ES).
		 * <p>
		 * The CMA-ES is a reliable stochastic optimization method which should be applied if derivative-based methods,
		 * e.g. conjugate gradient, fail due to a rugged search landscape. This is a bounded search.
		 */
		CMAES("CMAES"),
		/**
		 * Search using a non-linear conjugate gradient optimiser. Use the Fletcher-Reeves update formulas for the
		 * conjugate search directions.
		 * <p>
		 * This is a bounded search using simple truncation of coordinates at the bounds of the search space.
		 */
		CONJUGATE_GRADIENT_FR("Conjugate Gradient Fletcher-Reeves"),
		/**
		 * Search using a non-linear conjugate gradient optimiser. Use the Polak-Ribière update formulas for the
		 * conjugate search directions.
		 * <p>
		 * This is a bounded search using simple truncation of coordinates at the bounds of the search space.
		 */
		CONJUGATE_GRADIENT_PR("Conjugate Gradient Polak-Ribière"),
		/**
		 * Search using a Broyden-Fletcher-Goldfarb-Shanno (BFGS) gradient optimiser.
		 */
		BFGS("BFGS Gradient");

		private String name;

		private SearchMethod(String name)
		{
			this.name = name;
		}

		@Override
		public String toString()
		{
			return name;
		}
	}

	private SearchMethod searchMethod;
	private boolean gradientLineMinimisation = true;
	private double relativeThreshold = 1e-4, absoluteThreshold = 1e-10;
	private double[] lower, upper;
	private double[] lowerConstraint, upperConstraint;

	// The function to use for the Powell optimiser (which may have parameters mapped using the sqrt function) 
	private MultivariatePoisson powellFunction = null;
	private final boolean mapGaussian;

	/**
	 * Default constructor
	 * 
	 * @param f
	 *            The function
	 */
	public MaximumLikelihoodFitter(NonLinearFunction f)
	{
		super(f);
		mapGaussian = false;
	}

	/**
	 * Constructor for Gaussian2D functions. When using the Powell optimiser the background and signal parameters can be
	 * scaled using the sqrt() function. Parameters are reduced before passing to the Powell optimiser. The parameters
	 * are expanded before evaluation of the function. This allows faster exploration of the larger parameter range
	 * expected for the background and signal within the
	 * Powell optimiser.
	 * 
	 * @param f
	 *            The function
	 * @param mapGaussian
	 *            Set to true to map the background and signal parameters using sqrt() before passing to the Powell
	 *            optimiser.
	 */
	public MaximumLikelihoodFitter(Gaussian2DFunction f, boolean mapGaussian)
	{
		super(f);
		this.mapGaussian = mapGaussian;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#fit(int, double[], double[], double[], double[], double[], double)
	 */
	public FitStatus fit(int n, double[] y, double[] y_fit, double[] a, double[] a_dev, double[] error, double noise)
	{
		numberOfFittedPoints = n;

		PoissonLikelihoodWrapper maximumLikelihoodFunction = new PoissonLikelihoodWrapper(f, a, y, n);

		try
		{
			double[] startPoint = getInitialSolution(a);

			PointValuePair optimum = null;
			if (searchMethod == SearchMethod.POWELL || searchMethod == SearchMethod.POWELL_BOUNDED)
			{
				// Non-differentiable version using Powell Optimiser

				// This is as per the method in Numerical Recipes 10.5 (Direction Set (Powell's) method)
				// I could extend the optimiser and implement bounds on the directions moved. However the mapping
				// adapter seems to work OK.

				final boolean basisConvergence = false;
				CustomPowellOptimizer o = new CustomPowellOptimizer(relativeThreshold, absoluteThreshold, null,
						basisConvergence);

				OptimizationData maxIterationData = null;
				if (getMaxIterations() > 0)
					maxIterationData = new MaxIter(getMaxIterations());

				if (searchMethod == SearchMethod.POWELL)
				{
					if (powellFunction == null)
					{
						// We must map all the parameters into the same range. This is done in the Mortensen MLE 
						// Python code by using the sqrt of the number of photons and background.
						if (mapGaussian)
						{
							Gaussian2DFunction gf = (Gaussian2DFunction) f;
							// Re-map signal and background using the sqrt
							int[] indices = gf.gradientIndices();
							int[] map = new int[indices.length];
							int count = 0;
							// Background is always first
							if (indices[0] == Gaussian2DFunction.BACKGROUND)
							{
								map[count++] = 0;
							}
							// Look for the Signal in multiple peak 2D Gaussians
							for (int i = 1; i < indices.length; i++)
								if (indices[i] % 6 == Gaussian2DFunction.SIGNAL)
								{
									map[count++] = i;
								}
							if (count > 0)
							{
								powellFunction = new MappedMultivariatePoisson(maximumLikelihoodFunction,
										Arrays.copyOf(map, count));
							}
						}
						if (powellFunction == null)
						{
							powellFunction = new MultivariatePoisson(maximumLikelihoodFunction);
						}
					}

					// Update the maximum likelihood function in the Powell function wrapper
					powellFunction.fun = maximumLikelihoodFunction;

					OptimizationData positionChecker = null; //new PositionChecker(relativeThreshold, absoluteThreshold)
					if (powellFunction.isMapped())
					{
						MappedMultivariatePoisson adapter = (MappedMultivariatePoisson) powellFunction;
						optimum = o.optimize(maxIterationData, new MaxEval(getMaxEvaluations()), new ObjectiveFunction(
								powellFunction), GoalType.MINIMIZE, new InitialGuess(adapter.map(startPoint)),
								positionChecker);
						double[] solution = adapter.unmap(optimum.getPointRef());
						optimum = new PointValuePair(solution, optimum.getValue());
					}
					else
					{
						optimum = o.optimize(maxIterationData, new MaxEval(getMaxEvaluations()), new ObjectiveFunction(
								powellFunction), GoalType.MINIMIZE, new InitialGuess(startPoint), positionChecker);
					}
				}
				else
				{
					// Try using the mapping adapter for a bounded Powell search
					MultivariateFunctionMappingAdapter adapter = new MultivariateFunctionMappingAdapter(
							new MultivariatePoisson(maximumLikelihoodFunction), lower, upper);
					optimum = o.optimize(maxIterationData, new MaxEval(getMaxEvaluations()), new ObjectiveFunction(
							adapter), GoalType.MINIMIZE, new InitialGuess(adapter.boundedToUnbounded(startPoint)));
					double[] solution = adapter.unboundedToBounded(optimum.getPointRef());
					optimum = new PointValuePair(solution, optimum.getValue());
				}
				iterations = o.getIterations();
				evaluations = o.getEvaluations();
			}
			else if (searchMethod == SearchMethod.BOBYQA)
			{
				// Differentiable approximation using Powell's BOBYQA algorithm.
				// This is slower than the Powell optimiser and requires a high number of evaluations.
				int numberOfInterpolationPoints = this.getNumberOfFittedParameters() + 2;

				BOBYQAOptimizer o = new BOBYQAOptimizer(numberOfInterpolationPoints);
				optimum = o.optimize(new MaxEval(getMaxEvaluations()), new ObjectiveFunction(new MultivariatePoisson(
						maximumLikelihoodFunction)), GoalType.MINIMIZE, new InitialGuess(startPoint), new SimpleBounds(
						lower, upper));
				iterations = o.getIterations();
				evaluations = o.getEvaluations();
			}
			else if (searchMethod == SearchMethod.CMAES)
			{
				// TODO - Understand why the CMAES optimiser does not fit very well on test data. It appears 
				// to converge too early and the likelihood scores are not as low as the other optimisers.

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
				int popSize = (int) (4 + Math.floor(3 * Math.log(sigma.length)));

				// The CMAES optimiser is random and restarting can overcome problems with quick convergence.
				// The Apache commons documentations states that convergence should occur between 30N and 300N^2
				// function evaluations
				final int n30 = FastMath.min(sigma.length * sigma.length * 30, getMaxEvaluations() / 2);
				evaluations = 0;
				OptimizationData[] data = new OptimizationData[] { new InitialGuess(startPoint),
						new CMAESOptimizer.PopulationSize(popSize), new MaxEval(getMaxEvaluations()),
						new CMAESOptimizer.Sigma(sigma),
						new ObjectiveFunction(new MultivariatePoisson(maximumLikelihoodFunction)), GoalType.MINIMIZE,
						new SimpleBounds(lower, upper) };
				// Iterate to prevent early convergence
				int repeat = 0;
				while (evaluations < n30)
				{
					if (repeat++ > 1)
					{
						// Update the start point and population size
						data[0] = new InitialGuess(optimum.getPointRef());
						popSize *= 2;
						data[1] = new CMAESOptimizer.PopulationSize(popSize);
					}
					CMAESOptimizer o = new CMAESOptimizer(getMaxIterations(), stopFitness, isActiveCMA, diagonalOnly,
							checkFeasableCount, random, generateStatistics, new SimpleValueChecker(relativeThreshold,
									absoluteThreshold));
					PointValuePair result = o.optimize(data);
					iterations += o.getIterations();
					evaluations += o.getEvaluations();
					//System.out.printf("CMAES [%d] i=%d [%d], e=%d [%d]\n", repeat, o.getIterations(), iterations,
					//		o.getEvaluations(), totalEvaluations);
					if (optimum == null || result.getValue() < optimum.getValue())
					{
						optimum = result;
					}
				}
			}
			else if (searchMethod == SearchMethod.BFGS)
			{
				// BFGS can use an approximate line search minimisation where as Powell and conjugate gradient
				// methods require a more accurate line minimisation. The BFGS search does not do a full 
				// minimisation but takes appropriate steps in the direction of the current gradient.

				// Do not use the convergence checker on the value of the function. Use the convergence on the 
				// point coordinate and gradient
				//BFGSOptimizer o = new BFGSOptimizer(new SimpleValueChecker(rel, abs));
				BFGSOptimizer o = new BFGSOptimizer();

				// Configure maximum step length for each dimension using the bounds
				double[] stepLength = new double[lower.length];
				for (int i = 0; i < stepLength.length; i++)
					stepLength[i] = (upper[i] - lower[i]) * 0.3333333;

				// The GoalType is always minimise so no need to pass this in
				optimum = o.optimize(new MaxEval(getMaxEvaluations()), new ObjectiveFunctionGradient(
						new MultivariateVectorPoisson(maximumLikelihoodFunction)), new ObjectiveFunction(
						new MultivariatePoisson(maximumLikelihoodFunction)), new InitialGuess(startPoint),
						new SimpleBounds(lowerConstraint, upperConstraint), new GradientChecker(relativeThreshold,
								absoluteThreshold), new PositionChecker(relativeThreshold, absoluteThreshold),
						new BFGSOptimizer.StepLength(stepLength));
				iterations = o.getIterations();
				evaluations = o.getEvaluations();
			}
			else
			{
				// The line search algorithm often fails. This is due to searching into a region where the 
				// function evaluates to a negative so has been clipped. This means the upper bound of the line
				// cannot be found.
				// Note that running it on an easy problem (200 photons with fixed fitting (no background)) the algorithm
				// does sometimes produces results better than the Powell algorithm but it is slower.

				BoundedNonLinearConjugateGradientOptimizer o = new BoundedNonLinearConjugateGradientOptimizer(
						(searchMethod == SearchMethod.CONJUGATE_GRADIENT_FR) ? Formula.FLETCHER_REEVES
								: Formula.POLAK_RIBIERE, new SimpleValueChecker(relativeThreshold, absoluteThreshold));

				// Note: The gradients may become unstable at the edge of the bounds. Or they will not change 
				// direction if the true solution is on the bounds since the gradient will always continue 
				// towards the bounds. This is key to the conjugate gradient method. It searches along a vector 
				// until the direction of the gradient is in the opposite direction (using dot products, i.e. 
				// cosine of angle between them)

				// NR 10.7 states there is no advantage of the variable metric DFP or BFGS methods over
				// conjugate gradient methods. So I will try these first.

				// Try this:
				// Adapt the conjugate gradient optimiser to use the gradient to pick the search direction
				// and then for the line minimisation. However if the function is out of bounds then clip the 
				// variables at the bounds and continue. 
				// If the current point is at the bounds and the gradient is to continue out of bounds then 
				// clip the gradient too.

				// Or: just use the gradient for the search direction then use the line minimisation/rest
				// as per the Powell optimiser. The bounds should limit the search.

				// I tried a Bounded conjugate gradient optimiser with clipped variables:
				// This sometimes works. However when the variables go a long way out of the expected range the gradients
				// can have vastly different magnitudes. This results in the algorithm stalling since the gradients
				// can be close to zero and the some of the parameters are no longer adjusted.
				// Perhaps this can be looked for and the algorithm then gives up and resorts to a Powell optimiser from 
				// the current point.

				// Changed the bracketing step to very small (default is 1, changed to 0.001). This improves the 
				// performance. The gradient direction is very sensitive to small changes in the coordinates so a 
				// tighter bracketing of the line search helps.

				// Tried using a non-gradient method for the line search copied from the Powell optimiser:
				// This also works when the bracketing step is small but the number of iterations is higher.

				// 24.10.2014: I have tried to get conjugate gradient to work but the gradient function 
				// must not behave suitably for the optimiser. In the current state both methods of using a 
				// Bounded Conjugate Gradient Optimiser perform poorly relative to other optimisers:
				// Simulated : n=1000, signal=200, x=0.53, y=0.47
				// LVM : n=1000, signal=171, x=0.537, y=0.471 (1.003s)
				// Powell : n=1000, signal=187, x=0.537, y=0.48 (1.238s)
				// Gradient based PR (constrained): n=858, signal=161, x=0.533, y=0.474 (2.54s)
				// Gradient based PR (bounded): n=948, signal=161, x=0.533, y=0.473 (2.67s)
				// Non-gradient based : n=1000, signal=151.47, x=0.535, y=0.474 (1.626s)
				// The conjugate optimisers are slower, under predict the signal by the most and in the case of 
				// the gradient based optimiser, fail to converge on some problems. This is worse when constrained
				// fitting is used and not tightly bounded fitting.
				// I will leave the code in as an option but would not recommend using it. I may remove it in the 
				// future.

				// Note: It is strange that the non-gradient based line minimisation is more successful.
				// It may be that the gradient function is not accurate (due to round off error) or that it is
				// simple wrong when far from the optimum. My JUnit tests only evaluate the function within the 
				// expected range of the answer.

				// Note the default step size on the Powell optimiser is 1 but the initial directions are unit vectors.
				// So our bracketing step should be a minimum of 1 / average length of the first gradient vector to prevent
				// the first step being too large when bracketing.
				final double gradient[] = new double[startPoint.length];
				maximumLikelihoodFunction.value(startPoint, gradient);
				double l = 0;
				for (double d : gradient)
					l += d * d;
				final double bracketingStep = FastMath.min(0.001, ((l > 1) ? 1.0 / l : 1));
				//System.out.printf("Bracketing step = %f (length=%f)\n", bracketingStep, l);

				o.setUseGradientLineSearch(gradientLineMinimisation);

				optimum = o.optimize(new MaxEval(getMaxEvaluations()), new ObjectiveFunctionGradient(
						new MultivariateVectorPoisson(maximumLikelihoodFunction)), new ObjectiveFunction(
						new MultivariatePoisson(maximumLikelihoodFunction)), GoalType.MINIMIZE, new InitialGuess(
						startPoint), new SimpleBounds(lowerConstraint, upperConstraint),
						new BoundedNonLinearConjugateGradientOptimizer.BracketingStep(bracketingStep));
				iterations = o.getIterations();
				evaluations = o.getEvaluations();

				//maximumLikelihoodFunction.value(solution, gradient);
				//System.out.printf("Iter = %d, %g @ %s : %s\n", iterations, ll, Arrays.toString(solution),
				//		Arrays.toString(gradient));
			}

			final double[] solution = optimum.getPointRef();

			setSolution(a, solution);

			//System.out.printf("Iter = %d, Eval = %d, %g @ %s\n", iterations, evaluations, optimum.getValue(), 
			//	java.util.Arrays.toString(solution));

			// Compute residuals for the FunctionSolver interface
			if (y_fit == null || y_fit.length < n)
				y_fit = new double[n];
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
				// Assume the Maximum Likelihood estimator returns the optimum fit (achieves the Cramer Roa
				// lower bounds) and so the covariance can be obtained from the Fisher Information Matrix. 
				int[] gradientIndices = f.gradientIndices();
				final int nparams = gradientIndices.length;
				GradientCalculator calculator = GradientCalculatorFactory.newCalculator(nparams);
				double[][] I = calculator.fisherInformationMatrix(n, a, f);
				for (int i = 0; i < gradientIndices.length; i++)
					a_dev[gradientIndices[i]] = 1.0 / Math.sqrt(I[i][i]);
				// The following method just uses the sqrt but does not invert the covariance matrix
				//setDeviations(a_dev, covar);
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
			// TODO - Find out the other exceptions from the fitters and add return values to match.
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

	/**
	 * This setting applies to the conjugate gradient method of the Maximum Likelihood Estimator
	 * 
	 * @return the gradientLineMinimisation True if using the gradient for line minimisation
	 */
	public boolean isGradientLineMinimisation()
	{
		return gradientLineMinimisation;
	}

	/**
	 * This setting applies to the conjugate gradient method of the Maximum Likelihood Estimator
	 * 
	 * @param gradientLineMinimisation
	 *            Set to true to use the gradient for line minimisation
	 */
	public void setGradientLineMinimisation(boolean gradientLineMinimisation)
	{
		this.gradientLineMinimisation = gradientLineMinimisation;
	}

	/**
	 * @return the relative threshold for convergence in the Maximum Likelihood Estimator
	 */
	public double getRelativeThreshold()
	{
		return relativeThreshold;
	}

	/**
	 * @param relativeThreshold
	 *            the relative threshold for convergence in the Maximum Likelihood Estimator
	 */
	public void setRelativeThreshold(double relativeThreshold)
	{
		this.relativeThreshold = relativeThreshold;
	}

	/**
	 * @return the absolute threshold for convergence in the Maximum Likelihood Estimator
	 */
	public double getAbsoluteThreshold()
	{
		return absoluteThreshold;
	}

	/**
	 * @param absoluteThreshold
	 *            the absolute threshold for convergence in the Maximum Likelihood Estimator
	 */
	public void setAbsoluteThreshold(double absoluteThreshold)
	{
		this.absoluteThreshold = absoluteThreshold;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#isBounded()
	 */
	@Override
	public boolean isBounded()
	{
		switch (searchMethod)
		{
			case POWELL_BOUNDED:
			case BOBYQA:
			case CMAES:
			case BFGS:
				return true;
			default:
				return false;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#isConstrained()
	 */
	@Override
	public boolean isConstrained()
	{
		switch (searchMethod)
		{
			case CONJUGATE_GRADIENT_FR:
			case CONJUGATE_GRADIENT_PR:
			case BFGS:
				return true;
			default:
				return false;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#setBounds(double[], double[])
	 */
	@Override
	public void setBounds(double[] lowerB, double[] upperB)
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#setConstraints(double[], double[])
	 */
	@Override
	public void setConstraints(double[] lowerB, double[] upperB)
	{
		// Extract the bounds for the parameters we are fitting
		int[] indices = f.gradientIndices();

		lowerConstraint = new double[indices.length];
		upperConstraint = new double[indices.length];
		for (int i = 0; i < indices.length; i++)
		{
			lowerConstraint[i] = lowerB[indices[i]];
			upperConstraint[i] = upperB[indices[i]];
		}
	}
}
