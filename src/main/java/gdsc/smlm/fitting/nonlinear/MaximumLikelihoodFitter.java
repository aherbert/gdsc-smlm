package gdsc.smlm.fitting.nonlinear;

import gdsc.smlm.fitting.FisherInformationMatrix;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.function.LikelihoodWrapper;
import gdsc.smlm.function.NonLinearFunction;
import gdsc.smlm.function.PoissonGammaGaussianLikelihoodWrapper;
import gdsc.smlm.function.PoissonGaussianLikelihoodWrapper;
import gdsc.smlm.function.PoissonLikelihoodWrapper;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;

import java.util.Arrays;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.exception.TooManyIterationsException;
import org.apache.commons.math3.optim.BaseOptimizer;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.PointValuePair;
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
 * By default the probability mass function for observed value k is modelled as a Poisson process:<br/>
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
 * <p>
 * The probability mass function can be changed to a Poisson-Gaussian or Poisson-Gamma-Gaussian distribution in order to
 * model the counts from a CCD/EMCCD camera.
 */
public class MaximumLikelihoodFitter extends BaseFunctionSolver
{
	/**
	 * Wrap the LikelihoodFunction with classes that implement the required interfaces
	 */
	private class Likelihood
	{
		LikelihoodWrapper fun;

		public Likelihood(LikelihoodWrapper fun)
		{
			this.fun = fun;
		}
	}

	private class MultivariateLikelihood extends Likelihood implements MultivariateFunction
	{
		public MultivariateLikelihood(LikelihoodWrapper fun)
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
			return fun.likelihood(point);
		}

		public boolean isMapped()
		{
			return false;
		}
	}

	/**
	 * Map the specified indices using the sqrt function for use with the Powell optimiser
	 */
	private class MappedMultivariateLikelihood extends MultivariateLikelihood
	{
		final int[] map;

		public MappedMultivariateLikelihood(LikelihoodWrapper fun, int[] map)
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
			return fun.likelihood(unmap(point));
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

	private class MultivariateVectorLikelihood extends Likelihood implements MultivariateVectorFunction
	{
		public MultivariateVectorLikelihood(LikelihoodWrapper fun)
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
			fun.likelihood(point, gradient);
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
		POWELL("Powell", false),
		/**
		 * Search using Powell's conjugate direction method using hard limits to ensure a bounded search
		 */
		POWELL_BOUNDED("Powell (bounded)", false),
		/**
		 * Search using Powell's conjugate direction method using a mapping adapter to ensure a bounded search
		 */
		POWELL_ADAPTER("Powell (adapter)", false),
		/**
		 * Search using Powell's Bound Optimization BY Quadratic Approximation (BOBYQA) algorithm.
		 * <p>
		 * BOBYQA could also be considered as a replacement of any derivative-based optimizer when the derivatives are
		 * approximated by finite differences. This is a bounded search.
		 */
		BOBYQA("BOBYQA", false),
		/**
		 * Search using active Covariance Matrix Adaptation Evolution Strategy (CMA-ES).
		 * <p>
		 * The CMA-ES is a reliable stochastic optimization method which should be applied if derivative-based methods,
		 * e.g. conjugate gradient, fail due to a rugged search landscape. This is a bounded search.
		 */
		CMAES("CMAES", false),
		/**
		 * Search using a non-linear conjugate gradient optimiser. Use the Fletcher-Reeves update formulas for the
		 * conjugate search directions.
		 * <p>
		 * This is a bounded search using simple truncation of coordinates at the bounds of the search space.
		 */
		CONJUGATE_GRADIENT_FR("Conjugate Gradient Fletcher-Reeves", true),
		/**
		 * Search using a non-linear conjugate gradient optimiser. Use the Polak-Ribière update formulas for the
		 * conjugate search directions.
		 * <p>
		 * This is a bounded search using simple truncation of coordinates at the bounds of the search space.
		 */
		CONJUGATE_GRADIENT_PR("Conjugate Gradient Polak-Ribière", true),
		/**
		 * Search using a Broyden-Fletcher-Goldfarb-Shanno (BFGS) gradient optimiser.
		 */
		BFGS("BFGS", true);

		private final String name;
		private final boolean usesGradient;

		private SearchMethod(String name, boolean usesGradient)
		{
			this.name = name;
			this.usesGradient = usesGradient;
		}

		@Override
		public String toString()
		{
			return name;
		}

		/**
		 * @return True if the search method uses the gradient of the likelihood function
		 */
		public boolean usesGradients()
		{
			return usesGradient;
		}
	}

	public enum LikelihoodFunction
	{
		/**
		 * Use a Poisson likelihood model
		 */
		POISSON("Poisson"),
		/**
		 * Use a Poisson likelihood model
		 */
		POISSON_GAUSSIAN("Poisson+Gaussian"),
		/**
		 * Use a Poisson likelihood model
		 */
		POISSON_GAMMA_GAUSSIAN("Poisson+Gamma+Gaussian");

		private final String name;

		private LikelihoodFunction(String name)
		{
			this.name = name;
		}

		@Override
		public String toString()
		{
			return name;
		}
	}

	private SearchMethod searchMethod = SearchMethod.POWELL;
	private LikelihoodFunction likelihoodFunction = LikelihoodFunction.POISSON;
	private double alpha;
	private double sigma;

	private boolean gradientLineMinimisation = true;
	private double relativeThreshold = 1e-4, absoluteThreshold = 1e-10;
	private double[] lower, upper;
	private double[] lowerConstraint, upperConstraint;

	// The function to use for the Powell optimiser (which may have parameters mapped using the sqrt function) 
	private MultivariateLikelihood powellFunction = null;
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
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#computeFit(int, double[], double[], double[], double[],
	 * double[], double)
	 */
	public FitStatus computeFit(int n, double[] y, double[] y_fit, double[] a, double[] a_dev, double[] error,
			double noise)
	{
		LikelihoodWrapper maximumLikelihoodFunction = createLikelihoodWrapper(n, y, a);

		@SuppressWarnings("rawtypes")
		BaseOptimizer baseOptimiser = null;

		try
		{
			double[] startPoint = getInitialSolution(a);

			PointValuePair optimum = null;
			if (searchMethod == SearchMethod.POWELL || searchMethod == SearchMethod.POWELL_BOUNDED ||
					searchMethod == SearchMethod.POWELL_ADAPTER)
			{
				// Non-differentiable version using Powell Optimiser

				// This is as per the method in Numerical Recipes 10.5 (Direction Set (Powell's) method)
				// I could extend the optimiser and implement bounds on the directions moved. However the mapping
				// adapter seems to work OK.

				final boolean basisConvergence = false;

				// Perhaps these thresholds should be tighter?
				// The default is to use the sqrt() of the overall tolerance
				//final double lineRel = FastMath.sqrt(relativeThreshold);
				//final double lineAbs = FastMath.sqrt(absoluteThreshold);
				//final double lineRel = relativeThreshold * 1e2;
				//final double lineAbs = absoluteThreshold * 1e2;

				// Since we are fitting only a small number of parameters then just use the same tolerance 
				// for each search direction
				final double lineRel = relativeThreshold;
				final double lineAbs = absoluteThreshold;

				CustomPowellOptimizer o = new CustomPowellOptimizer(relativeThreshold, absoluteThreshold, lineRel,
						lineAbs, null, basisConvergence);
				baseOptimiser = o;

				OptimizationData maxIterationData = null;
				if (getMaxIterations() > 0)
					maxIterationData = new MaxIter(getMaxIterations());

				if (searchMethod == SearchMethod.POWELL_ADAPTER)
				{
					// Try using the mapping adapter for a bounded Powell search
					MultivariateFunctionMappingAdapter adapter = new MultivariateFunctionMappingAdapter(
							new MultivariateLikelihood(maximumLikelihoodFunction), lower, upper);
					optimum = o.optimize(maxIterationData, new MaxEval(getMaxEvaluations()),
							new ObjectiveFunction(adapter), GoalType.MINIMIZE,
							new InitialGuess(adapter.boundedToUnbounded(startPoint)));
					double[] solution = adapter.unboundedToBounded(optimum.getPointRef());
					optimum = new PointValuePair(solution, optimum.getValue());
				}
				else
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
								powellFunction = new MappedMultivariateLikelihood(maximumLikelihoodFunction,
										Arrays.copyOf(map, count));
							}
						}
						if (powellFunction == null)
						{
							powellFunction = new MultivariateLikelihood(maximumLikelihoodFunction);
						}
					}

					// Update the maximum likelihood function in the Powell function wrapper
					powellFunction.fun = maximumLikelihoodFunction;

					OptimizationData positionChecker = null;
					// new org.apache.commons.math3.optim.PositionChecker(relativeThreshold, absoluteThreshold);
					SimpleBounds simpleBounds = null;
					if (powellFunction.isMapped())
					{
						MappedMultivariateLikelihood adapter = (MappedMultivariateLikelihood) powellFunction;
						if (searchMethod == SearchMethod.POWELL_BOUNDED)
							simpleBounds = new SimpleBounds(adapter.map(lower), adapter.map(upper));
						optimum = o.optimize(maxIterationData, new MaxEval(getMaxEvaluations()),
								new ObjectiveFunction(powellFunction), GoalType.MINIMIZE,
								new InitialGuess(adapter.map(startPoint)), positionChecker, simpleBounds);
						double[] solution = adapter.unmap(optimum.getPointRef());
						optimum = new PointValuePair(solution, optimum.getValue());
					}
					else
					{
						if (searchMethod == SearchMethod.POWELL_BOUNDED)
							simpleBounds = new SimpleBounds(lower, upper);
						optimum = o.optimize(maxIterationData, new MaxEval(getMaxEvaluations()),
								new ObjectiveFunction(powellFunction), GoalType.MINIMIZE, new InitialGuess(startPoint),
								positionChecker, simpleBounds);
					}
				}
			}
			else if (searchMethod == SearchMethod.BOBYQA)
			{
				// Differentiable approximation using Powell's BOBYQA algorithm.
				// This is slower than the Powell optimiser and requires a high number of evaluations.
				int numberOfInterpolationPoints = this.getNumberOfFittedParameters() + 2;

				BOBYQAOptimizer o = new BOBYQAOptimizer(numberOfInterpolationPoints);
				baseOptimiser = o;
				optimum = o.optimize(new MaxEval(getMaxEvaluations()),
						new ObjectiveFunction(new MultivariateLikelihood(maximumLikelihoodFunction)), GoalType.MINIMIZE,
						new InitialGuess(startPoint), new SimpleBounds(lower, upper));
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
						new ObjectiveFunction(new MultivariateLikelihood(maximumLikelihoodFunction)), GoalType.MINIMIZE,
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
							checkFeasableCount, random, generateStatistics,
							new SimpleValueChecker(relativeThreshold, absoluteThreshold));
					baseOptimiser = o;
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
				// Prevent incrementing the iterations again
				baseOptimiser = null;
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
				baseOptimiser = o;

				// Configure maximum step length for each dimension using the bounds
				double[] stepLength = new double[lower.length];
				for (int i = 0; i < stepLength.length; i++)
				{
					stepLength[i] = (upper[i] - lower[i]) * 0.3333333;
					if (stepLength[i] <= 0)
						stepLength[i] = Double.POSITIVE_INFINITY;
				}

				// The GoalType is always minimise so no need to pass this in
				OptimizationData positionChecker = null;
				//new org.apache.commons.math3.optim.PositionChecker(relativeThreshold, absoluteThreshold);
				optimum = o.optimize(new MaxEval(getMaxEvaluations()),
						new ObjectiveFunctionGradient(new MultivariateVectorLikelihood(maximumLikelihoodFunction)),
						new ObjectiveFunction(new MultivariateLikelihood(maximumLikelihoodFunction)),
						new InitialGuess(startPoint), new SimpleBounds(lowerConstraint, upperConstraint),
						new BFGSOptimizer.GradientTolerance(relativeThreshold), positionChecker,
						new BFGSOptimizer.StepLength(stepLength));
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
								: Formula.POLAK_RIBIERE,
						new SimpleValueChecker(relativeThreshold, absoluteThreshold));
				baseOptimiser = o;

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
				// simply wrong when far from the optimum. My JUnit tests only evaluate the function within the 
				// expected range of the answer.

				// Note the default step size on the Powell optimiser is 1 but the initial directions are unit vectors.
				// So our bracketing step should be a minimum of 1 / average length of the first gradient vector to prevent
				// the first step being too large when bracketing.
				final double gradient[] = new double[startPoint.length];
				maximumLikelihoodFunction.likelihood(startPoint, gradient);
				double l = 0;
				for (double d : gradient)
					l += d * d;
				final double bracketingStep = FastMath.min(0.001, ((l > 1) ? 1.0 / l : 1));
				//System.out.printf("Bracketing step = %f (length=%f)\n", bracketingStep, l);

				o.setUseGradientLineSearch(gradientLineMinimisation);

				optimum = o.optimize(new MaxEval(getMaxEvaluations()),
						new ObjectiveFunctionGradient(new MultivariateVectorLikelihood(maximumLikelihoodFunction)),
						new ObjectiveFunction(new MultivariateLikelihood(maximumLikelihoodFunction)), GoalType.MINIMIZE,
						new InitialGuess(startPoint), new SimpleBounds(lowerConstraint, upperConstraint),
						new BoundedNonLinearConjugateGradientOptimizer.BracketingStep(bracketingStep));

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
				FisherInformationMatrix m = new FisherInformationMatrix(maximumLikelihoodFunction.fisherInformation(a));
				double[] crlb = m.crlb(true);
				Arrays.fill(a_dev, 0);
				final int[] gradientIndices = f.gradientIndices();
				for (int i = gradientIndices.length; i-- > 0;)
					a_dev[gradientIndices[i]] = crlb[i];
			}

			error[0] = NonLinearFit.getError(residualSumOfSquares, noise, n, f.gradientIndices().length);
			// Reverse negative log likelihood for maximum likelihood score
			value = -optimum.getValue();
		}
		catch (TooManyIterationsException e)
		{
			//System.out.printf("Too many iterations = %d\n", e.getMax());
			//e.printStackTrace();
			return FitStatus.TOO_MANY_ITERATIONS;
		}
		catch (TooManyEvaluationsException e)
		{
			//System.out.printf("Too many evaluations = %d\n", e.getMax());
			//e.printStackTrace();
			return FitStatus.TOO_MANY_EVALUATIONS;
		}
		catch (ConvergenceException e)
		{
			// Occurs when QR decomposition fails - mark as a singular non-linear model (no solution)
			//System.out.printf("Singular non linear model = %s\n", e.getMessage());
			return FitStatus.SINGULAR_NON_LINEAR_MODEL;
		}
		catch (BFGSOptimizer.LineSearchRoundoffException e)
		{
			//System.out.println("BFGS error: " + e.getMessage());
			//e.printStackTrace();
			return FitStatus.FAILED_TO_CONVERGE;
		}
		catch (Exception e)
		{
			//System.out.printf("Unknown error = %s\n", e.getMessage());
			e.printStackTrace();
			return FitStatus.UNKNOWN;
		}
		finally
		{
			if (baseOptimiser != null)
			{
				iterations += baseOptimiser.getIterations();
				evaluations += baseOptimiser.getEvaluations();
			}
		}

		// Check this as likelihood functions can go wrong
		if (Double.isInfinite(value) || Double.isNaN(value))
			return FitStatus.INVALID_LIKELIHOOD;

		return FitStatus.OK;
	}

	private LikelihoodWrapper createLikelihoodWrapper(int n, double[] y, double[] a)
	{
		LikelihoodWrapper maximumLikelihoodFunction = null;

		final double myAlpha = (this.alpha > 0) ? this.alpha : 1;

		// We can use different likelihood wrapper functions:
		switch (likelihoodFunction)
		{
			case POISSON_GAMMA_GAUSSIAN:
				// Poisson-Gamma-Gaussian - EM-CCD data
				maximumLikelihoodFunction = new PoissonGammaGaussianLikelihoodWrapper(f, a, y, n, myAlpha, sigma);
				break;

			case POISSON_GAUSSIAN:
				// Poisson-Gaussian - CCD data
				// Sigma must be positive, otherwise fall back to a Poisson likelihood function
				if (sigma > 0)
				{
					maximumLikelihoodFunction = new PoissonGaussianLikelihoodWrapper(f, a, y, n, myAlpha, sigma);
					break;
				}

			case POISSON:
			default:
				// Poisson - most counting data
		}

		// Check if the method requires the gradient but it cannot be computed
		if (maximumLikelihoodFunction == null ||
				(searchMethod.usesGradient && !maximumLikelihoodFunction.canComputeGradient()))
		{
			// Ensure no negative data for the Poisson likelihood method.
			// Just truncate the counts for now. These are from noise in the count estimates that we do not model.
			final double[] y2 = new double[n];
			for (int i = 0; i < n; i++)
			{
				if (y[i] < 0)
					y2[i] = 0;
				else
					y2[i] = y[i];
			}
			PoissonLikelihoodWrapper function = new PoissonLikelihoodWrapper(f, a, y2, n, myAlpha);
			// This will allow Powell searches. The effect on the gradient search algorithms may be weird so leave alone.
			if (!searchMethod.usesGradient)
				function.setAllowNegativeExpectedValues(true);
			maximumLikelihoodFunction = function;
		}
		return maximumLikelihoodFunction;
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
	 * @return the likelihood function to model the count
	 */
	public LikelihoodFunction getLikelihoodFunction()
	{
		return likelihoodFunction;
	}

	/**
	 * @param likelihoodFunction
	 *            the likelihood function to model the count
	 */
	public void setLikelihoodFunction(LikelihoodFunction likelihoodFunction)
	{
		this.likelihoodFunction = likelihoodFunction;
	}

	/**
	 * @return the alpha for the gamma component of the Poisson-Gamma-Gaussian likelihood function
	 */
	public double getAlpha()
	{
		return alpha;
	}

	/**
	 * @param alpha
	 *            the alpha for the gamma component of the Poisson-Gamma-Gaussian likelihood function
	 */
	public void setAlpha(double alpha)
	{
		this.alpha = alpha;
	}

	/**
	 * @return the sigma for the Gaussian component of the Poisson-Gaussian/Poisson-Gamma-Gaussian likelihood function
	 */
	public double getSigma()
	{
		return sigma;
	}

	/**
	 * @param sigma
	 *            the sigma for the Gaussian component of the Poisson-Gaussian/Poisson-Gamma-Gaussian likelihood
	 *            function
	 */
	public void setSigma(double sigma)
	{
		this.sigma = sigma;
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

	/**
	 * Note that certain search methods require bounds to function. A null pointer exception will be thrown by
	 * fitter if the bounds are not set for these methods.
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#isBounded()
	 */
	@Override
	public boolean isBounded()
	{
		switch (searchMethod)
		{
			case POWELL_ADAPTER:
			case POWELL_BOUNDED:
			case BOBYQA:
			case CMAES:
			case BFGS:
				return true;
			default:
				return false;
		}
	}

	/**
	 * Note that certain search methods require constraints to function. A null pointer exception will be thrown by
	 * fitter if the constraints are not set for these methods.
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

	@Override
	public boolean computeValue(int n, double[] y, double[] y_fit, double[] a)
	{
		LikelihoodWrapper maximumLikelihoodFunction = createLikelihoodWrapper(n, y, a);

		final double l = maximumLikelihoodFunction.likelihood(a);
		if (l == Double.POSITIVE_INFINITY)
			return false;

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

		// Reverse negative log likelihood for maximum likelihood score
		value = -l;

		return false;
	}
}
