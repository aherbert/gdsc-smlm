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
package uk.ac.sussex.gdsc.smlm.fitting;

import java.util.Arrays;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.MultivariateMatrixFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.exception.TooManyIterationsException;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer.Optimum;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.SimpleValueChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;

import uk.ac.sussex.gdsc.core.ij.Utils;
import uk.ac.sussex.gdsc.core.logging.Logger;
import uk.ac.sussex.gdsc.core.utils.Maths;

/**
 * Fit a binomial distribution to a histogram
 */
public class BinomialFitter
{
	private Logger logger = null;
	private boolean maximumLikelihood = true;
	private int fitRestarts = 5;

	/**
	 * Instantiates a new binomial fitter.
	 */
	public BinomialFitter()
	{
	}

	/**
	 * Instantiates a new binomial fitter.
	 *
	 * @param logger
	 *            Logging interface to report progress messages
	 */
	public BinomialFitter(Logger logger)
	{
		this.logger = logger;
	}

	/**
	 * Create a histogram from n=0 to n=N as a normalised probability.
	 * N = p.length - 1;
	 *
	 * @param data
	 *            the data
	 * @param cumulative
	 *            Build a cumulative histogram
	 * @return The cumulative histogram (p)
	 * @throws IllegalArgumentException
	 *             If any of the input data values are negative
	 */
	public static double[] getHistogram(int[] data, boolean cumulative)
	{
		final double[] newData = new double[data.length];
		for (int i = 0; i < data.length; i++)
		{
			if (data[i] < 0)
				throw new IllegalArgumentException("Input data must be positive");
			newData[i] = data[i];
		}
		return calculateHistogram(newData, cumulative);
	}

	/**
	 * Create a histogram from n=0 to n=N as a normalised probability.
	 * N = p.length - 1;
	 *
	 * @param data
	 *            the data
	 * @param cumulative
	 *            Build a cumulative histogram
	 * @return The cumulative histogram (p)
	 * @throws IllegalArgumentException
	 *             If any of the input data values are negative or non-integer
	 */
	public static double[] getHistogram(double[] data, boolean cumulative)
	{
		for (int i = 0; i < data.length; i++)
		{
			if (data[i] < 0)
				throw new IllegalArgumentException("Input data must be positive");
			if ((int) data[i] != data[i])
				throw new IllegalArgumentException("Input data must be integers");
		}
		return calculateHistogram(data, cumulative);
	}

	/**
	 * Create a histogram from n=0 to n=N as a normalised probability.
	 * N = p.length - 1;
	 *
	 * @param data
	 *            the data
	 * @param cumulative
	 *            Build a cumulative histogram
	 * @return The histogram (p)
	 */
	private static double[] calculateHistogram(double[] data, boolean cumulative)
	{
		final double[][] histogram = Maths.cumulativeHistogram(data, true);
		if (histogram[0].length == 0)
			return new double[] { 1 };
		// Pad to include all values
		final double[] nValues = histogram[0];
		final double[] pValues = histogram[1];
		final int N = (int) nValues[nValues.length - 1];
		final double[] p = new double[N + 1];

		// Pad the histogram out for any missing values between 0 and N
		for (int i = 1; i < nValues.length; i++)
		{
			final int j = (int) nValues[i - 1];
			final int k = (int) nValues[i];
			for (int ii = j; ii < k; ii++)
				p[ii] = pValues[i - 1];
		}
		p[N] = pValues[pValues.length - 1];

		// We need the original histogram, not the cumulative histogram
		if (!cumulative)
			for (int i = p.length; i-- > 1;)
				p[i] -= p[i - 1];

		return p;
	}

	/**
	 * Fit the binomial distribution (n,p) to the input data. Performs fitting assuming a fixed n value and attempts to
	 * optimise p. All n from minN to maxN are evaluated. If maxN is zero then all possible n from minN are evaluated
	 * until the fit is worse.
	 *
	 * @param data
	 *            The input data (all value must be positive)
	 * @param minN
	 *            The minimum n to evaluate
	 * @param maxN
	 *            The maximum n to evaluate. Set to zero to evaluate all possible values.
	 * @param zeroTruncated
	 *            True if the model should ignore n=0 (zero-truncated binomial)
	 * @return The best fit (n, p)
	 * @throws IllegalArgumentException
	 *             If any of the input data values are negative
	 */
	public double[] fitBinomial(int[] data, int minN, int maxN, boolean zeroTruncated)
	{
		double[] histogram = getHistogram(data, false);

		final double initialSS = Double.POSITIVE_INFINITY;
		double bestSS = initialSS;
		double[] parameters = null;
		int worse = 0;
		int N = histogram.length - 1;
		if (minN < 1)
			minN = 1;
		if (maxN > 0)
			if (N > maxN)
				// Limit the number fitted to maximum
				N = maxN;
			else if (N < maxN)
			{
				// Expand the histogram to the maximum
				histogram = Arrays.copyOf(histogram, maxN + 1);
				N = maxN;
			}
		if (minN > N)
			minN = N;

		final double mean = getMean(histogram);

		final String name = (zeroTruncated) ? "Zero-truncated Binomial distribution" : "Binomial distribution";

		log("Mean cluster size = %s", Utils.rounded(mean));
		log("Fitting cumulative " + name);

		// Since varying the N should be done in integer steps do this
		// for n=1,2,3,... until the SS peaks then falls off (is worse than the best
		// score several times in succession)
		for (int n = minN; n <= N; n++)
		{
			final PointValuePair solution = fitBinomial(histogram, mean, n, zeroTruncated);
			if (solution == null)
				continue;

			final double p = solution.getPointRef()[0];

			log("Fitted %s : N=%d, p=%s. SS=%g", name, n, Utils.rounded(p), solution.getValue());

			if (bestSS > solution.getValue())
			{
				bestSS = solution.getValue();
				parameters = new double[] { n, p };
				worse = 0;
			}
			else if (bestSS != initialSS)
				if (++worse >= 3)
					break;
		}

		return parameters;
	}

	/**
	 * Fit the binomial distribution (n,p) to the cumulative histogram. Performs fitting assuming a fixed n value and
	 * attempts to optimise p.
	 *
	 * @param histogram
	 *            The input histogram
	 * @param n
	 *            The n to evaluate
	 * @param zeroTruncated
	 *            True if the model should ignore n=0 (zero-truncated binomial)
	 * @return The best fit (n, p)
	 * @throws IllegalArgumentException
	 *             If any of the input data values are negative
	 */
	public PointValuePair fitBinomial(double[] histogram, int n, boolean zeroTruncated)
	{
		return fitBinomial(histogram, Double.NaN, n, zeroTruncated);
	}

	/**
	 * Fit the binomial distribution (n,p) to the cumulative histogram. Performs fitting assuming a fixed n value and
	 * attempts to optimise p.
	 *
	 * @param histogram
	 *            The input histogram
	 * @param mean
	 *            The histogram mean (used to estimate p). Calculated if NaN.
	 * @param n
	 *            The n to evaluate
	 * @param zeroTruncated
	 *            True if the model should ignore n=0 (zero-truncated binomial)
	 * @return The best fit (n, p)
	 * @throws IllegalArgumentException
	 *             If any of the input data values are negative
	 * @throws IllegalArgumentException
	 *             If any fitting a zero truncated binomial and there are no values above zero
	 */
	public PointValuePair fitBinomial(double[] histogram, double mean, int n, boolean zeroTruncated)
	{
		if (Double.isNaN(mean))
			mean = getMean(histogram);

		if (zeroTruncated && histogram[0] > 0)
		{
			log("Fitting zero-truncated histogram but there are zero values - Renormalising to ignore zero");
			double cumul = 0;
			for (int i = 1; i < histogram.length; i++)
				cumul += histogram[i];
			if (cumul == 0)
				throw new IllegalArgumentException("Fitting zero-truncated histogram but there are no non-zero values");
			histogram[0] = 0;
			for (int i = 1; i < histogram.length; i++)
				histogram[i] /= cumul;
		}

		final int nFittedPoints = Math.min(histogram.length, n + 1) - ((zeroTruncated) ? 1 : 0);
		if (nFittedPoints < 1)
		{
			log("No points to fit (%d): Histogram.length = %d, n = %d, zero-truncated = %b", nFittedPoints,
					histogram.length, n, zeroTruncated);
			return null;
		}

		// The model is only fitting the probability p
		// For a binomial n*p = mean => p = mean/n
		final double[] initialSolution = new double[] { FastMath.min(mean / n, 1) };

		// Create the function
		final BinomialModelFunction function = new BinomialModelFunction(histogram, n, zeroTruncated);

		final double[] lB = new double[1];
		final double[] uB = new double[] { 1 };
		final SimpleBounds bounds = new SimpleBounds(lB, uB);

		// Fit
		// CMAESOptimizer or BOBYQAOptimizer support bounds

		// CMAESOptimiser based on Matlab code:
		// https://www.lri.fr/~hansen/cmaes.m
		// Take the defaults from the Matlab documentation
		final int maxIterations = 2000;
		final double stopFitness = 0; //Double.NEGATIVE_INFINITY;
		final boolean isActiveCMA = true;
		final int diagonalOnly = 0;
		final int checkFeasableCount = 1;
		final RandomGenerator random = new Well19937c();
		final boolean generateStatistics = false;
		final ConvergenceChecker<PointValuePair> checker = new SimpleValueChecker(1e-6, 1e-10);
		// The sigma determines the search range for the variables. It should be 1/3 of the initial search region.
		final OptimizationData sigma = new CMAESOptimizer.Sigma(new double[] { (uB[0] - lB[0]) / 3 });
		final OptimizationData popSize = new CMAESOptimizer.PopulationSize((int) (4 + Math.floor(3 * Math.log(2))));

		try
		{
			PointValuePair solution = null;
			boolean noRefit = maximumLikelihood;
			if (n == 1 && zeroTruncated)
			{
				// No need to fit
				solution = new PointValuePair(new double[] { 1 }, 0);
				noRefit = true;
			}
			else
			{
				final GoalType goalType = (maximumLikelihood) ? GoalType.MAXIMIZE : GoalType.MINIMIZE;

				// Iteratively fit
				final CMAESOptimizer opt = new CMAESOptimizer(maxIterations, stopFitness, isActiveCMA, diagonalOnly,
						checkFeasableCount, random, generateStatistics, checker);
				for (int iteration = 0; iteration <= fitRestarts; iteration++)
				{
					try
					{
						// Start from the initial solution
						final PointValuePair result = opt.optimize(new InitialGuess(initialSolution),
								new ObjectiveFunction(function), goalType, bounds, sigma, popSize,
								new MaxIter(maxIterations), new MaxEval(maxIterations * 2));
						//System.out.printf("CMAES Iter %d initial = %g (%d)\n", iteration, result.getValue(),
						//		opt.getEvaluations());
						if (solution == null || result.getValue() < solution.getValue())
							solution = result;
					}
					catch (final TooManyEvaluationsException e)
					{
						// No solution
					}
					catch (final TooManyIterationsException e)
					{
						// No solution
					}
					if (solution == null)
						continue;
					try
					{
						// Also restart from the current optimum
						final PointValuePair result = opt.optimize(new InitialGuess(solution.getPointRef()),
								new ObjectiveFunction(function), goalType, bounds, sigma, popSize,
								new MaxIter(maxIterations), new MaxEval(maxIterations * 2));
						//System.out.printf("CMAES Iter %d restart = %g (%d)\n", iteration, result.getValue(),
						//		opt.getEvaluations());
						if (result.getValue() < solution.getValue())
							solution = result;
					}
					catch (final TooManyEvaluationsException e)
					{
						// No solution
					}
					catch (final TooManyIterationsException e)
					{
						// No solution
					}
				}
				if (solution == null)
					return null;
			}

			if (noRefit)
			{
				// Although we fit the log-likelihood, return the sum-of-squares to allow
				// comparison across different n
				final double p = solution.getPointRef()[0];
				double ss = 0;
				final double[] obs = function.p;
				final double[] exp = function.getP(p);
				for (int i = 0; i < obs.length; i++)
					ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);
				return new PointValuePair(solution.getPointRef(), ss);
			}
			// We can do a LVM refit if the number of fitted points is more than 1
			else if (nFittedPoints > 1)
			{
				// Improve SS fit with a gradient based LVM optimizer
				final LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();

				try
				{
					final BinomialModelFunctionGradient gradientFunction = new BinomialModelFunctionGradient(histogram,
							n, zeroTruncated);

					//@formatter:off
					final LeastSquaresProblem problem = new LeastSquaresBuilder()
							.maxEvaluations(Integer.MAX_VALUE)
							.maxIterations(3000)
							.start(solution.getPointRef())
							.target(gradientFunction.p)
							.weight(new DiagonalMatrix(gradientFunction.getWeights()))
							.model(gradientFunction, new MultivariateMatrixFunction() {
								@Override
								public double[][] value(double[] point) throws IllegalArgumentException
								{
									return gradientFunction.jacobian(point);
								}} )
							//.checker (checker)
							.build();
					//@formatter:on

					final Optimum lvmSolution = optimizer.optimize(problem);

					// Check the pValue is valid since the LVM is not bounded.
					final double p = lvmSolution.getPoint().getEntry(0);
					if (p <= 1 && p >= 0)
					{
						// True if the weights are 1
						final double ss = lvmSolution.getResiduals().dotProduct(lvmSolution.getResiduals());
						//double ss = 0;
						//double[] obs = gradientFunction.p;
						//double[] exp = gradientFunction.value(lvmSolution.getPoint().toArray());
						//for (int i = 0; i < obs.length; i++)
						//	ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);
						if (ss < solution.getValue())
							//log("Re-fitting improved the SS from %s to %s (-%s%%)",
							//		Utils.rounded(solution.getValue(), 4), Utils.rounded(ss, 4),
							//		Utils.rounded(100 * (solution.getValue() - ss) / solution.getValue(), 4));
							return new PointValuePair(lvmSolution.getPoint().toArray(), ss);
					}
				}
				catch (final TooManyIterationsException e)
				{
					log("Failed to re-fit: Too many iterations: %s", e.getMessage());
				}
				catch (final ConvergenceException e)
				{
					log("Failed to re-fit: %s", e.getMessage());
				}
				catch (final Exception e)
				{
					// Ignore this ...
				}
			}

			return solution;
		}
		catch (final Exception e)
		{
			log("Failed to fit Binomial distribution with N=%d : %s", n, e.getMessage());
		}
		return null;
	}

	/**
	 * Gets the mean.
	 *
	 * @param histogram
	 *            the histogram
	 * @return the mean
	 */
	private static double getMean(double[] histogram)
	{
		double sum = 0;
		double count = 0;
		for (int i = 0; i < histogram.length; i++)
		{
			sum += histogram[i] * i;
			count += histogram[i];
		}
		final double mean = sum / count;
		return mean;
	}

	/**
	 * Evaluates the cumulative binomial probability distribution. Assumes the
	 * input data is a cumulative histogram from 0 to N in integer increments.
	 */
	private class BinomialModel
	{
		int trials;
		double[] p;
		int startIndex;

		/**
		 * Create a new Binomial model using the input p-values
		 *
		 * @param p
		 *            The observed p-value
		 * @param trials
		 *            The number of trials
		 * @param zeroTruncated
		 *            Set to true to ignore the x=0 datapoint
		 */
		public BinomialModel(double[] p, int trials, boolean zeroTruncated)
		{
			this.trials = trials;
			startIndex = (zeroTruncated) ? 1 : 0;
			this.p = p;
		}

		/**
		 * Get the probability function for the input pValue.
		 *
		 * @param pValue
		 *            the value
		 * @return the probability function
		 */
		public double[] getP(double pValue)
		{
			final BinomialDistribution dist = new BinomialDistribution(trials, pValue);

			// Optionally ignore x=0 since we cannot see a zero size cluster.
			// This is done by re-normalising the cumulative probability excluding x=0
			// to match the input curve.
			//
			// See Zero-truncated (zt) binomial distribution:
			// http://www.vosesoftware.com/ModelRiskHelp/index.htm#Distributions/Discrete_distributions/Zero-truncated_binomial_distribution.htm
			// pi = 1 / ( 1 - f(0) )
			// Fzt(x) = pi . F(x)

			final double[] p2 = new double[p.length];
			for (int i = startIndex; i <= trials; i++)
				p2[i] = dist.probability(i);

			// Renormalise if necessary
			if (startIndex == 1)
			{
				final double pi = 1.0 / (1.0 - dist.probability(0));
				for (int i = 1; i <= trials; i++)
					p2[i] *= pi;
			}

			return p2;
		}
	}

	/**
	 * Allow optimisation using Apache Commons Math 3 MultivariateFunction optimisers
	 */
	public class BinomialModelFunction extends BinomialModel implements MultivariateFunction
	{
		/**
		 * Instantiates a new binomial model function.
		 *
		 * @param p
		 *            the p
		 * @param trials
		 *            the trials
		 * @param zeroTruncated
		 *            the zero truncated
		 */
		public BinomialModelFunction(double[] p, int trials, boolean zeroTruncated)
		{
			super(p, trials, zeroTruncated);
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see org.apache.commons.math3.analysis.MultivariateFunction#value(double[])
		 */
		@Override
		public double value(double[] parameters)
		{
			final double[] p2 = getP(parameters[0]);
			if (maximumLikelihood)
			{
				// Calculate the log-likelihood
				double ll = 0;
				// We cannot produce a likelihood for any n>N
				final int limit = trials + 1; // p.length
				for (int i = startIndex; i < limit; i++)
					// Sum for all observations the probability of the observation.
					// Use p[i] to indicate the frequency of this observation.
					ll += p[i] * Math.log(p2[i]);
				//System.out.printf("%f => %f\n", parameters[0], ll);
				return ll;
			}
			// Calculate the sum of squares
			double ss = 0;
			for (int i = startIndex; i < p.length; i++)
			{
				final double dx = p[i] - p2[i];
				ss += dx * dx;
			}
			return ss;
		}
	}

	/**
	 * Allow optimisation using Apache Commons Math 3 MultivariateFunction optimisers
	 */
	private class BinomialModelFunctionGradient extends BinomialModel implements MultivariateVectorFunction
	{
		long[] nC;

		/**
		 * Instantiates a new binomial model function gradient.
		 *
		 * @param histogram
		 *            the histogram
		 * @param trials
		 *            the trials
		 * @param zeroTruncated
		 *            the zero truncated
		 */
		public BinomialModelFunctionGradient(double[] histogram, int trials, boolean zeroTruncated)
		{
			super(histogram, trials, zeroTruncated);

			// We could ignore the first p value as it is always zero:
			//p = Arrays.copyOfRange(p, 1, p.length);
			// BUT then we would have to override the getP() method since this has
			// an offset of 1 and assumes the index of p is X.

			final int n = trials;
			nC = new long[n + 1];
			for (int k = 0; k <= n; k++)
				nC[k] = CombinatoricsUtils.binomialCoefficient(n, k);
		}

		private double[] getWeights()
		{
			final double[] w = new double[p.length];
			Arrays.fill(w, 1);
			return w;
		}

		/*
		 * (non-Javadoc)
		 *
		 * @see org.apache.commons.math3.analysis.MultivariateFunction#value(double[])
		 */
		@Override
		public double[] value(double[] point) throws IllegalArgumentException
		{
			return getP(point[0]);
		}

		double[][] jacobian(double[] variables)
		{
			// We could do analytical differentiation for the normal binomial:
			// pmf = nCk * p^k * (1-p)^(n-k)
			// pmf' = nCk * k*p^(k-1) * (1-p)^(n-k) +
			//        nCk * p^k * (n-k) * (1-p)^(n-k-1) * -1

			final double p = variables[0];
			final double[][] jacobian = new double[this.p.length][1];

			// Compute the gradient using analytical differentiation
			final int n = trials;

			// Note FastMath has specific case for integer power argument which is faster.

			if (startIndex == 0)
				for (int k = 0; k <= n; ++k)
				{
					//jacobian[k][0] = nC[k] * k * Math.pow(p, k - 1) * Math.pow(1 - p, n - k) +
					//		nC[k] * Math.pow(p, k) * (n - k) * Math.pow(1 - p, n - k - 1) * -1;

					// Optimise
					final double pk_1 = FastMath.pow(p, k - 1);
					final double pk = p * pk_1;
					final double q = 1 - p;
					final double q_n_k_1 = FastMath.pow(q, n - k - 1);
					final double q_n_k = q * q_n_k_1;
					jacobian[k][0] = nC[k] * (k * pk_1 * q_n_k - pk * (n - k) * q_n_k_1);
				}
			else
			{
				// Account for zero-truncated distribution
				jacobian[0][0] = 0;

				// In the zero-truncated Binomial all values are scaled by a factor
				// pi = 1.0 / (1.0 - dist.probability(0));

				// We must apply the product rule with pi as f
				// (f.g)' = f'.g +f.g'

				// So far we have only computed g' for the original Binomial

				//double pi = dist.probability(0);
				final double q = 1 - p;
				final double p_n = FastMath.pow(1 - p, n);
				final double f = 1.0 / (1.0 - nC[0] * p_n);
				final double ff = -1 / FastMath.pow(1.0 - nC[0] * p_n, 2) + n * FastMath.pow(q, n - 1);

				for (int k = 1; k <= n; ++k)
				{
					final double pk_1 = FastMath.pow(p, k - 1);
					final double pk = p * pk_1;
					final double q_n_k_1 = FastMath.pow(q, n - k - 1);
					final double q_n_k = q * q_n_k_1;

					final double g = nC[k] * pk * q_n_k;
					// Differentiate as above
					final double gg = nC[k] * (k * pk_1 * q_n_k - pk * (n - k) * q_n_k_1);
					jacobian[k][0] = ff * g + f * gg;
				}
			}

			//			// Compute the gradients using numerical differentiation.
			//			// Set the step h for computing the function around the desired point.
			//			// Set the delta using the desired fractional accuracy.
			//			// See Numerical Recipes, The Art of Scientific Computing (2nd edition) Chapter 5.7
			//			// on numerical derivatives
			//			final double delta = FastMath.cbrt(1e-6);
			//			final double h = delta * p;
			//
			//			// Ensure we stay within the 0-1 bounds
			//			final double upperP = Math.min(1, p + h);
			//			final double lowerP = Math.max(0, p - h);
			//			final double diff = upperP - lowerP;
			//			double[] upper = getP(upperP);
			//			double[] lower = getP(lowerP);
			//
			//			for (int i = startIndex; i <= trials; i++)
			//			{
			//				double g = (upper[i] - lower[i]) / diff;
			//				if (trials > 1)
			//					System.out.printf("(%d,%f)[%d] %f vs %f\n", trials, p, i, jacobian[i][0], g);
			//				jacobian[i][0] = g;
			//			}
			return jacobian;
		}
	}

	private void log(String format, Object... args)
	{
		if (logger != null)
			logger.info(format, args);
	}

	/**
	 * @return True if use maximum likelihood fitting
	 */
	public boolean isMaximumLikelihood()
	{
		return maximumLikelihood;
	}

	/**
	 * @param maximumLikelihood
	 *            True if use maximum likelihood fitting
	 */
	public void setMaximumLikelihood(boolean maximumLikelihood)
	{
		this.maximumLikelihood = maximumLikelihood;
	}

	/**
	 * @return the number of restarts for fitting
	 */
	public int getFitRestarts()
	{
		return fitRestarts;
	}

	/**
	 * Since fitting uses a bounded search seeded with random movements, restarting can improve the fit. Control the
	 * number of restarts used fot fitting.
	 *
	 * @param fitRestarts
	 *            the number of restarts for fitting
	 */
	public void setFitRestarts(int fitRestarts)
	{
		this.fitRestarts = Math.max(0, fitRestarts);
	}
}
