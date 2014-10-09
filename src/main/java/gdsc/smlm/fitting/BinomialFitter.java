package gdsc.smlm.fitting;

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

import gdsc.smlm.fitting.logging.Logger;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.utils.Maths;

import java.util.Arrays;

import org.apache.commons.math3.analysis.DifferentiableMultivariateVectorFunction;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.MultivariateMatrixFunction;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
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
import org.apache.commons.math3.optimization.PointVectorValuePair;
import org.apache.commons.math3.optimization.general.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

/**
 * Fit a binomial distribution to a histogram
 */
public class BinomialFitter
{
	private Logger logger = null;
	private boolean maximumLikelihood = true;

	public BinomialFitter()
	{

	}

	/**
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
	 * @param cumulative
	 *            Build a cumulative histogram
	 * @return The cumulative histogram (p)
	 * @throws IllegalArgumentException
	 *             If any of the input data values are negative
	 */
	public static double[] getHistogram(int[] data, boolean cumulative)
	{
		double[] newData = new double[data.length];
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
	 * @param cumulative
	 *            Build a cumulative histogram
	 * @return The histogram (p)
	 */
	private static double[] calculateHistogram(double[] data, boolean cumulative)
	{
		double[][] histogram = Maths.cumulativeHistogram(data, true);
		if (histogram[0].length == 0)
			return new double[] { 1 };
		// Pad to include all values
		double[] nValues = histogram[0];
		double[] pValues = histogram[1];
		int N = (int) nValues[nValues.length - 1];
		double[] p = new double[N + 1];

		// Pad the histogram out for any missing values between 0 and N
		for (int i = 1; i < nValues.length; i++)
		{
			int j = (int) nValues[i - 1];
			int k = (int) nValues[i];
			for (int ii = j; ii < k; ii++)
				p[ii] = pValues[i - 1];
		}
		p[N] = pValues[pValues.length - 1];

		// We need the original histogram, not the cumulative histogram
		if (!cumulative)
		{
			for (int i = p.length; i-- > 1;)
			{
				p[i] -= p[i - 1];
			}
		}

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
		int N = (int) histogram.length - 1;
		if (minN < 1)
			minN = 1;
		if (maxN > 0)
		{
			if (N > maxN)
			{
				// Limit the number fitted to maximum
				N = maxN;
			}
			else if (N < maxN)
			{
				// Expand the histogram to the maximum
				histogram = Arrays.copyOf(histogram, maxN + 1);
				N = maxN;
			}
		}
		if (minN > N)
			minN = N;

		final double mean = getMean(histogram);

		String name = (zeroTruncated) ? "Zero-truncated Binomial distribution" : "Binomial distribution";

		log("Mean cluster size = %s", Utils.rounded(mean));
		log("Fitting cumulative " + name);

		// Since varying the N should be done in integer steps do this
		// for n=1,2,3,... until the SS peaks then falls off (is worse than the best 
		// score several times in succession)
		for (int n = minN; n <= N; n++)
		{
			PointValuePair solution = fitBinomial(histogram, mean, n, zeroTruncated);
			if (solution == null)
				continue;

			double p = solution.getPointRef()[0];

			log("Fitted %s : N=%d, p=%s. SS=%g", name, n, Utils.rounded(p), solution.getValue());

			if (bestSS > solution.getValue())
			{
				bestSS = solution.getValue();
				parameters = new double[] { n, p };
				worse = 0;
			}
			else if (bestSS != initialSS)
			{
				if (++worse >= 3)
					break;
			}
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

		// The model is only fitting the probability p
		// For a binomial n*p = mean => p = mean/n
		double[] initialSolution = new double[] { Math.min(mean / n, 1) };

		// Create the function
		BinomialModelFunction function = new BinomialModelFunction(histogram, n, zeroTruncated);

		double[] lB = new double[1];
		double[] uB = new double[] { 1 };
		SimpleBounds bounds = new SimpleBounds(lB, uB);

		// Fit
		// CMAESOptimizer or BOBYQAOptimizer support bounds

		// CMAESOptimiser based on Matlab code:
		// https://www.lri.fr/~hansen/cmaes.m
		// Take the defaults from the Matlab documentation
		int maxIterations = 2000;
		double stopFitness = 0; //Double.NEGATIVE_INFINITY;
		boolean isActiveCMA = true;
		int diagonalOnly = 0;
		int checkFeasableCount = 1;
		RandomGenerator random = new Well19937c();
		boolean generateStatistics = false;
		ConvergenceChecker<PointValuePair> checker = new SimpleValueChecker(1e-6, 1e-10);
		// The sigma determines the search range for the variables. It should be 1/3 of the initial search region.
		OptimizationData sigma = new CMAESOptimizer.Sigma(new double[] { (uB[0] - lB[0]) / 3 });
		OptimizationData popSize = new CMAESOptimizer.PopulationSize((int) (4 + Math.floor(3 * Math.log(2))));

		try
		{
			GoalType goalType = (maximumLikelihood) ? GoalType.MAXIMIZE : GoalType.MINIMIZE;
			CMAESOptimizer opt = new CMAESOptimizer(maxIterations, stopFitness, isActiveCMA, diagonalOnly,
					checkFeasableCount, random, generateStatistics, checker);
			PointValuePair solution = opt.optimize(new InitialGuess(initialSolution), new ObjectiveFunction(function),
					goalType, bounds, sigma, popSize, new MaxIter(maxIterations), new MaxEval(maxIterations * 2));
			if (solution == null)
				return null;

			if (maximumLikelihood)
			{
				// Although we fit the log-likelihood, return the sum-of-squares to allow 
				// comparison across different n
				double p = solution.getPointRef()[0];
				double ss = 0;
				double[] obs = function.p;
				double[] exp = function.getP(p);
				for (int i = 0; i < obs.length; i++)
					ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);
				return new PointValuePair(solution.getPointRef(), ss);
			}
			else
			{
				// Improve SS fit with a gradient based LVM optimizer
				LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();
				try
				{
					BinomialModelFunctionGradient gradientFunction = new BinomialModelFunctionGradient(histogram, n,
							zeroTruncated);
					PointVectorValuePair lvmSolution = optimizer.optimize(3000, gradientFunction, gradientFunction.p,
							gradientFunction.getWeights(), solution.getPoint());

					double ss = 0;
					double[] obs = gradientFunction.p;
					double[] exp = lvmSolution.getValue();
					for (int i = 0; i < obs.length; i++)
						ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);
					// Check the pValue is valid since the LVM is not bounded.
					double p = lvmSolution.getPointRef()[0];
					if (ss < solution.getValue() && p <= 1 && p >= 0)
					{
						//log("Re-fitting improved the SS from %s to %s (-%s%%)",
						//		Utils.rounded(solution.getValue(), 4), Utils.rounded(ss, 4),
						//		Utils.rounded(100 * (solution.getValue() - ss) / solution.getValue(), 4));
						return new PointValuePair(lvmSolution.getPoint(), ss);
					}
				}
				catch (TooManyEvaluationsException e)
				{
					log("Failed to re-fit: Too many evaluations (%d)", optimizer.getEvaluations());
				}
				catch (ConvergenceException e)
				{
					log("Failed to re-fit: %s", e.getMessage());
				}
				catch (Exception e)
				{
					// Ignore this ...
				}
			}

			return solution;
		}
		catch (Exception e)
		{
			log("Failed to fit Binomial distribution with N=%d : %s", n, e.getMessage());
		}
		return null;
	}

	private double getMean(double[] histogram)
	{
		double sum = 0;
		double count = 0;
		for (int i = 0; i < histogram.length; i++)
		{
			sum += histogram[i] * i;
			count += histogram[i];
		}
		double mean = sum / count;
		return mean;
	}

	/**
	 * Evaluates the cumulative binomial probability distribution. Assumes the
	 * input data is a cumulative histogram from 0 to N in integer increments.
	 */
	public class BinomialModel
	{
		int trials;
		double[] p;
		int startIndex;

		public BinomialModel(double[] p, int trials, boolean zeroTruncated)
		{
			this.trials = trials;
			startIndex = (zeroTruncated) ? 1 : 0;
			this.p = p;
		}

		/**
		 * Get the cumulative probability function for the input pValue ignoring the X=0 data point
		 * 
		 * @param pValue
		 * @return
		 */
		public double[] getP(double pValue)
		{
			BinomialDistribution dist = new BinomialDistribution(trials, pValue);

			// Optionally ignore x=0 since we cannot see a zero size cluster.
			// This is done by re-normalising the cumulative probability excluding x=0 
			// to match the input curve.
			//
			// See Zero-truncated (zt) binomial distribution:
			// http://www.vosesoftware.com/ModelRiskHelp/index.htm#Distributions/Discrete_distributions/Zero-truncated_binomial_distribution.htm
			// pi = 1 / ( 1 - f(0) )
			// Fzt(x) = pi . F(x)
			//
			// This is equivalent to:
			// Fzt(x) = F(x) / ( 1 - f(0) ) = F(x) / cumul( F(x) )

			double cumul = 0;
			double[] p2 = new double[p.length];
			for (int i = startIndex; i < p.length; i++)
			{
				p2[i] = dist.probability(i);
				cumul += p2[i];
			}

			for (int i = 1; i < p.length; i++)
			{
				p2[i] /= cumul;
			}

			return p2;
		}
	}

	/**
	 * Allow optimisation using Apache Commons Math 3 MultivariateFunction optimisers
	 */
	public class BinomialModelFunction extends BinomialModel implements MultivariateFunction
	{
		public BinomialModelFunction(double[] p, int trials, boolean zeroTruncated)
		{
			super(p, trials, zeroTruncated);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see org.apache.commons.math3.analysis.MultivariateFunction#value(double[])
		 */
		public double value(double[] parameters)
		{
			double[] p2 = getP(parameters[0]);
			if (maximumLikelihood)
			{
				// Calculate the log-likelihood
				double ll = 0;
				// We cannot produce a likelihood for any n>N 
				int limit = trials + 1; // p.length
				for (int i = startIndex; i < limit; i++)
				{
					// Sum for all observations the probability of the observation.
					// Use p[i] to indicate the frequency of this observation. 
					ll += p[i] * Math.log(p2[i]);
				}
				//System.out.printf("%f => %f\n", parameters[0], ll);
				return ll;
			}
			else
			{
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
	}

	/**
	 * Allow optimisation using Apache Commons Math 3 MultivariateFunction optimisers
	 */
	public class BinomialModelFunctionGradient extends BinomialModel implements
			DifferentiableMultivariateVectorFunction
	{
		public BinomialModelFunctionGradient(double[] histogram, int trials, boolean zeroTruncated)
		{
			super(histogram, trials, zeroTruncated);

			// We could ignore the first p value as it is always zero:
			//p = Arrays.copyOfRange(p, 1, p.length);
			// BUT then we would have to override the getP() method since this has 
			// an offset of 1 and assumes the index of p is X.
		}

		public double[] getWeights()
		{
			double[] w = new double[p.length];
			Arrays.fill(w, 1);
			return w;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see org.apache.commons.math3.analysis.MultivariateFunction#value(double[])
		 */
		public double[] value(double[] point) throws IllegalArgumentException
		{
			return getP(point[0]);
		}

		public MultivariateMatrixFunction jacobian()
		{
			return new MultivariateMatrixFunction()
			{
				public double[][] value(double[] variables)
				{
					return jacobian(variables);
				}
			};
		}

		double[][] jacobian(double[] variables)
		{
			// Compute the gradients using numerical differentiation
			final double pValue = variables[0];
			double[][] jacobian = new double[p.length][1];

			final double delta = 0.001 * pValue;
			double[] p2 = getP(pValue);
			double[] p3 = getP(pValue + delta);

			for (int i = 0; i < jacobian.length; ++i)
			{
				jacobian[i][0] = (p3[i] - p2[i]) / delta;
			}
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
}
