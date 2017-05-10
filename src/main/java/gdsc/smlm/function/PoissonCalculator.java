package gdsc.smlm.function;

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Computes likelihood values for a Poisson function
 */
public class PoissonCalculator
{
	/**
	 * Get the Poisson log likelihood of value x given the mean. The mean must be strictly positive. x must be positive.
	 *
	 * @param u
	 *            the mean
	 * @param x
	 *            the x
	 * @return the log likelihood
	 */
	public static double logLikelihood(double u, double x)
	{
		// Unlikely so skip this ...
		//if (x == 0)
		//	return -u;
		return x * Math.log(u) - u - logFactorial(x);
	}

	private static double logFactorial(double k)
	{
		if (k <= 1)
			return 0.0;
		return Gamma.logGamma(k + 1);
	}

	/**
	 * Get the Poisson log likelihood of value x given the mean. The mean must be strictly positive. x must be positive.
	 *
	 * @param u
	 *            the mean
	 * @param x
	 *            the x
	 * @return the log likelihood
	 */
	public static double logLikelihood(double[] u, double[] x)
	{
		double ll = 0.0;
		for (int i = u.length; i-- > 0;)
			ll += logLikelihood(u[i], x[i]);
		return ll;
	}

	/**
	 * Get the Poisson likelihood of value x given the mean. The mean must be strictly positive. x must be positive.
	 *
	 * @param u
	 *            the mean
	 * @param x
	 *            the x
	 * @return the likelihood
	 */
	public static double likelihood(double u, double x)
	{
		//return Math.pow(u, x) * FastMath.exp(-u) / factorial(x);
		return FastMath.exp(logLikelihood(u, x));
	}

	/**
	 * Get the Poisson likelihood of value x given the mean. The mean must be strictly positive. x must be positive.
	 *
	 * @param u
	 *            the mean
	 * @param x
	 *            the x
	 * @return the likelihood
	 */
	public static double likelihood(double[] u, double[] x)
	{
		return FastMath.exp(logLikelihood(u, x));
	}

	/**
	 * Get the Poisson maximum log likelihood of value x given the mean is value x. x must be positive.
	 *
	 * @param x
	 *            the x
	 * @return the log likelihood
	 */
	public static double maximumLogLikelihood(double x)
	{
		return (x > 0.0) ? logLikelihood(x, x) : 0.0;
	}

	/**
	 * Get the Poisson maximum log likelihood of value x given the mean is value x. x must be positive.
	 *
	 * @param x
	 *            the x
	 * @return the log likelihood
	 */
	public static double maximumLogLikelihood(double[] x)
	{
		double ll = 0.0;
		for (int i = x.length; i-- > 0;)
			ll += maximumLogLikelihood(x[i]);
		return ll;
	}

	/**
	 * Get the Poisson maximum likelihood of value x given the mean is value x. x must be positive.
	 *
	 * @param x
	 *            the x
	 * @return the likelihood
	 */
	public static double maximumLikelihood(double x)
	{
		return (x > 0.0) ? likelihood(x, x) : 1;
	}

	/**
	 * Get the Poisson maximum likelihood of value x given the mean is value x. x must be positive.
	 *
	 * @param x
	 *            the x
	 * @return the likelihood
	 */
	public static double maximumLikelihood(double[] x)
	{
		return FastMath.exp(maximumLogLikelihood(x));
	}

	/**
	 * Get the Poisson log likelihood ratio of value x given the mean. The mean must be strictly positive. x must be
	 * positive.
	 *
	 * @param u
	 *            the mean
	 * @param x
	 *            the x
	 * @return the log likelihood ratio
	 */
	public static double logLikelihoodRatio(double[] u, double[] x)
	{
		// From https://en.wikipedia.org/wiki/Likelihood-ratio_test#Use:
		// LLR = -2 * [ ln(likelihood for alternative model) - ln(likelihood for null model)]
		// The model with more parameters (here alternative) will always fit at least as well—
		// i.e., have the same or greater log-likelihood—than the model with fewer parameters 
		// (here null)

		double ll = 0.0;
		for (int i = u.length; i-- > 0;)
		{
			//ll += logLikelihood(u[i], x[i]) - maximumLogLikelihood(x[i]);

			if (x[i] > 0.0)
			{
				//ll += (x[i] * Math.log(u[i]) - u[i]) - (x[i] * Math.log(x[i]) - x[i]);
				//ll += x[i] * Math.log(u[i]) - u[i] - x[i] * Math.log(x[i]) + x[i];
				//ll += x[i] * (Math.log(u[i]) - Math.log(x[i])) - u[i] + x[i];
				ll += x[i] * Math.log(u[i] / x[i]) - u[i] + x[i];
			}
			else
			{
				ll -= u[i];
			}
		}
		return -2.0 * ll;
	}

	/**
	 * Compute the p-value of the log-likelihood ratio using Wilk's theorem that the ratio asymptotically approaches a
	 * Chi-squared distribution with degrees of freedom equal to the difference in dimensionality of the two models
	 * (alternative and null).
	 *
	 * @param logLikelihoodRatio
	 *            the log-likelihood ratio
	 * @param degreesOfFreedom
	 *            the degrees of freedom
	 * @return the p-value
	 * @see https://en.wikipedia.org/wiki/Likelihood-ratio_test#Wilks.27_theorem
	 */
	public static double computePValue(double logLikelihoodRatio, int degreesOfFreedom)
	{
		// The ChiSquaredDistribution just wraps the gamma distribution for this function
		//return new ChiSquaredDistribution(degreesOfFreedom).cumulativeProbability(logLikelihoodRatio);
		//return new GammaDistribution(null, degreesOfFreedom / 2.0, 2.0).cumulativeProbability(logLikelihoodRatio);
		if (logLikelihoodRatio <= 0)
		{
			return 0;
		}
		else
		{
			return Gamma.regularizedGammaP(degreesOfFreedom / 2.0, logLikelihoodRatio / 2.0);
		}
	}
}