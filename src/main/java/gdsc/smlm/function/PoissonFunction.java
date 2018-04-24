package gdsc.smlm.function;

import org.apache.commons.math3.distribution.CustomPoissonDistribution;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;

/**
 * Implements the probability density function for a Poisson distribution.
 * <p>
 * This is a simple implementation of the LikelihoodFunction interface.
 * <p>
 * The likelihood function is designed to model on-chip amplification of a EMCCD/CCD/sCMOS camera which captures a
 * Poisson process of emitted light, converted to electrons on the camera chip, amplified by a gain and then read.
 */
public class PoissonFunction implements LikelihoodFunction, LogLikelihoodFunction
{
	private final CustomPoissonDistribution pd;
	
	/**
	 * The inverse of the on-chip gain multiplication factor
	 */
	final double alpha;

	/**
	 * The log of the inverse on-chip gain multiplication factor
	 */
	final double logAlpha;

	/**
	 * Allow non-integer observed values
	 */
	final boolean nonInteger;

	/**
	 * @param alpha
	 *            The inverse of the on-chip gain multiplication factor
	 * @param nonInteger
	 *            Allow non-integer observed values
	 */
	public PoissonFunction(double alpha, boolean nonInteger)
	{
		this.alpha = Math.abs(alpha);
		logAlpha = Math.log(alpha);
		this.nonInteger = nonInteger;
		pd = (nonInteger) ? new CustomPoissonDistribution(null, 1) : null;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.LikelihoodFunction#likelihood(double, double)
	 */
	public double likelihood(double o, double e)
	{
		if (o < 0 || e <= 0)
			return 0;

		// convert to photons
		o *= alpha;

		// Allow non-integer observed value using the gamma function to provide a factorial for non-integer values
		// PMF(l,k) = C * e^-l * l^k / gamma(k+1)
		// log(PMF) = -l + k * log(l) - logGamma(k+1)
		if (nonInteger)
		{
			//return (FastMath.exp(-e) * Math.pow(e, o) / factorial(o)) * alpha;

			final double ll = -e + o * Math.log(e) - logFactorial(o);
			return FastMath.exp(ll) * alpha;
		}

		pd.setMeanUnsafe(e);
		return pd.probability((int) o) * alpha;
	}

	/**
	 * Return the log of the factorial for the given real number, using the gamma function
	 * 
	 * @param k
	 * @return the log factorial
	 */
	public static double logFactorial(double k)
	{
		if (k <= 1)
			return 0;
		return Gamma.logGamma(k + 1);
	}

	/**
	 * Return the factorial for the given real number, using the gamma function
	 * 
	 * @param k
	 * @return the factorial
	 */
	public static double factorial(double k)
	{
		if (k <= 1)
			return 1;
		return Gamma.gamma(k + 1);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.LogLikelihoodFunction#logLikelihood(double, double)
	 */
	public double logLikelihood(double o, double e)
	{
		if (o < 0 || e <= 0)
			return Double.NEGATIVE_INFINITY;

		// convert to photons
		o *= alpha;

		// Allow non-integer observed value using the gamma function to provide a factorial for non-integer values
		// PMF(l,k) = C * e^-l * l^k / gamma(k+1)
		// log(PMF) = -l + k * log(l) - logGamma(k+1)
		if (nonInteger)
		{
			//return (FastMath.exp(-e) * Math.pow(e, o) / factorial(o)) * alpha;

			final double ll = -e + o * Math.log(e) - logFactorial(o);
			return ll + logAlpha;
		}

		pd.setMeanUnsafe(e);
		return pd.logProbability((int) o) + logAlpha;
	}
}