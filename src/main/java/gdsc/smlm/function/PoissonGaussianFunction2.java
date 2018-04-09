package gdsc.smlm.function;

import org.apache.commons.math3.util.FastMath;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Implements the probability density function for a Poisson-Gaussian Mixture. The Gaussian is assumed to have mean of
 * zero. If no mean (zero or below) is provided for the Poisson distribution then the probability density function
 * matches that of the Gaussian.
 * <p>
 * The implementation uses the saddle-point approximation described in Snyder, et al (1995) Compensation for readout
 * noise in CCD images. J.Opt. Soc. Am. 12, 272-283. The method is adapted from the C source code provided in the
 * appendix.
 * <p>
 * This is just a wrapper for the PoissonGaussianFunction that handles smart switching between a Gaussian likelihood
 * function and a Poisson-Gaussian likelihood function.
 * <p>
 * The likelihood function is designed to model on-chip amplification of a EMCCD/CCD/sCMOS camera which captures a
 * Poisson process of emitted light, converted to electrons on the camera chip, amplified by a gain and then read with
 * Gaussian noise.
 */
public class PoissonGaussianFunction2 implements LikelihoodFunction, LogLikelihoodFunction
{
	/**
	 * The inverse of the on-chip gain multiplication factor
	 */
	final double alpha;

	private boolean usePicardApproximation = false;
	private final double sigmasquared;

	private final double probabilityNormalisation;
	private final double logNormalisation;
	private final double probabilityNormalisationNoPoisson;
	private final double logNormalisationNoPoisson;

	/**
	 * @param alpha
	 *            The inverse of the on-chip gain multiplication factor
	 * @param sigmasquared
	 *            The variance of the Gaussian distribution at readout (must be positive)
	 */
	private PoissonGaussianFunction2(double alpha, double sigmasquared)
	{
		if (sigmasquared <= 0)
			throw new IllegalArgumentException("Gaussian variance must be strictly positive");
		alpha = Math.abs(alpha);
		
		// Apply gain to the readout standard deviation. 
		// This compresses the probability distribution by alpha. Thus we can compute the
		// probability using a Poisson or Poisson-Gaussian mixture and then compress the
		// output probability so the cumulative probability is 1 over the uncompressed range.
		sigmasquared *= (alpha * alpha); 
		
		this.alpha = alpha;
		this.sigmasquared = sigmasquared;

		// As per PoissonGaussianFunction
		probabilityNormalisation = PoissonGaussianFunction.NORMALISATION * alpha;
		logNormalisation = PoissonGaussianFunction.LOG_NORMALISATION + Math.log(alpha);
		probabilityNormalisationNoPoisson = PoissonGaussianFunction.getProbabilityNormalisation(sigmasquared) * alpha;
		logNormalisationNoPoisson = PoissonGaussianFunction.getLogNormalisation(sigmasquared) + Math.log(alpha);
	}

	/**
	 * Creates the with standard deviation.
	 *
	 * @param alpha
	 *            The inverse of the on-chip gain multiplication factor
	 * @param s
	 *            The standard deviation of the Gaussian distribution at readout
	 * @return the poisson gaussian function 2
	 * @throws IllegalArgumentException
	 *             if the mean or variance is zero or below
	 */
	public static PoissonGaussianFunction2 createWithStandardDeviation(final double alpha, final double s)
	{
		return new PoissonGaussianFunction2(alpha, s * s);
	}

	/**
	 * Creates the with variance.
	 *
	 * @param alpha
	 *            The inverse of the on-chip gain multiplication factor
	 * @param var
	 *            The variance of the Gaussian distribution at readout (must be positive)
	 * @return the poisson gaussian function 2
	 * @throws IllegalArgumentException
	 *             if the mean or variance is zero or below
	 */
	public static PoissonGaussianFunction2 createWithVariance(final double alpha, final double var)
	{
		return new PoissonGaussianFunction2(alpha, var);
	}

	/**
	 * Return if using the Picard approximation for the initial saddle point
	 * 
	 * @return True if using the Picard approximation
	 */
	public boolean isUsePicardApproximation()
	{
		return usePicardApproximation;
	}

	/**
	 * Specify whether to use the Picard approximation for the initial saddle point. The alternative is Pade.
	 * 
	 * @param usePicardApproximation
	 *            True to use the Picard approximation
	 */
	public void setUsePicardApproximation(boolean usePicardApproximation)
	{
		this.usePicardApproximation = usePicardApproximation;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.LikelihoodFunction#likelihood(double, double)
	 */
	public double likelihood(double o, double e)
	{
		// convert to photons
		o *= alpha;
		if (e <= 0)
		{
			// If no Poisson mean then just use the Gaussian
			return FastMath.exp(-0.5 * o * o / sigmasquared) * probabilityNormalisationNoPoisson;
		}
		else
		{
			e *= alpha;
			double saddlepoint = (usePicardApproximation) ? PoissonGaussianFunction.picard(o, e, sigmasquared)
					: PoissonGaussianFunction.pade(o, e, sigmasquared);
			saddlepoint = PoissonGaussianFunction.newton_iteration(o, e, sigmasquared, saddlepoint);
			final double logP = PoissonGaussianFunction.sp_approx(o, e, sigmasquared, saddlepoint);
			return FastMath.exp(logP) * probabilityNormalisation;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.LogLikelihoodFunction#logLikelihood(double, double)
	 */
	public double logLikelihood(double o, double e)
	{
		// convert to photons
		o *= alpha;
		if (e <= 0)
		{
			// If no Poisson mean then just use the Gaussian
			return (-0.5 * o * o / sigmasquared) + logNormalisationNoPoisson;
		}
		else
		{
			e *= alpha;
			double saddlepoint = (usePicardApproximation) ? PoissonGaussianFunction.picard(o, e, sigmasquared)
					: PoissonGaussianFunction.pade(o, e, sigmasquared);
			saddlepoint = PoissonGaussianFunction.newton_iteration(o, e, sigmasquared, saddlepoint);
			final double logP = PoissonGaussianFunction.sp_approx(o, e, sigmasquared, saddlepoint);
			return logP + logNormalisation;
		}
	}
}