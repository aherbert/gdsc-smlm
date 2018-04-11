package gdsc.smlm.function;

import org.apache.commons.math3.util.FastMath;

import gdsc.core.utils.Maths;

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
 * The implementation uses full convolution described from Huang, et al (2013), Supplementary Notes Eq 1.1:<br/>
 * P(D) = A Sum_q e^-u * u^q / q! * 1/sqrt(2pi var) * e ^ -((D-q*g)^2 / 2*var)<br/>
 * Where:<br/>
 * A = normalisation constant
 * var = the variance of the pixel <br/>
 * g = the gain of the pixel <br/>
 * u = the function value (expected number of photons) <br/>
 * D = the observed value at the pixel
 * <p>
 * The likelihood function is designed to model on-chip amplification of a EMCCD/CCD/sCMOS camera which captures a
 * Poisson process of emitted light, converted to electrons on the camera chip, amplified by a gain and then read with
 * Gaussian noise.
 */
public class PoissonGaussianConvolutionFunction implements LikelihoodFunction, LogLikelihoodFunction
{
	private static final LogFactorial logFactorial = new LogFactorial();

	/**
	 * The on-chip gain multiplication factor
	 */
	final double g;

	private final double var;
	private final double s;
	private final double var_by_2;

	private final double logNormalisationGaussian;

	/**
	 * @param alpha
	 *            The inverse of the on-chip gain multiplication factor
	 * @param sigmasquared
	 *            The variance of the Gaussian distribution at readout (must be positive)
	 */
	private PoissonGaussianConvolutionFunction(double alpha, double sigmasquared)
	{
		if (sigmasquared <= 0)
			throw new IllegalArgumentException("Gaussian variance must be strictly positive");
		alpha = Math.abs(alpha);

		this.g = 1.0 / alpha;
		this.var = sigmasquared;
		s = Math.sqrt(var);
		var_by_2 = var * 2;

		// Determine the normalisation factor A in the event that the probability 
		// distribution is being used as a discrete distribution.
		logNormalisationGaussian = PoissonGaussianFunction.getLogNormalisation(var);
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
	 *             if the variance is zero or below
	 */
	public static PoissonGaussianConvolutionFunction createWithStandardDeviation(final double alpha, final double s)
	{
		return new PoissonGaussianConvolutionFunction(alpha, s * s);
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
	 *             if the variance is zero or below
	 */
	public static PoissonGaussianConvolutionFunction createWithVariance(final double alpha, final double var)
	{
		return new PoissonGaussianConvolutionFunction(alpha, var);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.LikelihoodFunction#likelihood(double, double)
	 */
	public double likelihood(double o, double e)
	{
		if (e <= 0)
		{
			// If no Poisson mean then just use the Gaussian
			return FastMath.exp((-0.5 * o * o / var) + logNormalisationGaussian);
		}
		else
		{
			// Use same nomenclature as Huang et al

			final double u = e / g; // expected photoelectrons
			final double D = o; // Camera counts
			// g == gain
			// var = readout variance

			// This is the probability of a Poisson convolved with a Gaussian.
			// Evaluate the Poisson only in the range where the Gaussian is significant.
			// I.e. when the Gaussian probability is zero then the Poisson is not relevant.
			// Use +/- 5 SD
			// x = D - q*g => q = (D-x) / g
			int qmax = (int) Math.ceil((D + 5 * s) / g);
			if (qmax < 0)
				return 0;
			int qmin = (int) Math.floor((D - 5 * s) / g);
			if (qmin < 0)
				qmin = 0;

			// Note: If D is camera counts then it will likely be limited to a 16-bit range
			// Assuming the gain is at least 1 then the max q is:
			// 65536 + 5 * s => This is an acceptable table size to pre-compute the log 
			// factorial if s is reasonable. 

			logFactorial.ensureRange(qmin, qmax);

			final double logu = Math.log(u);
			double p = 0;
			for (int q = qmin; q <= qmax; q++)
			{
				// P(D|q) = e^-u * u^q / q! * 1/sqrt(2pi var) * e ^ -((D-q*g)^2 / 2*var)
				// log(P(D|q) = -u + q * log(u) - log(q!) - (D-q*g)^2/2*var - log(sqrt(2pi var))

				//// Poisson
				//double pp = q * logu - u - LogFactorial.logF(q);
				//
				//// Gaussian
				//double gp = -(Maths.pow2(D - q * g) / var_by_2) + logNormalisationGaussian;
				//
				////System.out.printf("D=%f,q=%d,pp=%g,gp=%g  %g\n", D, q, FastMath.exp(pp), FastMath.exp(gp)
				////		, FastMath.exp(-(Maths.pow2(D - q * g) / var_by_2)) / Math.sqrt(Math.PI*var_by_2));
				//
				//// Combine
				//p += FastMath.exp(pp + gp);

				p += FastMath.exp(
						// Poisson
						q * logu - u - LogFactorial.logF(q)
						// Gaussian
								- (Maths.pow2(D - q * g) / var_by_2) + logNormalisationGaussian);
			}

			// Determine normalisation
			// Note: This is needed when using this as a discrete probability distribution, 
			// e.g. input observed count is integer

			return p;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.LogLikelihoodFunction#logLikelihood(double, double)
	 */
	public double logLikelihood(double o, double e)
	{
		if (e <= 0)
		{
			// If no Poisson mean then just use the Gaussian
			return (-0.5 * o * o / var) + logNormalisationGaussian;
		}
		else
		{
			// Use same nomenclature as Huang et al

			final double u = e / g; // expected photoelectrons
			final double D = o; // Camera counts
			// g == gain
			// var = readout variance

			// This is the probability of a Poisson convolved with a Gaussian.
			// Evaluate the Poisson only in the range where the Gaussian is significant.
			// I.e. when the Gaussian probability is zero then the Poisson is not relevant.
			// Use +/- 5 SD
			// x = D - q*g => q = (D-x) / g
			int qmax = (int) Math.ceil((D + 5 * s) / g);
			if (qmax < 0)
				return Double.NEGATIVE_INFINITY;
			int qmin = (int) Math.floor((D - 5 * s) / g);
			if (qmin < 0)
				qmin = 0;

			logFactorial.ensureRange(qmin, qmax);

			final double logu = Math.log(u);
			double p = 0;
			for (int q = qmin; q <= qmax; q++)
			{
				// P(D|q) = e^-u * u^q / q! * 1/sqrt(2pi var) * e ^ -((D-q*g)^2 / 2*var)
				// log(P(D|q) = -u + q * log(u) - log(q!) - (D-q*g)^2/2*var - log(sqrt(2pi var))

				p += FastMath.exp(
						// Poisson
						q * logu - u - LogFactorial.logF(q)
						// Gaussian
								- (Maths.pow2(D - q * g) / var_by_2) + logNormalisationGaussian);
			}

			// Determine normalisation
			// Note: This is needed when using this as a discrete probability distribution, 
			// e.g. input observed count is integer

			return Math.log(p);
		}
	}
}