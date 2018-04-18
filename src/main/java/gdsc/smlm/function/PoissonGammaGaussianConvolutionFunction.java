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
 * Implements the probability density function for a Poisson-Gamma-Gaussian Mixture. The Gaussian is assumed to have
 * mean of
 * zero. If no mean (zero or below) is provided for the Poisson distribution then the probability density function
 * matches that of the Gaussian.
 * <p>
 * The implementation uses full convolution described from Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI
 * equation 3:<br/>
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
public class PoissonGammaGaussianConvolutionFunction implements LikelihoodFunction, LogLikelihoodFunction
{
	/**
	 * The on-chip gain multiplication factor
	 */
	final double g;

	private final double var;
	private final double s;
	private final double range;
	private final double var_by_2;

	private final double logNormalisationGaussian;

	/**
	 * Instantiates a new poisson gaussian convolution function.
	 *
	 * @param alpha
	 *            The inverse of the on-chip gain multiplication factor
	 * @param variance
	 *            The variance of the Gaussian distribution at readout (must be positive)
	 * @param isVariance
	 *            Set to true if the input parameter is variance; otherwise is is standard deviation
	 */
	private PoissonGammaGaussianConvolutionFunction(double alpha, double variance, boolean isVariance)
	{
		if (variance <= 0)
			throw new IllegalArgumentException("Gaussian variance must be strictly positive");
		alpha = Math.abs(alpha);

		this.g = 1.0 / alpha;
		if (isVariance)
		{
			s = Math.sqrt(variance);
			this.var = variance;
		}
		else
		{
			s = variance;
			this.var = s * s;
		}
		var_by_2 = var * 2;

		// Use a range to cover the Gaussian convolution
		range = 5 * this.s;

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
	public static PoissonGammaGaussianConvolutionFunction createWithStandardDeviation(final double alpha,
			final double s)
	{
		return new PoissonGammaGaussianConvolutionFunction(alpha, s, false);
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
	public static PoissonGammaGaussianConvolutionFunction createWithVariance(final double alpha, final double var)
	{
		return new PoissonGammaGaussianConvolutionFunction(alpha, var, true);
	}

	private static final double twoPi = 2 * Math.PI;

	/**
	 * Calculate the probability density function for a Poisson-Gamma distribution model of EM-gain.
	 * <p>
	 * See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
	 * <p>
	 * Note: This implementation will underestimate the cumulative probability (sum<1) when the mean is close to 1 and
	 * the gain is low (<10).
	 * 
	 * @param c
	 *            The count to evaluate
	 * @param p
	 *            The average number of photons per pixel input to the EM-camera (must be positive)
	 * @param m
	 *            The multiplication factor (gain)
	 * @return The probability
	 */
	public static double poissonGamma(double c, double p, double m)
	{
		// The default evaluation is:
		//if (true)
		//return Math.sqrt(p / (c * m)) * FastMath.exp(-c / m - p) * Bessel.I1(2 * Math.sqrt(c * p / m));

		// However 
		// Bessel.I1(x) -> Infinity
		// The current implementation of Bessel.I1(x) is Infinity at x==710 so we switch 
		// to an approximation...

		// Any observed count above zero
		if (c > 0.0)
		{
			// The observed count converted to photons
			final double c_m = c / m;
			final double cp_m = p * c_m;

			// The current implementation of Bessel.I1(x) is Infinity at x==710
			// The limit on (c * p / m) is therefore (709/2)^2 = 125670.25
			final double x = 2 * Math.sqrt(cp_m);
			if (cp_m > 10000)
			{
				// Approximate Bessel function i1(x) when using large x:
				// i1(x) ~ exp(x)/sqrt(2*pi*x)
				// However the entire equation is logged (creating transform),
				// evaluated then raised to e to prevent overflow error on 
				// large exp(x)

				// p = sqrt(p / (c * m)) * exp(-c_m - p) * exp(2 * sqrt(cp_m)) / sqrt(2*pi*2*sqrt(cp_m))
				// p = sqrt(p / (c * m)) * exp(-c_m - p) * exp(x) / sqrt(2*pi*x)
				// log(p) = 0.5 * log(p / (c * m)) - c_m - p + x - 0.5 * log(2*pi*x)

				// This is the transform from the Python source code within the supplementary information of 
				// the paper Mortensen, et al (2010) Nature Methods 7, 377-383.
				// p = sqrt(p / (c * m)) * exp(-c_m - p) * exp(2 * sqrt(cp_m)) / (sqrt(2*pi)*sqrt(2*sqrt(cp_m)))
				// log(p) = 0.5 * log(p / (c * m)) - c_m - p + 2 * sqrt(cp_m) - log(sqrt(2*pi)*sqrt(2*sqrt(cp_m)))
				// log(p) = 0.5 * log(p / (c * m)) - c_m - p + 2 * sqrt(cp_m) - log(sqrt(2)*sqrt(pi)*sqrt(2)*sqrt(sqrt(cp_m)))
				// log(p) = 0.5 * log(p / (c * m)) - c_m - p + 2 * sqrt(cp_m) - log(2*sqrt(pi)*sqrt(sqrt(cp_m)))

				// This avoids a call to Math.pow 
				final double transform = 0.5 * Math.log(p / (c * m)) - c_m - p + x - 0.5 * Math.log(twoPi * x);

				//final double transform2 = 0.5 * Math.log(p / (c * m)) - c_m - p + x -
				//		Math.log(2 * Math.sqrt(Math.PI) * Math.pow(p * c_m, 0.25));
				//System.out.printf("t1=%g, t2=%g error=%g\n", transform, transform2,
				//		gdsc.core.utils.DoubleEquality.relativeError(transform, transform2));

				return FastMath.exp(transform);
			}
			else
			{
				return Math.sqrt(p / (c * m)) * FastMath.exp(-c_m - p) * Bessel.I1(x);
			}
		}
		else if (c == 0.0)
		{
			return FastMath.exp(-p);
		}
		else
		{
			return 0;
		}
	}

	/**
	 * Calculate the log probability density function for a Poisson-Gamma distribution model of EM-gain.
	 * <p>
	 * See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
	 * <p>
	 * Note: This implementation will underestimate the cumulative probability (sum<1) when the mean is close to 1 and
	 * the gain is low (<10).
	 * 
	 * @param c
	 *            The count to evaluate
	 * @param p
	 *            The average number of photons per pixel input to the EM-camera (must be positive)
	 * @param m
	 *            The multiplication factor (gain)
	 * @return The log probability
	 */
	public static double logPoissonGamma(double c, double p, double m)
	{
		// As above without final exp
		if (c > 0.0)
		{
			final double c_m = c / m;
			final double cp_m = p * c_m;
			final double x = 2 * Math.sqrt(cp_m);
			if (cp_m > 10000)
			{
				return 0.5 * Math.log(p / (c * m)) - c_m - p + x - 0.5 * Math.log(twoPi * x);
			}
			else
			{
				return 0.5 * Math.log(p / (c * m)) - c_m - p + Math.log(Bessel.I1(x));
			}
		}
		else if (c == 0.0)
		{
			return -p;
		}
		else
		{
			return Double.NEGATIVE_INFINITY;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.LikelihoodFunction#likelihood(double, double)
	 */
	public double likelihood(final double o, final double e)
	{
		if (e <= 0)
		{
			// If no Poisson mean then just use the Gaussian
			return FastMath.exp((-o * o / var_by_2) + logNormalisationGaussian);
		}
		else
		{
			// Note:
			// This is a convolution of two continuous probability distributions.
			// It does not compute a good estimate when the variance is small since the 
			// single point approximation to the gaussian is not valid. This is not
			// relevant for a EM-CCD since the variance is likely to be above 10 counts.
			// It also underestimates the cumulative distribution (sum < 1) when the Poisson 
			// mean is close to 1 or the gain is small (<4) due to underestimation in the 
			// Poisson-Gamma distribution. 

			// Use a range to cover the Gaussian convolution
			double max = o + range;
			if (max < 0)
				return 0;
			double min = o - range;
			if (min < 0)
				min = 0;

			return computeP(o, e, max, min);
		}
	}

	private double computeP(final double o, final double e, double max, double min)
	{
		double p = 0;

		// Overcome the problem with small variance using a set number of steps to 
		// cover the range. This effectively makes the Poisson-Gamma a continuous
		// probability distribution.
		// Note:
		// This does seem to be valid. The Poisson-Gamma is a discrete PMF.
		// The CameraModelAnalysis plugin works with full integration if this function
		// is computed using integer steps.

		//		if (s < 0)
		//		{
		//			double step = (max - min) / 10;
		//
		//			for (int i = 0; i <= 10; i++)
		//			{
		//				double c = min + i * step;
		//				p += FastMath.exp(
		//						// Poisson-Gamma
		//						logPoissonGamma(c, e, g)
		//						// Gaussian
		//								- (Maths.pow2(c - o) / var_by_2) + logNormalisationGaussian);
		//			}
		//			p *= step;
		//		}
		//		else
		//		{

		int cmax = (int) Math.ceil(max);
		int cmin = (int) Math.floor(min);

		for (int c = cmin; c <= cmax; c++)
		{
			p += FastMath.exp(
					// Poisson-Gamma
					logPoissonGamma(c, e, g)
					// Gaussian
							- (Maths.pow2(c - o) / var_by_2) + logNormalisationGaussian);
		}

		//		}

		return p;
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
			return (-o * o / var_by_2) + logNormalisationGaussian;
		}
		else
		{
			double max = o + range;
			if (max < 0)
				return Double.NEGATIVE_INFINITY;
			double min = o - range;
			if (min < 0)
				min = 0;

			return Math.log(computeP(o, e, max, min));
		}
	}
}