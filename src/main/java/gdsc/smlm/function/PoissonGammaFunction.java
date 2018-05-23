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
 * Implements the probability density function for a Poisson-Gamma Mixture.
 * <p>
 * The implementation uses the Poisson-Gamma mixture described from Ulbrich & Isacoff (2007). Nature Methods 4, 319-321,
 * SI equation 3:<br/>
 * Gp,m(0) = e^-p<br/>
 * Gp,m(c|c>0) = sqrt(p/(c*m)) * e^(-c/m -p) * I1(2*sqrt(c*p/m))<br/>
 * Where:<br/>
 * c = the observed value at the pixel <br/>
 * p = the function value (expected number of photons) <br/>
 * m = the gain of the pixel <br/>
 * I1 = Modified Bessel function of the first kind <br/>
 * <p>
 * The likelihood function is designed to model on-chip amplification of a EMCCD/CCD/sCMOS camera which captures a
 * Poisson process of emitted light, converted to electrons on the camera chip, amplified by a gain and then read with
 * Gaussian noise.
 */
public class PoissonGammaFunction implements LikelihoodFunction, LogLikelihoodFunction, GradientLikelihoodFunction
{
	/**
	 * The on-chip gain multiplication factor
	 */
	final double m;

	/**
	 * Instantiates a new poisson gamma function.
	 *
	 * @param m
	 *            The on-chip gain multiplication factor
	 * @throws IllegalArgumentException
	 *             if the gain is zero or below
	 */
	public PoissonGammaFunction(double m)
	{
		if (!(m > 0))
			throw new IllegalArgumentException("Gain must be strictly positive");
		this.m = m;
	}

	/**
	 * Creates the with standard deviation.
	 *
	 * @param alpha
	 *            The inverse of the on-chip gain multiplication factor
	 * @return the poisson gamma function
	 * @throws IllegalArgumentException
	 *             if the gain is zero or below
	 */
	public static PoissonGammaFunction createWithAlpha(final double alpha)
	{
		return new PoissonGammaFunction(1.0 / alpha);
	}

	private static final double twoPi = 2 * Math.PI;

	/**
	 * Calculate the probability density function for a Poisson-Gamma distribution model of EM-gain.
	 * <p>
	 * See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
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
		// Any observed count above zero
		if (c > 0.0)
		{
			// The observed count converted to photons
			final double c_m = c / m;

			// The current implementation of Bessel.II(x) is Infinity at x==710. 
			// Also exp(-c/m -p) will be sub-normal at < -709. 
			// Switch to an approximation.
			final double x = 2 * Math.sqrt(p * c_m);
			final double _c_m_p = -c_m - p;
			if (x > 709 || _c_m_p < -709)
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
				//final double transform = 0.5 * Math.log(p / (c * m)) + _c_m_p + x - 0.5 * Math.log(twoPi * x);

				//final double transform2 = 0.5 * Math.log(p / (c * m)) - c_m - p + x -
				//		Math.log(2 * Math.sqrt(Math.PI) * Math.pow(p * c_m, 0.25));
				//System.out.printf("t1=%g, t2=%g error=%g\n", transform, transform2,
				//		gdsc.core.utils.DoubleEquality.relativeError(transform, transform2));

				//return FastMath.exp(0.5 * Math.log(p / (c * m)) + _c_m_p + x - 0.5 * Math.log(twoPi * x));
				return (x / (2 * c)) * FastMath.exp(_c_m_p + x - 0.5 * Math.log(twoPi * x));
			}
			else
			{
				//return Math.sqrt(p / (c * m)) * FastMath.exp(_c_m_p) * Bessel.I1(x);
				return (x / (2 * c)) * FastMath.exp(_c_m_p) * Bessel.I1(x);
			}
		}
		else if (c == 0.0)
		{
			// This is the Dirac delta function plus the probability of 
			// the Poisson-Gamma distribution with shape n=1 at c=0 (reduced to an exponential):
			// Dirac = exp^-p

			// Note:
			// Poisson:
			// 1/n! p^n*e^-p
			// Gamma:
			// 1/((n-1)!m^n) c^(n-1) * e^-c/m

			// If the Gamma takes positive integer arguments it is an Erlang distribution,
			// i.e., the sum of n independent exponentially distributed random variables, 
			// each of which has a mean of p.
			// The Gamma is only non-zero at c==0 when n=1. 
			// Then it is just an exponential distribution.

			// Poisson probability of n=1: FastMath.exp(-p) * p 
			// Gamma probability of c=0 given n=1, Gamma(shape=1,scale=m) = 1 / m

			//System.out.printf("p=%g, m=%g gamma=%g  pp=%g\n", p, m, 
			//		new CustomGammaDistribution(null, 1, m).density(0), 1/m);

			return FastMath.exp(-p) * (1 + p / m);
		}
		else
		{
			return 0;
		}
	}

	/**
	 * Calculate the probability density function for a Poisson-Gamma distribution model of EM-gain for observed Poisson
	 * counts. This avoids the computation of the Dirac delta function at c=0.
	 * <p>
	 * This method is suitable for use in integration routines.
	 * <p>
	 * If c==0 then the true probability is obtained by adding Math.exp(-p).
	 * <p>
	 * See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
	 * 
	 * @param c
	 *            The count to evaluate
	 * @param p
	 *            The average number of photons per pixel input to the EM-camera (must be positive)
	 * @param m
	 *            The multiplication factor (gain)
	 * @return The probability function for observed Poisson counts
	 * @see #poissonGamma(double, double, double)
	 * @see #dirac(double)
	 */
	public static double poissonGammaN(double c, double p, double m)
	{
		// As above with no Dirac delta function at c=0

		if (c > 0.0)
		{
			final double c_m = c / m;
			final double x = 2 * Math.sqrt(p * c_m);
			final double _c_m_p = -c_m - p;
			if (x > 709 || _c_m_p < -709)
			{
				//return FastMath.exp(0.5 * Math.log(p / (c * m)) + _c_m_p + x - 0.5 * Math.log(twoPi * x));
				return (x / (2 * c)) * FastMath.exp(_c_m_p + x - 0.5 * Math.log(twoPi * x));
			}
			else
			{
				//return Math.sqrt(p / (c * m)) * FastMath.exp(_c_m_p) * Bessel.I1(x);
				return (x / (2 * c)) * FastMath.exp(_c_m_p) * Bessel.I1(x);
			}
		}
		else if (c == 0.0)
		{
			// No Dirac delta function
			return FastMath.exp(-p) * p / m;
		}
		else
		{
			return 0;
		}
	}

	/**
	 * Calculate the probability density function for a Poisson-Gamma distribution model of EM-gain for no observed
	 * Poisson counts. This is the Dirac delta function at c=0.
	 * <p>
	 * See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
	 *
	 * @param p
	 *            The average number of photons per pixel input to the EM-camera (must be positive)
	 * @return The probability function for observed Poisson counts
	 * @see #poissonGamma(double, double, double)
	 */
	public static double dirac(double p)
	{
		return FastMath.exp(-p);
	}

	/**
	 * Calculate the probability density function for a Poisson-Gamma distribution model of EM-gain for no observed
	 * Poisson counts. This is the Dirac delta function at c=0.
	 * <p>
	 * See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
	 *
	 * @param p
	 *            The average number of photons per pixel input to the EM-camera (must be positive)
	 * @param dG_dp
	 *            the gradient of the function G(c) with respect to parameter p
	 * @return The probability function for observed Poisson counts
	 * @see #poissonGamma(double, double, double)
	 */
	public static double dirac(double p, double[] dG_dp)
	{
		final double exp_p = FastMath.exp(-p);
		dG_dp[0] = -exp_p;
		return exp_p;
	}

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
	 * @param dG_dp
	 *            the gradient of the function G(c) with respect to parameter p
	 * @return The probability
	 */
	public static double poissonGamma(double c, double p, double m, double[] dG_dp)
	{
		// Any observed count above zero
		if (c > 0.0)
		{
			// The observed count converted to photons
			final double c_m = c / m;
			final double cp_m = p * c_m;

			// The current implementation of Bessel.II(x) is Infinity at x==710
			// due to the use of Math.exp(x). Switch to an approximation.
			final double x = 2 * Math.sqrt(cp_m);
			final double _c_m_p = -c_m - p;
			if (x > 709 || _c_m_p < -709)
			{
				// Approximate Bessel function i0(x)/i1(x) when using large x:
				// In(x) ~ exp(x)/sqrt(2*pi*x)
				double exp_transform = FastMath.exp(_c_m_p + x - 0.5 * Math.log(twoPi * x));
				double G = (x / (2 * c)) * exp_transform;
				dG_dp[0] = exp_transform / m - G;
				return G;
			}
			else
			{
				// G(c) = e^-p . e^-c/m . sum n=1 to inf { 1/(n!(n-1)!) . p^n c^(n-1) / m^n }
				// dG(c)/dp = e^-p . e^-c/m . sum n=1 to inf { 1/(n!(n-1)!) . n * p^(n-1) c^(n-1) / m^n } - e^-p . e^-c/m . sum n=1 to inf { 1/(n!(n-1)!) . p^n c^(n-1) / m^n } 
				// dG(c)/dp = e^-p . e^-c/m . 1/m . sum n=1 to inf { 1/((n-1)!(n-1)!) . p^(n-1) c^(n-1) / m^(n-1) } - G(c) 
				// dG(c)/dp = e^-p . e^-c/m . 1/m . sum n=0 to inf { 1/(n!^2) . (pc/m)^n } - G(c)

				// Bessel I0 = sum n=0 to inf { 1/(n!^2) . ((x/2)^2)^n }
				// x = 2 * sqrt(cp/m)

				// dG(c)/dp = e^-p . e^-c/m . 1/m . I0(2*sqrt(cp/m)) - G(c)
				// dG(c)/dp = e^(-c/m -p) . I0(2*sqrt(cp/m))/m - G(c)

				double exp_c_m_p = FastMath.exp(_c_m_p);
				//double G = Math.sqrt(p / (c * m)) * exp_c_m_p * Bessel.I1(x);
				double G = (x / (2 * c)) * exp_c_m_p * Bessel.I1(x);
				dG_dp[0] = exp_c_m_p * Bessel.I0(x) / m - G;
				return G;
			}
		}
		else if (c == 0.0)
		{
			// f(p) = exp(-p) * (1 + p / m)
			// df/dp = (-exp(-p) * (1 + p / m)) + (exp(-p) / m)   
			double exp_p = FastMath.exp(-p);
			double G = exp_p * (1 + p / m);
			dG_dp[0] = exp_p / m - G;
			return G;
		}
		else
		{
			dG_dp[0] = 0;
			return 0;
		}
	}


	/**
	 * Calculate the probability density function for a Poisson-Gamma distribution model of EM-gain for observed Poisson
	 * counts. This avoids the computation of the Dirac delta function at c=0.
	 * <p>
	 * This method is suitable for use in integration routines.
	 * <p>
	 * If c==0 then the true probability is obtained by adding Math.exp(-p).
	 * <p>
	 * See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
	 * 
	 * @param c
	 *            The count to evaluate
	 * @param p
	 *            The average number of photons per pixel input to the EM-camera (must be positive)
	 * @param m
	 *            The multiplication factor (gain)
	 * @param dG_dp
	 *            the gradient of the function G(c) with respect to parameter p
	 * @return The probability function for observed Poisson counts
	 * @see #poissonGamma(double, double, double)
	 * @see #dirac(double)
	 */
	public static double poissonGammaN(double c, double p, double m, double[] dG_dp)
	{
		// As above with no Dirac delta function at c=0

		if (c > 0.0)
		{
			final double c_m = c / m;
			final double cp_m = p * c_m;
			final double x = 2 * Math.sqrt(cp_m);
			final double _c_m_p = -c_m - p;
			if (x > 709 || _c_m_p < -709)
			{
				double exp_transform = FastMath.exp(_c_m_p + x - 0.5 * Math.log(twoPi * x));
				double G = (x / (2 * c)) * exp_transform;
				dG_dp[0] = exp_transform / m - G;
				return G;
			}
			else
			{
				double exp_c_m_p = FastMath.exp(_c_m_p);
				double G = (x / (2 * c)) * exp_c_m_p * Bessel.I1(x);
				dG_dp[0] = exp_c_m_p * Bessel.I0(x) / m - G;
				return G;
			}
		}
		else if (c == 0.0)
		{
			// No Dirac delta function
			double exp_p_m = FastMath.exp(-p) / m;
			double G = exp_p_m * p;
			dG_dp[0] = exp_p_m - G;
			return G;
		}
		else
		{
			dG_dp[0] = 0;
			return 0;
		}
	}
	
	/**
	 * Calculate the probability density function for a Poisson-Gamma distribution model of EM-gain.
	 * <p>
	 * See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
	 * <p>
	 * This is a special version which computes only part of the gradient.
	 * The partial gradient is equal to the actual gradient plus the value of the function.
	 *
	 * @param c
	 *            The count to evaluate
	 * @param p
	 *            The average number of photons per pixel input to the EM-camera (must be positive)
	 * @param m
	 *            The multiplication factor (gain)
	 * @param dG_dp
	 *            the partial gradient of the function G(c) with respect to parameter p
	 * @return The probability
	 */
	static double poissonGammaPartial(double c, double p, double m, double[] dG_dp)
	{
		// As above but do not subtract the G from the gradient. 
		if (c > 0.0)
		{
			final double c_m = c / m;
			final double cp_m = p * c_m;
			final double x = 2 * Math.sqrt(cp_m);
			final double _c_m_p = -c_m - p;
			if (x > 709 || _c_m_p < -709)
			{
				double exp_transform = FastMath.exp(_c_m_p + x - 0.5 * Math.log(twoPi * x));
				double G = (x / (2 * c)) * exp_transform;
				dG_dp[0] = exp_transform / m;
				return G;
			}
			else
			{
				double exp_c_m_p = FastMath.exp(_c_m_p);
				//double G = Math.sqrt(p / (c * m)) * exp_c_m_p * Bessel.I1(x);
				double G = (x / (2 * c)) * exp_c_m_p * Bessel.I1(x);
				dG_dp[0] = exp_c_m_p * Bessel.I0(x) / m;
				return G;
			}
		}
		else if (c == 0.0)
		{
			double exp_p = FastMath.exp(-p);
			double G = exp_p * (1 + p / m);
			dG_dp[0] = exp_p / m;
			return G;
		}
		else
		{
			dG_dp[0] = 0;
			return 0;
		}
	}

	/**
	 * Calculate the log probability density function for a Poisson-Gamma distribution model of EM-gain.
	 * <p>
	 * See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
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
			if (x > 709)
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
			// log (FastMath.exp(-p) * (1 + p / m))
			return -p + Math.log(1 + p / m);
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
		return poissonGamma(o, e, m);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.LogLikelihoodFunction#logLikelihood(double, double)
	 */
	public double logLikelihood(double o, double e)
	{
		return logPoissonGamma(o, e, m);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.GradientLikelihoodFunction#likelihood(double, double, double[])
	 */
	public double likelihood(double o, double t, double[] dp_dt)
	{
		return poissonGamma(o, t, m, dp_dt);
	}
}