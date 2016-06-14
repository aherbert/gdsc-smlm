package gdsc.smlm.function;

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

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.util.FastMath;

/**
 * This is a function to compute the likelihood assuming a Poisson-Gamma-Gaussian distribution.
 * <p>
 * For each observed value the log-likelihood is computed from the Poisson-Gamma-Gaussian distribution (a Poisson
 * convolved with a Gamma distribution convolved with a Gaussian). The Poisson-Gamma distribution is derived
 * analytically in the paper Maximilian H Ulbrich & Ehud Y Isacoff, Nature Methods - 4, 319 - 321 (2007) to explain the
 * probability distribution of ADUs given a fixed photon level per pixel and set gain in an EM-CCD camera (The Poisson
 * distribution models the photon count and the Gamma distribution models the EM-gain). This is then numerically
 * convolved with a Gaussian distribution to model the read noise of the camera.
 * <p>
 * The distribution of Ulbrich & Isacoff has no analytical solution to the convolution with a Gaussian. However the
 * convolution with a Gaussian has the most effect when the counts are low. The Poisson-Gamma-Gaussian can be
 * approximated using the Poisson-Gamma and a partial convolution with a Gaussian at low counts. This method is provided
 * as Python source code within the supplementary information of the paper Mortensen, et al (2010) Nature Methods 7,
 * 377-383. This Java implementation is based on the Python code.
 * <P>
 * The mean of the Poisson distribution is set using the expected value. The scale (EM-gain) for the Gamma distribution
 * and standard deviation of the Gaussian is fixed and set in the constructor. The mean of the Gaussian is assumed to be
 * zero.
 */
public class PoissonGammaGaussianFunction implements LikelihoodFunction
{
	/**
	 * The inverse scale of the Gamma distribution (e.g. the inverse of the EM-gain multiplication factor)
	 */
	final private double alpha;
	/**
	 * The standard deviation of the Gaussian (e.g. Width of the noise distribution in the EMCCD output)
	 */
	final private double sigma;

	private final double twoSigma2;
	private final double sqrt2sigma2;
	private final double sqrt2piSigma2;

	private boolean useApproximation = true;
	private boolean useSimpleIntegration = true;
	private double minimumProbability = Double.MIN_VALUE;

	/**
	 * Initialise the function.
	 * <p>
	 * The input parameters must be the full parameters for the non-linear function. Only those parameters with gradient
	 * indices should be passed in to the functions to obtain the value (and gradient).
	 * <p>
	 * Note: Negative parameters are made absolute
	 * 
	 * @param alpha
	 *            Inverse gain of the EMCCD chip
	 * @param s
	 *            The Gaussian standard deviation
	 */
	public PoissonGammaGaussianFunction(double alpha, double s)
	{
		this.alpha = Math.abs(alpha);
		this.sigma = Math.abs(s);
		twoSigma2 = 2 * s * s;
		sqrt2sigma2 = Math.sqrt(2 * s * s);
		sqrt2piSigma2 = Math.sqrt(2 * Math.PI * s * s);
	}

	private static final double twoSqrtPi = 2 * Math.sqrt(Math.PI);
	private static final double sqrt2pi = Math.sqrt(2 * Math.PI);

	/**
	 * Compute the likelihood
	 * 
	 * @param o
	 *            The observed count
	 * @param e
	 *            The expected count
	 * @return The likelihood
	 */
	public double likelihood(final double o, final double e)
	{
		// Use the same variables as the Python code
		final double cij = o;
		final double eta = alpha * e; // convert to photons

		if (sigma == 0)
		{
			// No convolution with a Gaussian. Simply evaluate for a Poisson-Gamma distribution.
			final double p;

			// Any observed count above zero
			if (cij > 0.0)
			{
				// The observed count converted to photons
				final double nij = alpha * cij;

				// The current implementation of Bessel.I1(x) is Infinity at x==710
				// The limit on eta * nij is therefore (709/2)^2 = 125670.25
				if (eta * nij > 10000)
				{
					// Approximate Bessel function i1(x) when using large x:
					// i1(x) ~ exp(x)/sqrt(2*pi*x)
					// However the entire equation is logged (creating transform),
					// evaluated then raised to e to prevent overflow error on 
					// large exp(x)

					final double transform = 0.5 * Math.log(alpha * eta / cij) - nij - eta + 2 * Math.sqrt(eta * nij) -
							Math.log(twoSqrtPi * Math.pow(eta * nij, 0.25));
					p = FastMath.exp(transform);
				}
				else
				{
					// Second part of equation 135
					p = Math.sqrt(alpha * eta / cij) * FastMath.exp(-nij - eta) *
							Bessel.I1(2 * Math.sqrt(eta * nij));
				}
			}
			else if (cij == 0.0)
			{
				p = FastMath.exp(-eta);
			}
			else
			{
				p = 0;
			}
			
			return (p > minimumProbability) ? p : minimumProbability;			
		}
		else if (useApproximation)
		{
			return mortensenApproximation(cij, eta);
		}
		else
		{
			// This code is the full evaluation of equation 7 from the supplementary information  
			// of the paper Chao, et al (2013) Nature Methods 10, 335-338.
			// It is the full evaluation of a Poisson-Gamma-Gaussian convolution PMF. 

			final double sk = sigma; // Read noise
			final double g = 1.0 / alpha; // Gain
			final double z = o; // Observed pixel value
			final double vk = eta; // Expected number of photons

			// Compute the integral to infinity of:
			// exp( -((z-u)/(sqrt(2)*s)).^2 - u/g ) * sqrt(vk*u/g) .* besseli(1, 2 * sqrt(vk*u/g)) ./ u;

			final double vk_g = vk * alpha; // vk / g
			final double sqrt2sigma = Math.sqrt(2) * sk;

			// Specify the function to integrate
			UnivariateFunction f = new UnivariateFunction()
			{
				public double value(double u)
				{
					return eval(sqrt2sigma, z, vk_g, g, u);
				}
			};

			// Integrate to infinity is not necessary. The convolution of the function with the 
			// Gaussian should be adequately sampled using a nxSD around the maximum.
			// Find a bracket containing the maximum
			double lower, upper;
			double maxU = Math.max(1, o);
			double rLower = maxU;
			double rUpper = maxU + 1;
			double f1 = f.value(rLower);
			double f2 = f.value(rUpper);

			// Calculate the simple integral and the range
			double sum = f1 + f2;
			boolean searchUp = f2 > f1;

			if (searchUp)
			{
				while (f2 > f1)
				{
					f1 = f2;
					rUpper += 1;
					f2 = f.value(rUpper);
					sum += f2;
				}
				maxU = rUpper - 1;
			}
			else
			{
				// Ensure that u stays above zero
				while (f1 > f2 && rLower > 1)
				{
					f2 = f1;
					rLower -= 1;
					f1 = f.value(rLower);
					sum += f1;
				}
				maxU = (rLower > 1) ? rLower + 1 : rLower;
			}

			lower = Math.max(0, maxU - 5 * sk);
			upper = maxU + 5 * sk;

			if (useSimpleIntegration && lower > 0)
			{
				// If we are not at the zero boundary then we can use a simple integration by adding the 
				// remaining points in the range
				for (double u = rLower - 1; u >= lower; u -= 1)
				{
					sum += f.value(u);
				}
				for (double u = rUpper + 1; u <= upper; u += 1)
				{
					sum += f.value(u);
				}
			}
			else
			{
				// Use Legendre-Gauss integrator
				try
				{
					final double relativeAccuracy = 1e-4;
					final double absoluteAccuracy = 1e-8;
					final int minimalIterationCount = 3;
					final int maximalIterationCount = 32;
					final int integrationPoints = 16;

					// Use an integrator that does not use the boundary since u=0 is undefined (divide by zero)
					UnivariateIntegrator i = new IterativeLegendreGaussIntegrator(integrationPoints, relativeAccuracy,
							absoluteAccuracy, minimalIterationCount, maximalIterationCount);

					sum = i.integrate(2000, f, lower, upper);
				}
				catch (TooManyEvaluationsException ex)
				{
					return mortensenApproximation(cij, eta);
				}
			}

			// Compute the final probability
			//final double 
			f1 = z / sqrt2sigma;
			final double p = (FastMath.exp(-vk) / (sqrt2pi * sk)) * (FastMath.exp(-(f1 * f1)) + sum);
			return (p > minimumProbability) ? p : minimumProbability;			
		}
	}

	//private static double pMinObserved = 1;

	private double mortensenApproximation(final double cij, final double eta)
	{
		// This code is adapted from the Python source code within the supplementary information of 
		// the paper Mortensen, et al (2010) Nature Methods 7, 377-383.

		// [Poisson PMF] multiplied by the [value at zero]:
		// [(eta^0 / 0!) * FastMath.exp(-eta)] * [eta * alpha]
		// FastMath.exp(-eta) * [eta * alpha]
		final double f0 = alpha * FastMath.exp(-eta) * eta;

		// ?
		final double fp0 = f0 * 0.5 * alpha * (eta - 2);

		// The cumulative normal distribution of the read noise
		// at the observed count
		final double conv0 = 0.5 * (1 + Erf.erf(cij / (sqrt2sigma2)));

		// [Noise * Gaussian PMF at observed count] + 
		//  [observed count * cumulative distribution of read noise at observed count]
		// [sigma*FastMath.exp(-cij**2/(twoSigma2))/Math.sqrt(2*pi)] + [cij*conv0]
		final double conv1 = sigma * FastMath.exp(-(cij * cij) / twoSigma2) / sqrt2pi + cij * conv0;

		// ? 
		double temp = (f0 * conv0 + fp0 * conv1 + FastMath.exp(-eta) * gauss(cij));

		if (cij > 0.0)
		{
			// The observed count converted to photons
			final double nij = alpha * cij;

			if (eta * nij > 10000)
			{
				// Approximate Bessel function i1(x) when using large x:
				// i1(x) ~ exp(x)/sqrt(2*pi*x)
				// However the entire equation is logged (creating transform),
				// evaluated then raised to e to prevent overflow error on 
				// large exp(x)

				final double transform = 0.5 * Math.log(alpha * eta / cij) - nij - eta + 2 * Math.sqrt(eta * nij) -
						Math.log(twoSqrtPi * Math.pow(eta * nij, 0.25));
				temp += (FastMath.exp(transform) - f0 - fp0 * cij);
			}
			else
			{
				// Second part of equation 135 but not sure what the 
				// -f0-fp0*cij term is.
				// This indicates that temp should already be the first
				// part of eq.135: exp(-eta)*delta(cij)
				temp += (Math.sqrt(alpha * eta / cij) * FastMath.exp(-nij - eta) * Bessel.I1(2 * Math.sqrt(eta * nij)) -
						f0 - fp0 * cij);
			}
		}

		// XXX : Debugging: Store the smallest likelihood we ever see. 
		// This can be used to set a limit for the likelihood
		//if (pMinObserved > temp && temp > 0)
		//{
		//	pMinObserved = temp;
		//}

		return (temp > minimumProbability) ? temp : minimumProbability;
	}

	private double eval(double sqrt2sigma, double z, double vk_g, double g, double u)
	{
		final double f1 = (z - u) / sqrt2sigma;
		final double f2 = Math.sqrt(vk_g * u);
		return FastMath.exp(-(f1 * f1) - u / g) * f2 * Bessel.I1(2 * f2) / u;
	}

	/**
	 * This code is adapted from the Python source code within the supplementary information of the paper Mortensen, et
	 * al (2010) Nature Methods 7, 377-383.
	 * <p>
	 * Note this will return Double.NEGATIVE_INFINITY if there is no Gaussian standard deviation and the observed count
	 * is below zero (since the likelihood is zero).
	 * 
	 * @param o
	 *            The observed count
	 * @param e
	 *            The expected count
	 * @return The log-likelihood
	 */
	public double logLikelihood(final double o, final double e)
	{
		return Math.log(likelihood(o, e));
	}

	private double gauss(final double x)
	{
		return FastMath.exp(-(x * x) / twoSigma2) / sqrt2piSigma2;
	}

	/**
	 * @return the useApproximation
	 */
	public boolean isUseApproximation()
	{
		return useApproximation;
	}

	/**
	 * @param useApproximation
	 *            the useApproximation to set
	 */
	public void setUseApproximation(boolean useApproximation)
	{
		this.useApproximation = useApproximation;
	}

	/**
	 * @return the useSimpleIntegration
	 */
	public boolean isUseSimpleIntegration()
	{
		return useSimpleIntegration;
	}

	/**
	 * @param useSimpleIntegration
	 *            the useSimpleIntegration to set
	 */
	public void setUseSimpleIntegration(boolean useSimpleIntegration)
	{
		this.useSimpleIntegration = useSimpleIntegration;
	}

	/**
	 * @return the alpha
	 */
	public double getAlpha()
	{
		return alpha;
	}

	/**
	 * @return the sigma
	 */
	public double getSigma()
	{
		return sigma;
	}

	/**
	 * Gets the minimum probability that will ever be returned. Setting this above zero allows the use of Math.log() on
	 * the likelihood value.
	 *
	 * @return the minimum probability
	 */
	public double getMinimumProbability()
	{
		return minimumProbability;
	}

	/**
	 * Sets the minimum probability that will ever be returned. Setting this above zero allows the use of Math.log() on
	 * the likelihood value.
	 *
	 * @param p
	 *            the new minimum probability
	 */
	public void setMinimumProbability(double p)
	{
		this.minimumProbability = p;
	}
}