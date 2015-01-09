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
public class PoissonGammaGaussianFunction
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
	 * This code is adapted from the Python source code within the supplementary information of the paper Mortensen, et
	 * al (2010) Nature Methods 7, 377-383.
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

			// Any observed count above zero
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
					return FastMath.exp(transform);
				}
				else
				{
					// Second part of equation 135
					return Math.sqrt(alpha * eta / cij) * FastMath.exp(-nij - eta) *
							Bessel.I1(2 * Math.sqrt(eta * nij));
				}
			}
			else if (cij == 0.0)
			{
				return FastMath.exp(-eta);
			}
			else
			{
				return 0;
			}
		}
		else
		{
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
					temp += (Math.sqrt(alpha * eta / cij) * FastMath.exp(-nij - eta) *
							Bessel.I1(2 * Math.sqrt(eta * nij)) - f0 - fp0 * cij);
				}
			}

			return temp;
		}
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
}