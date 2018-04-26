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
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
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
 * <p>
 * The likelihood function is designed to model on-chip amplification of a EMCCD camera which captures a
 * Poisson process of emitted light, converted to electrons on the camera chip, amplified by gain modelled by the Gamma
 * distribution and then read with Gaussian noise.
 */
public class PoissonGammaGaussianFunction implements LikelihoodFunction, LogLikelihoodFunction
{
	/** The minimum read-noise. Values below this will be set to zerro and read-noise is not modelled. */
	public static final double MIN_READ_NOISE = 1e-3;

	/** Math.sqrt(2 * Math.PI). */
	private static final double sqrt2pi = Math.sqrt(2 * Math.PI);

	/**
	 * Define the convolution mode for combining the Poisson-Gamma PMF with the Gaussian PDF
	 */
	public enum ConvolutionMode
	{
		/**
		 * Use the approximation described in the Python source code within the supplementary information of the paper
		 * Mortensen, et al (2010) Nature Methods 7, 377-383.
		 */
		APPROXIMATION,

		/**
		 * Convolve the Poisson-Gamma on discrete (integer) intervals with the Gaussian PDF. This method is accurate
		 * when the read noise is above 1.
		 */
		DISCRETE_PDF,

		/**
		 * Convolve the Poisson-Gamma on discrete (integer) intervals with the Gaussian cumulative probability density
		 * function. This method is accurate for all read noise. This computes a PMF(X=x) on the integer scale.
		 */
		DISCRETE_PMF,

		//@formatter:off
		/**
		 * Convolve the Poisson-Gamma as a continuous PDF with the Gaussian using Simpson's Rule. If integration fails
		 * then the approximation will be used.
		 */
		SIMPSON_PDF { @Override public boolean validAtBoundary() { return false; }},

		/**
		 * Convolve the Poisson-Gamma as a continuous PDF with the Gaussian using a
		 * <a href="http://mathworld.wolfram.com/Legendre-GaussQuadrature.html">
		 * Legendre-Gauss</a> quadrature. If integration fails then the approximation will be used.
		 */
		LEGENDRE_GAUSS_PDF { @Override public boolean validAtBoundary() { return false; }};
		//@formatter:on

		/**
		 * Specify if this method is valid at the x=0 boundary.
		 *
		 * @return true, if valid when the Gaussian kernal can cross the x=0 boundary
		 */
		public boolean validAtBoundary()
		{
			return true;
		}
	}

	private class PGGFunction implements UnivariateFunction
	{
		int i = 0;
		final double o, e;

		public PGGFunction(double o, double e)
		{
			this.o = o;
			this.e = e;
		}

		public double value(double u)
		{
			i++;
			double pg = PoissonGammaFunction.poissonGammaNonZero(u, e, m);
			return (pg == 0) ? 0 : pg * gaussianPDF(u - o);
		}
	}

	/** The convolution mode. */
	private ConvolutionMode convolutionMode = ConvolutionMode.APPROXIMATION;

	/** The boundary convolution mode. */
	private ConvolutionMode boundaryConvolutionMode = ConvolutionMode.APPROXIMATION;

	/**
	 * The scale of the Gamma distribution (e.g. the on-chip gain multiplication factor)
	 */
	final private double m;
	/**
	 * The standard deviation of the Gaussian (e.g. Width of the noise distribution in the EMCCD output)
	 */
	final private double sigma;

	/** The two sigma 2. */
	private final double twoSigma2;

	/** The sqrt 2 sigma 2. */
	private final double sqrt2sigma2;

	/** The sqrt 2 pi sigma 2. */
	private final double sqrt2piSigma2;

	/** The minimum probability. */
	private double minimumProbability = Double.MIN_VALUE;

	/** The integrator. */
	private UnivariateIntegrator integrator = null;

	/**
	 * Instantiates a new poisson gamma gaussian function.
	 *
	 * @param alpha
	 *            Inverse gain of the EMCCD chip
	 * @param s
	 *            The Gaussian standard deviation at readout
	 * @throws IllegalArgumentException
	 *             If the gain is below 1 (the Gamma distribution is not modelled for a scale below 1).
	 */
	public PoissonGammaGaussianFunction(double alpha, double s) throws IllegalArgumentException
	{
		if (!(alpha > 0 && alpha <= 1))
			throw new IllegalArgumentException("Gain must be above 1");
		this.m = 1.0 / alpha;
		s = Math.abs(s);
		// Ignore tiny read noise
		if (s < MIN_READ_NOISE)
		{
			sigma = twoSigma2 = sqrt2sigma2 = sqrt2piSigma2 = 0;
		}
		else
		{
			sigma = s;
			twoSigma2 = 2 * s * s;
			sqrt2sigma2 = Math.sqrt(2 * s * s);
			sqrt2piSigma2 = Math.sqrt(2 * Math.PI * s * s);
		}
	}

	/**
	 * {@inheritDoc}
	 * <p>
	 * This code is adapted from the Python source code within the supplementary information of the paper Mortensen, et
	 * al (2010) Nature Methods 7, 377-383.
	 * <p>
	 * The output is a PMF. Ideally the input x should be discrete but this is not a requirement.
	 * 
	 * @see gdsc.smlm.function.LikelihoodFunction#likelihood(double, double)
	 */
	public double likelihood(final double o, final double e)
	{
		// This did not speed up MLE fitting so has been commented out.
		//		// When the observed ADUs and expected ADUs are much higher than the sigma then 
		//		// there is no point in convolving with a Gaussian
		//
		//		double mySigma = sigma;
		//		final double sLimit = sigma * 10;
		//		if (o > sLimit && e > sLimit)
		//		{
		//			//System.out.println("Skipping convolution");
		//			//mySigma = sigma;
		//		}
		//		
		//		if (mySigma == 0)

		if (sigma == 0)
		{
			// No convolution with a Gaussian. Simply evaluate for a Poisson-Gamma distribution.
			// This can handle e<=0.
			return checkMinProbability(PoissonGammaFunction.poissonGamma(o, e, m));
		}

		// If no Poisson mean then just use the Gaussian (Poisson-Gamma p=1 at x=0, p=0 otherwise)
		if (e <= 0)
		{
			// Output as PDF
			//return checkMinProbability(gaussianPDF(o));

			// Output as PMF
			return checkMinProbability(gaussianCDF(o - 0.5, o + 0.5));
		}

		if (convolutionMode == ConvolutionMode.APPROXIMATION)
		{
			return mortensenApproximation(o, e);
		}
		else
		{
			// Full evaluation of a Poisson-Gamma-Gaussian convolution PMF.
			ConvolutionMode mode = convolutionMode;

			// Integrate to infinity is not necessary. The convolution of the function with the 
			// Gaussian should be adequately sampled using a nxSD around the function value.
			// Find a bracket around the value.
			double range = 5 * sigma;
			int upper = (int) Math.ceil(o + range);
			if (upper < 0)
				return 0;
			int lower = (int) Math.floor(o - range);
			if (lower < 0)
				lower = 0;

			if (lower == 0 && !mode.validAtBoundary())
			{
				mode = boundaryConvolutionMode;
				if (mode == ConvolutionMode.APPROXIMATION)
					return mortensenApproximation(o, e);
			}

			double p = 0;

			// Special case:
			// Due to the Poisson-Gamma function delta function at u=0
			// The probability at u=0 may be very large compared to the rest of 
			// the Poisson-Gamma when e is low. To compensate always compute the 
			// full CDF integration at u=0.
			double erf2 = gaussianErf(-o + 0.5);
			if (erf2 != -1)
			{
				// Assume u==0
				p = FastMath.exp(-e) * (erf2 - gaussianErf(-o - 0.5)) * 0.5;
			}
			// Mark that u==0 has been computed
			if (lower == 0)
			{
				lower++;
				if (lower >= upper)
					// Nothing more to compute
					return checkMinProbability(p);
			}

			if (mode == ConvolutionMode.DISCRETE_PDF)
			{
				// Use a simple integration by adding the points in the range.
				if (lower == upper)
				{
					// Avoid double computation when lower==upper
					p += PoissonGammaFunction.poissonGammaNonZero(lower, e, m) * gaussianPDF(lower - o);
				}
				else
				{
					for (int u = lower; u <= upper; u++)
					{
						p += PoissonGammaFunction.poissonGammaNonZero(lower, e, m) * gaussianPDF(u - o);
					}
				}
			}
			else if (mode == ConvolutionMode.DISCRETE_PMF)
			{
				// Use a simple integration by adding the points in the range.
				// Use the error function to obtain the integral of the Gaussian
				double u_o = lower - o - 0.5;
				double erf = gaussianErf(u_o);
				if (lower == upper)
				{
					// Avoid double computation when lower==upper
					p += PoissonGammaFunction.poissonGammaNonZero(lower, e, m) * (gaussianErf(u_o + 1) - erf) * 0.5;
				}
				else
				{
					for (int u = lower; u <= upper; u++)
					{
						final double prevErf = erf;
						u_o += 1.0;
						erf = gaussianErf(u_o);
						p += PoissonGammaFunction.poissonGammaNonZero(lower, e, m) * (erf - prevErf) * 0.5;
					}
				}
			}
			else
			{
				// Use integrator

				// Note that the Poisson-Gamma function has a delta function at u=0. 
				// This prevents integration close to the zero boundary as the function 
				// is not smooth.

				// However the integrator may be faster when the range (upper-lower) 
				// is large as it uses fewer points.

				// Specify the function to integrate.
				PGGFunction f = new PGGFunction(o, e);

				try
				{
					p += createIntegrator().integrate(2000, f, lower, upper);
				}
				catch (TooManyEvaluationsException ex)
				{
					System.out.printf("Integration failed: o=%g, e=%g, eval=%d\n", o, e, f.i);
					return mortensenApproximation(o, e);
				}
				//System.out.printf("Integration eval=%d\n", f.i);
			}

			return checkMinProbability(p);
		}
	}

	private double checkMinProbability(double p)
	{
		return (p > minimumProbability) ? p : minimumProbability;
	}

	/**
	 * Mortensen approximation.
	 *
	 * @param cij
	 *            the cij
	 * @param eta
	 *            the eta
	 * @return the double
	 */
	private double mortensenApproximation(final double cij, final double eta)
	{
		// This code is adapted from the Python source code within the supplementary information of 
		// the paper Mortensen, et al (2010) Nature Methods 7, 377-383.

		// The implementation of the approximation is not documented.
		// This is meant to be convolving a PMF of a Poisson-Gamma mixture with the PDF of a Gaussian.

		// See Ulbrich & Isacoff (2007). Nature Methods 4, 319-321, SI equation 3.
		// G n>0 (c) = sum n {  (1 / n!) p^n e^-p (1 / ((n-1!)m^n)) c^n-1 e^-c/m }
		// G n>0 (c) = sqrt(p/(c*m)) * exp(-c/m - p) * Bessel.I1(2 * sqrt(c*p/m))
		// G n>0 (c) = sqrt(p/(c*m)) * exp(-c/m - p) * exp(2 * sqrt(c*p/m)) / sqrt(2*pi*sqrt(c*p/m))
		// G n=0 (c) = exp(-p)

		// p = eta
		// m = 1/alpha
		// c = cij

		// [Poisson PMF] multiplied by the [value at zero]:
		// [(eta^0 / 0!) * FastMath.exp(-eta)] * [eta * alpha]
		// FastMath.exp(-eta) * [eta * alpha]
		final double exp_eta = FastMath.exp(-eta);
		double f0 = exp_eta * eta / m;

		// ?
		double fp0 = f0 * 0.5 * (eta - 2) / m;

		// The cumulative normal distribution of the read noise
		// at the observed count
		final double conv0 = gaussianCDF(cij);

		// [Noise * Gaussian PMF at observed count] + 
		//  [observed count * cumulative distribution of read noise at observed count]
		// [sigma*FastMath.exp(-cij**2/(twoSigma2))/Math.sqrt(2*pi)] + [cij*conv0]
		final double conv1 = sigma * FastMath.exp(-(cij * cij) / twoSigma2) / sqrt2pi + cij * conv0;

		// ? 
		double temp = f0 * conv0 + fp0 * conv1 + exp_eta * gaussianPDF(cij);

		//		// TESTING
		//		// Simple method:
		//		temp = FastMath.exp(-eta) * gauss(cij); // G(c==0) * Gaussian;
		//		if (cij <= 0)
		//			return temp;
		//		// Reset. The remaining will be the Poisson-Gamma and no convolution
		//		f0 = fp0 = 0;
		//
		//		// Q. How to normalise so that at low cij there is a mixture and at high cij there is no mixture
		//		// and the result is the Poisson-Gamma. Perhaps this is what the above code is doing.

		if (cij > 0.0)
		{
			temp += PoissonGammaFunction.poissonGammaNonZero(cij, eta, m) - f0 - fp0 * cij;
		}

		// XXX : Debugging: Store the smallest likelihood we ever see. 
		// This can be used to set a limit for the likelihood
		//if (pMinObserved > temp && temp > 0)
		//{
		//	pMinObserved = temp;
		//}

		return checkMinProbability(temp);
	}

	/**
	 * {@inheritDoc}
	 * <p>
	 * This computes the log of {@link #likelihood(double, double)}.
	 * <p>
	 * The output is a PMF. Ideally the input x should be discrete but this is not a requirement.
	 * 
	 * @param o
	 *            The observed count
	 * @param e
	 *            The expected count
	 * @return The log-likelihood
	 * 
	 * @see #likelihood(double, double)
	 * @see gdsc.smlm.function.LogLikelihoodFunction#logLikelihood(double, double)
	 */
	public double logLikelihood(final double o, final double e)
	{
		return Math.log(likelihood(o, e));
	}

	/**
	 * Gaussian PDF.
	 *
	 * @param x
	 *            the x
	 * @return the density
	 */
	double gaussianPDF(final double x)
	{
		return FastMath.exp(-(x * x) / twoSigma2) / sqrt2piSigma2;
	}

	/**
	 * Gaussian CDF.
	 *
	 * @param x
	 *            the x
	 * @return the cumulative density
	 */
	double gaussianCDF(final double x)
	{
		//return 0.5 * (1 + org.apache.commons.math3.special.Erf.erf(x / sqrt2sigma2));
		// This may not be precise enough. 
		// Absolute error is <3e-7. Not sure what relative error is.
		// The standard Erf is much slower.
		return 0.5 * (1 + Erf.erf(x / sqrt2sigma2));
	}

	/**
	 * Gaussian CDF.
	 *
	 * @param x
	 *            the x
	 * @param x2
	 *            the x 2
	 * @return the cumulative density
	 */
	double gaussianCDF(final double x, final double x2)
	{
		//return 0.5 * (1 + org.apache.commons.math3.special.Erf.erf(x / sqrt2sigma2, x2 / sqrt2sigma2));
		// This may not be precise enough. 
		// Absolute error is <3e-7. Not sure what relative error is.
		// The standard Erf is much slower.
		return 0.5 * Erf.erf(x / sqrt2sigma2, x2 / sqrt2sigma2);
	}

	/**
	 * Gaussian CDF.
	 *
	 * @param x
	 *            the x
	 * @return the cumulative density
	 */
	double gaussianErf(final double x)
	{
		//return org.apache.commons.math3.special.Erf.erf(x / sqrt2sigma2);
		// This may not be precise enough. 
		// Absolute error is <3e-7. Not sure what relative error is.
		// The standard Erf is much slower.
		return Erf.erf(x / sqrt2sigma2);
	}

	/**
	 * Gets the alpha.
	 *
	 * @return the alpha
	 */
	public double getAlpha()
	{
		return 1 / m;
	}

	/**
	 * Gets the sigma.
	 *
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

	/**
	 * Gets the convolution mode.
	 *
	 * @return the convolution mode
	 */
	public ConvolutionMode getConvolutionMode()
	{
		return convolutionMode;
	}

	/**
	 * Sets the convolution mode.
	 * <p>
	 * WARNING:
	 * Currently the integration modes do not work when the Gaussian kernel will sample
	 * across x=0. The delta contribution of the Poisson-Gamma at x=0 only
	 * makes sense as part of a PMF. Using the Poisson-Gamma approximation
	 * as a PDF is incorrect when the Gaussian will cross the x=0 boundary.
	 * A different mode can be specified for use at the boundary.
	 *
	 * @param convolutionMode
	 *            the new convolution mode
	 * @see #setBoundaryConvolutionMode(ConvolutionMode)
	 */
	public void setConvolutionMode(ConvolutionMode convolutionMode)
	{
		this.convolutionMode = convolutionMode;
		integrator = null;
	}

	/**
	 * Gets the boundary convolution mode.
	 *
	 * @return the boundary convolution mode
	 */
	public ConvolutionMode getBoundaryConvolutionMode()
	{
		return boundaryConvolutionMode;
	}

	/**
	 * Sets the boundary convolution mode.
	 * WARNING:
	 * Currently the integration modes do not work when the Gaussian kernel will sample
	 * across x=0. The delta contribution of the Poisson-Gamma at x=0 only
	 * makes sense as part of a PMF. Using the Poisson-Gamma approximation
	 * as a PDF is incorrect when the Gaussian will cross the x=0 boundary.
	 * Using an invalid mode will cause an exception.
	 *
	 * @param boundaryConvolutionMode
	 *            the new boundary convolution mode
	 * @throws IllegalArgumentException
	 *             If the mode is not valid at the boundary
	 */
	public void setBoundaryConvolutionMode(ConvolutionMode boundaryConvolutionMode) throws IllegalArgumentException
	{
		if (!boundaryConvolutionMode.validAtBoundary())
			throw new IllegalArgumentException("Not valid at the boundary");
		this.boundaryConvolutionMode = boundaryConvolutionMode;
	}

	private UnivariateIntegrator createIntegrator()
	{
		UnivariateIntegrator i = integrator;
		if (i == null)
		{
			// This is the integrator for the Poisson-Gamma when observed count x>=1
			// i.e. not at the boundary x=0.

			final double relativeAccuracy = 1e-4;
			final double absoluteAccuracy = 1e-8;
			int minimalIterationCount;

			switch (convolutionMode)
			{
				case SIMPSON_PDF:
					// Number of function evaluations = 2^iteration + 1 
					// => 5 for 2 iterations
					// => 9 for 3 iterations
					minimalIterationCount = 2;
					i = new SimpsonIntegrator(relativeAccuracy, absoluteAccuracy, minimalIterationCount,
							SimpsonIntegrator.SIMPSON_MAX_ITERATIONS_COUNT);
					break;

				case LEGENDRE_GAUSS_PDF:
					// Not sure how to configure this.
					// The integration points are used for each sub-interval.
					// Function evaluations = integrationPoints * intervals.
					// The intervals start at 1,2 and increase by at least 4 at each stage after that.
					// At least 1 stage is done thus 3 * integrationPoints functions evaluations 
					// will be done for minimalIterationCount=1.
					minimalIterationCount = 1;
					final int maximalIterationCount = 32;
					final int integrationPoints = 8;
					i = new IterativeLegendreGaussIntegrator(integrationPoints, relativeAccuracy, absoluteAccuracy,
							minimalIterationCount, maximalIterationCount);
					break;

				default:
					throw new IllegalStateException();
			}

			integrator = i;
		}
		return i;
	}
}