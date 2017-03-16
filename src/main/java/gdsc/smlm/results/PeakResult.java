/*
 * 
 */
package gdsc.smlm.results;

import gdsc.smlm.function.gaussian.Gaussian2DFunction;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.exception.TooManyEvaluationsException;

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

/**
 * Specifies a peak fitting result
 */
public class PeakResult implements Comparable<PeakResult>
{
	public int peak;
	public int origX;
	public int origY;
	public float origValue;
	public double error;
	public float noise;
	public float[] params;
	public float[] paramsStdDev;

	public PeakResult(int peak, int origX, int origY, float origValue, double error, float noise, float[] params,
			float[] paramsStdDev)
	{
		this.peak = peak;
		this.origX = origX;
		this.origY = origY;
		this.origValue = origValue;
		this.error = error;
		this.noise = noise;
		this.params = params;
		this.paramsStdDev = paramsStdDev;
	}

	/**
	 * Simple constructor to create a result with location, width and strength
	 * 
	 * @param x
	 * @param y
	 * @param sd
	 * @param signal
	 */
	public PeakResult(float x, float y, float sd, float signal)
	{
		origX = (int) Math.round(x);
		origY = (int) Math.round(y);
		params = new float[7];
		params[Gaussian2DFunction.X_POSITION] = x;
		params[Gaussian2DFunction.Y_POSITION] = y;
		params[Gaussian2DFunction.X_SD] = params[Gaussian2DFunction.Y_SD] = sd;
		params[Gaussian2DFunction.SIGNAL] = signal;
	}

	/**
	 * Get the signal strength (i.e. the volume under the Gaussian peak, amplitude * 2 * pi * sx * sy)
	 * 
	 * @return The signal of the first peak
	 */
	public float getSignal()
	{
		return params[Gaussian2DFunction.SIGNAL];
	}

	/**
	 * Get the signal strength for the nth peak (i.e. the volume under the Gaussian peak, amplitude * 2 * pi * sx * sy)
	 * 
	 * @param peakId
	 *            The peak number
	 * @return The signal of the nth peak
	 */
	public float getSignal(int peakId)
	{
		final int index = peakId * 6 + Gaussian2DFunction.SIGNAL;
		return (index < params.length) ? params[index] : 0;
	}

	/**
	 * Calculate the localisation precision. Uses the Mortensen formula for an EMCCD camera
	 * (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 6.
	 * <p>
	 * This method will use the background noise of the peak to approximate the expected background value of each pixel.
	 * 
	 * @param a
	 *            The size of the pixels in nm
	 * @param gain
	 *            The gain used to convert ADUs to photons
	 * @param emCCD
	 *            True if an emCCD camera
	 * @return The location precision in nm of the first peak
	 */
	public double getPrecision(final double a, final double gain, boolean emCCD)
	{
		// Get peak standard deviation in nm. Just use the average of the X & Y.
		final double s = a * getSD();
		final double N = getSignal();
		return getPrecision(a, s, N / gain, noise / gain, emCCD);
	}

	/**
	 * Calculate the localisation variance. Uses the Mortensen formula for an EMCCD camera
	 * (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 6.
	 * <p>
	 * This method will use the background noise of the peak to approximate the expected background value of each pixel.
	 * 
	 * @param a
	 *            The size of the pixels in nm
	 * @param gain
	 *            The gain used to convert ADUs to photons
	 * @param emCCD
	 *            True if an emCCD camera
	 * @return The location variance in nm of the first peak
	 */
	public double getVariance(final double a, final double gain, boolean emCCD)
	{
		// Get peak standard deviation in nm. Just use the average of the X & Y.
		final double s = a * getSD();
		final double N = getSignal();
		return getVariance(a, s, N / gain, noise / gain, emCCD);
	}

	/**
	 * Calculate the localisation precision for Maximum Likelihood estimation. Uses the Mortensen formula for an EMCCD
	 * camera (Mortensen, et al (2010) Nature Methods 7, 377-383), SI equation 54.
	 * <p>
	 * This method will use the background noise of the peak to approximate the expected background value of each pixel.
	 * 
	 * @param a
	 *            The size of the pixels in nm
	 * @param gain
	 *            The gain used to convert ADUs to photons
	 * @param emCCD
	 *            True if an emCCD camera
	 * @return The location precision in nm of the first peak
	 */
	public double getMLPrecision(final double a, final double gain, boolean emCCD)
	{
		// Get peak standard deviation in nm. Just use the average of the X & Y.
		final double s = a * getSD();
		final double N = getSignal();
		return getMLPrecision(a, s, N / gain, noise / gain, emCCD);
	}

	/**
	 * Calculate the localisation variance for Maximum Likelihood estimation. Uses the Mortensen formula for an EMCCD
	 * camera (Mortensen, et al (2010) Nature Methods 7, 377-383), SI equation 54.
	 * <p>
	 * This method will use the background noise of the peak to approximate the expected background value of each pixel.
	 * 
	 * @param a
	 *            The size of the pixels in nm
	 * @param gain
	 *            The gain used to convert ADUs to photons
	 * @param emCCD
	 *            True if an emCCD camera
	 * @return The location variance in nm of the first peak
	 */
	public double getMLVariance(final double a, final double gain, boolean emCCD)
	{
		// Get peak standard deviation in nm. Just use the average of the X & Y.
		final double s = a * getSD();
		final double N = getSignal();
		return getMLVariance(a, s, N / gain, noise / gain, emCCD);
	}

	/**
	 * Convert the local background to an estimate of noise. Local background and noise are in ADU count units.
	 * <p>
	 * This assumes the local background is photon shot noise. The background is first converted to photons using the
	 * gain. The shot noise is taken assuming a Poisson distribution (thus the variance equals the number of photons).
	 * This is amplified by 2 if the data was taken on an EM-CCD camera. The square root is the noise in photons. This
	 * is converted back to ADUs using the gain. E.G.
	 * 
	 * <pre>
	 * return Math.sqrt((background / gain) * ((emCCD) ? 2 : 1)) * gain;
	 * </pre>
	 *
	 * @param background
	 *            the background
	 * @param gain
	 *            the gain
	 * @param emCCD
	 *            True if an emCCD camera
	 * @return the noise estimate
	 */
	public static double localBackgroundToNoise(double background, double gain, boolean emCCD)
	{
		if (background <= 0)
			return 0;
		return Math.sqrt((background / gain) * ((emCCD) ? 2 : 1)) * gain;
	}

	/**
	 * Convert the noise to local background. Local background and noise are in ADU count units.
	 * <p>
	 * This assumes the local background is photon shot noise. This is the opposite conversion to
	 * {@link #localBackgroundToNoise(double, double, boolean)}.
	 *
	 * @param noise
	 *            the noise
	 * @param gain
	 *            the gain
	 * @param emCCD
	 *            True if an emCCD camera
	 * @return the local background estimate
	 */
	public static double noiseToLocalBackground(double noise, double gain, boolean emCCD)
	{
		if (noise <= 0)
			return 0;
		noise /= gain;
		noise *= noise;
		if (emCCD)
			noise /= 2;
		return noise * gain;
	}

	/**
	 * Calculate the localisation precision. Uses the Mortensen formula for an EMCCD camera
	 * (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 6.
	 * <p>
	 * This method will use the background noise to approximate the expected background value of each pixel.
	 * 
	 * @param a
	 *            The size of the pixels in nm
	 * @param s
	 *            The peak standard deviation in nm
	 * @param N
	 *            The peak signal in photons
	 * @param b
	 *            The background noise standard deviation in photons
	 * @param emCCD
	 *            True if an emCCD camera
	 * @return The location precision in nm in each dimension (X/Y)
	 */
	public static double getPrecision(double a, double s, double N, double b, boolean emCCD)
	{
		// Get background expected value. This is what is actually used in the Mortensen method
		final double b2 = b * b;

		if (emCCD)
		{
			// If an emCCD camera was used then the input standard deviation will already be amplified 
			// by the EM-gain sqrt(2) factor. To prevent double counting this factor we must divide by it.
			// Since this has been squared then divide by 2.
			return getPrecisionX(a, s, N, b2 / 2.0, 2);
		}
		else
		{
			return getPrecisionX(a, s, N, b2, 1);
		}
	}

	/**
	 * Calculate the localisation variance. Uses the Mortensen formula for an EMCCD camera
	 * (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 6.
	 * <p>
	 * This method will use the background noise to approximate the expected background value of each pixel.
	 * 
	 * @param a
	 *            The size of the pixels in nm
	 * @param s
	 *            The peak standard deviation in nm
	 * @param N
	 *            The peak signal in photons
	 * @param b
	 *            The background noise standard deviation in photons
	 * @param emCCD
	 *            True if an emCCD camera
	 * @return The location variance in nm in each dimension (X/Y)
	 */
	public static double getVariance(double a, double s, double N, double b, boolean emCCD)
	{
		// Get background expected value. This is what is actually used in the Mortensen method
		final double b2 = b * b;

		if (emCCD)
		{
			// If an emCCD camera was used then the input standard deviation will already be amplified 
			// by the EM-gain sqrt(2) factor. To prevent double counting this factor we must divide by it.
			// Since this has been squared then divide by 2.
			return getVarianceX(a, s, N, b2 / 2.0, 2);
		}
		else
		{
			return getVarianceX(a, s, N, b2, 1);
		}
	}

	/**
	 * Calculate the localisation precision for Maximum Likelihood estimation. Uses the Mortensen formula for an EMCCD
	 * camera (Mortensen, et al (2010) Nature Methods 7, 377-383), SI equation 54.
	 * <p>
	 * This method will use the background noise to approximate the expected background value of each pixel.
	 * 
	 * @param a
	 *            The size of the pixels in nm
	 * @param s
	 *            The peak standard deviation in nm
	 * @param N
	 *            The peak signal in photons
	 * @param b
	 *            The background noise standard deviation in photons
	 * @param emCCD
	 *            True if an emCCD camera
	 * @return The location precision in nm in each dimension (X/Y)
	 */
	public static double getMLPrecision(double a, double s, double N, double b, boolean emCCD)
	{
		// Get background expected value. This is what is actually used in the Mortensen method
		final double b2 = b * b;

		if (emCCD)
		{
			// If an emCCD camera was used then the input standard deviation will already be amplified 
			// by the EM-gain sqrt(2) factor. To prevent double counting this factor we must divide by it.
			// Since this has been squared then divide by 2.
			return getMLPrecisionX(a, s, N, b2 / 2.0, true);
		}
		else
		{
			return getMLPrecisionX(a, s, N, b2, false);
		}
	}

	/**
	 * Calculate the localisation variance for Maximum Likelihood estimation. Uses the Mortensen formula for an EMCCD
	 * camera (Mortensen, et al (2010) Nature Methods 7, 377-383), SI equation 54.
	 * <p>
	 * This method will use the background noise to approximate the expected background value of each pixel.
	 * 
	 * @param a
	 *            The size of the pixels in nm
	 * @param s
	 *            The peak standard deviation in nm
	 * @param N
	 *            The peak signal in photons
	 * @param b
	 *            The background noise standard deviation in photons
	 * @param emCCD
	 *            True if an emCCD camera
	 * @return The location variance in nm in each dimension (X/Y)
	 */
	public static double getMLVariance(double a, double s, double N, double b, boolean emCCD)
	{
		// Get background expected value. This is what is actually used in the Mortensen method
		final double b2 = b * b;

		if (emCCD)
		{
			// If an emCCD camera was used then the input standard deviation will already be amplified 
			// by the EM-gain sqrt(2) factor. To prevent double counting this factor we must divide by it.
			// Since this has been squared then divide by 2.
			return getMLVarianceX(a, s, N, b2 / 2.0, true);
		}
		else
		{
			return getMLVarianceX(a, s, N, b2, false);
		}
	}

	/**
	 * Calculate the localisation precision for least squares estimation. Uses the Mortensen formula for an EMCCD camera
	 * (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 6.
	 * <p>
	 * If the expected photons per pixel is unknown then use the standard deviation across the image and the method
	 * {@link #getPrecision(double, double, double, double, boolean)}.
	 * 
	 * @param a
	 *            The size of the pixels in nm
	 * @param s
	 *            The peak standard deviation in nm
	 * @param N
	 *            The peak signal in photons
	 * @param b2
	 *            The expected number of photons per pixel from a background with spatially constant
	 *            expectation value across the image (Note that this is b^2 not b, which could be the standard deviation
	 *            of the image pixels)
	 * @param emCCD
	 *            True if an emCCD camera
	 * @return The location precision in nm in each dimension (X/Y)
	 */
	public static double getPrecisionX(final double a, final double s, final double N, final double b2, boolean emCCD)
	{
		// EM-CCD noise factor
		final double F = (emCCD) ? 2 : 1;
		return getPrecisionX(a, s, N, b2, F);
	}

	/**
	 * Calculate the localisation precision for least squares estimation. Uses the Mortensen formula for an EMCCD camera
	 * (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 6.
	 * 
	 * @param a
	 *            The size of the pixels in nm
	 * @param s
	 *            The peak standard deviation in nm
	 * @param N
	 *            The peak signal in photons
	 * @param b2
	 *            The expected number of photons per pixel from a background with spatially constant
	 *            expectation value across the image (Note that this is b^2 not b, which could be the standard deviation
	 *            of the image pixels)
	 * @param F
	 *            EM-CCD noise factor (usually 2 for an EM-CCD camera, else 1)
	 * @return The location precision in nm in each dimension (X/Y)
	 */
	public static double getPrecisionX(final double a, final double s, final double N, final double b2, final double F)
	{
		return Math.sqrt(getVarianceX(a, s, N, b2, F));
	}

	/**
	 * Calculate the localisation variance for least squares estimation. Uses the Mortensen formula for an EMCCD camera
	 * (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 6.
	 * <p>
	 * If the expected photons per pixel is unknown then use the standard deviation across the image and the method
	 * {@link #getPrecision(double, double, double, double, boolean)}.
	 * 
	 * @param a
	 *            The size of the pixels in nm
	 * @param s
	 *            The peak standard deviation in nm
	 * @param N
	 *            The peak signal in photons
	 * @param b2
	 *            The expected number of photons per pixel from a background with spatially constant
	 *            expectation value across the image (Note that this is b^2 not b, which could be the standard deviation
	 *            of the image pixels)
	 * @param emCCD
	 *            True if an emCCD camera
	 * @return The location variance in nm in each dimension (X/Y)
	 */
	public static double getVarianceX(final double a, final double s, final double N, final double b2, boolean emCCD)
	{
		// EM-CCD noise factor
		final double F = (emCCD) ? 2 : 1;
		return getVarianceX(a, s, N, b2, F);
	}

	/**
	 * Calculate the localisation variance for least squares estimation. Uses the Mortensen formula for an EMCCD camera
	 * (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 6.
	 * 
	 * @param a
	 *            The size of the pixels in nm
	 * @param s
	 *            The peak standard deviation in nm
	 * @param N
	 *            The peak signal in photons
	 * @param b2
	 *            The expected number of photons per pixel from a background with spatially constant
	 *            expectation value across the image (Note that this is b^2 not b, which could be the standard deviation
	 *            of the image pixels)
	 * @param F
	 *            EM-CCD noise factor (usually 2 for an EM-CCD camera, else 1)
	 * @return The location variance in nm in each dimension (X/Y)
	 */
	public static double getVarianceX(final double a, final double s, final double N, final double b2, final double F)
	{
		if (N <= 0)
			return Double.POSITIVE_INFINITY;

		// Note that we input b^2 directly to this equation. This is the expected value of the pixel background. 
		// If the background is X then the variance of a Poisson distribution will be X 
		// and the standard deviation at each pixel will be sqrt(X). Thus the Mortensen formula
		// can be used without knowing the background explicitly by using the variance of the pixels.

		final double a2 = a * a;
		// Adjustment for square pixels
		final double sa2 = s * s + a2 / 12.0;
		// 16 / 9 = 1.7777777778
		// 8 * pi = 25.13274123
		return F * (sa2 / N) * (1.7777777778 + (25.13274123 * sa2 * b2) / (N * a2));
	}

	/**
	 * Calculate the localisation precision for Maximum Likelihood estimation. Uses the Mortensen formula for an EMCCD
	 * camera (Mortensen, et al (2010) Nature Methods 7, 377-383), SI equation 54.
	 * 
	 * @param a
	 *            The size of the pixels in nm
	 * @param s
	 *            The peak standard deviation in nm
	 * @param N
	 *            The peak signal in photons
	 * @param b2
	 *            The expected number of photons per pixel from a background with spatially constant
	 *            expectation value across the image (Note that this is b^2 not b, which could be the standard deviation
	 *            of the image pixels)
	 * @param emCCD
	 *            True if an emCCD camera
	 * @return The location precision in nm in each dimension (X/Y)
	 */
	public static double getMLPrecisionX(double a, double s, double N, double b2, boolean emCCD)
	{
		return Math.sqrt(getMLVarianceX(a, s, N, b2, emCCD));
	}

	/**
	 * Calculate the localisation precision for Maximum Likelihood estimation. Uses the Mortensen formula for an EMCCD
	 * camera (Mortensen, et al (2010) Nature Methods 7, 377-383), SI equation 54.
	 * 
	 * @param a
	 *            The size of the pixels in nm
	 * @param s
	 *            The peak standard deviation in nm
	 * @param N
	 *            The peak signal in photons
	 * @param b2
	 *            The expected number of photons per pixel from a background with spatially constant
	 *            expectation value across the image (Note that this is b^2 not b, which could be the standard deviation
	 *            of the image pixels)
	 * @param emCCD
	 *            True if an emCCD camera
	 * @param integrationPoints
	 *            the number of integration points for the LegendreGaussIntegrator
	 * @return The location precision in nm in each dimension (X/Y)
	 */
	public static double getMLPrecisionX(double a, double s, double N, double b2, boolean emCCD, int integrationPoints)
	{
		return Math.sqrt(getMLVarianceX(a, s, N, b2, emCCD, integrationPoints));
	}

	/**
	 * Calculate the localisation variance for Maximum Likelihood estimation. Uses the Mortensen formula for an EMCCD
	 * camera (Mortensen, et al (2010) Nature Methods 7, 377-383), SI equation 54.
	 * 
	 * @param a
	 *            The size of the pixels in nm
	 * @param s
	 *            The peak standard deviation in nm
	 * @param N
	 *            The peak signal in photons
	 * @param b2
	 *            The expected number of photons per pixel from a background with spatially constant
	 *            expectation value across the image (Note that this is b^2 not b, which could be the standard deviation
	 *            of the image pixels)
	 * @param emCCD
	 *            True if an emCCD camera
	 * @return The location variance in nm in each dimension (X/Y)
	 */
	public static double getMLVarianceX(double a, double s, double N, double b2, boolean emCCD)
	{
		// The class JUnit test shows that 10 integration points is the fastest for realistic input parameters.
		return getMLVarianceX(a, s, N, b2, emCCD, 10);
	}

	/**
	 * Calculate the localisation variance for Maximum Likelihood estimation. Uses the Mortensen formula for an EMCCD
	 * camera (Mortensen, et al (2010) Nature Methods 7, 377-383), SI equation 54.
	 * <p>
	 * In the event of failure to integrate the formula the variance for Least Squares Estimation is returned.
	 * 
	 * @param a
	 *            The size of the pixels in nm
	 * @param s
	 *            The peak standard deviation in nm
	 * @param N
	 *            The peak signal in photons
	 * @param b2
	 *            The expected number of photons per pixel from a background with spatially constant
	 *            expectation value across the image (Note that this is b^2 not b, which could be the standard deviation
	 *            of the image pixels)
	 * @param emCCD
	 *            True if an emCCD camera
	 * @param integrationPoints
	 *            the number of integration points for the LegendreGaussIntegrator
	 * @return The location variance in nm in each dimension (X/Y)
	 */
	public static double getMLVarianceX(double a, double s, double N, double b2, boolean emCCD, int integrationPoints)
	{
		if (N <= 0)
			return Double.POSITIVE_INFINITY;

		// EM-CCD noise factor
		final double F = (emCCD) ? 2 : 1;
		final double a2 = a * a;
		// Adjustment for square pixels
		final double sa2 = s * s + a2 / 12.0;

		final double rho = 2 * Math.PI * sa2 * b2 / (N * a2);
		try
		{
			final double I1 = computeI1(rho, integrationPoints);
			if (I1 > 0)
				return F * (sa2 / N) * (1 / I1);
			//else 
			//	System.out.printf("Invalid I1 = %f\n", I1);
		}
		catch (TooManyEvaluationsException e)
		{

		}
		return getVarianceX(a, s, N, b2, emCCD);
	}

	/**
	 * Compute the function I1 using numerical integration. See Mortensen, et al (2010) Nature Methods 7, 377-383), SI
	 * equation 43.
	 * 
	 * <pre>
	 * I1 = 1 + sum [ ln(t) / (1 + t/rho) ] dt
	 *    = - sum [ t * ln(t) / (t + rho) ] dt
	 * </pre>
	 * 
	 * Where sum is the integral between 0 and 1. In the case of rho=0 the function returns 1;
	 * 
	 * @param rho
	 * @param integrationPoints
	 *            the number of integration points for the LegendreGaussIntegrator
	 * @return the I1 value
	 */
	private static double computeI1(final double rho, int integrationPoints)
	{
		if (rho == 0)
			return 1;

		final double relativeAccuracy = 1e-4;
		final double absoluteAccuracy = 1e-8;
		final int minimalIterationCount = 3;
		final int maximalIterationCount = 32;

		// Use an integrator that does not use the boundary since log(0) is undefined.
		UnivariateIntegrator i = new IterativeLegendreGaussIntegrator(integrationPoints, relativeAccuracy,
				absoluteAccuracy, minimalIterationCount, maximalIterationCount);

		// Specify the function to integrate
		UnivariateFunction f = new UnivariateFunction()
		{
			public double value(double x)
			{
				return x * Math.log(x) / (x + rho);
			}
		};
		final double i1 = -i.integrate(2000, f, 0, 1);
		//System.out.printf("I1 = %f (%d)\n", i1, i.getEvaluations());

		// The function requires more evaluations and sometimes does not converge,
		// presumably because log(x) significantly changes as x -> 0 where as x log(x) in the function above 
		// is more stable

		//		UnivariateFunction f2 = new UnivariateFunction()
		//		{
		//			@Override
		//			public double value(double x)
		//			{
		//				return Math.log(x) / ( 1 + x / rho);
		//			}
		//		};
		//		double i2 = 1 + i.integrate(2000, f2, 0, 1);
		//		System.out.printf("I1 (B) = %f (%d)\n", i2, i.getEvaluations());

		return i1;
	}

	public int compareTo(PeakResult o)
	{
		// Sort by peak number: Ascending
		if (peak == o.peak)
		{
			// Sort by peak height: Descending
			if (params[Gaussian2DFunction.SIGNAL] > o.params[Gaussian2DFunction.SIGNAL])
				return -1;
			if (params[Gaussian2DFunction.SIGNAL] < o.params[Gaussian2DFunction.SIGNAL])
				return 1;
			return 0;
		}
		return peak - o.peak;
	}

	/**
	 * Calculate the combined peak standard deviation. This is equal to the square root of the product of the width in X
	 * and Y.
	 * 
	 * @return The combined peak standard deviation in the X and Y dimension
	 */
	public float getSD()
	{
		return getSD(params[Gaussian2DFunction.X_SD], params[Gaussian2DFunction.Y_SD]);
	}

	/**
	 * Calculate the combined peak standard deviation. This is equal to the square root of the product of the width in X
	 * and Y.
	 * 
	 * @param xsd
	 * @param ysd
	 * @return The combined peak standard deviation in the X and Y dimension
	 */
	public static float getSD(float xsd, float ysd)
	{
		if (xsd == ysd)
			return xsd;
		return (float) Math.sqrt(Math.abs(xsd * ysd));
	}

	/**
	 * Calculate the combined peak standard deviation. This is equal to the square root of the product of the width in X
	 * and Y.
	 * 
	 * @param xsd
	 * @param ysd
	 * @return The combined peak standard deviation in the X and Y dimension
	 */
	public static double getSD(double xsd, double ysd)
	{
		if (xsd == ysd)
			return xsd;
		return Math.sqrt(Math.abs(xsd * ysd));
	}

	/**
	 * @return The background for the first peak
	 */
	public float getBackground()
	{
		return params[Gaussian2DFunction.BACKGROUND];
	}

	/**
	 * Get the amplitude for the first peak. Amplitude = Signal / (2*pi*sd0*sd1).
	 * 
	 * @return The amplitude for the first peak
	 */
	public float getAmplitude()
	{
		return getAmplitude(params);
	}

	/**
	 * Get the amplitude for the first peak. Amplitude = Signal / (2*pi*sd0*sd1).
	 * 
	 * @return The amplitude for the first peak
	 */
	public static float getAmplitude(float[] params)
	{
		return (float) (params[Gaussian2DFunction.SIGNAL] /
				(2 * Math.PI * params[Gaussian2DFunction.X_SD] * params[Gaussian2DFunction.Y_SD]));
	}

	/**
	 * @return The angle for the first peak
	 */
	public float getAngle()
	{
		return params[Gaussian2DFunction.ANGLE];
	}

	/**
	 * @return The x position for the first peak
	 */
	public float getXPosition()
	{
		return params[Gaussian2DFunction.X_POSITION];
	}

	/**
	 * @return The y position for the first peak
	 */
	public float getYPosition()
	{
		return params[Gaussian2DFunction.Y_POSITION];
	}

	/**
	 * @return The x-dimension standard deviation for the first peak
	 */
	public float getXSD()
	{
		return params[Gaussian2DFunction.X_SD];
	}

	/**
	 * @return The y-dimension standard deviation for the first peak
	 */
	public float getYSD()
	{
		return params[Gaussian2DFunction.Y_SD];
	}

	/**
	 * Gets the frame.
	 *
	 * @return The time frame that this result corresponds to
	 */
	public int getFrame()
	{
		return peak;
	}

	/**
	 * Checks for end frame. Derived classes can override this.
	 *
	 * @return true, if successful
	 */
	public boolean hasEndFrame()
	{
		return false;
	}

	/**
	 * Gets the end frame. Default = {@link #getFrame()}. Derived classes can override this.
	 *
	 * @return The last time frame that this result corresponds to
	 */
	public int getEndFrame()
	{
		return peak;
	}

	/**
	 * Checks for id. Derived classes can override this.
	 *
	 * @return true, if successful
	 */
	public boolean hasId()
	{
		return false;
	}

	/**
	 * Gets the id. Default = 0. Derived classes can override this.
	 *
	 * @return The results identifier
	 */
	public int getId()
	{
		return 0;
	}

	/**
	 * Checks for precision. Derived classes can override this.
	 *
	 * @return true, if successful
	 */
	public boolean hasPrecision()
	{
		return false;
	}

	/**
	 * Gets the localisation precision. Default = Double.NaN. Derived classes can override this.
	 * <p>
	 * This is provided so that results can be loaded from external sources.
	 *
	 * @return the precision (in nm)
	 */
	public double getPrecision()
	{
		return Double.NaN;
	}

	/**
	 * Return the true positive score for use in classification analysis
	 * 
	 * @return The true positive score
	 */
	public double getTruePositiveScore()
	{
		return (origValue != 0) ? 1 : 0;
	}

	/**
	 * Return the false positive score for use in classification analysis
	 * 
	 * @return The false positive score
	 */
	public double getFalsePositiveScore()
	{
		return 1 - getTruePositiveScore();
	}

	/**
	 * Return the true negative score for use in classification analysis
	 * 
	 * @return The true negative score
	 */
	public double getTrueNegativeScore()
	{
		return (origValue != 0) ? 0 : 1;
	}

	/**
	 * Return the false negative score for use in classification analysis
	 * 
	 * @return The false negative score
	 */
	public double getFalseNegativeScore()
	{
		return 1 - getTrueNegativeScore();
	}

	/**
	 * Return the squared distance to the other peak result
	 * 
	 * @param r
	 *            The result
	 * @return The squared distance
	 */
	public double distance2(PeakResult r)
	{
		final double dx = getXPosition() - r.getXPosition();
		final double dy = getYPosition() - r.getYPosition();
		return dx * dx + dy * dy;
	}

	/**
	 * Return the distance to the other peak result
	 * 
	 * @param r
	 *            The result
	 * @return The distance
	 */
	public double distance(PeakResult r)
	{
		return Math.sqrt(distance2(r));
	}

	/**
	 * Return the squared distance to the other coordinate.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @return The squared distance
	 */
	public double distance2(final double x, final double y)
	{
		final double dx = getXPosition() - x;
		final double dy = getYPosition() - y;
		return dx * dx + dy * dy;
	}

	/**
	 * Return the distance to the other coordinate.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @return The distance
	 */
	public double distance(final double x, final double y)
	{
		return Math.sqrt(distance2(x, y));
	}

	/**
	 * This methods return the x-position. To allow filters to use the actual shift requires either off-setting the
	 * position with the initial fit position, or extending this class so the shift can be stored.
	 */
	public float getXShift()
	{
		return getXPosition();
	}

	/**
	 * This methods return the y-position. To allow filters to use the actual shift requires either off-setting the
	 * position with the initial fit position, or extending this class so the shift can be stored.
	 */
	public float getYShift()
	{
		return getYPosition();
	}

	/**
	 * Gets the noise.
	 *
	 * @return the noise
	 */
	public float getNoise()
	{
		return noise;
	}
}
