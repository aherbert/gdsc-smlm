/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm.results;

/**
 * Contains calculator functions for working with peak results
 */
public interface Gaussian2DPeakResultCalculator
{
	/**
	 * Gets the single Gaussian 2D standard deviation from independent x and y standard deviations. s =
	 * sqrt(abs(sx*sy)).
	 *
	 * @param params
	 *            the params
	 * @return the single Gaussian 2D standard deviation
	 */
	public float getStandardDeviation(float[] params);

	/**
	 * Gets the single Gaussian 2D standard deviation squared from independent x and y standard deviations. s2 =
	 * abs(sx*sy).
	 *
	 * @param params
	 *            the params
	 * @return the single Gaussian 2D standard deviation squared
	 */
	public float getStandardDeviation2(float[] params);

	/**
	 * Get the amplitude of a Gaussian 2D PSF. Amplitude = intensity / (2*pi*sx*sy).
	 *
	 * @param params
	 *            the params
	 * @return the amplitude
	 */
	public float getAmplitude(float[] params);

	/**
	 * Get the height of the central pixel of a Gaussian 2D PSF. The integral of the pixel containing the
	 * centre of the Gaussian is computed.
	 *
	 * @param params
	 *            the params
	 * @return the pixel amplitude
	 */
	public float getPixelAmplitude(float[] params);

	/**
	 * Calculate the localisation precision for least squares fitting a Gaussian2D PSF to a Gaussian2D PSF. This is
	 * an approximation of the precision of fitting to an optical PSF. Uses the Mortensen formula for an EMCCD camera
	 * (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 6.
	 * <p>
	 * This method will use the background noise to approximate the expected background value of each pixel.
	 *
	 * @param params
	 *            the params
	 * @param noise
	 *            the noise
	 * @return The location precision in nm in each dimension (X/Y)
	 */
	public double getLSEPrecision(float[] params, float noise);

	/**
	 * Calculate the localisation precision for least squares fitting a Gaussian2D PSF to a Gaussian2D PSF. This is
	 * an approximation of the precision of fitting to an optical PSF. Uses the Mortensen formula for an EMCCD camera
	 * (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 6.
	 *
	 * @param params
	 *            the params
	 * @return The location precision in nm in each dimension (X/Y)
	 */
	public double getLSEPrecision(float[] params);

	/**
	 * Calculate the localisation variance for least squares fitting a Gaussian2D PSF to a Gaussian2D PSF. This is
	 * an approximation of the variance of fitting to an optical PSF. Uses the Mortensen formula for an EMCCD camera
	 * (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 6.
	 * <p>
	 * This method will use the background noise to approximate the expected background value of each pixel.
	 *
	 * @param params
	 *            the params
	 * @param noise
	 *            the noise
	 * @return The location variance in nm in each dimension (X/Y)
	 */
	public double getLSEVariance(float[] params, float noise);

	/**
	 * Calculate the localisation variance for least squares fitting a Gaussian2D PSF to a Gaussian2D PSF. This is
	 * an approximation of the variance of fitting to an optical PSF. Uses the Mortensen formula for an EMCCD camera
	 * (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 6.
	 *
	 * @param params
	 *            the params
	 * @return The location variance in nm in each dimension (X/Y)
	 */
	public double getLSEVariance(float[] params);

	/**
	 * Calculate the localisation precision for maximum likelihood fitting a Gaussian2D PSF to a Gaussian2D PSF. This is
	 * an approximation of the precision of fitting to an optical PSF. Uses the Mortensen formula for an EMCCD camera
	 * (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 6.
	 * <p>
	 * This method will use the background noise to approximate the expected background value of each pixel.
	 *
	 * @param params
	 *            the params
	 * @param noise
	 *            the noise
	 * @return The location precision in nm in each dimension (X/Y)
	 */
	public double getMLEPrecision(float[] params, float noise);

	/**
	 * Calculate the localisation precision for maximum likelihood fitting a Gaussian2D PSF to a Gaussian2D PSF. This is
	 * an approximation of the precision of fitting to an optical PSF. Uses the Mortensen formula for an EMCCD camera
	 * (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 6.
	 *
	 * @param params
	 *            the params
	 * @return The location precision in nm in each dimension (X/Y)
	 */
	public double getMLEPrecision(float[] params);

	/**
	 * Calculate the localisation variance for maximum likelihood fitting a Gaussian2D PSF to a Gaussian2D PSF. This is
	 * an approximation of the variance of fitting to an optical PSF. Uses the Mortensen formula for an EMCCD camera
	 * (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 6.
	 * <p>
	 * This method will use the background noise to approximate the expected background value of each pixel.
	 *
	 * @param params
	 *            the params
	 * @param noise
	 *            the noise
	 * @return The location variance in nm in each dimension (X/Y)
	 */
	public double getMLEVariance(float[] params, float noise);

	/**
	 * Calculate the localisation variance for maximum likelihood fitting a Gaussian2D PSF to a Gaussian2D PSF. This is
	 * an approximation of the variance of fitting to an optical PSF. Uses the Mortensen formula for an EMCCD camera
	 * (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 6.
	 *
	 * @param params
	 *            the params
	 * @return The location variance in nm in each dimension (X/Y)
	 */
	public double getMLEVariance(float[] params);
}
