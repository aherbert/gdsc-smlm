package gdsc.smlm.results;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.exception.TooManyEvaluationsException;

import gdsc.core.data.utils.TypeConverter;
import gdsc.core.utils.BitFlags;
import gdsc.smlm.data.config.CalibrationHelper;
import gdsc.smlm.data.config.ConfigurationException;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.data.config.SMLMSettings.CameraType;
import gdsc.smlm.data.config.SMLMSettings.DistanceUnit;
import gdsc.smlm.data.config.SMLMSettings.IntensityUnit;
import gdsc.smlm.data.config.SMLMSettings.PSF;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Contains helper functions for working with Gaussian 2D peak results
 */
public class Gaussian2DPeakResultHelper
{
	private static class BaseGaussian2DPeakResultCalculator implements Gaussian2DPeakResultCalculator
	{
		final static double twoPi = 2 * Math.PI;

		final CalibrationHelper calibration;
		final int isx, isy;
		boolean oneAxisSD;

		// Set dynamically when needed
		double nmPerPixel;
		TypeConverter<IntensityUnit> toPhoton;
		TypeConverter<DistanceUnit> toPixel;
		TypeConverter<DistanceUnit> toNM;
		boolean emCCD;

		/**
		 * Instantiates a new Gaussian 2D peak result helper.
		 * <p>
		 * The calibration need only contain the information required for the given helper function. It is suggested to
		 * create this helper instance and then call the helper function once to determine if the helper is valid,
		 * catching
		 * and handling the configuration exception as appropriate.
		 * <p>
		 * Note: A factory method is provided to simplify creation for a specific helper function.
		 *
		 * @param psf
		 *            the psf
		 * @param calibration
		 *            the calibration (used for converting the parameters)
		 * @see #create(PSF, CalibrationHelper, int)
		 * @throws ConfigurationException
		 *             If not a Gaussian 2D PSF
		 */
		public BaseGaussian2DPeakResultCalculator(PSF psf, CalibrationHelper calibration) throws ConfigurationException
		{
			int[] indices = PSFHelper.getGaussian2DWxWyIndices(psf);
			isx = indices[0];
			isy = indices[1];
			oneAxisSD = isx == isy;

			// This may not be needed
			this.calibration = calibration;
		}

		private BaseGaussian2DPeakResultCalculator(BaseGaussian2DPeakResultCalculator helper)
		{
			this.isx = helper.isx;
			this.isy = helper.isy;
			this.calibration = helper.calibration;
			this.nmPerPixel = helper.nmPerPixel;
			this.toPhoton = helper.toPhoton;
			this.toPixel = helper.toPixel;
			this.toNM = helper.toNM;
			this.emCCD = helper.emCCD;
		}

		public float getStandardDeviation(float[] params)
		{
			return (oneAxisSD) ? params[isx]
					: (float) Gaussian2DPeakResultHelper.getStandardDeviation(params[isx], params[isy]);
		}

		public float getAmplitude(float[] params) throws ConfigurationException
		{
			// Try to create the converter
			if (toPixel == null)
			{
				if (calibration == null)
					throw new ConfigurationException("Not a valid calibration");
				toPixel = calibration.getDistanceConverter(DistanceUnit.PIXEL);
			}

			return (float) (params[PeakResult.INTENSITY] /
					(twoPi * toPixel.convert(params[isx]) * toPixel.convert(params[isy])));
		}

		private static boolean isCCD(CalibrationHelper calibration)
		{
			return calibration.isCCDCamera();
		}

		public double getPrecision(float[] params, float noise) throws ConfigurationException
		{
			// Try to create the converter
			if (toPhoton == null)
			{
				if (calibration == null || !calibration.hasNmPerPixel() || !isCCD(calibration))
					throw new ConfigurationException("Not a valid calibration");
				nmPerPixel = calibration.getNmPerPixel();
				emCCD = calibration.getCameraType() == CameraType.EMCCD;
				toPhoton = calibration.getIntensityConverter(IntensityUnit.PHOTON);
				toNM = calibration.getDistanceConverter(DistanceUnit.NM);
			}

			return Gaussian2DPeakResultHelper.getPrecision(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(noise), emCCD);
		}

		public double getPrecisionX(float[] params) throws ConfigurationException
		{
			// Try to create the converter
			if (toPhoton == null)
			{
				if (calibration == null || !calibration.hasNmPerPixel() || !isCCD(calibration))
					throw new ConfigurationException("Not a valid calibration");
				nmPerPixel = calibration.getNmPerPixel();
				emCCD = calibration.getCameraType() == CameraType.EMCCD;
				toPhoton = calibration.getIntensityConverter(IntensityUnit.PHOTON);
				toNM = calibration.getDistanceConverter(DistanceUnit.NM);
			}

			return Gaussian2DPeakResultHelper.getPrecisionX(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(params[PeakResult.BACKGROUND]),
					emCCD);
		}

		public double getVariance(float[] params, float noise) throws ConfigurationException
		{
			// Try to create the converter
			if (toPhoton == null)
			{
				if (calibration == null || !calibration.hasNmPerPixel() || !isCCD(calibration))
					throw new ConfigurationException("Not a valid calibration");
				nmPerPixel = calibration.getNmPerPixel();
				emCCD = calibration.getCameraType() == CameraType.EMCCD;
				toPhoton = calibration.getIntensityConverter(IntensityUnit.PHOTON);
				toNM = calibration.getDistanceConverter(DistanceUnit.NM);
			}

			return Gaussian2DPeakResultHelper.getVariance(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(noise), emCCD);
		}

		public double getVarianceX(float[] params) throws ConfigurationException
		{
			// Try to create the converter
			if (toPhoton == null)
			{
				if (calibration == null || !calibration.hasNmPerPixel() || !isCCD(calibration))
					throw new ConfigurationException("Not a valid calibration");
				nmPerPixel = calibration.getNmPerPixel();
				emCCD = calibration.getCameraType() == CameraType.EMCCD;
				toPhoton = calibration.getIntensityConverter(IntensityUnit.PHOTON);
				toNM = calibration.getDistanceConverter(DistanceUnit.NM);
			}

			return Gaussian2DPeakResultHelper.getVarianceX(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(params[PeakResult.BACKGROUND]),
					emCCD);
		}
	}

	private static abstract class BaseFastGaussian2DPeakResultCalculator extends BaseGaussian2DPeakResultCalculator
	{
		public BaseFastGaussian2DPeakResultCalculator(BaseGaussian2DPeakResultCalculator helper)
				throws ConfigurationException
		{
			super(helper);
		}

		@Override
		public abstract float getStandardDeviation(float[] params);

		@Override
		public float getAmplitude(float[] params) throws ConfigurationException
		{
			return (float) (params[PeakResult.INTENSITY] /
					(twoPi * toPixel.convert(params[isx]) * toPixel.convert(params[isy])));
		}

		@Override
		public double getPrecision(float[] params, float noise) throws ConfigurationException
		{
			return Gaussian2DPeakResultHelper.getPrecision(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(noise), emCCD);
		}

		@Override
		public double getPrecisionX(float[] params) throws ConfigurationException
		{
			return Gaussian2DPeakResultHelper.getPrecisionX(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(params[PeakResult.BACKGROUND]),
					emCCD);
		}

		@Override
		public double getVariance(float[] params, float noise) throws ConfigurationException
		{
			return Gaussian2DPeakResultHelper.getVariance(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(noise), emCCD);
		}

		@Override
		public double getVarianceX(float[] params) throws ConfigurationException
		{
			return Gaussian2DPeakResultHelper.getVarianceX(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(params[PeakResult.BACKGROUND]),
					emCCD);
		}
	}

	/**
	 * Private class to allow caching the converters for two-axis Gaussian 2D
	 */
	private static class TwoAxisFastGaussian2DPeakResultCalculator extends BaseFastGaussian2DPeakResultCalculator
	{
		public TwoAxisFastGaussian2DPeakResultCalculator(BaseGaussian2DPeakResultCalculator helper)
				throws ConfigurationException
		{
			super(helper);
		}

		@Override
		public float getStandardDeviation(float[] params)
		{
			return (float) Gaussian2DPeakResultHelper.getStandardDeviation(params[isx], params[isy]);
		}
	}

	/**
	 * Private class to allow caching the converters for one-axis Gaussian 2D
	 */
	private static class OneAxisFastGaussian2DPeakResultCalculator extends BaseFastGaussian2DPeakResultCalculator
	{
		public OneAxisFastGaussian2DPeakResultCalculator(BaseGaussian2DPeakResultCalculator helper)
				throws ConfigurationException
		{
			super(helper);
		}

		@Override
		public float getStandardDeviation(float[] params)
		{
			return params[isx];
		}
	}

	/** Flag for the {@link Gaussian2DPeakResultCalculator#getAmplitude(float[])} function. */
	public static final int AMPLITUDE = 0x00000001;
	/** Flag for the {@link Gaussian2DPeakResultCalculator#getPrecision(float[], float)} function. */
	public static final int PRECISION = 0x00000002;
	/** Flag for the {@link Gaussian2DPeakResultCalculator#getPrecisionX(float[])} function. */
	public static final int PRECISION_X = 0x00000004;

	/** Dummy Gaussian 2D parameters */
	private static final float[] PARAMS = new float[7];

	/** The index of the Sx parameter in the PeakResult parameters array. */
	public static final int INDEX_SX = PeakResult.STANDARD_PARAMETERS;
	/** The index of the Sy parameter in the PeakResult parameters array. */
	public static final int INDEX_SY = PeakResult.STANDARD_PARAMETERS + 1;
	/** The index of the Angle parameter in the PeakResult parameters array. */
	public static final int INDEX_A = PeakResult.STANDARD_PARAMETERS + 2;

	/**
	 * Creates a new Gaussian 2D peak result calculator.
	 * <p>
	 * The calibration need only contain the information required for the specified functions.
	 *
	 * @param psf
	 *            the psf
	 * @param calibration
	 *            the calibration (used for converting the parameters)
	 * @param flags
	 *            the flags specifying the helper functions to support
	 * @return the gaussian 2 D peak result helper
	 * @throws ConfigurationException
	 *             If not a Gaussian 2D PSF or the calibration is invalid
	 */
	public static Gaussian2DPeakResultCalculator create(PSF psf, CalibrationHelper calibration, int flags)
			throws ConfigurationException
	{
		BaseGaussian2DPeakResultCalculator helper = new BaseGaussian2DPeakResultCalculator(psf, calibration);

		// Try the desired methods
		if (BitFlags.anySet(flags, AMPLITUDE))
			helper.getAmplitude(PARAMS);
		if (BitFlags.anySet(flags, PRECISION))
			helper.getPrecision(PARAMS, 0);
		if (BitFlags.anySet(flags, PRECISION_X))
			helper.getPrecisionX(PARAMS);

		// Get a fast implementation
		return (helper.oneAxisSD) ? new OneAxisFastGaussian2DPeakResultCalculator(helper)
				: new TwoAxisFastGaussian2DPeakResultCalculator(helper);
	}

	/**
	 * Calculate the localisation precision for least squares fitting a Gaussian2D PSF to a Gaussian2D PSF. This is an
	 * approximation of the precision of fitting to an optical PSF. Uses the Mortensen formula for an EMCCD camera
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
	 * Calculate the localisation variance for least squares fitting a Gaussian2D PSF to a Gaussian2D PSF. This is an
	 * approximation of the precision of fitting to an optical PSF. Uses the Mortensen formula for an EMCCD camera
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
	 * Calculate the localisation precision for maximum likelihood fitting a Gaussian2D PSF to a Gaussian2D PSF. This is
	 * an approximation of the precision of fitting to an optical PSF. Uses the Mortensen formula for an EMCCD
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
	 * Calculate the localisation variance for maximum likelihood fitting a Gaussian2D PSF to a Gaussian2D PSF. This is
	 * an approximation of the precision of fitting to an optical PSF. Uses the Mortensen formula for an EMCCD
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
	 * Calculate the localisation precision for least squares fitting a Gaussian2D PSF to a Gaussian2D PSF. This is an
	 * approximation of the precision of fitting to an optical PSF for least squares estimation. Uses the Mortensen
	 * formula for an EMCCD camera
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
	 * Calculate the localisation precision for least squares fitting a Gaussian2D PSF to a Gaussian2D PSF. This is an
	 * approximation of the precision of fitting to an optical PSF for least squares estimation. Uses the Mortensen
	 * formula for an EMCCD camera
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
	 * Calculate the localisation variance for least squares fitting a Gaussian2D PSF to a Gaussian2D PSF. This is an
	 * approximation of the precision of fitting to an optical PSF for least squares estimation. Uses the Mortensen
	 * formula for an EMCCD camera
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
	 * Calculate the localisation variance for least squares fitting a Gaussian2D PSF to a Gaussian2D PSF. This is an
	 * approximation of the precision of fitting to an optical PSF for least squares estimation. Uses the Mortensen
	 * formula for an EMCCD camera
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
	 * Calculate the localisation precision for maximum likelihood fitting a Gaussian2D PSF to a Gaussian2D PSF. This is
	 * an approximation of the precision of fitting to an optical PSF. Uses the Mortensen formula for an EMCCD
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
	 * Calculate the localisation precision for maximum likelihood fitting a Gaussian2D PSF to a Gaussian2D PSF. This is
	 * an approximation of the precision of fitting to an optical PSF. Uses the Mortensen formula for an EMCCD
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
	 * Calculate the localisation variance for maximum likelihood fitting a Gaussian2D PSF to a Gaussian2D PSF. This is
	 * an approximation of the precision of fitting to an optical PSF. Uses the Mortensen formula for an EMCCD
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
	 * Calculate the localisation variance for maximum likelihood fitting a Gaussian2D PSF to a Gaussian2D PSF. This is
	 * an approximation of the precision of fitting to an optical PSF. Uses the Mortensen formula for an EMCCD
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

	/**
	 * Get the amplitude of a Gaussian 2D PSF. Amplitude = intensity / (2*pi*sx*sy).
	 *
	 * @param intensity
	 *            the intensity
	 * @param sx
	 *            the sx
	 * @param sy
	 *            the sy
	 */
	public static double getAmplitude(double intensity, double sx, double sy)
	{
		return (intensity / (2 * Math.PI * sx * sy));
	}

	/**
	 * Gets the single Gaussian 2D standard deviation from independent x and y standard deviations. s =
	 * sqrt(abs(sx*sy)).
	 *
	 * @param sx
	 *            the sx
	 * @param sy
	 *            the sy
	 * @return the single Gaussian 2D standard deviation
	 */
	public static double getStandardDeviation(double sx, double sy)
	{
		return Math.sqrt(Math.abs(sx * sy));
	}

	/**
	 * Creates the params array for a Gaussian 2D peak result.
	 * <p>
	 * The psf parameters are used to determine if the PSF is one axis, two axis or two axis and angle. If no PSF
	 * parameters are given then a 1 axis Gaussian is created with a standard deviation of 1.
	 *
	 * @param background
	 *            the background
	 * @param intensity
	 *            the intensity
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param z
	 *            the z
	 * @param psfParameters
	 *            the psf parameters
	 * @return the params
	 * @throws IllegalArgumentException
	 *             If the number of PSF parameters is invalid
	 */
	public static float[] createParams(float background, float intensity, float x, float y, float z,
			float... psfParameters) throws IllegalArgumentException
	{
		if (psfParameters == null)
			return createOneAxisParams(background, intensity, x, y, z, 1);
		switch (psfParameters.length)
		{
			case 1:
				return createOneAxisParams(background, intensity, x, y, z, psfParameters[0]);
			case 2:
				return createTwoAxisParams(background, intensity, x, y, z, psfParameters[0], psfParameters[1]);
			case 3:
				return createTwoAxisAndAngleParams(background, intensity, x, y, z, psfParameters[0], psfParameters[1],
						psfParameters[2]);
			default:
				throw new IllegalArgumentException("Invalid number of Gaussian 2D PSF parameters");
		}
	}

	/**
	 * Creates the params array for a one axis Gaussian 2D peak result.
	 *
	 * @param background
	 *            the background
	 * @param intensity
	 *            the intensity
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param z
	 *            the z
	 * @param s
	 *            the s
	 * @return the params
	 */
	public static float[] createOneAxisParams(float background, float intensity, float x, float y, float z, float s)
	{
		float[] params = new float[PeakResult.STANDARD_PARAMETERS + 1];
		params[PeakResult.BACKGROUND] = background;
		params[PeakResult.INTENSITY] = intensity;
		params[PeakResult.X] = x;
		params[PeakResult.Y] = y;
		params[PeakResult.Z] = z;
		params[INDEX_SX] = s;
		return params;
	}

	/**
	 * Creates the params array for a two axis Gaussian 2D peak result.
	 *
	 * @param background
	 *            the background
	 * @param intensity
	 *            the intensity
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param z
	 *            the z
	 * @param sx
	 *            the sx
	 * @param sy
	 *            the sy
	 * @return the params
	 */
	public static float[] createTwoAxisParams(float background, float intensity, float x, float y, float z, float sx,
			float sy)
	{
		float[] params = new float[PeakResult.STANDARD_PARAMETERS + 2];
		params[PeakResult.BACKGROUND] = background;
		params[PeakResult.INTENSITY] = intensity;
		params[PeakResult.X] = x;
		params[PeakResult.Y] = y;
		params[PeakResult.Z] = z;
		params[INDEX_SX] = sx;
		params[INDEX_SY] = sy;
		return params;
	}

	/**
	 * Creates the params array for a two axis and angle Gaussian 2D peak result.
	 *
	 * @param background
	 *            the background
	 * @param intensity
	 *            the intensity
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param z
	 *            the z
	 * @param sx
	 *            the sx
	 * @param sy
	 *            the sy
	 * @param a
	 *            the a
	 * @return the params
	 */
	public static float[] createTwoAxisAndAngleParams(float background, float intensity, float x, float y, float z,
			float sx, float sy, float a)
	{
		float[] params = new float[PeakResult.STANDARD_PARAMETERS + 3];
		params[PeakResult.BACKGROUND] = background;
		params[PeakResult.INTENSITY] = intensity;
		params[PeakResult.X] = x;
		params[PeakResult.Y] = y;
		params[PeakResult.Z] = z;
		params[INDEX_SX] = sx;
		params[INDEX_SY] = sy;
		params[INDEX_A] = a;
		return params;
	}
}
