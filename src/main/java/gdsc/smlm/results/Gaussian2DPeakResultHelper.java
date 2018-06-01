package gdsc.smlm.results;

import java.util.Arrays;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.exception.TooManyEvaluationsException;

import gdsc.core.data.utils.ConversionException;
import gdsc.core.data.utils.TypeConverter;
import gdsc.core.utils.BitFlags;
import gdsc.smlm.data.config.CalibrationProtos.Calibration;
import gdsc.smlm.data.config.CalibrationProtos.CameraType;
import gdsc.smlm.data.config.CalibrationReader;
import gdsc.smlm.data.config.ConfigurationException;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.PSFProtos.PSFType;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import gdsc.smlm.function.Erf;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;

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
	private static final double ROOT2 = Math.sqrt(2);
	private static final double F1 = cumulative(1);
	private static final double F2 = cumulative(2) / 4;

	private static class BaseGaussian2DPeakResultCalculator implements Gaussian2DPeakResultCalculator
	{
		final static double twoPi = 2 * Math.PI;
		final static double ONE_OVER_ROOT2 = 1.0 / ROOT2;

		final CalibrationReader calibration;
		final int isx, isy;
		private boolean oneAxisSD;

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
		 * @throws ConfigurationException
		 *             If not a Gaussian 2D PSF
		 * @throws ConversionException
		 *             If unit conversion fails
		 * @see #create(PSF, CalibrationReader, int)
		 */
		public BaseGaussian2DPeakResultCalculator(PSF psf, CalibrationReader calibration)
				throws ConfigurationException, ConversionException
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

		public float getStandardDeviation2(float[] params)
		{
			return (oneAxisSD) ? params[isx] * params[isx] : Math.abs(params[isx] * params[isy]);
		}

		public float getAmplitude(float[] params) throws ConfigurationException, ConversionException
		{
			// Try to create the converter
			if (toPixel == null)
			{
				if (calibration == null)
					throw new ConfigurationException("No calibration");
				toPixel = calibration.getDistanceConverter(DistanceUnit.PIXEL);
			}

			return (float) (params[PeakResult.INTENSITY] /
					(twoPi * toPixel.convert(params[isx]) * toPixel.convert(params[isy])));
		}

		public float getPixelAmplitude(float[] params) throws ConfigurationException, ConversionException
		{
			// Try to create the converter
			if (toPixel == null)
			{
				if (calibration == null)
					throw new ConfigurationException("No calibration");
				toPixel = calibration.getDistanceConverter(DistanceUnit.PIXEL);
			}

			return getPixelAmplitudeImpl(params);
		}

		protected float getPixelAmplitudeImpl(float[] params)
		{
			// Get the Gaussian parameters in pixels
			double x = toPixel.convert(params[PeakResult.X]);
			double y = toPixel.convert(params[PeakResult.Y]);
			double sx = toPixel.convert(params[isx]);
			double sy = toPixel.convert(params[isy]);

			return (float) (params[PeakResult.INTENSITY] * 0.25 * gaussianPixelIntegral(x, sx) *
					gaussianPixelIntegral(y, sy));
		}

		/**
		 * Compute the integral of the pixel using the error function
		 *
		 * @param x
		 *            the x
		 * @param s
		 *            the s
		 * @return the double
		 */
		double gaussianPixelIntegral(double x, double s)
		{
			// Find the pixel boundary. Assume 0.5 is the centre of the pixel, we round down 
			// to find the distance to the lower pixel boundary (lx):
			// lx = x - l = x - floor(x)
			//
			// l          x             u
			// | <- lx -> |  <- ux ->   |
			//
			// We compute the integral over the pixel from the lower (l) to the upper (u).
			// erf(-lx, ux)
			double lx = Math.floor(x) - x; // Reversed for convenience, i.e. compute -lx
			return Erf.erf(lx * (ONE_OVER_ROOT2 / s), (lx + 1) * (ONE_OVER_ROOT2 / s));
		}

		public double getLSEPrecision(float[] params, float noise) throws ConfigurationException, ConversionException
		{
			// Try to create the converter
			if (toPhoton == null)
			{
				checkPrecisionCalibration();
				nmPerPixel = calibration.getNmPerPixel();
				emCCD = calibration.getCameraType() == CameraType.EMCCD;
				toPhoton = calibration.getIntensityConverter(IntensityUnit.PHOTON);
				toNM = calibration.getDistanceConverter(DistanceUnit.NM);
			}

			return Gaussian2DPeakResultHelper.getPrecision(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(noise), emCCD);
		}

		private void checkPrecisionCalibration()
		{
			if (calibration == null)
				throw new ConfigurationException("No calibration");
			if (!calibration.hasNmPerPixel())
				throw new ConfigurationException("Not a valid calibration: nm/pixel is required");
			// Note: The Mortensen formula can be used for a sCMOS since that is like a CCD just with 
			// per pixel read noise. The noise component then should represent an average across the
			// region used to fit the data.
			//if (!calibration.isCCDCamera())
			//	throw new ConfigurationException("Not a valid calibration: CCD/EM-CCD camera type is required");
		}

		public double getLSEPrecision(float[] params) throws ConfigurationException, ConversionException
		{
			// Try to create the converter
			if (toPhoton == null)
			{
				checkPrecisionCalibration();
				nmPerPixel = calibration.getNmPerPixel();
				emCCD = calibration.getCameraType() == CameraType.EMCCD;
				toPhoton = calibration.getIntensityConverter(IntensityUnit.PHOTON);
				toNM = calibration.getDistanceConverter(DistanceUnit.NM);
			}

			return Gaussian2DPeakResultHelper.getPrecisionX(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(params[PeakResult.BACKGROUND]),
					emCCD);
		}

		public double getLSEVariance(float[] params, float noise) throws ConfigurationException, ConversionException
		{
			// Try to create the converter
			if (toPhoton == null)
			{
				checkPrecisionCalibration();
				nmPerPixel = calibration.getNmPerPixel();
				emCCD = calibration.getCameraType() == CameraType.EMCCD;
				toPhoton = calibration.getIntensityConverter(IntensityUnit.PHOTON);
				toNM = calibration.getDistanceConverter(DistanceUnit.NM);
			}

			return Gaussian2DPeakResultHelper.getVariance(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(noise), emCCD);
		}

		public double getLSEVariance(float[] params) throws ConfigurationException, ConversionException
		{
			// Try to create the converter
			if (toPhoton == null)
			{
				checkPrecisionCalibration();
				nmPerPixel = calibration.getNmPerPixel();
				emCCD = calibration.getCameraType() == CameraType.EMCCD;
				toPhoton = calibration.getIntensityConverter(IntensityUnit.PHOTON);
				toNM = calibration.getDistanceConverter(DistanceUnit.NM);
			}

			return Gaussian2DPeakResultHelper.getVarianceX(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(params[PeakResult.BACKGROUND]),
					emCCD);
		}

		public double getMLEPrecision(float[] params, float noise) throws ConfigurationException, ConversionException
		{
			// Try to create the converter
			if (toPhoton == null)
			{
				checkPrecisionCalibration();
				nmPerPixel = calibration.getNmPerPixel();
				emCCD = calibration.getCameraType() == CameraType.EMCCD;
				toPhoton = calibration.getIntensityConverter(IntensityUnit.PHOTON);
				toNM = calibration.getDistanceConverter(DistanceUnit.NM);
			}

			return Gaussian2DPeakResultHelper.getMLPrecision(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(noise), emCCD);
		}

		public double getMLEPrecision(float[] params)
				throws ConfigurationException, ConversionException, ConversionException
		{
			// Try to create the converter
			if (toPhoton == null)
			{
				checkPrecisionCalibration();
				nmPerPixel = calibration.getNmPerPixel();
				emCCD = calibration.getCameraType() == CameraType.EMCCD;
				toPhoton = calibration.getIntensityConverter(IntensityUnit.PHOTON);
				toNM = calibration.getDistanceConverter(DistanceUnit.NM);
			}

			return Gaussian2DPeakResultHelper.getMLPrecisionX(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(params[PeakResult.BACKGROUND]),
					emCCD);
		}

		public double getMLEVariance(float[] params, float noise)
				throws ConfigurationException, ConversionException, ConversionException
		{
			// Try to create the converter
			if (toPhoton == null)
			{
				checkPrecisionCalibration();
				nmPerPixel = calibration.getNmPerPixel();
				emCCD = calibration.getCameraType() == CameraType.EMCCD;
				toPhoton = calibration.getIntensityConverter(IntensityUnit.PHOTON);
				toNM = calibration.getDistanceConverter(DistanceUnit.NM);
			}

			return Gaussian2DPeakResultHelper.getMLVariance(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(noise), emCCD);
		}

		public double getMLEVariance(float[] params)
				throws ConfigurationException, ConversionException, ConversionException
		{
			// Try to create the converter
			if (toPhoton == null)
			{
				checkPrecisionCalibration();
				nmPerPixel = calibration.getNmPerPixel();
				emCCD = calibration.getCameraType() == CameraType.EMCCD;
				toPhoton = calibration.getIntensityConverter(IntensityUnit.PHOTON);
				toNM = calibration.getDistanceConverter(DistanceUnit.NM);
			}

			return Gaussian2DPeakResultHelper.getMLVarianceX(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(params[PeakResult.BACKGROUND]),
					emCCD);
		}
	}

	private static abstract class BaseFastGaussian2DPeakResultCalculator extends BaseGaussian2DPeakResultCalculator
	{
		public BaseFastGaussian2DPeakResultCalculator(BaseGaussian2DPeakResultCalculator helper)
				throws ConfigurationException, ConversionException
		{
			super(helper);
		}

		@Override
		public abstract float getStandardDeviation(float[] params);

		@Override
		public float getAmplitude(float[] params) throws ConfigurationException, ConversionException
		{
			return (float) (params[PeakResult.INTENSITY] /
					(twoPi * toPixel.convert(params[isx]) * toPixel.convert(params[isy])));
		}

		public float getPixelAmplitude(float[] params) throws ConfigurationException, ConversionException
		{
			return getPixelAmplitudeImpl(params);
		}

		@Override
		public double getLSEPrecision(float[] params, float noise) throws ConfigurationException, ConversionException
		{
			return Gaussian2DPeakResultHelper.getPrecision(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(noise), emCCD);
		}

		@Override
		public double getLSEPrecision(float[] params) throws ConfigurationException, ConversionException
		{
			return Gaussian2DPeakResultHelper.getPrecisionX(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(params[PeakResult.BACKGROUND]),
					emCCD);
		}

		@Override
		public double getLSEVariance(float[] params, float noise) throws ConfigurationException, ConversionException
		{
			return Gaussian2DPeakResultHelper.getVariance(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(noise), emCCD);
		}

		@Override
		public double getLSEVariance(float[] params) throws ConfigurationException, ConversionException
		{
			return Gaussian2DPeakResultHelper.getVarianceX(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(params[PeakResult.BACKGROUND]),
					emCCD);
		}

		@Override
		public double getMLEPrecision(float[] params, float noise) throws ConfigurationException, ConversionException
		{
			return Gaussian2DPeakResultHelper.getMLPrecision(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(noise), emCCD);
		}

		@Override
		public double getMLEPrecision(float[] params) throws ConfigurationException, ConversionException
		{
			return Gaussian2DPeakResultHelper.getMLPrecisionX(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(params[PeakResult.BACKGROUND]),
					emCCD);
		}

		@Override
		public double getMLEVariance(float[] params, float noise) throws ConfigurationException, ConversionException
		{
			return Gaussian2DPeakResultHelper.getMLVariance(nmPerPixel, toNM.convert(getStandardDeviation(params)),
					toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(noise), emCCD);
		}

		@Override
		public double getMLEVariance(float[] params) throws ConfigurationException, ConversionException
		{
			return Gaussian2DPeakResultHelper.getMLVarianceX(nmPerPixel, toNM.convert(getStandardDeviation(params)),
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
				throws ConfigurationException, ConversionException
		{
			super(helper);
		}

		@Override
		public float getStandardDeviation(float[] params)
		{
			return (float) Gaussian2DPeakResultHelper.getStandardDeviation(params[isx], params[isy]);
		}

		@Override
		public float getStandardDeviation2(float[] params)
		{
			return Math.abs(params[isx] * params[isy]);
		}
	}

	/**
	 * Private class to allow caching the converters for one-axis Gaussian 2D
	 */
	private static class OneAxisFastGaussian2DPeakResultCalculator extends BaseFastGaussian2DPeakResultCalculator
	{
		public OneAxisFastGaussian2DPeakResultCalculator(BaseGaussian2DPeakResultCalculator helper)
				throws ConfigurationException, ConversionException
		{
			super(helper);
		}

		@Override
		public float getStandardDeviation(float[] params)
		{
			return params[isx];
		}

		@Override
		public float getStandardDeviation2(float[] params)
		{
			return params[isx] * params[isx];
		}
	}

	/** Flag for the {@link Gaussian2DPeakResultCalculator#getAmplitude(float[])} function. */
	public static final int AMPLITUDE = 0x00000001;
	/** Flag for the {@link Gaussian2DPeakResultCalculator#getLSEPrecision(float[], float)} function. */
	public static final int LSE_PRECISION = 0x00000002;
	/** Flag for the {@link Gaussian2DPeakResultCalculator#getLSEPrecision(float[])} function. */
	public static final int LSE_PRECISION_X = 0x00000004;
	/** Flag for the {@link Gaussian2DPeakResultCalculator#getMLEPrecision(float[], float)} function. */
	public static final int MLE_PRECISION = 0x00000008;
	/** Flag for the {@link Gaussian2DPeakResultCalculator#getMLEPrecision(float[])} function. */
	public static final int MLE_PRECISION_X = 0x00000010;
	/** Flag for the {@link Gaussian2DPeakResultCalculator#getPixelAmplitude(float[])} function. */
	public static final int PIXEL_AMPLITUDE = 0x00000020;

	/** Dummy Gaussian 2D parameters */
	private static final float[] PARAMS = new float[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];

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
	 * @throws ConversionException
	 *             If unit conversion fails
	 */
	public static Gaussian2DPeakResultCalculator create(PSF psf, Calibration calibration, int flags)
			throws ConfigurationException, ConversionException
	{
		CalibrationReader helper = (calibration != null) ? new CalibrationReader(calibration) : null;
		return create(psf, helper, flags);
	}

	/**
	 * Creates a new Gaussian 2D peak result calculator.
	 * <p>
	 * The calibration need only contain the information required for the specified functions.
	 *
	 * @param psf
	 *            the psf
	 * @param calibrationReader
	 *            the calibration reader (used for converting the parameters)
	 * @param flags
	 *            the flags specifying the helper functions to support
	 * @return the gaussian 2 D peak result helper
	 * @throws ConfigurationException
	 *             If not a Gaussian 2D PSF or the calibration is invalid
	 * @throws ConversionException
	 *             If unit conversion fails
	 */
	public static Gaussian2DPeakResultCalculator create(PSF psf, CalibrationReader calibrationReader, int flags)
			throws ConfigurationException, ConversionException
	{
		BaseGaussian2DPeakResultCalculator helper = new BaseGaussian2DPeakResultCalculator(psf, calibrationReader);

		// Try the desired methods
		if (BitFlags.anySet(flags, AMPLITUDE))
			helper.getAmplitude(PARAMS);
		if (BitFlags.anySet(flags, LSE_PRECISION))
			helper.getLSEPrecision(PARAMS, 0);
		if (BitFlags.anySet(flags, LSE_PRECISION_X))
			helper.getLSEPrecision(PARAMS);
		if (BitFlags.anySet(flags, MLE_PRECISION))
			helper.getMLEPrecision(PARAMS, 0);
		if (BitFlags.anySet(flags, MLE_PRECISION_X))
			helper.getMLEPrecision(PARAMS);
		if (BitFlags.anySet(flags, PIXEL_AMPLITUDE))
			helper.getPixelAmplitude(PARAMS);

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
	 * Gets the single Gaussian 2D standard deviation squared from independent x and y standard deviations. s2 =
	 * abs(sx*sy). This avoids a sqrt in the computation of {@link #getStandardDeviation(double, double)}.
	 *
	 * @param sx
	 *            the sx
	 * @param sy
	 *            the sy
	 * @return the single Gaussian 2D standard deviation squared
	 */
	public static double getStandardDeviation2(double sx, double sy)
	{
		return Math.abs(sx * sy);
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

	/**
	 * Creates the params array for a Gaussian 2D peak result using the PSF type.
	 *
	 * @param psfType
	 *            the psf type
	 * @param gaussian2DParams
	 *            the full gaussian 2D function parameters
	 * @return the PeakResult params
	 * @throws IllegalArgumentException
	 *             if the psf is not a Gaussian 2D or the input parameters are not the correct length for a Gaussian 2D
	 *             function
	 */
	public static float[] createParams(PSFType psfType, float[] gaussian2DParams) throws IllegalArgumentException
	{
		if (gaussian2DParams.length != PeakResult.STANDARD_PARAMETERS + 3)
			throw new IllegalArgumentException("Parameters must be a full Gaussian 2D parameters array of length " +
					(PeakResult.STANDARD_PARAMETERS + 3));

		switch (psfType)
		{
			case ONE_AXIS_GAUSSIAN_2D:
				return Arrays.copyOf(gaussian2DParams, PeakResult.STANDARD_PARAMETERS + 1);

			case ASTIGMATIC_GAUSSIAN_2D:
			case TWO_AXIS_GAUSSIAN_2D:
				return Arrays.copyOf(gaussian2DParams, PeakResult.STANDARD_PARAMETERS + 2);

			case TWO_AXIS_AND_THETA_GAUSSIAN_2D:
				return gaussian2DParams;

			default:
				throw new IllegalArgumentException("PSF type must be a Gaussian 2D PSF");
		}
	}

	/**
	 * Gets the Signal-to-Noise Ratio (SNR). This is ratio of the average signal value to the standard deviation of the
	 * background. Ideally the standard deviation of the background is computed in the region around the Gaussian.
	 * <p>
	 * The average signal value is taken using the expected sum of the Gaussian within the range divided by the
	 * elliptial area of the same range. For example a range of 1 would expect the Gaussian to be 0.6827 of the total
	 * intensity; the area would be pi * sx * sx; and the SNR = intensity * 0.6827 / (pi * sx * sy * noise).
	 * <p>
	 * As an alternative definition, the standard deviation of the background can be computed using the standard
	 * deviation of the signal in the region around the Gaussian.
	 * <p>
	 * Note: Arguments are not checked
	 *
	 * @param intensity
	 *            the total intensity of the Gaussian
	 * @param sx
	 *            the Gaussian standard deviation in the X dimension
	 * @param sy
	 *            the Gaussian standard deviation in the Y dimension
	 * @param range
	 *            the range over which to compute the average signal
	 * @param noise
	 *            the noise (standard deviation of the background)
	 * @return the snr
	 * @see <a href=
	 *      "https://en.wikipedia.org/wiki/Signal-to-noise_ratio_(imaging)">https://en.wikipedia.org/wiki/Signal-to-noise_ratio_(imaging)</a>
	 */
	public static double getSNR(double intensity, double sx, double sy, double range, double noise)
	{
		// Get the expected total signal using the Error function assuming a standard deviation of 1
		double exp = cumulative(range);
		return intensity * exp / (Math.PI * sx * sy * range * range * noise);
	}

	/**
	 * Compute the cumulative normal distribution within the range -x to x:
	 * 
	 * <pre>
	 * return erf(x / sqrt(2));
	 * </pre>
	 * 
	 * This uses a fast approximation to the Error function.
	 *
	 * @param x
	 *            the x
	 * @return the cumulative normal distribution within the range -x to x
	 */
	public static double cumulative(double x)
	{
		return Erf.erf(x / ROOT2);
	}

	/**
	 * Gets the Signal-to-Noise Ratio (SNR). This is ratio of the average signal value to the standard deviation of the
	 * background. Ideally the standard deviation of the background is computed in the region around the Gaussian.
	 * <p>
	 * The average signal value is taken using the expected sum of the Gaussian within a range of 1 SD divided by the
	 * elliptial area of the same range. For example a range of 1 would expect the Gaussian to be 0.6827 of the total
	 * intensity; the area would be pi * sx * sx; and the SNR = intensity * 0.6827 / (pi * sx * sy * noise).
	 * <p>
	 * As an alternative definition, the standard deviation of the background can be computed using the standard
	 * deviation of the signal in the region around the Gaussian.
	 * <p>
	 * Note: Arguments are not checked
	 *
	 * @param intensity
	 *            the total intensity of the Gaussian
	 * @param sx
	 *            the Gaussian standard deviation in the X dimension
	 * @param sy
	 *            the Gaussian standard deviation in the Y dimension
	 * @param range
	 *            the range over which to compute the average signal
	 * @param noise
	 *            the noise (standard deviation of the background)
	 * @return the snr
	 * @see <a href=
	 *      "https://en.wikipedia.org/wiki/Signal-to-noise_ratio_(imaging)">https://en.wikipedia.org/wiki/Signal-to-noise_ratio_(imaging)</a>
	 */
	public static double getSNR1(double intensity, double sx, double sy, double noise)
	{
		return intensity * F1 / (Math.PI * sx * sy * noise);
	}

	/**
	 * Gets the Signal-to-Noise Ratio (SNR). This is ratio of the average signal value to the standard deviation of the
	 * background. Ideally the standard deviation of the background is computed in the region around the Gaussian.
	 * <p>
	 * The average signal value is taken using the expected sum of the Gaussian within a range of 2 SD divided by the
	 * elliptial area of the same range. For example a range of 2 would expect the Gaussian to be 0.9545 of the total
	 * intensity; the area would be pi * sx * 2 * sx * 2; and the SNR = intensity * 0.9545 / (pi * sx * sy * 4 * noise).
	 * <p>
	 * As an alternative definition, the standard deviation of the background can be computed using the standard
	 * deviation of the signal in the region around the Gaussian.
	 * <p>
	 * Note: Arguments are not checked
	 *
	 * @param intensity
	 *            the total intensity of the Gaussian
	 * @param sx
	 *            the Gaussian standard deviation in the X dimension
	 * @param sy
	 *            the Gaussian standard deviation in the Y dimension
	 * @param range
	 *            the range over which to compute the average signal
	 * @param noise
	 *            the noise (standard deviation of the background)
	 * @return the snr
	 * @see <a href=
	 *      "https://en.wikipedia.org/wiki/Signal-to-noise_ratio_(imaging)">https://en.wikipedia.org/wiki/Signal-to-noise_ratio_(imaging)</a>
	 */
	public static double getSNR2(double intensity, double sx, double sy, double noise)
	{
		return intensity * F2 / (Math.PI * sx * sy * noise);
	}
}
