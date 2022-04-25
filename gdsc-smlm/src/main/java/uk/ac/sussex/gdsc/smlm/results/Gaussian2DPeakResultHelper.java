/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.results;

import java.util.Arrays;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import uk.ac.sussex.gdsc.core.data.utils.ConversionException;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.utils.BitFlagUtils;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.ConfigurationException;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFOrBuilder;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.data.config.PsfHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.function.Erf;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;

/**
 * Contains helper functions for working with Gaussian 2D peak results.
 */
public final class Gaussian2DPeakResultHelper {
  // Allow parameters to be name following standard conventions:
  // s = standard deviation
  // N = photons
  // b = background
  // CHECKSTYLE.OFF: ParameterName

  /** sqrt(0.5), or 1 / sqrt(2). */
  static final double ONE_OVER_ROOT2 = 0.7071067811865475244008443621048490392;
  private static final double R1 = cumulative2D(1) / Math.PI;
  private static final double R2 = cumulative2D(2) / (Math.PI * 4);
  /**
   * The Mahalanobis distance r for a 2D Gaussian that contains 50 percent of the integral.
   */
  public static final double R_2D_50 = inverseCumulative2D(0.5);
  private static final double P05 = 0.5 / (Math.PI * MathUtils.pow2(R_2D_50));

  /**
   * Flag for the {@link Gaussian2DPeakResultCalculator#getAmplitude(float[])} function.
   */
  public static final int AMPLITUDE = 0x00000001;
  /**
   * Flag for the {@link Gaussian2DPeakResultCalculator#getLsePrecision(float[], float)} function.
   */
  public static final int LSE_PRECISION = 0x00000002;
  /**
   * Flag for the {@link Gaussian2DPeakResultCalculator#getLsePrecision(float[])} function.
   */
  public static final int LSE_PRECISION_X = 0x00000004;
  /**
   * Flag for the {@link Gaussian2DPeakResultCalculator#getMlePrecision(float[], float)} function.
   */
  public static final int MLE_PRECISION = 0x00000008;
  /**
   * Flag for the {@link Gaussian2DPeakResultCalculator#getMlePrecision(float[])} function.
   */
  public static final int MLE_PRECISION_X = 0x00000010;
  /**
   * Flag for the {@link Gaussian2DPeakResultCalculator#getPixelAmplitude(float[])} function.
   */
  public static final int PIXEL_AMPLITUDE = 0x00000020;

  /** Dummy Gaussian 2D parameters. */
  private static final float[] PARAMS = new float[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];

  /** The index of the Sx parameter in the PeakResult parameters array. */
  public static final int INDEX_SX = PeakResult.STANDARD_PARAMETERS;
  /** The index of the Sy parameter in the PeakResult parameters array. */
  public static final int INDEX_SY = PeakResult.STANDARD_PARAMETERS + 1;
  /** The index of the Angle parameter in the PeakResult parameters array. */
  public static final int INDEX_A = PeakResult.STANDARD_PARAMETERS + 2;

  /**
   * The default points to use for maximum likelihood precision computation, see
   * {@link #getMLVarianceX(double, double, double, double, boolean, int)}
   *
   * <p>Testing shows that 10 integration points is the fastest for realistic input parameters.
   */
  public static final int POINTS = 10;

  /**
   * Base class for computing Gaussian 2D peak result data.
   */
  private static class BaseGaussian2DPeakResultCalculator
      implements Gaussian2DPeakResultCalculator {
    private static final String NO_CALIBRATION = "No calibration";

    static final double TWO_PI = 2 * Math.PI;

    final CalibrationReader calibration;
    final int isx;
    final int isy;
    boolean oneAxisSD;

    // Set dynamically when needed
    double nmPerPixel;
    TypeConverter<IntensityUnit> toPhoton;
    TypeConverter<DistanceUnit> toPixel;
    TypeConverter<DistanceUnit> toNM;
    boolean emCcd;

    /**
     * Instantiates a new Gaussian 2D peak result helper.
     *
     * <p>The calibration need only contain the information required for the given helper function.
     * It is suggested to create this helper instance and then call the helper function once to
     * determine if the helper is valid, catching and handling the configuration exception as
     * appropriate.
     *
     * <p>Note: A factory method is provided to simplify creation for a specific helper function.
     *
     * @param psf the psf
     * @param calibration the calibration (used for converting the parameters)
     * @throws ConfigurationException If not a Gaussian 2D PSF
     * @throws ConversionException If unit conversion fails
     * @see #create(PSF, CalibrationReader, int)
     */
    BaseGaussian2DPeakResultCalculator(PSF psf, CalibrationReader calibration) {
      final int[] indices = PsfHelper.getGaussian2DWxWyIndices(psf);
      isx = indices[0];
      isy = indices[1];
      oneAxisSD = isx == isy;

      // This may not be needed
      this.calibration = calibration;
    }

    private BaseGaussian2DPeakResultCalculator(BaseGaussian2DPeakResultCalculator helper) {
      this.isx = helper.isx;
      this.isy = helper.isy;
      this.calibration = helper.calibration;
      this.nmPerPixel = helper.nmPerPixel;
      this.toPhoton = helper.toPhoton;
      this.toPixel = helper.toPixel;
      this.toNM = helper.toNM;
      this.emCcd = helper.emCcd;
    }

    @Override
    public float getStandardDeviation(float[] params) {
      return (oneAxisSD) ? params[isx]
          : (float) Gaussian2DPeakResultHelper.getStandardDeviation(params[isx], params[isy]);
    }

    @Override
    public float getStandardDeviation2(float[] params) {
      return (oneAxisSD) ? params[isx] * params[isx] : Math.abs(params[isx] * params[isy]);
    }

    @Override
    public float getAmplitude(float[] params) {
      // Try to create the converter
      if (toPixel == null) {
        if (calibration == null) {
          throw new ConfigurationException(NO_CALIBRATION);
        }
        toPixel = calibration.getDistanceConverter(DistanceUnit.PIXEL);
      }

      return (float) (params[PeakResult.INTENSITY]
          / (TWO_PI * toPixel.convert(params[isx]) * toPixel.convert(params[isy])));
    }

    @Override
    public float getPixelAmplitude(float[] params) {
      // Try to create the converter
      if (toPixel == null) {
        if (calibration == null) {
          throw new ConfigurationException(NO_CALIBRATION);
        }
        toPixel = calibration.getDistanceConverter(DistanceUnit.PIXEL);
      }

      return getPixelAmplitudeImpl(params);
    }

    /**
     * Gets the pixel amplitude.
     *
     * @param params the params
     * @return the pixel amplitude
     */
    float getPixelAmplitudeImpl(float[] params) {
      // Get the Gaussian parameters in pixels
      final double x = toPixel.convert(params[PeakResult.X]);
      final double y = toPixel.convert(params[PeakResult.Y]);
      final double sx = toPixel.convert(params[isx]);
      final double sy = toPixel.convert(params[isy]);

      return (float) (params[PeakResult.INTENSITY] * 0.25 * gaussianPixelIntegral(x, sx)
          * gaussianPixelIntegral(y, sy));
    }

    /**
     * Compute the integral of the pixel using the error function.
     *
     * @param x the x
     * @param s the standard deviation
     * @return the integral
     */
    double gaussianPixelIntegral(double x, double s) {
      // Find the pixel boundary. Assume 0.5 is the centre of the pixel, we round down
      // to find the distance to the lower pixel boundary (lx):
      // lx = x - l = x - floor(x)
      //
      // l x u
      // | <- lx -> | <- ux -> |
      //
      // We compute the integral over the pixel from the lower (l) to the upper (u).
      // erf(-lx, ux)
      final double lx = Math.floor(x) - x; // Reversed for convenience, i.e. compute -lx
      return Erf.erf(lx * (ONE_OVER_ROOT2 / s), (lx + 1) * (ONE_OVER_ROOT2 / s));
    }

    @Override
    public double getLsePrecision(float[] params, float noise) {
      // Try to create the converter
      if (toPhoton == null) {
        checkPrecisionCalibration();
        nmPerPixel = calibration.getNmPerPixel();
        emCcd = calibration.getCameraType() == CameraType.EMCCD;
        toPhoton = calibration.getIntensityConverter(IntensityUnit.PHOTON);
        toNM = calibration.getDistanceConverter(DistanceUnit.NM);
      }

      return Gaussian2DPeakResultHelper.getPrecision(nmPerPixel,
          toNM.convert(getStandardDeviation(params)),
          toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(noise), emCcd);
    }

    @Override
    public double getLsePrecision(float[] params) {
      // Try to create the converter
      if (toPhoton == null) {
        checkPrecisionCalibration();
        nmPerPixel = calibration.getNmPerPixel();
        emCcd = calibration.getCameraType() == CameraType.EMCCD;
        toPhoton = calibration.getIntensityConverter(IntensityUnit.PHOTON);
        toNM = calibration.getDistanceConverter(DistanceUnit.NM);
      }

      return Gaussian2DPeakResultHelper.getPrecisionX(nmPerPixel,
          toNM.convert(getStandardDeviation(params)),
          toPhoton.convert(params[PeakResult.INTENSITY]),
          toPhoton.convert(params[PeakResult.BACKGROUND]), emCcd);
    }

    @Override
    public double getLseVariance(float[] params, float noise) {
      // Try to create the converter
      if (toPhoton == null) {
        checkPrecisionCalibration();
        nmPerPixel = calibration.getNmPerPixel();
        emCcd = calibration.getCameraType() == CameraType.EMCCD;
        toPhoton = calibration.getIntensityConverter(IntensityUnit.PHOTON);
        toNM = calibration.getDistanceConverter(DistanceUnit.NM);
      }

      return Gaussian2DPeakResultHelper.getVariance(nmPerPixel,
          toNM.convert(getStandardDeviation(params)),
          toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(noise), emCcd);
    }

    @Override
    public double getLseVariance(float[] params) {
      // Try to create the converter
      if (toPhoton == null) {
        checkPrecisionCalibration();
        nmPerPixel = calibration.getNmPerPixel();
        emCcd = calibration.getCameraType() == CameraType.EMCCD;
        toPhoton = calibration.getIntensityConverter(IntensityUnit.PHOTON);
        toNM = calibration.getDistanceConverter(DistanceUnit.NM);
      }

      return Gaussian2DPeakResultHelper.getVarianceX(nmPerPixel,
          toNM.convert(getStandardDeviation(params)),
          toPhoton.convert(params[PeakResult.INTENSITY]),
          toPhoton.convert(params[PeakResult.BACKGROUND]), emCcd);
    }

    @Override
    public double getMlePrecision(float[] params, float noise) {
      // Try to create the converter
      if (toPhoton == null) {
        checkPrecisionCalibration();
        nmPerPixel = calibration.getNmPerPixel();
        emCcd = calibration.getCameraType() == CameraType.EMCCD;
        toPhoton = calibration.getIntensityConverter(IntensityUnit.PHOTON);
        toNM = calibration.getDistanceConverter(DistanceUnit.NM);
      }

      return Gaussian2DPeakResultHelper.getMLPrecision(nmPerPixel,
          toNM.convert(getStandardDeviation(params)),
          toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(noise), emCcd);
    }

    @Override
    public double getMlePrecision(float[] params) {
      // Try to create the converter
      if (toPhoton == null) {
        checkPrecisionCalibration();
        nmPerPixel = calibration.getNmPerPixel();
        emCcd = calibration.getCameraType() == CameraType.EMCCD;
        toPhoton = calibration.getIntensityConverter(IntensityUnit.PHOTON);
        toNM = calibration.getDistanceConverter(DistanceUnit.NM);
      }

      return Gaussian2DPeakResultHelper.getMLPrecisionX(nmPerPixel,
          toNM.convert(getStandardDeviation(params)),
          toPhoton.convert(params[PeakResult.INTENSITY]),
          toPhoton.convert(params[PeakResult.BACKGROUND]), emCcd);
    }

    @Override
    public double getMleVariance(float[] params, float noise) {
      // Try to create the converter
      if (toPhoton == null) {
        checkPrecisionCalibration();
        nmPerPixel = calibration.getNmPerPixel();
        emCcd = calibration.getCameraType() == CameraType.EMCCD;
        toPhoton = calibration.getIntensityConverter(IntensityUnit.PHOTON);
        toNM = calibration.getDistanceConverter(DistanceUnit.NM);
      }

      return Gaussian2DPeakResultHelper.getMLVariance(nmPerPixel,
          toNM.convert(getStandardDeviation(params)),
          toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(noise), emCcd);
    }

    @Override
    public double getMleVariance(float[] params) {
      // Try to create the converter
      if (toPhoton == null) {
        checkPrecisionCalibration();
        nmPerPixel = calibration.getNmPerPixel();
        emCcd = calibration.getCameraType() == CameraType.EMCCD;
        toPhoton = calibration.getIntensityConverter(IntensityUnit.PHOTON);
        toNM = calibration.getDistanceConverter(DistanceUnit.NM);
      }

      return Gaussian2DPeakResultHelper.getMLVarianceX(nmPerPixel,
          toNM.convert(getStandardDeviation(params)),
          toPhoton.convert(params[PeakResult.INTENSITY]),
          toPhoton.convert(params[PeakResult.BACKGROUND]), emCcd);
    }

    private void checkPrecisionCalibration() {
      if (calibration == null) {
        throw new ConfigurationException(NO_CALIBRATION);
      }
      if (!calibration.hasNmPerPixel()) {
        throw new ConfigurationException("Not a valid calibration: nm/pixel is required");
        // Note: The Mortensen formula can be used for a sCMOS since that is like a CCD just with
        // per pixel read noise. The noise component then should represent an average across the
        // region used to fit the data.
        // if (!calibration.isCCDCamera())
        // throw new ConfigurationException("Not a valid calibration: CCD/EM-CCD camera type is
        // required");
      }
    }
  }

  /**
   * Private class to allow caching the converters from an instance of
   * BaseGaussian2DPeakResultCalculator. The instance must have been used for each target
   * method to ensure the converter(s) are not null.
   */
  private abstract static class BaseFastGaussian2DPeakResultCalculator
      extends BaseGaussian2DPeakResultCalculator {
    BaseFastGaussian2DPeakResultCalculator(BaseGaussian2DPeakResultCalculator helper) {
      super(helper);
    }

    @Override
    public abstract float getStandardDeviation(float[] params);

    @Override
    public float getAmplitude(float[] params) {
      return (float) (params[PeakResult.INTENSITY]
          / (TWO_PI * toPixel.convert(params[isx]) * toPixel.convert(params[isy])));
    }

    @Override
    public float getPixelAmplitude(float[] params) {
      return getPixelAmplitudeImpl(params);
    }

    @Override
    public double getLsePrecision(float[] params, float noise) {
      return Gaussian2DPeakResultHelper.getPrecision(nmPerPixel,
          toNM.convert(getStandardDeviation(params)),
          toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(noise), emCcd);
    }

    @Override
    public double getLsePrecision(float[] params) {
      return Gaussian2DPeakResultHelper.getPrecisionX(nmPerPixel,
          toNM.convert(getStandardDeviation(params)),
          toPhoton.convert(params[PeakResult.INTENSITY]),
          toPhoton.convert(params[PeakResult.BACKGROUND]), emCcd);
    }

    @Override
    public double getLseVariance(float[] params, float noise) {
      return Gaussian2DPeakResultHelper.getVariance(nmPerPixel,
          toNM.convert(getStandardDeviation(params)),
          toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(noise), emCcd);
    }

    @Override
    public double getLseVariance(float[] params) {
      return Gaussian2DPeakResultHelper.getVarianceX(nmPerPixel,
          toNM.convert(getStandardDeviation(params)),
          toPhoton.convert(params[PeakResult.INTENSITY]),
          toPhoton.convert(params[PeakResult.BACKGROUND]), emCcd);
    }

    @Override
    public double getMlePrecision(float[] params, float noise) {
      return Gaussian2DPeakResultHelper.getMLPrecision(nmPerPixel,
          toNM.convert(getStandardDeviation(params)),
          toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(noise), emCcd);
    }

    @Override
    public double getMlePrecision(float[] params) {
      return Gaussian2DPeakResultHelper.getMLPrecisionX(nmPerPixel,
          toNM.convert(getStandardDeviation(params)),
          toPhoton.convert(params[PeakResult.INTENSITY]),
          toPhoton.convert(params[PeakResult.BACKGROUND]), emCcd);
    }

    @Override
    public double getMleVariance(float[] params, float noise) {
      return Gaussian2DPeakResultHelper.getMLVariance(nmPerPixel,
          toNM.convert(getStandardDeviation(params)),
          toPhoton.convert(params[PeakResult.INTENSITY]), toPhoton.convert(noise), emCcd);
    }

    @Override
    public double getMleVariance(float[] params) {
      return Gaussian2DPeakResultHelper.getMLVarianceX(nmPerPixel,
          toNM.convert(getStandardDeviation(params)),
          toPhoton.convert(params[PeakResult.INTENSITY]),
          toPhoton.convert(params[PeakResult.BACKGROUND]), emCcd);
    }
  }

  /**
   * Private class to allow caching the converters for two-axis Gaussian 2D.
   */
  private static class TwoAxisFastGaussian2DPeakResultCalculator
      extends BaseFastGaussian2DPeakResultCalculator {
    TwoAxisFastGaussian2DPeakResultCalculator(BaseGaussian2DPeakResultCalculator helper) {
      super(helper);
    }

    @Override
    public float getStandardDeviation(float[] params) {
      return (float) Gaussian2DPeakResultHelper.getStandardDeviation(params[isx], params[isy]);
    }

    @Override
    public float getStandardDeviation2(float[] params) {
      return Math.abs(params[isx] * params[isy]);
    }
  }

  /**
   * Private class to allow caching the converters for one-axis Gaussian 2D.
   */
  private static class OneAxisFastGaussian2DPeakResultCalculator
      extends BaseFastGaussian2DPeakResultCalculator {
    OneAxisFastGaussian2DPeakResultCalculator(BaseGaussian2DPeakResultCalculator helper) {
      super(helper);
    }

    @Override
    public float getStandardDeviation(float[] params) {
      return params[isx];
    }

    @Override
    public float getStandardDeviation2(float[] params) {
      return params[isx] * params[isx];
    }
  }

  /**
   * No public constructor.
   */
  private Gaussian2DPeakResultHelper() {}

  /**
   * Creates a new Gaussian 2D peak result calculator.
   *
   * <p>The calibration need only contain the information required for the specified functions.
   *
   * @param psf the psf
   * @param calibration the calibration (used for converting the parameters)
   * @param flags the flags specifying the helper functions to support
   * @return the gaussian 2 D peak result helper
   * @throws ConfigurationException If not a Gaussian 2D PSF or the calibration is invalid
   * @throws ConversionException If unit conversion fails
   */
  public static Gaussian2DPeakResultCalculator create(PSF psf, Calibration calibration, int flags) {
    final CalibrationReader helper =
        (calibration != null) ? new CalibrationReader(calibration) : null;
    return create(psf, helper, flags);
  }

  /**
   * Creates a new Gaussian 2D peak result calculator.
   *
   * <p>The calibration need only contain the information required for the specified functions.
   *
   * @param psf the psf
   * @param calibrationReader the calibration reader (used for converting the parameters)
   * @param flags the flags specifying the helper functions to support
   * @return the gaussian 2 D peak result helper
   * @throws ConfigurationException If not a Gaussian 2D PSF or the calibration is invalid
   * @throws ConversionException If unit conversion fails
   */
  public static Gaussian2DPeakResultCalculator create(PSF psf, CalibrationReader calibrationReader,
      int flags) {
    final BaseGaussian2DPeakResultCalculator helper =
        new BaseGaussian2DPeakResultCalculator(psf, calibrationReader);

    // Try the desired methods
    if (BitFlagUtils.anySet(flags, AMPLITUDE)) {
      helper.getAmplitude(PARAMS);
    }
    if (BitFlagUtils.anySet(flags, LSE_PRECISION)) {
      helper.getLsePrecision(PARAMS, 0);
    }
    if (BitFlagUtils.anySet(flags, LSE_PRECISION_X)) {
      helper.getLsePrecision(PARAMS);
    }
    if (BitFlagUtils.anySet(flags, MLE_PRECISION)) {
      helper.getMlePrecision(PARAMS, 0);
    }
    if (BitFlagUtils.anySet(flags, MLE_PRECISION_X)) {
      helper.getMlePrecision(PARAMS);
    }
    if (BitFlagUtils.anySet(flags, PIXEL_AMPLITUDE)) {
      helper.getPixelAmplitude(PARAMS);
    }

    // Get a fast implementation
    return (helper.oneAxisSD) ? new OneAxisFastGaussian2DPeakResultCalculator(helper)
        : new TwoAxisFastGaussian2DPeakResultCalculator(helper);
  }

  /**
   * Calculate the localisation precision for least squares fitting a Gaussian2D PSF to a Gaussian2D
   * PSF. This is an approximation of the precision of fitting to an optical PSF. Uses the Mortensen
   * formula for an EMCCD camera (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 6.
   *
   * <p>This method will use the background noise to approximate the expected background value of
   * each pixel.
   *
   * @param a The size of the pixels in nm
   * @param s The peak standard deviation in nm
   * @param N The peak signal in photons
   * @param b The background noise standard deviation in photons
   * @param emCcd True if an emCcd camera
   * @return The location precision in nm in each dimension (X/Y)
   */
  public static double getPrecision(double a, double s, double N, double b, boolean emCcd) {
    // Get background expected value. This is what is actually used in the Mortensen method
    final double b2 = b * b;

    if (emCcd) {
      // If an emCcd camera was used then the input standard deviation will already be amplified
      // by the EM-gain sqrt(2) factor. To prevent double counting this factor we must divide by it.
      // Since this has been squared then divide by 2.
      return getPrecisionX(a, s, N, b2 / 2.0, 2);
    }
    return getPrecisionX(a, s, N, b2, 1);
  }

  /**
   * Calculate the localisation variance for least squares fitting a Gaussian2D PSF to a Gaussian2D
   * PSF. This is an approximation of the precision of fitting to an optical PSF. Uses the Mortensen
   * formula for an EMCCD camera (Mortensen, et al (2010) Nature Methods 7, 377-383), equation 6.
   *
   * <p>This method will use the background noise to approximate the expected background value of
   * each pixel.
   *
   * @param a The size of the pixels in nm
   * @param s The peak standard deviation in nm
   * @param N The peak signal in photons
   * @param b The background noise standard deviation in photons
   * @param emCcd True if an emCcd camera
   * @return The location variance in nm in each dimension (X/Y)
   */
  public static double getVariance(double a, double s, double N, double b, boolean emCcd) {
    // Get background expected value. This is what is actually used in the Mortensen method
    final double b2 = b * b;

    if (emCcd) {
      // If an emCcd camera was used then the input standard deviation will already be amplified
      // by the EM-gain sqrt(2) factor. To prevent double counting this factor we must divide by it.
      // Since this has been squared then divide by 2.
      return getVarianceX(a, s, N, b2 / 2.0, 2);
    }
    return getVarianceX(a, s, N, b2, 1);
  }

  /**
   * Calculate the localisation precision for maximum likelihood fitting a Gaussian2D PSF to a
   * Gaussian2D PSF. This is an approximation of the precision of fitting to an optical PSF. Uses
   * the Mortensen formula for an EMCCD camera (Mortensen, et al (2010) Nature Methods 7, 377-383),
   * SI equation 54.
   *
   * <p>This method will use the background noise to approximate the expected background value of
   * each pixel.
   *
   * @param a The size of the pixels in nm
   * @param s The peak standard deviation in nm
   * @param N The peak signal in photons
   * @param b The background noise standard deviation in photons
   * @param emCcd True if an emCcd camera
   * @return The location precision in nm in each dimension (X/Y)
   */
  public static double getMLPrecision(double a, double s, double N, double b, boolean emCcd) {
    // Get background expected value. This is what is actually used in the Mortensen method
    final double b2 = b * b;

    if (emCcd) {
      // If an emCcd camera was used then the input standard deviation will already be amplified
      // by the EM-gain sqrt(2) factor. To prevent double counting this factor we must divide by it.
      // Since this has been squared then divide by 2.
      return getMLPrecisionX(a, s, N, b2 / 2.0, true);
    }
    return getMLPrecisionX(a, s, N, b2, false);
  }

  /**
   * Calculate the localisation variance for maximum likelihood fitting a Gaussian2D PSF to a
   * Gaussian2D PSF. This is an approximation of the precision of fitting to an optical PSF. Uses
   * the Mortensen formula for an EMCCD camera (Mortensen, et al (2010) Nature Methods 7, 377-383),
   * SI equation 54.
   *
   * <p>This method will use the background noise to approximate the expected background value of
   * each pixel.
   *
   * @param a The size of the pixels in nm
   * @param s The peak standard deviation in nm
   * @param N The peak signal in photons
   * @param b The background noise standard deviation in photons
   * @param emCcd True if an emCcd camera
   * @return The location variance in nm in each dimension (X/Y)
   */
  public static double getMLVariance(double a, double s, double N, double b, boolean emCcd) {
    // Get background expected value. This is what is actually used in the Mortensen method
    final double b2 = b * b;

    if (emCcd) {
      // If an emCcd camera was used then the input standard deviation will already be amplified
      // by the EM-gain sqrt(2) factor. To prevent double counting this factor we must divide by it.
      // Since this has been squared then divide by 2.
      return getMLVarianceX(a, s, N, b2 / 2.0, true);
    }
    return getMLVarianceX(a, s, N, b2, false);
  }

  /**
   * Calculate the localisation precision for least squares fitting a Gaussian2D PSF to a Gaussian2D
   * PSF. This is an approximation of the precision of fitting to an optical PSF for least squares
   * estimation. Uses the Mortensen formula for an EMCCD camera (Mortensen, et al (2010) Nature
   * Methods 7, 377-383), equation 6.
   *
   * <p>If the expected photons per pixel is unknown then use the standard deviation across the
   * image and the method {@link #getPrecision(double, double, double, double, boolean)}.
   *
   * @param a The size of the pixels in nm
   * @param s The peak standard deviation in nm
   * @param N The peak signal in photons
   * @param b2 The expected number of photons per pixel from a background with spatially constant
   *        expectation value across the image (Note that this is b^2 not b, which could be the
   *        standard deviation of the image pixels)
   * @param emCcd True if an emCcd camera
   * @return The location precision in nm in each dimension (X/Y)
   */
  public static double getPrecisionX(final double a, final double s, final double N,
      final double b2, boolean emCcd) {
    // EM-CCD noise factor
    final double F = (emCcd) ? 2 : 1;
    return getPrecisionX(a, s, N, b2, F);
  }

  /**
   * Calculate the localisation precision for least squares fitting a Gaussian2D PSF to a Gaussian2D
   * PSF. This is an approximation of the precision of fitting to an optical PSF for least squares
   * estimation. Uses the Mortensen formula for an EMCCD camera (Mortensen, et al (2010) Nature
   * Methods 7, 377-383), equation 6.
   *
   * @param a The size of the pixels in nm
   * @param s The peak standard deviation in nm
   * @param N The peak signal in photons
   * @param b2 The expected number of photons per pixel from a background with spatially constant
   *        expectation value across the image (Note that this is b^2 not b, which could be the
   *        standard deviation of the image pixels)
   * @param F EM-CCD noise factor (usually 2 for an EM-CCD camera, else 1)
   * @return The location precision in nm in each dimension (X/Y)
   */
  public static double getPrecisionX(final double a, final double s, final double N,
      final double b2, final double F) {
    return Math.sqrt(getVarianceX(a, s, N, b2, F));
  }

  /**
   * Calculate the localisation variance for least squares fitting a Gaussian2D PSF to a Gaussian2D
   * PSF. This is an approximation of the precision of fitting to an optical PSF for least squares
   * estimation. Uses the Mortensen formula for an EMCCD camera (Mortensen, et al (2010) Nature
   * Methods 7, 377-383), equation 6.
   *
   * <p>If the expected photons per pixel is unknown then use the standard deviation across the
   * image and the method {@link #getPrecision(double, double, double, double, boolean)}.
   *
   * @param a The size of the pixels in nm
   * @param s The peak standard deviation in nm
   * @param N The peak signal in photons
   * @param b2 The expected number of photons per pixel from a background with spatially constant
   *        expectation value across the image (Note that this is b^2 not b, which could be the
   *        standard deviation of the image pixels)
   * @param emCcd True if an emCcd camera
   * @return The location variance in nm in each dimension (X/Y)
   */
  public static double getVarianceX(final double a, final double s, final double N, final double b2,
      boolean emCcd) {
    // EM-CCD noise factor
    final double F = (emCcd) ? 2 : 1;
    return getVarianceX(a, s, N, b2, F);
  }

  /**
   * Calculate the localisation variance for least squares fitting a Gaussian2D PSF to a Gaussian2D
   * PSF. This is an approximation of the precision of fitting to an optical PSF for least squares
   * estimation. Uses the Mortensen formula for an EMCCD camera (Mortensen, et al (2010) Nature
   * Methods 7, 377-383), equation 6.
   *
   * @param a The size of the pixels in nm
   * @param s The peak standard deviation in nm
   * @param N The peak signal in photons
   * @param b2 The expected number of photons per pixel from a background with spatially constant
   *        expectation value across the image (Note that this is b^2 not b, which could be the
   *        standard deviation of the image pixels)
   * @param F EM-CCD noise factor (usually 2 for an EM-CCD camera, else 1)
   * @return The location variance in nm in each dimension (X/Y)
   */
  public static double getVarianceX(final double a, final double s, final double N, final double b2,
      final double F) {
    if (N <= 0) {
      return Double.POSITIVE_INFINITY;
    }

    // Note that we input b^2 directly to this equation. This is the expected value of the pixel
    // background.
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
   * Calculate the localisation precision for maximum likelihood fitting a Gaussian2D PSF to a
   * Gaussian2D PSF. This is an approximation of the precision of fitting to an optical PSF. Uses
   * the Mortensen formula for an EMCCD camera (Mortensen, et al (2010) Nature Methods 7, 377-383),
   * SI equation 54.
   *
   * @param a The size of the pixels in nm
   * @param s The peak standard deviation in nm
   * @param N The peak signal in photons
   * @param b2 The expected number of photons per pixel from a background with spatially constant
   *        expectation value across the image (Note that this is b^2 not b, which could be the
   *        standard deviation of the image pixels)
   * @param emCcd True if an emCcd camera
   * @return The location precision in nm in each dimension (X/Y)
   */
  public static double getMLPrecisionX(double a, double s, double N, double b2, boolean emCcd) {
    return Math.sqrt(getMLVarianceX(a, s, N, b2, emCcd));
  }

  /**
   * Calculate the localisation precision for maximum likelihood fitting a Gaussian2D PSF to a
   * Gaussian2D PSF. This is an approximation of the precision of fitting to an optical PSF. Uses
   * the Mortensen formula for an EMCCD camera (Mortensen, et al (2010) Nature Methods 7, 377-383),
   * SI equation 54.
   *
   * @param a The size of the pixels in nm
   * @param s The peak standard deviation in nm
   * @param N The peak signal in photons
   * @param b2 The expected number of photons per pixel from a background with spatially constant
   *        expectation value across the image (Note that this is b^2 not b, which could be the
   *        standard deviation of the image pixels)
   * @param emCcd True if an emCcd camera
   * @param integrationpoints the number of integration points for the LegendreGaussIntegrator
   * @return The location precision in nm in each dimension (X/Y)
   */
  public static double getMLPrecisionX(double a, double s, double N, double b2, boolean emCcd,
      int integrationpoints) {
    return Math.sqrt(getMLVarianceX(a, s, N, b2, emCcd, integrationpoints));
  }

  /**
   * Calculate the localisation variance for maximum likelihood fitting a Gaussian2D PSF to a
   * Gaussian2D PSF. This is an approximation of the precision of fitting to an optical PSF. Uses
   * the Mortensen formula for an EMCCD camera (Mortensen, et al (2010) Nature Methods 7, 377-383),
   * SI equation 54.
   *
   * @param a The size of the pixels in nm
   * @param s The peak standard deviation in nm
   * @param N The peak signal in photons
   * @param b2 The expected number of photons per pixel from a background with spatially constant
   *        expectation value across the image (Note that this is b^2 not b, which could be the
   *        standard deviation of the image pixels)
   * @param emCcd True if an emCcd camera
   * @return The location variance in nm in each dimension (X/Y)
   */
  public static double getMLVarianceX(double a, double s, double N, double b2, boolean emCcd) {
    return getMLVarianceX(a, s, N, b2, emCcd, POINTS);
  }

  /**
   * Calculate the localisation variance for maximum likelihood fitting a Gaussian2D PSF to a
   * Gaussian2D PSF. This is an approximation of the precision of fitting to an optical PSF. Uses
   * the Mortensen formula for an EMCCD camera (Mortensen, et al (2010) Nature Methods 7, 377-383),
   * SI equation 54.
   *
   * <p>In the event of failure to integrate the formula the variance for Least Squares Estimation
   * is returned.
   *
   * @param a The size of the pixels in nm
   * @param s The peak standard deviation in nm
   * @param N The peak signal in photons
   * @param b2 The expected number of photons per pixel from a background with spatially constant
   *        expectation value across the image (Note that this is b^2 not b, which could be the
   *        standard deviation of the image pixels)
   * @param emCcd True if an emCcd camera
   * @param integrationpoints the number of integration points for the LegendreGaussIntegrator
   * @return The location variance in nm in each dimension (X/Y)
   */
  public static double getMLVarianceX(double a, double s, double N, double b2, boolean emCcd,
      int integrationpoints) {
    if (N <= 0) {
      return Double.POSITIVE_INFINITY;
    }

    // EM-CCD noise factor
    final double F = (emCcd) ? 2 : 1;
    final double a2 = a * a;
    // Adjustment for square pixels
    final double sa2 = s * s + a2 / 12.0;

    final double rho = 2 * Math.PI * sa2 * b2 / (N * a2);
    try {
      final double i1 = computeI1(rho, integrationpoints);
      if (i1 > 0) {
        return F * (sa2 / N) * (1 / i1);
      }
    } catch (final TooManyEvaluationsException ignored) {
      // Ignore
    }
    return getVarianceX(a, s, N, b2, emCcd);
  }

  /**
   * Compute the function I1 using numerical integration. See Mortensen, et al (2010) Nature Methods
   * 7, 377-383), SI equation 43.
   *
   * <pre>
   * I1 = 1 + sum [ ln(t) / (1 + t/rho) ] dt
   *    = - sum [ t * ln(t) / (t + rho) ] dt
   * </pre>
   *
   * <p>Where sum is the integral between 0 and 1. In the case of rho=0 the function returns 1;
   *
   * @param rho the rho
   * @param integrationpoints the number of integration points for the LegendreGaussIntegrator
   * @return the I1 value
   */
  private static double computeI1(final double rho, int integrationpoints) {
    if (rho == 0) {
      return 1;
    }

    final double relativeAccuracy = 1e-4;
    final double absoluteAccuracy = 1e-8;
    final int minimalIterationCount = 3;
    final int maximalIterationCount = 32;

    // Use an integrator that does not use the boundary since log(0) is undefined.
    final UnivariateIntegrator i = new IterativeLegendreGaussIntegrator(integrationpoints,
        relativeAccuracy, absoluteAccuracy, minimalIterationCount, maximalIterationCount);

    // Specify the function to integrate.
    // Note:
    // The alternative function Math.log(x) / ( 1 + x / rho) requires more evaluations and
    // sometimes does not converge, presumably because log(x) significantly changes as x -> 0
    // where as x log(x) in the following function is more stable.
    final UnivariateFunction f = x -> x * Math.log(x) / (x + rho);
    return -i.integrate(2000, f, 0, 1);
  }

  /**
   * Get the amplitude of a Gaussian 2D PSF. Amplitude = intensity / (2*pi*sx*sy).
   *
   * @param intensity the intensity
   * @param sx the sx
   * @param sy the sy
   * @return the amplitude
   */
  public static double getAmplitude(double intensity, double sx, double sy) {
    return (intensity / (2 * Math.PI * sx * sy));
  }

  /**
   * Gets the single Gaussian 2D standard deviation from independent x and y standard deviations. s
   * = sqrt(abs(sx*sy)).
   *
   * @param sx the sx
   * @param sy the sy
   * @return the single Gaussian 2D standard deviation
   */
  public static double getStandardDeviation(double sx, double sy) {
    return Math.sqrt(Math.abs(sx * sy));
  }

  /**
   * Gets the single Gaussian 2D standard deviation squared from independent x and y standard
   * deviations. s2 = abs(sx*sy). This avoids a sqrt in the computation of
   * {@link #getStandardDeviation(double, double)}.
   *
   * @param sx the sx
   * @param sy the sy
   * @return the single Gaussian 2D standard deviation squared
   */
  public static double getStandardDeviation2(double sx, double sy) {
    return Math.abs(sx * sy);
  }

  /**
   * Creates the params array for a Gaussian 2D peak result using the PSF type.
   *
   * @param psfType the psf type
   * @param gaussian2DParams the full gaussian 2D function parameters
   * @return the PeakResult params
   * @throws IllegalArgumentException if the psf is not a Gaussian 2D or the input parameters are
   *         not the correct length for a Gaussian 2D function
   */
  public static float[] createParams(PSFType psfType, float[] gaussian2DParams) {
    if (gaussian2DParams.length != PeakResult.STANDARD_PARAMETERS + 3) {
      throw new IllegalArgumentException(
          "Parameters must be a full Gaussian 2D parameters array of length "
              + (PeakResult.STANDARD_PARAMETERS + 3));
    }

    switch (psfType) {
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
   * Creates the params array for a Gaussian 2D peak result.
   *
   * <p>The psf parameters are used to determine if the PSF is one axis, two axis or two axis and
   * angle. If no PSF parameters are given then a 1 axis Gaussian is created with a standard
   * deviation of 1.
   *
   * @param background the background
   * @param intensity the intensity
   * @param x the x
   * @param y the y
   * @param z the z
   * @param psfParameters the psf parameters
   * @return the params
   * @throws IllegalArgumentException If the number of PSF parameters is invalid
   */
  public static float[] createParams(float background, float intensity, float x, float y, float z,
      float... psfParameters) {
    if (psfParameters == null) {
      return createOneAxisParams(background, intensity, x, y, z, 1);
    }
    switch (psfParameters.length) {
      case 1:
        return createOneAxisParams(background, intensity, x, y, z, psfParameters[0]);
      case 2:
        return createTwoAxisParams(background, intensity, x, y, z, psfParameters[0],
            psfParameters[1]);
      case 3:
        return createTwoAxisAndAngleParams(background, intensity, x, y, z, psfParameters[0],
            psfParameters[1], psfParameters[2]);
      default:
        throw new IllegalArgumentException("Invalid number of Gaussian 2D PSF parameters");
    }
  }

  /**
   * Creates the params array for a one axis Gaussian 2D peak result.
   *
   * @param background the background
   * @param intensity the intensity
   * @param x the x
   * @param y the y
   * @param z the z
   * @param s the s
   * @return the params
   */
  public static float[] createOneAxisParams(float background, float intensity, float x, float y,
      float z, float s) {
    final float[] params = new float[PeakResult.STANDARD_PARAMETERS + 1];
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
   * @param background the background
   * @param intensity the intensity
   * @param x the x
   * @param y the y
   * @param z the z
   * @param sx the sx
   * @param sy the sy
   * @return the params
   */
  public static float[] createTwoAxisParams(float background, float intensity, float x, float y,
      float z, float sx, float sy) {
    final float[] params = new float[PeakResult.STANDARD_PARAMETERS + 2];
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
   * @param background the background
   * @param intensity the intensity
   * @param x the x
   * @param y the y
   * @param z the z
   * @param sx the sx
   * @param sy the sy
   * @param a the a
   * @return the params
   */
  public static float[] createTwoAxisAndAngleParams(float background, float intensity, float x,
      float y, float z, float sx, float sy, float a) {
    final float[] params = new float[PeakResult.STANDARD_PARAMETERS + 3];
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
   * Compute the cumulative normal distribution within the range -x to x:
   *
   * <pre>
   * return erf(x / sqrt(2))
   * </pre>
   *
   * <p>This uses a fast approximation to the Error function.
   *
   * @param x the x
   * @return the cumulative normal distribution {@code CDF(X<x)}
   */
  public static double cumulative(double x) {
    return Erf.erf(x * ONE_OVER_ROOT2);
  }

  /**
   * Compute the cumulative 2D normal distribution as the probability that a sample lies inside the
   * ellipsoid determined by its Mahalanobis distance r from the Gaussian, a direct generalization
   * of the standard deviation.
   *
   * <pre>
   * return 1 - exp(-r * r / 2)
   * </pre>
   *
   * <p>This formula is provided in <a href=
   * "https://en.wikipedia.org/wiki/Mahalanobis_distance#Normal_distributions"> Wikipedia:
   * Mahalanobis distance - Normal distributions</a>
   *
   * @param r the distance r
   * @return the cumulative 2D normal distribution {@code F(r)}
   */
  public static double cumulative2D(double r) {
    // https://en.wikipedia.org/wiki/Multivariate_normal_distribution#Cumulative_distribution_function
    // https://upload.wikimedia.org/wikipedia/commons/a/a2/Cumulative_function_n_dimensional_Gaussians_12.2013.pdf
    // 1 - exp(x) == -(exp(x) - 1)
    return -Math.expm1(-r * r * 0.5);
  }

  /**
   * Compute the inverse cumulative 2D normal distribution as the Mahalanobis distance r from the
   * Gaussian given the probability that a sample lies inside the ellipsoid determined by its
   * distance r.
   *
   * <pre>
   * return sqrt(-2 ln (1-p) )
   * </pre>
   *
   * <p>This formula is provided in <a href=
   * "https://en.wikipedia.org/wiki/Mahalanobis_distance#Normal_distributions"> Wikipedia:
   * Mahalanobis distance - Normal distributions</a>
   *
   * @param p the cumulative 2D normal distribution {@code F(r)}
   * @return Mahalanobis distance r from the Gaussian
   * @throws IllegalArgumentException If p is not in the range 0-1
   */
  public static double inverseCumulative2D(double p) {
    if (p <= 0 || p > 1) {
      if (p == 0) {
        return 0; // Avoid returning -0 (the result of Math.sqrt(-0))
      }
      throw new IllegalArgumentException("P must be in the range 0 - 1");
    }
    return Math.sqrt(-2 * Math.log1p(-p));
  }

  /**
   * Gets the average signal value within a range r standard deviations of the centre.
   *
   * <p>The average signal value is taken using the expected sum of the Gaussian within the range
   * divided by the elliptical area of the same range. The expected sum is computed using
   * {@link #cumulative2D(double)} and the area would be pi * sx * sy * r * r.
   *
   * <p>Note: Argument sx and sy are not checked that they are positive.
   *
   * @param intensity the total intensity of the Gaussian
   * @param sx the Gaussian standard deviation in the X dimension
   * @param sy the Gaussian standard deviation in the Y dimension
   * @param r the range over which to compute the average signal (Mahalanobis distance r)
   * @return the mean
   */
  public static double getMeanSignalUsingR(double intensity, double sx, double sy, double r) {
    return intensity * cumulative2D(r) / (Math.PI * sx * sy * r * r);
  }

  /**
   * Gets the average signal value within 1 standard deviations of the centre.
   *
   * <p>The average signal value is taken using the expected sum of the Gaussian within the range
   * divided by the elliptical area of the same range. The expected sum is computed using
   * {@link #cumulative2D(double)} and the area would be pi * sx * sy.
   *
   * <p>Note: Argument sx and sy are not checked that they are positive.
   *
   * @param intensity the total intensity of the Gaussian
   * @param sx the Gaussian standard deviation in the X dimension
   * @param sy the Gaussian standard deviation in the Y dimension
   * @return the mean
   */
  public static double getMeanSignalUsingR1(double intensity, double sx, double sy) {
    return intensity * R1 / (sx * sy);
  }

  /**
   * Gets the average signal value within 2 standard deviations of the centre.
   *
   * <p>The average signal value is taken using the expected sum of the Gaussian within the range
   * divided by the elliptical area of the same range. The expected sum is computed using
   * {@link #cumulative2D(double)} and the area would be pi * sx * sy * 2 * 2.
   *
   * <p>As an alternative definition, the standard deviation of the background can be computed using
   * the standard deviation of the signal in the region around the Gaussian.
   *
   * <p>Note: Argument sx and sy are not checked that they are positive.
   *
   * @param intensity the total intensity of the Gaussian
   * @param sx the Gaussian standard deviation in the X dimension
   * @param sy the Gaussian standard deviation in the Y dimension
   * @return the mean
   */
  public static double getMeanSignalUsingR2(double intensity, double sx, double sy) {
    return intensity * R2 / (sx * sy);
  }

  /**
   * Gets the average signal value using the range r standard deviations of the centre defined by
   * the given the cumulative 2D normal distribution {@code F(r)}.
   *
   * <p>The average signal value is taken using the expected sum of the Gaussian within the range
   * divided by the elliptical area of the same range. The expected range is computed using
   * {@link #inverseCumulative2D(double)} and the area would be pi * sx * sy * r * r.
   *
   * <p>Note: Argument sx and sy are not checked that they are positive.
   *
   * @param intensity the total intensity of the Gaussian
   * @param sx the Gaussian standard deviation in the X dimension
   * @param sy the Gaussian standard deviation in the Y dimension
   * @param p the cumulative 2D normal distribution {@code F(r)}
   * @return the mean
   * @throws IllegalArgumentException If p is not in the range 0-1
   */
  public static double getMeanSignalUsingP(double intensity, double sx, double sy, double p) {
    final double r = inverseCumulative2D(p);
    if (r == 0) {
      return 0;
    }
    return intensity * p / (Math.PI * sx * sy * r * r);
  }

  /**
   * Gets the average signal value using the range r standard deviations of the centre that covers
   * half of the total 2D Gaussian (cumulative 2D normal distribution {@code F(r)=0.5}).
   *
   * <p>The average signal value is taken using the expected sum of the Gaussian within the range
   * divided by the elliptical area of the same range. The expected range is computed using
   * {@link #inverseCumulative2D(double)} and the area would be pi * sx * sy * r * r.
   *
   * <p>Note: When F(r)=0.5 then the inverseCumulative2D function computes the factor to convert a
   * standard deviation of a 1D Gaussian to a Half-Width at Half Maxima (HWHM):
   *
   * <pre>
   * HWHM = sqrt(2 * log(2))
   * (F(r)= 0.5
   *      = sqrt(-2 * log(0.5)
   *      = sqrt(-2 * log(1 / 2))
   *      = sqrt(-2 * (log(1) - log(2)))
   *      = sqrt(-2 * -log(2))
   * </pre>
   *
   * <p>Thus this computes the mean signal within the HWHM of a 2D Gaussian.
   *
   * <p>Note: Argument sx and sy are not checked that they are positive.
   *
   * @param intensity the total intensity of the Gaussian
   * @param sx the Gaussian standard deviation in the X dimension
   * @param sy the Gaussian standard deviation in the Y dimension
   * @return the mean
   */
  public static double getMeanSignalUsingP05(double intensity, double sx, double sy) {
    return intensity * P05 / (sx * sy);
  }

  /**
   * Add mean intensity if the PSF is a Gaussian 2D function and no mean intensity is present in the
   * results. The mean intensity is the mean signal within the elliptical area of the half-width
   * half-maxima of the Gaussian.
   *
   * @param psf the psf
   * @param results the results
   * @see MemoryPeakResults#hasMeanIntensity()
   * @see #getMeanSignalUsingP05(double, double, double)
   */
  public static void addMeanIntensity(PSFOrBuilder psf, MemoryPeakResults results) {
    // Some formats already read the mean intensity, e.g. SMLM, TSF
    if (PsfHelper.isGaussian2D(psf) && !results.hasMeanIntensity()) {
      final int[] indices = PsfHelper.getGaussian2DWxWyIndices(psf);
      final int isx = indices[0];
      final int isy = indices[1];
      results.forEach((PeakResultProcedure) peakResult -> {
        final float[] p = peakResult.getParameters();
        final float u = (float) Gaussian2DPeakResultHelper
            .getMeanSignalUsingP05(p[PeakResult.INTENSITY], p[isx], p[isy]);
        peakResult.setMeanIntensity(u);
      });
    }
  }
}
