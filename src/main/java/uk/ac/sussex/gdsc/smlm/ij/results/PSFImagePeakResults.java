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
package uk.ac.sussex.gdsc.smlm.ij.results;

import java.awt.Rectangle;

import org.apache.commons.math3.util.FastMath;

import uk.ac.sussex.gdsc.core.data.NotImplementedException;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.smlm.data.config.ConfigurationException;
import uk.ac.sussex.gdsc.smlm.data.config.PSFHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitConverterFactory;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Draws the fit results using the Gaussian PSF to an ImageJ image.
 */
public class PSFImagePeakResults extends IJImagePeakResults {
  private boolean fixedWidth = false;
  private float psfWidth = 0f;
  private boolean calculatedPrecision = false;

  private Gaussian2DPeakResultCalculator calculator;
  private TypeConverter<DistanceUnit> dc;
  private int isx, isy, ia;

  private boolean requirePSFParameters;

  // Multiplication factors and variables for plotting the fixed Gaussian
  private double[] fixedParams = null;

  /**
   * Instantiates a new PSF image peak results.
   *
   * @param title Title of the image (appended with a suffix)
   * @param bounds Define the bounding rectangle of the image coordinates. Any results outside this
   *        will not be displayed.
   * @param scale the scale
   */
  public PSFImagePeakResults(String title, Rectangle bounds, float scale) {
    super(title, bounds, scale);
  }

  /** {@inheritDoc} */
  @Override
  protected void preBegin() {
    // this.displayFlags should be OK so don't call super.preBegin()

    requirePSFParameters = true;

    int flags = 0;
    if ((displayFlags & DISPLAY_SIGNAL) != 0) {
      flags |= Gaussian2DPeakResultHelper.AMPLITUDE;
    }

    // Note: we do not catch configuration exceptions if required Gaussian 2D configuration is
    // missing

    if (fixedWidth) {
      // Check if we need the amplitude for the fixed width PSF
      requirePSFParameters = flags != 0;
    } else if (calculatedPrecision) {
      flags |= Gaussian2DPeakResultHelper.LSE_PRECISION;

      // To convert the precision to pixels
      if (!hasCalibration()) {
        throw new ConfigurationException("nm/pixel is required when drawing using the precision");
      }

      dc = UnitConverterFactory.createConverter(DistanceUnit.NM, DistanceUnit.PIXEL,
          getCalibrationReader().getNmPerPixel());
    } else {
      // We need to know the parameters for the Gaussian 2D PSF
      final int[] indices = PSFHelper.getGaussian2DWxWyIndices(getPSF());
      isx = indices[0];
      isy = indices[1];
      try {
        ia = PSFHelper.getGaussian2DAngleIndex(getPSF());
      } catch (final ConfigurationException e) {
        // No rotation angle
        ia = 0;
      }
    }

    if (flags != 0) {
      calculator = Gaussian2DPeakResultHelper.create(getPSF(), getCalibrationReader(), flags);
    }
  }

  /** {@inheritDoc} */
  @Override
  public boolean isUncalibrated() {
    // Not supported
    return false;
  }

  /** {@inheritDoc} */
  @Override
  public void setUncalibrated(boolean uncalibrated) {
    throw new NotImplementedException(
        "This method is not supported. The PSF assumes the units are in pixels.");
  }

  /** {@inheritDoc} */
  @Override
  public void add(int peak, float x, float y, float v) {
    if (requirePSFParameters) {
      throw new IllegalStateException("The PSF image requires missing Gaussian 2D parameters");
    }
    add(new PeakResult(peak, x, y, v));
  }

  /** {@inheritDoc} */
  @Override
  public void add(float x, float y, float v) {
    if (requirePSFParameters) {
      throw new IllegalStateException("The PSF image requires missing Gaussian 2D parameters");
    }
    add(new PeakResult(x, y, v));
  }

  /** {@inheritDoc} */
  @Override
  public void add(int[] allpeak, float[] allx, float[] ally, float[] allv) {
    if (requirePSFParameters) {
      throw new IllegalStateException("The PSF image requires missing Gaussian 2D parameters");
    }
    for (int i = 0; i < allx.length; i++) {
      add(new PeakResult(allpeak[i], allx[i], ally[i], allv[i]));
    }
  }

  /** {@inheritDoc} */
  @Override
  public void add(float[] allx, float[] ally, float[] allv) {
    if (requirePSFParameters) {
      throw new IllegalStateException("The PSF image requires missing Gaussian 2D parameters");
    }
    for (int i = 0; i < allx.length; i++) {
      add(new PeakResult(allx[i], ally[i], allv[i]));
    }
  }

  /** {@inheritDoc} */
  @Override
  public void add(int peak, int origX, int origY, float origValue, double error, float noise,
      float meanIntensity, float[] params, float[] paramsDev) {
    if (!imageActive) {
      return;
    }
    addPeak(peak, noise, params);
    updateImage();
  }

  private void addPeak(int peak, float noise, float[] params) {
    float x = mapX(params[PeakResult.X]);
    float y = mapY(params[PeakResult.Y]);

    // Check bounds
    if (x < 0 || x >= imageWidth || y < 0 || y >= imageHeight) {
      return;
    }

    checkAndUpdateToFrame(peak);

    // Initialise for a free Gaussian function:
    // f(x,y) = A exp(-(a(x-x0)(x-x0) + 2b(x-x0)(y-y0) + c(y-y0)(y-y0)))
    // See: http://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function

    final float amplitude =
        ((displayFlags & DISPLAY_SIGNAL) != 0) ? calculator.getAmplitude(params) : 1;

    final double[] psfParams;
    if (fixedWidth) {
      psfParams = fixedParams;
    } else {
      // Precalculate multiplication factors
      final double t, sx, sy;
      if (calculatedPrecision) {
        t = 0.0;
        final double precision = calculator.getLSEPrecision(params, noise);
        sx = sy = dc.convert(precision);
      } else {
        sx = params[isx];
        sy = params[isy];
        t = (ia != 0) ? params[ia] : 0;
      }
      psfParams = getPSFParameters(t, sx, sy);
    }

    final double a = psfParams[0];
    final double b = psfParams[1];
    final double c = psfParams[2];
    final double width = psfParams[3];
    final double height = psfParams[4];

    // Use 0.5 offset to centre the value in the middle of each pixel
    x -= 0.5 / scale;
    y -= 0.5 / scale;

    int xmin = (int) Math.floor(x - width * scale);
    int xmax = (int) Math.ceil(x + width * scale);
    int ymin = (int) Math.floor(y - height * scale);
    int ymax = (int) Math.ceil(y + height * scale);

    // Clip range
    xmin = FastMath.max(xmin, 0);
    xmax = (int) FastMath.min(xmax, xlimit);
    ymin = FastMath.max(ymin, 0);
    ymax = (int) FastMath.min(ymax, ylimit);

    // Compute Gaussian PSF
    final int[] index = new int[(xmax - xmin + 1) * (ymax - ymin + 1)];
    final float[] value = new float[index.length];
    int i = 0;
    for (int y0 = ymin; y0 <= ymax; y0++) {
      for (int x0 = xmin; x0 <= xmax; x0++) {
        final int ii = y0 * imageWidth + x0;
        index[i] = ii;
        final float dx = (x0 - x) / scale;
        final float dy = (y0 - y) / scale;
        value[i] = (float) (amplitude * FastMath.exp(a * dx * dx + b * dx * dy + c * dy * dy));
        i++;
      }
    }

    // Now add the values to the configured indices
    synchronized (data) {
      size++;
      while (i-- > 0) {
        data[index[i]] += value[i];
      }
    }
  }

  /**
   * Gets the PSF parameters.
   *
   * @param t the rotation angle
   * @param sx the standard deviation in the x dimension
   * @param sy the standard deviation in the y dimension
   * @return the parameters
   */
  private static double[] getPSFParameters(double t, double sx, double sy) {
    double a, b, c;
    double height, width;

    if (t == 0) {
      // sin(0) == 0
      // cos(0) == 1
      a = (-1.0 / (2.0 * sx * sx));
      b = 0.0;
      c = (-1.0 / (2.0 * sy * sy));

      // Calculate the range for the PSF as 3 sigma.
      width = 3.0 * sx;
      height = 3.0 * sy;
    } else {
      a = -(Math.cos(t) * Math.cos(t) / (2.0 * sx * sx)
          + Math.sin(t) * Math.sin(t) / (2.0 * sy * sy));
      b = -(-Math.sin(2.0 * t) / (2.0 * sx * sx) + Math.sin(2.0 * t) / (2.0 * sy * sy));
      c = -(Math.sin(t) * Math.sin(t) / (2.0 * sx * sx)
          + Math.cos(t) * Math.cos(t) / (2.0 * sy * sy));

      // Note that the Gaussian2DFitter returns the angle of the major axis (sx) relative to the
      // x-axis.
      // The angle is in the range -pi/2 to pi/2

      // The width and height for the range to be plotted can be derived from the general parametric
      // form of the ellipse.
      // See: http://en.wikipedia.org/wiki/Ellipse#General_parametric_form

      // Ensure the angle is the correct range (0 to pi)
      if (t < 0) {
        t += Math.PI;
      }

      final double phi = t; // Angle between x-axis and major axis of ellipse
      final double t1 = -t; // Angle around the ellipse from 0 to 2pi for the x-axis
      final double t2 = t1 + Math.PI / 2; // Angle around the ellipse from 0 to 2pi for the y-axis

      // Calculate the size of the ellipse at 3 sigma
      width =
          Math.abs(3.0 * (sx * Math.cos(t1) * Math.cos(phi) - sy * Math.sin(t1) * Math.sin(phi)));
      height =
          Math.abs(3.0 * (sx * Math.cos(t2) * Math.sin(phi) + sy * Math.sin(t2) * Math.cos(phi)));
    }

    final double[] params = new double[5];
    params[0] = a;
    params[1] = b;
    params[2] = c;
    params[3] = width;
    params[4] = height;
    return params;
  }

  /** {@inheritDoc} */
  @Override
  public void addAll(PeakResult[] results) {
    if (!imageActive) {
      return;
    }

    // TODO - Make this more efficient. It could use worker threads to increase speed.
    int i = 0;
    for (final PeakResult result : results) {
      addPeak(result.getFrame(), result.getNoise(), result.getParameters());
      if (++i % 64 == 0) {
        updateImage();
        if (!imageActive) {
          return;
        }
      }
    }
    updateImage();
  }

  /**
   * @return the width for a fixed-width Gaussian.
   */
  public float getWidth() {
    return psfWidth;
  }

  /**
   * Set the width of a fixed-width PSF. This is before scale adjustment so is provided in terms of
   * the original fitted data.
   *
   * @param width the width to set for a fixed-width Gaussian
   */
  public void setWidth(float width) {
    this.psfWidth = width;

    if (width > 0) {
      fixedWidth = true;
      fixedParams = getPSFParameters(0.0, width, width);
    } else {
      fixedWidth = false;
      fixedParams = null;
    }
  }

  /**
   * Checks if is calculated precision.
   *
   * @return if true plot the width of the PSF using the calculated precision
   */
  public boolean isCalculatedPrecision() {
    return calculatedPrecision;
  }

  /**
   * Set to true to plot the width of the PSF using the calculated precision.
   *
   * @param calculatedPrecision the new calculated precision
   */
  public void setCalculatedPrecision(boolean calculatedPrecision) {
    this.calculatedPrecision = calculatedPrecision;
  }
}
