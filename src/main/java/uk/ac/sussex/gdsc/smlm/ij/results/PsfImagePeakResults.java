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

package uk.ac.sussex.gdsc.smlm.ij.results;

import java.awt.Rectangle;
import org.apache.commons.math3.util.FastMath;
import uk.ac.sussex.gdsc.core.data.NotImplementedException;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.smlm.data.config.ConfigurationException;
import uk.ac.sussex.gdsc.smlm.data.config.PsfHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitConverterUtils;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Draws the fit results using the Gaussian PSF to an ImageJ image.
 */
public class PsfImagePeakResults extends ImageJImagePeakResults {
  private boolean fixedWidth;
  private float psfWidth;
  private boolean calculatedPrecision;

  private Gaussian2DPeakResultCalculator calculator;
  private TypeConverter<DistanceUnit> dc;
  private int isx;
  private int isy;
  private int ia;

  private boolean requirePsfParameters;

  // Multiplication factors and variables for plotting the fixed Gaussian
  private double[] fixedParams;

  /**
   * Instantiates a new PSF image peak results.
   *
   * @param title Title of the image (appended with a suffix)
   * @param bounds Define the bounding rectangle of the image coordinates. Any results outside this
   *        will not be displayed.
   * @param scale the scale
   */
  public PsfImagePeakResults(String title, Rectangle bounds, float scale) {
    super(title, bounds, scale);
  }

  @Override
  protected void preBegin() {
    // this.displayFlags should be OK so don't call super.preBegin()

    requirePsfParameters = true;

    int flags = 0;
    if ((displayFlags & DISPLAY_SIGNAL) != 0) {
      flags |= Gaussian2DPeakResultHelper.AMPLITUDE;
    }

    // Note: we do not catch configuration exceptions if required Gaussian 2D configuration is
    // missing

    if (fixedWidth) {
      // Check if we need the amplitude for the fixed width PSF
      requirePsfParameters = flags != 0;
    } else if (calculatedPrecision) {
      flags |= Gaussian2DPeakResultHelper.LSE_PRECISION;

      // To convert the precision to pixels
      if (!hasCalibration()) {
        throw new ConfigurationException("nm/pixel is required when drawing using the precision");
      }

      dc = UnitConverterUtils.createConverter(DistanceUnit.NM, DistanceUnit.PIXEL,
          getCalibrationReader().getNmPerPixel());
    } else {
      // We need to know the parameters for the Gaussian 2D PSF
      final int[] indices = PsfHelper.getGaussian2DWxWyIndices(getPsf());
      isx = indices[0];
      isy = indices[1];
      try {
        ia = PsfHelper.getGaussian2DAngleIndex(getPsf());
      } catch (final ConfigurationException ex) {
        // No rotation angle
        ia = 0;
      }
    }

    if (flags != 0) {
      calculator = Gaussian2DPeakResultHelper.create(getPsf(), getCalibrationReader(), flags);
    }
  }

  @Override
  public boolean isUncalibrated() {
    // Not supported
    return false;
  }

  @Override
  public void setUncalibrated(boolean uncalibrated) {
    throw new NotImplementedException(
        "This method is not supported. The PSF assumes the units are in pixels.");
  }

  @Override
  public void add(int peak, float x, float y, float value) {
    if (requirePsfParameters) {
      throw new IllegalStateException("The PSF image requires missing Gaussian 2D parameters");
    }
    add(new PeakResult(peak, x, y, value));
  }

  @Override
  public void add(float x, float y, float value) {
    if (requirePsfParameters) {
      throw new IllegalStateException("The PSF image requires missing Gaussian 2D parameters");
    }
    add(new PeakResult(x, y, value));
  }

  @Override
  public void add(int[] allpeak, float[] allx, float[] ally, float[] allv) {
    if (requirePsfParameters) {
      throw new IllegalStateException("The PSF image requires missing Gaussian 2D parameters");
    }
    for (int i = 0; i < allx.length; i++) {
      add(new PeakResult(allpeak[i], allx[i], ally[i], allv[i]));
    }
  }

  @Override
  public void add(float[] allx, float[] ally, float[] allv) {
    if (requirePsfParameters) {
      throw new IllegalStateException("The PSF image requires missing Gaussian 2D parameters");
    }
    for (int i = 0; i < allx.length; i++) {
      add(new PeakResult(allx[i], ally[i], allv[i]));
    }
  }

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
      final double t;
      final double sx;
      final double sy;
      if (calculatedPrecision) {
        t = 0.0;
        final double precision = calculator.getLsePrecision(params, noise);
        sx = sy = dc.convert(precision);
      } else {
        sx = params[isx];
        sy = params[isy];
        t = (ia != 0) ? params[ia] : 0;
      }
      psfParams = getPsfParameters(t, sx, sy);
    }

    final double a = psfParams[0];
    final double b = psfParams[1];
    final double c = psfParams[2];
    final double width = psfParams[3];
    final double height = psfParams[4];

    // Use 0.5 offset to centre the value in the middle of each pixel
    x -= 0.5 / scale;
    y -= 0.5 / scale;

    final int xmin = Math.max(0, (int) Math.floor(x - width * scale));
    int xmax = (int) Math.ceil(x + width * scale);
    final int ymin = Math.max(0, (int) Math.floor(y - height * scale));
    int ymax = (int) Math.ceil(y + height * scale);

    // Clip range
    xmax = (int) Math.min(xmax, xlimit);
    ymax = (int) Math.min(ymax, ylimit);

    // Compute Gaussian PSF
    final int[] index = new int[(xmax - xmin + 1) * (ymax - ymin + 1)];
    final float[] value = new float[index.length];
    int i1 = 0;
    for (int y0 = ymin; y0 <= ymax; y0++) {
      for (int x0 = xmin; x0 <= xmax; x0++) {
        final int i2 = y0 * imageWidth + x0;
        index[i1] = i2;
        final float dx = (x0 - x) / scale;
        final float dy = (y0 - y) / scale;
        value[i1] = (float) (amplitude * FastMath.exp(a * dx * dx + b * dx * dy + c * dy * dy));
        i1++;
      }
    }

    // Now add the values to the configured indices
    synchronized (data) {
      size++;
      while (i1-- > 0) {
        data[index[i1]] += value[i1];
      }
    }
  }

  /**
   * Gets the PSF parameters.
   *
   * @param angle the rotation angle
   * @param sx the standard deviation in the x dimension
   * @param sy the standard deviation in the y dimension
   * @return the parameters
   */
  private static double[] getPsfParameters(double angle, double sx, double sy) {
    // abc are defined factors for xx, xy and yy:
    // https://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function
    // CHECKSTYLE.OFF: LocalVariableName
    double a;
    double b;
    double c;
    // CHECKSTYLE.ON: LocalVariableName
    double height;
    double width;

    if (angle == 0) {
      // sin(0) == 0
      // cos(0) == 1
      a = (-1.0 / (2.0 * sx * sx));
      b = 0.0;
      c = (-1.0 / (2.0 * sy * sy));

      // Calculate the range for the PSF as 3 sigma.
      width = 3.0 * sx;
      height = 3.0 * sy;
    } else {
      a = -(Math.cos(angle) * Math.cos(angle) / (2.0 * sx * sx)
          + Math.sin(angle) * Math.sin(angle) / (2.0 * sy * sy));
      b = -(-Math.sin(2.0 * angle) / (2.0 * sx * sx) + Math.sin(2.0 * angle) / (2.0 * sy * sy));
      c = -(Math.sin(angle) * Math.sin(angle) / (2.0 * sx * sx)
          + Math.cos(angle) * Math.cos(angle) / (2.0 * sy * sy));

      // Note that the Gaussian2DFitter returns the angle of the major axis (sx) relative to the
      // x-axis. The angle is in the range -pi/2 to pi/2

      // The width and height for the range to be plotted can be derived from the general parametric
      // form of the ellipse.
      // See: http://en.wikipedia.org/wiki/Ellipse#General_parametric_form

      // Ensure the angle is the correct range (0 to pi)
      if (angle < 0) {
        angle += Math.PI;
      }

      final double phi = angle; // Angle between x-axis and major axis of ellipse
      final double t1 = -angle; // Angle around the ellipse from 0 to 2pi for the x-axis
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

  @Override
  public void addAll(PeakResult[] results) {
    if (!imageActive) {
      return;
    }

    // TODO - Make this more efficient. It could use worker threads to increase speed.
    int counter = 0;
    for (final PeakResult result : results) {
      addPeak(result.getFrame(), result.getNoise(), result.getParameters());
      if (++counter % 64 == 0) {
        updateImage();
        if (!imageActive) {
          return;
        }
      }
    }
    updateImage();
  }

  /**
   * Gets the width for a fixed-width Gaussian..
   *
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
      fixedParams = getPsfParameters(0.0, width, width);
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
