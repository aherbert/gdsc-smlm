/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.function.gaussian;

import org.apache.commons.math3.util.FastMath;

/**
 * Evaluates a 2-dimensional Gaussian function for a single peak.
 *
 * <p>The single parameter x in the {@link #eval(int, double[])} function is assumed to be a linear
 * index into 2-dimensional data. The dimensions of the data must be specified to allow unpacking to
 * coordinates.
 *
 * <p>Data should be packed in descending dimension order, e.g. Y,X : Index for [x,y] = MaxX*y + x.
 */
public class SingleFreeCircularGaussian2DFunction extends Gaussian2DFunction {
  private static final int[] gradientIndices;

  static {
    gradientIndices = createGradientIndices(1, new SingleFreeCircularGaussian2DFunction(1, 1));
  }

  /** The background. */
  protected double background;
  /** The x0 position. */
  protected double x0pos;
  /** The x1 position. */
  protected double x1pos;

  /** Set to true if the Gaussian has no rotation angle. */
  protected boolean zeroAngle;
  /** The amplitude./height normalisation: 1/(2*pi*sx*sy). */
  protected double norm;
  /** The amplitude./height. */
  protected double height;
  /** x0 position pre-factor. */
  protected double aa;
  /** x0*x1 position pre-factor (for rotation). */
  protected double bb;
  /** x1 position pre-factor. */
  protected double cc;
  /** x width gradient pre-factor. */
  protected double nx;
  /** x width gradient pre-factor. */
  protected double ax;
  /** x width gradient pre-factor. */
  protected double bx;
  /** x width gradient pre-factor. */
  protected double cx;
  /** y width gradient pre-factor. */
  protected double ny;
  /** y width gradient pre-factor. */
  protected double ay;
  /** y width gradient pre-factor. */
  protected double by;
  /** y width gradient pre-factor. */
  protected double cy;

  /**
   * Constructor.
   *
   * @param maxx The maximum x value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   * @param maxy The maximum y value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   */
  public SingleFreeCircularGaussian2DFunction(int maxx, int maxy) {
    super(maxx, maxy);
  }

  @Override
  public Gaussian2DFunction copy() {
    return new SingleFreeCircularGaussian2DFunction(maxx, maxy);
  }

  @Override
  public void initialise(double[] a) {
    background = a[BACKGROUND];
    x0pos = a[X_POSITION];
    x1pos = a[Y_POSITION];

    // Precalculate multiplication factors
    final double theta = a[ANGLE];
    final double sx = a[X_SD];
    final double sy = a[Y_SD];
    final double sx2 = sx * sx;
    final double sy2 = sy * sy;
    final double sx3 = sx2 * sx;
    final double sy3 = sy2 * sy;

    norm = ONE_OVER_TWO_PI / (sx * sy);
    height = a[SIGNAL] * norm;

    // All prefactors are negated since the Gaussian uses the exponential to the negative:
    // (A/2*pi*sx*sy) * exp( -( a(x-x0)^2 + 2b(x-x0)(y-y0) + c(y-y0)^2 ) )

    if (theta == 0) {
      zeroAngle = true;

      // cosSqt = 1
      // sinSqt = 0
      // sin2t = 0

      aa = -0.5 / sx2;
      bb = 0;
      cc = -0.5 / sy2;

      // For the x-width gradient
      nx = -1.0 / sx;
      ax = 1.0 / sx3;
      bx = 0;
      cx = 0;

      // For the y-width gradient
      ny = -1.0 / sy;
      ay = 0;
      by = 0;
      cy = 1.0 / sy3;
    } else {
      zeroAngle = false;

      final double cosSqt = Math.cos(theta) * Math.cos(theta);
      final double sinSqt = Math.sin(theta) * Math.sin(theta);
      final double sin2t = Math.sin(2 * theta);

      aa = -0.5 * (cosSqt / sx2 + sinSqt / sy2);
      bb = -0.25 * (-sin2t / sx2 + sin2t / sy2);
      cc = -0.5 * (sinSqt / sx2 + cosSqt / sy2);

      // For the x-width gradient
      nx = -1.0 / sx;
      ax = cosSqt / sx3;
      bx = -0.5 * sin2t / sx3;
      cx = sinSqt / sx3;

      // For the y-width gradient
      ny = -1.0 / sy;
      ay = sinSqt / sy3;
      by = 0.5 * sin2t / sy3;
      cy = cosSqt / sy3;
    }
  }

  /**
   * Evaluates a 2-dimensional elliptical Gaussian function for a single peak.
   *
   * <p>{@inheritDoc}
   */
  @Override
  public double eval(final int x, final double[] dyda) {
    // First parameter is the background level
    dyda[0] = 1.0; // Gradient for a constant background is 1

    // Unpack the predictor into the dimensions
    final int x1 = x / maxx;
    final int x0 = x % maxx;

    return background + gaussian(x0, x1, dyda);
  }

  /**
   * Evaluates a 2-dimensional elliptical Gaussian function for a single peak.
   *
   * <p>{@inheritDoc}
   */
  @Override
  public double eval(final int x) {
    // Unpack the predictor into the dimensions
    final int x1 = x / maxx;
    final int x0 = x % maxx;

    final double dx = x0 - x0pos;
    final double dy = x1 - x1pos;

    if (zeroAngle) {
      return background + height * FastMath.exp(aa * dx * dx + cc * dy * dy);
    }
    return background + height * FastMath.exp(aa * dx * dx + bb * dx * dy + cc * dy * dy);
  }

  private double gaussian(final int x0, final int x1, final double[] dyDa) {
    final double dx = x0 - x0pos;
    final double dy = x1 - x1pos;
    final double dx2 = dx * dx;
    final double dxy = dx * dy;
    final double dy2 = dy * dy;

    // Calculate gradients
    if (zeroAngle) {
      final double exp = FastMath.exp(aa * dx2 + cc * dy2);
      dyDa[1] = norm * exp;
      final double y = height * exp;

      dyDa[2] = y * (-2.0 * aa * dx);
      dyDa[3] = y * (-2.0 * cc * dy);

      dyDa[4] = y * (nx + ax * dx2);
      dyDa[5] = y * (ny + cy * dy2);
      return y;
    }

    final double exp = FastMath.exp(aa * dx2 + bb * dxy + cc * dy2);
    dyDa[1] = norm * exp;
    final double y = height * exp;

    dyDa[2] = y * (-2.0 * aa * dx - bb * dy);
    dyDa[3] = y * (-2.0 * cc * dy - bb * dx);

    dyDa[4] = y * (nx + ax * dx2 + bx * dxy + cx * dy2);
    dyDa[5] = y * (ny + ay * dx2 + by * dxy + cy * dy2);
    return y;
  }

  @Override
  public int getNPeaks() {
    return 1;
  }

  @Override
  public boolean evaluatesBackground() {
    return true;
  }

  @Override
  public boolean evaluatesSignal() {
    return true;
  }

  @Override
  public boolean evaluatesAngle() {
    return false;
  }

  @Override
  public boolean evaluatesPosition() {
    return true;
  }

  @Override
  public boolean evaluatesSD0() {
    return true;
  }

  @Override
  public boolean evaluatesSD1() {
    return true;
  }

  @Override
  public int getGradientParametersPerPeak() {
    return 5;
  }

  @Override
  public int[] gradientIndices() {
    return gradientIndices;
  }
}
