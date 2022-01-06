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

package uk.ac.sussex.gdsc.smlm.function.gaussian;

import org.apache.commons.math3.util.FastMath;

/**
 * Evaluates a 2-dimensional elliptical Gaussian function for a configured number of peaks.
 *
 * <p>The single parameter x in the {@link #eval(int, double[])} function is assumed to be a linear
 * index into 2-dimensional data. The dimensions of the data must be specified to allow unpacking to
 * coordinates.
 *
 * <p>Data should be packed in descending dimension order, e.g. Y,X : Index for [x,y] = MaxX*y + x.
 */
public class EllipticalGaussian2DFunction extends MultiPeakGaussian2DFunction {
  /** The number of gradient parameters for each Gaussian. */
  protected static final int GRADIENT_PARAMETERS_PER_PEAK = 6;

  /** set to true if the Gaussian has no rotation angle. */
  protected boolean[] zeroAngle;
  /** The pre-computed function factors for each Gaussian. */
  protected final double[][] peakFactors;
  /** The Gaussian parameters (a). */
  protected double[] params;

  /**
   * Constructor.
   *
   * @param numberOfPeaks The number of peaks
   * @param maxx The maximum x value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   * @param maxy The maximum y value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   */
  public EllipticalGaussian2DFunction(int numberOfPeaks, int maxx, int maxy) {
    super(numberOfPeaks, maxx, maxy);
    zeroAngle = new boolean[numberOfPeaks];
    peakFactors = new double[numberOfPeaks][16];
  }

  @Override
  public Gaussian2DFunction copy() {
    return new EllipticalGaussian2DFunction(numberOfPeaks, maxx, maxy);
  }

  /** The index for the The amplitude./height normalisation: 1/(2*pi*sx*sy). */
  protected static final int N = 0;
  /** The index for the The amplitude./height. */
  protected static final int HEIGHT = 1;
  /** The index for the x0 position pre-factor. */
  protected static final int AA = 2;
  /** The index for the x0*x1 position pre-factor (for rotation). */
  protected static final int BB = 3;
  /** The index for the x1 position pre-factor. */
  protected static final int CC = 4;
  /** The index for the x0 position gradient pre-factor. */
  protected static final int AA2 = 5;
  /** The index for the x0*x1 position gradient pre-factor. */
  protected static final int BB2 = 6;
  /** The index for the x1 position gradient pre-factor. */
  protected static final int CC2 = 7;
  /** The index for the x width gradient pre-factor. */
  protected static final int NX = 8;
  /** The index for the x width gradient pre-factor. */
  protected static final int AX = 9;
  /** The index for the x width gradient pre-factor. */
  protected static final int BX = 10;
  /** The index for the x width gradient pre-factor. */
  protected static final int CX = 11;
  /** The index for the y width gradient pre-factor. */
  protected static final int NY = 12;
  /** The index for the y width gradient pre-factor. */
  protected static final int AY = 13;
  /** The index for the y width gradient pre-factor. */
  protected static final int BY = 14;
  /** The index for the y width gradient pre-factor. */
  protected static final int CY = 15;

  @Override
  public void initialise(double[] a) {
    this.params = a;
    // Precalculate multiplication factors
    for (int j = 0; j < numberOfPeaks; j++) {
      final double theta = a[j * PARAMETERS_PER_PEAK + ANGLE];
      final double sx = a[j * PARAMETERS_PER_PEAK + X_SD];
      final double sy = a[j * PARAMETERS_PER_PEAK + Y_SD];
      final double sx2 = sx * sx;
      final double sy2 = sy * sy;
      final double sx3 = sx2 * sx;
      final double sy3 = sy2 * sy;

      peakFactors[j][N] = ONE_OVER_TWO_PI / (sx * sy);
      peakFactors[j][HEIGHT] = a[j * PARAMETERS_PER_PEAK + SIGNAL] * peakFactors[j][N];

      // All prefactors are negated since the Gaussian uses the exponential to the negative:
      // (A/2*pi*sx*sy) * exp( -( a(x-x0)^2 + 2b(x-x0)(y-y0) + c(y-y0)^2 ) )

      if (theta == 0) {
        // cosSqt = 1
        // sinSqt = 0
        // sincost = 0
        // sin2t = 0
        // cos2t = 1

        zeroAngle[j] = true;

        peakFactors[j][AA] = -0.5 / sx2;
        peakFactors[j][BB] = 0;
        peakFactors[j][CC] = -0.5 / sy2;

        // For the angle gradient
        peakFactors[j][AA2] = 0;
        peakFactors[j][BB2] = -0.5 * (-1.0 / sx2 + 1.0 / sy2);
        peakFactors[j][CC2] = 0;

        // For the x-width gradient
        peakFactors[j][NX] = -1.0 / sx;
        peakFactors[j][AX] = 1.0 / sx3;
        peakFactors[j][BX] = 0;
        peakFactors[j][CX] = 0;

        // For the y-width gradient
        peakFactors[j][NY] = -1.0 / sy;
        peakFactors[j][AY] = 0;
        peakFactors[j][BY] = 0;
        peakFactors[j][CY] = 1.0 / sy3;
      } else {
        zeroAngle[j] = false;

        final double cosSqt = Math.cos(theta) * Math.cos(theta);
        final double sinSqt = Math.sin(theta) * Math.sin(theta);
        final double sincost = Math.sin(theta) * Math.cos(theta);
        final double sin2t = Math.sin(2 * theta);
        final double cos2t = Math.cos(2 * theta);

        peakFactors[j][AA] = -0.5 * (cosSqt / sx2 + sinSqt / sy2);
        peakFactors[j][BB] = -0.25 * (-sin2t / sx2 + sin2t / sy2);
        peakFactors[j][CC] = -0.5 * (sinSqt / sx2 + cosSqt / sy2);

        // For the angle gradient
        peakFactors[j][AA2] = -(-sincost / sx2 + sincost / sy2);
        peakFactors[j][BB2] = -0.5 * (-cos2t / sx2 + cos2t / sy2);
        peakFactors[j][CC2] = -(sincost / sx2 - sincost / sy2);

        // For the x-width gradient
        peakFactors[j][NX] = -1.0 / sx;
        peakFactors[j][AX] = cosSqt / sx3;
        peakFactors[j][BX] = -0.5 * sin2t / sx3;
        peakFactors[j][CX] = sinSqt / sx3;

        // For the y-width gradient
        peakFactors[j][NY] = -1.0 / sy;
        peakFactors[j][AY] = sinSqt / sy3;
        peakFactors[j][BY] = 0.5 * sin2t / sy3;
        peakFactors[j][CY] = cosSqt / sy3;
      }
    }
  }

  /**
   * Evaluates a 2-dimensional elliptical Gaussian function for multiple peaks.
   *
   * <p>{@inheritDoc}
   */
  @Override
  public double eval(final int x, final double[] dyda) {
    // Track the position of the parameters
    int apos = 0;
    int dydapos = 0;

    // First parameter is the background level
    double y = params[BACKGROUND];
    dyda[dydapos++] = 1.0; // Gradient for a constant background is 1

    // Unpack the predictor into the dimensions
    final int x1 = x / maxx;
    final int x0 = x % maxx;

    for (int j = 0; j < numberOfPeaks; j++) {
      y += gaussian(x0, x1, dyda, apos, dydapos, zeroAngle[j], peakFactors[j]);
      apos += PARAMETERS_PER_PEAK;
      dydapos += GRADIENT_PARAMETERS_PER_PEAK;
    }

    return y;
  }

  /**
   * Evaluates a 2-dimensional elliptical Gaussian function for multiple peaks.
   *
   * <p>{@inheritDoc}
   */
  @Override
  public double eval(final int x) {
    // Track the position of the parameters
    int apos = 0;

    // First parameter is the background level
    double y = params[BACKGROUND];

    // Unpack the predictor into the dimensions
    final int x1 = x / maxx;
    final int x0 = x % maxx;

    for (int j = 0; j < numberOfPeaks; j++, apos += PARAMETERS_PER_PEAK) {
      y += gaussian(x0, x1, apos, zeroAngle[j], peakFactors[j]);
    }

    return y;
  }

  /**
   * Compute the Gaussian at a set offset from the centre.
   *
   * @param x0 the x0 offset
   * @param x1 the x1 offset
   * @param dyDa the first-order gradient
   * @param apos the parameter position for the current peak
   * @param dydapos the gradient position for the current peak
   * @param zeroAngle set to true if the Gaussian has no rotation angle
   * @param factors the factors
   * @return the Gaussian value
   */
  protected double gaussian(final int x0, final int x1, final double[] dyDa, final int apos,
      final int dydapos, final boolean zeroAngle, final double[] factors) {
    final double dx = x0 - params[apos + X_POSITION];
    final double dy = x1 - params[apos + Y_POSITION];
    final double dx2 = dx * dx;
    final double dxy = dx * dy;
    final double dy2 = dy * dy;

    final double aa = factors[AA];
    final double cc = factors[CC];

    // Calculate gradients
    if (zeroAngle) {
      final double exp = FastMath.exp(aa * dx2 + cc * dy2);
      dyDa[dydapos] = factors[N] * exp;
      final double y = factors[HEIGHT] * exp;

      dyDa[dydapos + 1] = y * (-2.0 * aa * dx);
      dyDa[dydapos + 2] = y * (-2.0 * cc * dy);

      dyDa[dydapos + 3] = y * (factors[NX] + factors[AX] * dx2);
      dyDa[dydapos + 4] = y * (factors[NY] + factors[CY] * dy2);

      dyDa[dydapos + 5] = y * (factors[BB2] * dxy);

      return y;
    }
    final double bb = factors[BB];

    final double exp = FastMath.exp(aa * dx2 + bb * dxy + cc * dy2);
    dyDa[dydapos] = factors[N] * exp;
    final double y = factors[HEIGHT] * exp;

    dyDa[dydapos + 1] = y * (-2.0 * aa * dx - bb * dy);
    dyDa[dydapos + 2] = y * (-2.0 * cc * dy - bb * dx);

    dyDa[dydapos + 3] =
        y * (factors[NX] + factors[AX] * dx2 + factors[BX] * dxy + factors[CX] * dy2);
    dyDa[dydapos + 4] =
        y * (factors[NY] + factors[AY] * dx2 + factors[BY] * dxy + factors[CY] * dy2);

    dyDa[dydapos + 5] = y * (factors[AA2] * dx2 + factors[BB2] * dxy + factors[CC2] * dy2);

    return y;
  }

  /**
   * Compute the Gaussian at a set offset from the centre.
   *
   * @param x0 the x0 offset
   * @param x1 the x1 offset
   * @param apos the parameter position for the current peak
   * @param zeroAngle set to true if the Gaussian has no rotation angle
   * @param factors the factors
   * @return the Gaussian value
   */
  protected double gaussian(final int x0, final int x1, final int apos, final boolean zeroAngle,
      final double[] factors) {
    final double dx = x0 - params[apos + X_POSITION];
    final double dy = x1 - params[apos + Y_POSITION];

    if (zeroAngle) {
      return factors[HEIGHT] * FastMath.exp(factors[AA] * dx * dx + factors[CC] * dy * dy);
    }

    return factors[HEIGHT]
        * FastMath.exp(factors[AA] * dx * dx + factors[BB] * dx * dy + factors[CC] * dy * dy);
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
    return true;
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
    return GRADIENT_PARAMETERS_PER_PEAK;
  }
}
