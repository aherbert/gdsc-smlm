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

package uk.ac.sussex.gdsc.smlm.function.gaussian;

import org.apache.commons.math3.util.FastMath;

/**
 * Evaluates an 2-dimensional Gaussian function for a configured number of peaks.
 *
 * <p>The single parameter x in the {@link #eval(int, double[])} function is assumed to be a linear
 * index into 2-dimensional data. The dimensions of the data must be specified to allow unpacking to
 * coordinates.
 *
 * <p>Data should be packed in descending dimension order, e.g. Y,X : Index for [x,y] = MaxX*y + x.
 */
public class FixedGaussian2DFunction extends MultiPeakGaussian2DFunction {
  /** The number of gradient parameters for each Gaussian. */
  protected static final int GRADIENT_PARAMETERS_PER_PEAK = 3;

  /** The pre-computed function factors for each Gaussian. */
  protected final double[][] peakFactors;
  /** The Gaussian parameters (a). */
  protected double[] a;

  /**
   * Constructor.
   *
   * @param npeaks The number of peaks
   * @param maxx The maximum x value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   * @param maxy The maximum y value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   */
  public FixedGaussian2DFunction(int npeaks, int maxx, int maxy) {
    super(npeaks, maxx, maxy);
    peakFactors = new double[npeaks][4];
  }

  /** {@inheritDoc} */
  @Override
  public Gaussian2DFunction copy() {
    return new FixedGaussian2DFunction(npeaks, maxx, maxy);
  }

  /** The index for the The amplitude./height normalisation: 1/(2*pi*sx*sy). */
  protected static final int N = 0;
  /** The index for the The amplitude./height. */
  protected static final int HEIGHT = 1;
  /** The index for the x0 position pre-factor. */
  protected static final int AA = 2;
  /** The index for the x0 position gradient pre-factor. */
  protected static final int AA2 = 3;

  /** {@inheritDoc} */
  @Override
  public void initialise(double[] a) {
    this.a = a;
    // Precalculate multiplication factors
    for (int j = 0; j < npeaks; j++) {
      final double sx = a[j * PARAMETERS_PER_PEAK + X_SD];
      final double sx2 = sx * sx;

      peakFactors[j][N] = ONE_OVER_TWO_PI / sx2;
      peakFactors[j][HEIGHT] = a[j * PARAMETERS_PER_PEAK + SIGNAL] * peakFactors[j][N];

      // All prefactors are negated since the Gaussian uses the exponential to the negative:
      // (A/2*pi*sx*sy) * exp( -( a(x-x0)^2 + 2b(x-x0)(y-y0) + c(y-y0)^2 ) )

      peakFactors[j][AA] = -0.5 / sx2;
      peakFactors[j][AA2] = -2.0 * peakFactors[j][AA];
    }
  }

  /**
   * Evaluates an 2-dimensional fixed circular Gaussian function for multiple peaks.
   *
   * <p>{@inheritDoc}
   */
  @Override
  public double eval(final int x, final double[] dyda) {
    // Track the position of the parameters
    int apos = 0;
    int dydapos = 0;

    // First parameter is the background level
    double y = a[BACKGROUND];
    dyda[dydapos++] = 1.0; // Gradient for a constant background is 1

    // Unpack the predictor into the dimensions
    final int x1 = x / maxx;
    final int x0 = x % maxx;

    for (int j = 0; j < npeaks; j++) {
      y += gaussian(x0, x1, dyda, apos, dydapos, peakFactors[j]);
      apos += PARAMETERS_PER_PEAK;
      dydapos += GRADIENT_PARAMETERS_PER_PEAK;
    }

    return y;
  }

  /**
   * Compute the Gaussian at a set offset from the centre.
   *
   * @param x0 the x0 offset
   * @param x1 the x1 offset
   * @param dy_da the first-order gradient
   * @param apos the parameter position for the current peak
   * @param dydapos the gradient position for the current peak
   * @param factors the factors
   * @return the Gaussian value
   */
  protected double gaussian(final int x0, final int x1, final double[] dy_da, final int apos,
      final int dydapos, final double[] factors) {
    final double dx = x0 - a[apos + X_POSITION];
    final double dy = x1 - a[apos + Y_POSITION];

    // Calculate gradients

    final double aadx2dy2 = factors[AA] * (dx * dx + dy * dy);
    final double exp = FastMath.exp(aadx2dy2);
    dy_da[dydapos] = factors[N] * exp;
    final double y = factors[HEIGHT] * exp;
    final double yaa2 = y * factors[AA2];
    dy_da[dydapos + 1] = yaa2 * dx;
    dy_da[dydapos + 2] = yaa2 * dy;

    return y;
  }

  /**
   * Evaluates an 2-dimensional fixed circular Gaussian function for multiple peaks.
   *
   * <p>{@inheritDoc}
   */
  @Override
  public double eval(final int x) {
    // Track the position of the parameters
    int apos = 0;

    // First parameter is the background level
    double y = a[BACKGROUND];

    // Unpack the predictor into the dimensions
    final int x1 = x / maxx;
    final int x0 = x % maxx;

    for (int j = 0; j < npeaks; j++, apos += PARAMETERS_PER_PEAK) {
      y += gaussian(x0, x1, apos, peakFactors[j]);
    }

    return y;
  }

  /**
   * Compute the Gaussian at a set offset from the centre.
   *
   * @param x0 the x0 offset
   * @param x1 the x1 offset
   * @param apos the parameter position for the current peak
   * @param factors the factors
   * @return the Gaussian value
   */
  protected double gaussian(final int x0, final int x1, final int apos, final double[] factors) {
    final double dx = x0 - a[apos + X_POSITION];
    final double dy = x1 - a[apos + Y_POSITION];

    return factors[HEIGHT] * FastMath.exp(factors[AA] * (dx * dx + dy * dy));
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
    return false;
  }

  @Override
  public boolean evaluatesSD1() {
    return false;
  }

  @Override
  public int getGradientParametersPerPeak() {
    return GRADIENT_PARAMETERS_PER_PEAK;
  }
}
