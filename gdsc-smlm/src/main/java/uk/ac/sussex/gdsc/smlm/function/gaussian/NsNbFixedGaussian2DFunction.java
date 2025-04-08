/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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

import uk.ac.sussex.gdsc.smlm.utils.StdMath;

/**
 * Evaluates a 2-dimensional Gaussian function for a configured number of peaks.
 *
 * <p>The single parameter x in the {@link #eval(int, double[])} function is assumed to be a linear
 * index into 2-dimensional data. The dimensions of the data must be specified to allow unpacking to
 * coordinates.
 *
 * <p>Data should be packed in descending dimension order, e.g. Y,X : Index for [x,y] = MaxX*y + x.
 */
public class NsNbFixedGaussian2DFunction extends MultiPeakGaussian2DFunction {
  /** The number of gradient parameters for each Gaussian. */
  protected static final int GRADIENT_PARAMETERS_PER_PEAK = 2;
  /** The index for the The amplitude./height normalisation: 1/(2*pi*sx*sy). */
  protected static final int N = 0;
  /** The index for the The amplitude./height. */
  protected static final int HEIGHT = 1;
  /** The index for the x0 position pre-factor. */
  protected static final int AA = 2;
  /** The index for the x0 position gradient pre-factor. */
  protected static final int AA2 = 3;

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
  public NsNbFixedGaussian2DFunction(int numberOfPeaks, int maxx, int maxy) {
    super(numberOfPeaks, maxx, maxy);
    peakFactors = new double[numberOfPeaks][4];
  }

  @Override
  public Gaussian2DFunction copy() {
    return new NsNbFixedGaussian2DFunction(numberOfPeaks, maxx, maxy);
  }

  @Override
  public void initialise(double[] a) {
    this.params = a;
    // Precalculate multiplication factors
    for (int j = 0; j < numberOfPeaks; j++) {
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
   * Evaluates a 2-dimensional fixed circular Gaussian function for multiple peaks.
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

    // Unpack the predictor into the dimensions
    final int x1 = x / maxx;
    final int x0 = x % maxx;

    for (int j = 0; j < numberOfPeaks; j++) {
      y += gaussian(x0, x1, dyda, apos, dydapos, peakFactors[j]);
      apos += PARAMETERS_PER_PEAK;
      dydapos += GRADIENT_PARAMETERS_PER_PEAK;
    }

    return y;
  }

  /**
   * Evaluates a 2-dimensional fixed circular Gaussian function for multiple peaks.
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
      y += gaussian(x0, x1, apos, peakFactors[j]);
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
   * @param factors the factors
   * @return the Gaussian value
   */
  protected double gaussian(final int x0, final int x1, final double[] dyDa, final int apos,
      final int dydapos, final double[] factors) {
    final double dx = x0 - params[apos + X_POSITION];
    final double dy = x1 - params[apos + Y_POSITION];

    // Calculate gradients

    final double aadx2dy2 = factors[AA] * (dx * dx + dy * dy);
    final double y = factors[HEIGHT] * StdMath.exp(aadx2dy2);
    final double yaa2 = y * factors[AA2];
    dyDa[dydapos] = yaa2 * dx;
    dyDa[dydapos + 1] = yaa2 * dy;

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
    final double dx = x0 - params[apos + X_POSITION];
    final double dy = x1 - params[apos + Y_POSITION];

    return factors[HEIGHT] * StdMath.exp(factors[AA] * (dx * dx + dy * dy));
  }

  @Override
  public boolean evaluatesBackground() {
    return false;
  }

  @Override
  public boolean evaluatesSignal() {
    return false;
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
