/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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

import uk.ac.sussex.gdsc.core.data.NotImplementedException;
import uk.ac.sussex.gdsc.smlm.utils.StdMath;

/**
 * Evaluates a 2-dimensional Gaussian function for a single peak.
 *
 * <p>The single parameter x in the {@link #eval(int, double[])} function is assumed to be a linear
 * index into 2-dimensional data. The dimensions of the data must be specified to allow unpacking to
 * coordinates.
 *
 * <p>Data should be packed in descending dimension order, e.g. Y,X : Index for [x,y] = MaxX*y + x.
 *
 * <p>This function uses the signal, position and x-sd parameters. It ignores the background
 * parameter, y-sd and angle parameter. It can be used for fast evaluation of the sum of a Gaussian
 * within a region defined by maxx. maxy is defined by the number of calls to the eval() functions.
 * It does not support gradient evaluation.
 */
public class SingleSimpleGaussian2DFunction extends Gaussian2DFunction {
  private static final int[] gradientIndices = new int[0];

  /** The x0 position. */
  protected double x0pos;
  /** The x1 position. */
  protected double x1pos;

  /** The amplitude./height. */
  protected double height;
  /** position pre-factor. */
  protected double aa;

  /**
   * Constructor.
   *
   * @param maxx The maximum x value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   * @param maxy The maximum y value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   */
  public SingleSimpleGaussian2DFunction(int maxx, int maxy) {
    super(maxx, maxy);
  }

  @Override
  public Gaussian2DFunction copy() {
    return new SingleSimpleGaussian2DFunction(maxx, maxy);
  }

  @Override
  public void initialise(double[] a) {
    x0pos = a[X_POSITION];
    x1pos = a[Y_POSITION];

    final double sx = a[X_SD];
    final double sx2 = sx * sx;

    final double n = ONE_OVER_TWO_PI / sx2;
    height = a[SIGNAL] * n;

    // All prefactors are negated since the Gaussian uses the exponential to the negative:
    // A * exp( -( a(x-x0)^2 + 2b(x-x0)(y-y0) + c(y-y0)^2 ) )

    aa = -0.5 / sx2;
  }

  /**
   * Not implemented.
   *
   * <p>{@inheritDoc}
   */
  @Override
  public double eval(final int x, final double[] dyda) {
    throw new NotImplementedException();
  }

  /**
   * Evaluates a 2-dimensional circular Gaussian function for a single peak.
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

    return height * StdMath.exp(aa * (dx * dx + dy * dy));
  }

  @Override
  public int getNPeaks() {
    return 1;
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
    return false;
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
    return 0;
  }

  @Override
  public int[] gradientIndices() {
    return gradientIndices;
  }
}
