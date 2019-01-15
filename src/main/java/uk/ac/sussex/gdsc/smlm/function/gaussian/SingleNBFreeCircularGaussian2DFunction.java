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
 * Evaluates an 2-dimensional Gaussian function for a single peak.
 *
 * <p>The single parameter x in the {@link #eval(int, double[])} function is assumed to be a linear
 * index into 2-dimensional data. The dimensions of the data must be specified to allow unpacking to
 * coordinates.
 *
 * <p>Data should be packed in descending dimension order, e.g. Y,X : Index for [x,y] = MaxX*y + x.
 */
public class SingleNBFreeCircularGaussian2DFunction extends SingleFreeCircularGaussian2DFunction {
  private static final int[] gradientIndices;

  static {
    gradientIndices = createGradientIndices(1, new SingleNBFreeCircularGaussian2DFunction(1, 1));
  }

  /**
   * Constructor.
   *
   * @param maxx The maximum x value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   * @param maxy The maximum y value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   */
  public SingleNBFreeCircularGaussian2DFunction(int maxx, int maxy) {
    super(maxx, maxy);
  }

  @Override
  public Gaussian2DFunction copy() {
    return new SingleNBFreeCircularGaussian2DFunction(maxx, maxy);
  }

  /**
   * Evaluates an 2-dimensional elliptical Gaussian function for a single peak.
   *
   * <p>{@inheritDoc}
   */
  @Override
  public double eval(final int x, final double[] dyda) {
    // Unpack the predictor into the dimensions
    final int x1 = x / maxx;
    final int x0 = x % maxx;

    return background + gaussian(x0, x1, dyda);
  }

  private double gaussian(final int x0, final int x1, final double[] dy_da) {
    final double dx = x0 - x0pos;
    final double dy = x1 - x1pos;
    final double dx2 = dx * dx;
    final double dxy = dx * dy;
    final double dy2 = dy * dy;

    // Calculate gradients

    if (zeroAngle) {
      final double exp = FastMath.exp(aa * dx2 + cc * dy2);
      dy_da[0] = n * exp;
      final double y = height * exp;

      dy_da[1] = y * (-2.0 * aa * dx);
      dy_da[2] = y * (-2.0 * cc * dy);

      dy_da[3] = y * (nx + ax * dx2);
      dy_da[4] = y * (ny + cy * dy2);
      return y;
    }

    final double exp = FastMath.exp(aa * dx2 + bb * dxy + cc * dy2);
    dy_da[0] = n * exp;
    final double y = height * exp;

    dy_da[1] = y * (-2.0 * aa * dx - bb * dy);
    dy_da[2] = y * (-2.0 * cc * dy - bb * dx);

    dy_da[3] = y * (nx + ax * dx2 + bx * dxy + cx * dy2);
    dy_da[4] = y * (ny + ay * dx2 + by * dxy + cy * dy2);

    return y;
  }

  @Override
  public boolean evaluatesBackground() {
    return false;
  }

  @Override
  public int[] gradientIndices() {
    return gradientIndices;
  }
}
