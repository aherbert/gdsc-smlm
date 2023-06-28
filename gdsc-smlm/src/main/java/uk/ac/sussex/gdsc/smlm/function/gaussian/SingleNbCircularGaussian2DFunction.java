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

import uk.ac.sussex.gdsc.smlm.utils.StdMath;

/**
 * Evaluates a 2-dimensional Gaussian function for a single peak.
 *
 * <p>The single parameter x in the {@link #eval(int, double[])} function is assumed to be a linear
 * index into 2-dimensional data. The dimensions of the data must be specified to allow unpacking to
 * coordinates.
 *
 * <p>Data should be packed in descending dimension order, e.g. Y,X : Index for [x,y] = MaxX*y + x.
 */
public class SingleNbCircularGaussian2DFunction extends SingleCircularGaussian2DFunction {
  private static final int[] gradientIndices;

  static {
    gradientIndices = createGradientIndices(1, new SingleNbCircularGaussian2DFunction(1, 1));
  }

  /**
   * Constructor.
   *
   * @param maxx The maximum x value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   * @param maxy The maximum y value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   */
  public SingleNbCircularGaussian2DFunction(int maxx, int maxy) {
    super(maxx, maxy);
  }

  @Override
  public Gaussian2DFunction copy() {
    return new SingleNbCircularGaussian2DFunction(maxx, maxy);
  }

  /**
   * Evaluates a 2-dimensional circular Gaussian function for a single peak.
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

  private double gaussian(final int x0, final int x1, final double[] dyDa) {
    final double dx = x0 - x0pos;
    final double dy = x1 - x1pos;

    // Calculate gradients

    final double aadx2dy2 = aa * (dx * dx + dy * dy);
    final double exp = StdMath.exp(aadx2dy2);
    dyDa[0] = norm * exp;
    final double y = height * exp;
    final double yaa2 = y * aa2;
    dyDa[1] = yaa2 * dx;
    dyDa[2] = yaa2 * dy;

    dyDa[3] = ax * y * (1 + aadx2dy2);

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
