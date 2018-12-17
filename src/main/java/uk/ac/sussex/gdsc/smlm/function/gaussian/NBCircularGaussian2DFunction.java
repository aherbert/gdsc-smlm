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

/**
 * Evaluates an 2-dimensional Gaussian function for a configured number of peaks.
 *
 * <p>The single parameter x in the {@link #eval(int, double[])} function is assumed to be a linear
 * index into 2-dimensional data. The dimensions of the data must be specified to allow unpacking to
 * coordinates.
 *
 * <p>Data should be packed in descending dimension order, e.g. Y,X : Index for [x,y] = MaxX*y + x.
 */
public class NBCircularGaussian2DFunction extends CircularGaussian2DFunction {
  /**
   * Constructor.
   *
   * @param numberOfPeaks The number of peaks
   * @param maxx The maximum x value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   * @param maxy The maximum y value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   */
  public NBCircularGaussian2DFunction(int numberOfPeaks, int maxx, int maxy) {
    super(numberOfPeaks, maxx, maxy);
  }

  /** {@inheritDoc} */
  @Override
  public Gaussian2DFunction copy() {
    return new NBCircularGaussian2DFunction(numberOfPeaks, maxx, maxy);
  }

  /**
   * Evaluates an 2-dimensional circular Gaussian function for multiple peaks.
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

  @Override
  public boolean evaluatesBackground() {
    return false;
  }
}
