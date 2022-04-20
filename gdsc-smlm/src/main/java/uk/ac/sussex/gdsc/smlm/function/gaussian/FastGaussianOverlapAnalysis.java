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

/**
 * Compute the overlap between 2D Gaussian functions.
 *
 * <p>Given an input 2D Gaussian a region is created that covers a range of the function. The
 * overlap of other functions within this region can be computed.
 */
public class FastGaussianOverlapAnalysis {
  /**
   * A constant holding the maximum value an {@code int} can have, 2<sup>31</sup>-1.
   */
  private static final long MAX_VALUE = Integer.MAX_VALUE;

  private final int flags;
  private final AstigmatismZModel zModel;

  private final int maxx;
  private final int maxy;
  private final double centrex;
  private final double centrey;

  private double overlap;

  /**
   * Create a new overlap analysis object.
   *
   * @param flags The flags describing the Gaussian2DFunction function (see GaussianFunctionFactory)
   * @param zModel the z model (used to create the widths from the z position)
   * @param params The parameters for the Gaussian (assumes a single peak)
   * @param maxx The x-range over which to compute the function (assumed to be strictly positive)
   * @param maxy The y-range over which to compute the function (assumed to be strictly positive)
   * @throws IllegalArgumentException if {@code maxx * maxy} overflows an int
   */
  public FastGaussianOverlapAnalysis(int flags, AstigmatismZModel zModel, double[] params, int maxx,
      int maxy) {
    this.flags = flags;
    this.zModel = zModel;

    this.maxx = Math.max(1, maxx);
    this.maxy = Math.max(1, maxy);
    final long r = (long) maxx * (long) maxy;
    if ((int) r != r) {
      throw new IllegalArgumentException("Input range is too large: maxx * maxy = " + r);
    }
    // We will sample the Gaussian at integer intervals, i.e. on a pixel grid.
    // Pixels centres should be at 0.5,0.5. So if we want to draw a Gauss
    // centred in the middle of a pixel we need to adjust each centre.
    // Then we subtract the position of the target gaussian.
    centrex = maxx * 0.5 - 0.5 - params[Gaussian2DFunction.X_POSITION];
    centrey = maxy * 0.5 - 0.5 - params[Gaussian2DFunction.Y_POSITION];
  }

  /**
   * Gets the range over which to evaluate a Gaussian using a factor of the standard deviation.
   *
   * <p>The range is clipped to 1 to Integer.MAX_VALUE.
   *
   * @param sd the standard deviation
   * @param range the range factor
   * @return the range
   */
  public static int getRange(double sd, double range) {
    final long l = (long) Math.ceil(2 * sd * range);
    if (l < 1) {
      return 1;
    }
    if (l >= MAX_VALUE) {
      return Integer.MAX_VALUE;
    }
    return (int) l + 1;
  }

  /**
   * Gets the range over which to evaluate a Gaussian using a factor of the standard deviation.
   *
   * <p>The range is clipped to 1 to max.
   *
   * @param sd the standard deviation
   * @param range the range factor
   * @param max the max value to return
   * @return the range
   */
  public static int getRange(double sd, double range, int max) {
    final long l = (long) Math.ceil(2 * sd * range);
    if (l < 1) {
      return 1;
    }
    if (l >= max) {
      return max;
    }
    return (int) l + 1;
  }

  /**
   * Add the Gaussian function data to the overlap region. This is the region that contains the
   * input function within the range defined in the constructor.
   *
   * <p>The square region is masked using the expected sum of the function within the range. The
   * overlap of other functions within this masked region can be computed, or within the square
   * region.
   *
   * @param params The parameters for the Gaussian (can be multiple peaks)
   */
  public void add(double[] params) {
    // Add the function to the overlap
    final int numberOfPeaks = params.length / Gaussian2DFunction.PARAMETERS_PER_PEAK;
    final Gaussian2DFunction f =
        GaussianFunctionFactory.create2D(numberOfPeaks, maxx, maxy, flags, zModel);
    params = params.clone();
    for (int n = 0; n < numberOfPeaks; n++) {
      params[n * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_POSITION] += centrex;
      params[n * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.Y_POSITION] += centrey;
    }
    overlap += f.integral(params);
  }

  /**
   * Get the overlap. This is the sum of all the functions added to the analysis using
   * {@link #add(double[])} within the region defined in the constructor.
   *
   * @return The overlap
   */
  public double getOverlap() {
    return overlap;
  }

  /**
   * Gets the size of the overlap region. Typically this is the product of the x and y range used in
   * the constructor.
   *
   * @return the size
   */
  public int getSize() {
    return maxx * maxy;
  }
}
