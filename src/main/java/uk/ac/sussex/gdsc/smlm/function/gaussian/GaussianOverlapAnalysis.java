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

import uk.ac.sussex.gdsc.core.utils.SortUtils;
import uk.ac.sussex.gdsc.smlm.function.Erf;

import org.apache.commons.math3.util.FastMath;

/**
 * Compute the overlap between 2D Gaussian functions.
 *
 * <p>Given an input 2D Gaussian a region is created that covers a range of the function. The square
 * region is masked using a fraction of the expected sum of the function within the range. The
 * overlap of other functions within this region can be computed.
 */
public class GaussianOverlapAnalysis {
  /**
   * A constant holding the maximum value an {@code int} can have, 2<sup>31</sup>-1.
   */
  private static final long MAX_VALUE = Integer.MAX_VALUE;
  private static final double SQRT2 = FastMath.sqrt(2.0);

  private final int flags;
  private final AstigmatismZModel zModel;
  private final double[] params0;

  private final int maxx;
  private final int maxy;
  private final int size;
  private final double centrex;
  private final double centrey;

  private double[] data;
  private double[] overlap;
  private boolean[] mask;

  private double fraction = 0.95;

  /**
   * Create a new overlap analysis object.
   *
   * @param flags The flags describing the Gaussian2DFunction function (see GaussianFunctionFactory)
   * @param zModel the z model (used to create the widths from the z position)
   * @param params The parameters for the Gaussian (assumes a single peak)
   * @param maxx The x-range over which to compute the function (assumed to be strictly positive)
   * @param maxy The y-range over which to compute the function (assumed to be strictly positive)
   */
  public GaussianOverlapAnalysis(int flags, AstigmatismZModel zModel, double[] params, int maxx,
      int maxy) {
    this.flags = flags;
    this.zModel = zModel;
    this.params0 = params.clone();

    this.maxx = Math.max(1, maxx);
    this.maxy = Math.max(1, maxy);
    size = maxx * maxy;
    if (size < 0) {
      throw new IllegalArgumentException(
          "Input range is too large: maxx * maxy = " + ((long) maxx) * maxy);
    }
    // We will sample the Gaussian at integer intervals, i.e. on a pixel grid.
    // Pixels centres should be at 0.5,0.5. So if we want to draw a Gauss
    // centred in the middle of a pixel we need to adjust each centre
    centrex = maxx * 0.5 - 0.5;
    centrey = maxy * 0.5 - 0.5;
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
    if (l < 1L) {
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
    if (l < 1L) {
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
   * @param params The parameters for the Gaussian (can be multiple peaks)
   */
  public void add(double[] params) {
    add(params, false);
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
   * @param withinMask Set to true to only compute the overlap within the mask. This effects the
   *        computation of the weighted background (see {@link #getWeightedbackground()}.
   */
  public void add(double[] params, boolean withinMask) {
    // Note: When computing the overlap no adjustment is made for sampling on a pixel grid.
    // This is OK since the method will be used to evaluate the overlap between Gaussians that have
    // been fit using the functions.

    if (data == null) {
      // Initialise the input function
      data = new double[size];
      overlap = new double[size];
      final Gaussian2DFunction f =
          GaussianFunctionFactory.create2D(1, maxx, maxy, this.flags, zModel);

      // Note that the position should be in the centre of the sample region
      final double cx = params0[Gaussian2DFunction.X_POSITION];
      final double cy = params0[Gaussian2DFunction.Y_POSITION];
      params0[Gaussian2DFunction.X_POSITION] = centrex;
      params0[Gaussian2DFunction.Y_POSITION] = centrey;

      f.initialise(params0);
      for (int k = 0; k < size; k++) {
        data[k] = f.eval(k);
      }
      // Reset
      params0[Gaussian2DFunction.X_POSITION] = cx;
      params0[Gaussian2DFunction.Y_POSITION] = cy;

      if (fraction < 1) {
        // Create a mask with a fraction of the function value
        double sum = 0;
        final int[] indices = new int[size];
        for (int k = 0; k < size; k++) {
          sum += data[k];
          indices[k] = k;
        }
        SortUtils.sortIndices(indices, data, true);
        final double expected = sum * fraction;
        double last = 0;
        boolean useMask = false;
        mask = new boolean[data.length];
        for (int i = 0; i < data.length; i++) {
          final double v = data[indices[i]];
          // Note: We track the value since the Gaussian is symmetric and we want to include
          // all pixels with the same value
          final double newSum = sum + v;
          if (newSum >= expected && last != v) {
            // This is a new value that takes us over the expected signal
            useMask = true;
            break;
          }
          sum = newSum;
          mask[indices[i]] = true;
          last = v;
        }
        if (!useMask) {
          mask = null;
        }
      }
    }

    // Add the function to the overlap
    final int numberOfPeaks = params.length / Gaussian2DFunction.PARAMETERS_PER_PEAK;
    final Gaussian2DFunction f =
        GaussianFunctionFactory.create2D(numberOfPeaks, maxx, maxy, flags, zModel);
    params = params.clone();
    for (int n = 0; n < numberOfPeaks; n++) {
      params[n * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_POSITION] +=
          centrex - params0[Gaussian2DFunction.X_POSITION];
      params[n * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.Y_POSITION] +=
          centrey - params0[Gaussian2DFunction.Y_POSITION];
    }
    f.initialise(params);
    if (mask == null || !withinMask) {
      for (int k = 0; k < size; k++) {
        overlap[k] += f.eval(k);
      }
    } else {
      for (int k = 0; k < size; k++) {
        if (mask[k]) {
          overlap[k] += f.eval(k);
        }
      }
    }

    // // Debug
    // double[] combined = data.clone();
    // for (int k = 0; k < size; k++)
    // combined[k] += overlap[k];
    // uk.ac.sussex.gdsc.core.ij.Utils.display("Spot i", data, maxx, maxx);
    // uk.ac.sussex.gdsc.core.ij.Utils.display("Overlap", overlap, maxx, maxx);
    // uk.ac.sussex.gdsc.core.ij.Utils.display("Combined", combined, maxx, maxx);
    // System.out.printf("Signal %.2f, sd %.2f, overlap = %s\n", params0[Gaussian2DFunction.SIGNAL],
    // params0[Gaussian2DFunction.X_SD], Arrays.toString(params));
  }

  /**
   * Get the overlap data.
   *
   * <p>Computes the sum of the central function, and the sum of the overlap
   *
   * @return The data [sum f1, sum overlap]
   */
  public double[] getOverlapData() {
    final double[] result = new double[2];
    if (overlap != null) {
      double sumF = 0;
      double sumO = 0;
      if (mask == null) {
        for (int k = 0; k < size; k++) {
          sumF += data[k];
          sumO += overlap[k];
        }
      } else {
        for (int k = 0; k < size; k++) {
          if (mask[k]) {
            sumF += data[k];
            sumO += overlap[k];
          }
        }
      }

      result[0] = sumF;
      result[1] = sumO;
    }
    return result;
  }

  /**
   * Get the weighted background.
   *
   * <p>Computes the convolution of the central function and the overlap for the central pixel of
   * the region. This is an estimate of the background contributed to the region by overlapping
   * functions.
   *
   * <p>The result of this function is effected by how the overlap was computed, either within the
   * mask or within the entire square region (see {@link #add(double[], boolean)})
   *
   * @return The weighted background
   */
  public double getWeightedbackground() {
    if (overlap != null) {
      double sum = 0;
      double sumw = 0;
      final double norm = 1.0 / params0[Gaussian2DFunction.SIGNAL];
      // double[] combined = new double[size];
      for (int k = 0; k < size; k++) {
        final double v = data[k] * norm;
        sum += overlap[k] * v;
        sumw += v;
      }
      // uk.ac.sussex.gdsc.core.ij.Utils.display("O Spot i", data, maxx, maxx);
      // uk.ac.sussex.gdsc.core.ij.Utils.display("O Spot overlap", overlap, maxx, maxx);
      return sum / sumw;
    }
    return 0;
  }

  /**
   * Get the probability of a standard Gaussian within the given range x, i.e.
   * {@code P(-x < X <= x)}.
   *
   * @param x (must be positive)
   * @return The probability (0-1)
   */
  public static double getArea(double x) {
    // final NormalDistribution d = new NormalDistribution();
    // return d.probability(-x, x);

    // Since this is symmetrical then we only evaluate half of the error function.
    // return 0.5 * Erf.erf(-x / SQRT2, x / SQRT2);
    // return 0.5 * (Erf.erf(x / SQRT2) - Erf.erf(-x / SQRT2));

    // This only need to be approximate so use a fast error function
    return Erf.erf(x / SQRT2);
  }

  /**
   * Gets the fraction of the function value to use to create the mask region.
   *
   * @return the fraction
   */
  public double getFraction() {
    return fraction;
  }

  /**
   * Sets the fraction of the function value to use to create the mask region. Since a Gaussian 2D
   * function can be circular this can help ignore the corners of the function as it is evaluated on
   * a rectangular region.
   *
   * @param fraction the new fraction
   */
  public void setFraction(double fraction) {
    if (fraction > 1 || fraction < 0) {
      throw new IllegalArgumentException("Fraction must be in the range 0-1");
    }
    this.fraction = fraction;
  }
}
