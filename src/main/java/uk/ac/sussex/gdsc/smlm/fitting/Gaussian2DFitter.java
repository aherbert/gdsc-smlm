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

package uk.ac.sussex.gdsc.smlm.fitting;

import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;

import org.apache.commons.math3.util.FastMath;

import java.util.Arrays;

/**
 * Fits a 2-dimensional Gaussian function for the specified peak. Can optionally fit an elliptical
 * Gaussian function.
 *
 * <p>Performs fitting using the configured algorithm.
 */
public class Gaussian2DFitter {
  /** The fit configuration. */
  protected Gaussian2DFitConfiguration fitConfiguration;
  /** The solver. */
  protected FunctionSolver solver;
  /** The last successful fit. Used to compute the residuals. */
  protected double[] residuals = null;
  /** Allow calculation of residuals to be turned off (overwrite constructor fit configuration). */
  protected boolean computeResiduals = true;
  /** The lower bounds for function solvers. */
  protected double[] lower;
  /** The upper bounds for function solvers. */
  protected double[] upper;

  /**
   * Instantiates a new gaussian 2D fitter.
   *
   * @param fitConfiguration the fit configuration
   */
  public Gaussian2DFitter(Gaussian2DFitConfiguration fitConfiguration) {
    if (fitConfiguration == null) {
      throw new NullPointerException("No fit configuration");
    }
    this.fitConfiguration = fitConfiguration;
    computeResiduals = fitConfiguration.isComputeResiduals();
  }

  private static double half_max_position(double[] data, int index, int[] point, int[] dim,
      int dimension, int[] cumul_region, int dirn, double background) {
    int i;
    int i_start;
    int i_end;
    int i_step;
    final double v = data[index];
    final double v_half = 0.5f * (v + background);
    double v_prev = v;
    double v_this;
    int jump;

    if (dirn == 1) {
      i_start = point[dimension] + 1;
      i_end = dim[dimension];
      i_step = 1;
    } else {
      i_start = point[dimension] - 1;
      i_end = -1;
      i_step = -1;
    }

    jump = i_step * cumul_region[dimension];

    for (i = i_start; i != i_end; i += i_step) {
      index += jump;
      v_this = data[index];

      if (v_this < v_half) {
        return i - i_step * (v_half - v_this) / (v_prev - v_this);
      }

      v_prev = v_this;
    }

    // Not reached the half-max point. Return the dimension limit.
    if (dirn == 1) {
      return dim[dimension];
    }
    return 0f;
  }

  /**
   * Compute the full-width at half-maximum (FWHM) using a line profile through 2D data.
   *
   * @param data the data
   * @param index the index
   * @param point the point
   * @param dim the maximum size of each dimension
   * @param dimension the dimension (0 or 1)
   * @param cumul_region the cumulative region sizes (1, dim[0], dim[0] * dim[1])
   * @param background the background
   * @return the FWHM
   */
  public static double half_max_linewidth(double[] data, int index, int[] point, int[] dim,
      int dimension, int[] cumul_region, double background) {
    double linewidth;
    double a;
    double b;

    a = half_max_position(data, index, point, dim, dimension, cumul_region, 1, background);
    b = half_max_position(data, index, point, dim, dimension, cumul_region, -1, background);

    linewidth = a - b;

    return linewidth;
  }

  /**
   * Accepts a single array containing 2-dimensional data and a list of the peaks to fit. Data
   * should be packed in descending dimension order, e.g. Y,X : Index for [y,z] = MaxX*y + x.
   *
   * <p>Performs fitting using the specified method with a Levenberg-Marquardt algorithm.
   *
   * <p>Adapted from the CCPN fit_peaks routine for Python.
   *
   * @param data The data to fit
   * @param maxx The data size in the x dimension
   * @param maxy The data size in the y dimension
   * @param peaks The index of the peaks
   * @return The fit result
   */
  public FitResult fit(final double[] data, final int maxx, final int maxy, final int[] peaks) {
    return fit(data, maxx, maxy, peaks, null);
  }

  /**
   * Accepts a single array containing 2-dimensional data and a list of the peaks to fit. Data
   * should be packed in descending dimension order, e.g. Y,X : Index for [y,z] = MaxX*y + x.
   *
   * <p>Performs fitting using the specified method with a Levenberg-Marquardt algorithm.
   *
   * @param data The data to fit
   * @param maxx The data size in the x dimension
   * @param maxy The data size in the y dimension
   * @param peaks The index of the peaks (must be within the data bounds if the heights are null)
   * @param heights An initial estimate of the peak heights (can be null)
   * @return The fit result
   */
  public FitResult fit(final double[] data, final int maxx, final int maxy, final int[] peaks,
      double[] heights) {
    final int npeaks = peaks.length;

    final int paramsPerPeak = Gaussian2DFunction.PARAMETERS_PER_PEAK;

    final double[] params = new double[1 + paramsPerPeak * npeaks];

    // Get peak heights (if multiple peaks)
    if (npeaks > 1 && (heights == null || heights.length != peaks.length)) {
      heights = new double[peaks.length];
      for (int i = 0; i < peaks.length; i++) {
        heights[i] = data[peaks[i]];
      }
    }

    final double background = getBackground(data, maxx, maxy, npeaks);

    // Set the initial parameters
    params[Gaussian2DFunction.BACKGROUND] = background;

    final boolean[] amplitudeEstimate = new boolean[npeaks];
    if (npeaks == 1) {
      double sum = 0;
      final int size = maxx * maxy;
      for (int i = size; i-- > 0;) {
        sum += data[i];
      }
      params[Gaussian2DFunction.SIGNAL] = sum - background * size;
      params[Gaussian2DFunction.X_POSITION] = peaks[0] % maxx;
      params[Gaussian2DFunction.Y_POSITION] = peaks[0] / maxx;
    } else {
      for (int i = 0, j = 0; i < peaks.length; i++, j += paramsPerPeak) {
        final int index = peaks[i];
        params[j + Gaussian2DFunction.SIGNAL] = heights[i] - background;
        params[j + Gaussian2DFunction.X_POSITION] = index % maxx;
        params[j + Gaussian2DFunction.Y_POSITION] = index / maxx;
        amplitudeEstimate[i] = true;
      }
    }

    // We have estimated the background already
    return fit(data, maxx, maxy, npeaks, params, amplitudeEstimate, true);
  }

  /**
   * Guess the background from the data given the estimated peak heights.
   *
   * <p>For a single peak the method assumes the peak is in the centre. In this case the minimum
   * average value of the four edges in used.
   *
   * <p>If multiple peaks heights are provided then always use the minimum value in the data since
   * it cannot be assumed that all peaks are away from the edge of the data.
   *
   * @param data the data
   * @param maxx the maxx
   * @param maxy the maxy
   * @param npeaks the number of peaks
   * @return The background estimate
   */
  public static double getBackground(final double[] data, final int maxx, final int maxy,
      final int npeaks) {
    // TODO - What is the best method for setting the background?
    // 1. Min in data
    // 2. Average of border?
    // 3. Evaluate total volume under the initial Gaussian params and subtract that from the sum of
    // the image?
    // (note initial gaussian requires guess of the amplitude which needs background (or bias))

    // -----
    // Noted that if the peak height is negative then fitting becomes unstable.
    // This is likely when fitting multiple peaks since the initial edge guess
    // for the background may be wrong.
    // Use the edge value for single peaks but the minimum value in the data if fitting
    // multiple peaks.
    // -----

    double background = 0;

    // Utils.display("Spot", data, maxx, maxy);

    if (npeaks == 1) {
      //// Set background using the average value of the edge in the data
      // final int s2 = getIndex(0, maxy - 1, maxx);
      // for (int xi = 0, xi2 = s2; xi < maxx; xi++, xi2++)
      // background += data[xi] + data[xi2];
      // for (int yi = maxx, yi2 = getIndex(maxx - 1, 1, maxx); yi < s2; yi += maxx, yi2 += maxx)
      // background += data[yi] + data[yi2];
      // background /= 2 * (maxx + maxy - 2);

      // Liu, et al (2013), Optics Express 21, 29462-87:
      // Get background using the lowest mean of the four edges.
      final int corner2 = (maxy - 1) * maxx;
      double b1 = 0;
      double b2 = 0;
      double b3 = 0;
      double b4 = 0;
      for (int xi = 0, xi2 = corner2; xi < maxx; xi++, xi2++) {
        b1 += data[xi];
        b2 += data[xi2];
      }
      for (int yi = 0, yi2 = maxx - 1; yi <= corner2; yi += maxx, yi2 += maxx) {
        b3 += data[yi];
        b4 += data[yi2];
      }
      //// Debugging which is best
      // System.out.printf("Background %f vs %f (%f)\n", background, min(min(b1, b2) / maxx, min(b3,
      //// b4) / maxy),
      // uk.ac.sussex.gdsc.core.utils.DoubleEquality.relativeError(background,
      // min(min(b1, b2) / maxx, min(b3, b4) / maxy)));
      background = min(min(b1, b2) / maxx, min(b3, b4) / maxy);
    } else {
      // Set background using the minimum value in the data
      background = data[0];
      for (int i = maxx * maxy; --i > 0;) {
        if (background > data[i]) {
          background = data[i];
        }
      }
    }

    return background;
  }

  @SuppressWarnings("unused")
  private static int getIndex(int x, int y, int maxx) {
    return y * maxx + x;
  }

  private static double min(double a, double b) {
    return (a < b) ? a : b;
  }

  /**
   * Accepts a single array containing 2-dimensional data and a list of the peaks to fit. Data
   * should be packed in descending dimension order, e.g. Y,X : Index for [y,z] = MaxX*y + x.
   *
   * <p>Performs fitting using the specified method with a Levenberg-Marquardt algorithm.
   *
   * <p>The input parameter can estimate the signal (the total volume of the Gaussian) or the
   * amplitude (the height of the Gaussian). The signal = amplitude * 2 * pi * sd0 * sd1. The
   * amplitude is the recommended method to estimate parameters for multiple peaks. The signal can
   * be estimated for a single peak by summing all the pixels (minus the background).
   *
   * <p>Note that if the background parameter is zero it will be assumed that the input amplitude is
   * the total height of the peak. A new background will be estimated and the heights lowered by
   * this estimated background to create the amplitude estimate.
   *
   * <p>If a peak location is outside the region bounds and has no input width parameters set or
   * from the fit configuration then fitting will fail (this is because they cannot be estimated).
   *
   * @param data The data to fit
   * @param maxx The data size in the x dimension
   * @param maxy The data size in the y dimension
   * @param npeaks The number of peaks
   * @param params The parameters of the peaks. Must have the Signal/Amplitude,Xpos,Ypos set. Other
   *        parameters that are zero will be estimated.
   * @param amplitudeEstimate Set to true if the peak has amplitude estimated in the
   *        {@link Gaussian2DFunction#SIGNAL} field. The default is signal.
   * @return The fit result
   */
  public FitResult fit(final double[] data, final int maxx, final int maxy, final int npeaks,
      final double[] params, final boolean[] amplitudeEstimate) {
    return fit(data, maxx, maxy, npeaks, params, amplitudeEstimate, false);
  }

  /**
   * Accepts a single array containing 2-dimensional data and a list of the peaks to fit. Data
   * should be packed in descending dimension order, e.g. Y,X : Index for [y,z] = MaxX*y + x.
   *
   * <p>Performs fitting using the specified method with a Levenberg-Marquardt algorithm.
   *
   * <p>The input parameter can estimate the signal (the total volume of the Gaussian) or the
   * amplitude (the height of the Gaussian). The signal = amplitude * 2 * pi * sd0 * sd1. The
   * amplitude is the recommended method to estimate parameters for multiple peaks. The signal can
   * be estimated for a single peak by summing all the pixels (minus the background).
   *
   * <p>Note that if the background parameter is zero it will be assumed that the input
   * signal/amplitude is the total volume/height of the peak. A new background will be estimated and
   * the volume/heights lowered by this estimated background to create the new parameters.
   *
   * <p>If a peak location is outside the region bounds and has no input width parameters set or
   * from the fit configuration then fitting will fail (this is because they cannot be estimated).
   *
   * @param data The data to fit
   * @param maxx The data size in the x dimension
   * @param maxy The data size in the y dimension
   * @param npeaks The number of peaks
   * @param params The parameters of the peaks. Must have the Signal/Amplitude,Xpos,Ypos set. Other
   *        parameters that are zero will be estimated. This can optionally be ignored for the
   *        background parameter which is valid if zero.
   * @param zeroBackground Set to true if a zero value for the background parameter is the estimate
   * @param amplitudeEstimate Set to true if the peak has amplitude estimated in the
   *        {@link Gaussian2DFunction#SIGNAL} field. The default is signal.
   * @return The fit result
   */
  public FitResult fit(final double[] data, final int maxx, final int maxy, final int npeaks,
      double[] params, final boolean[] amplitudeEstimate, final boolean zeroBackground) {
    // Reset
    residuals = null;
    solver = null;

    // Fitting variables
    final int ySize = maxx * maxy;
    // Y must be the correct size
    final double[] y = (data.length == ySize) ? data : Arrays.copyOf(data, ySize); // Value at index

    final int paramsPerPeak = Gaussian2DFunction.PARAMETERS_PER_PEAK; // Fixed for a
                                                                      // Gaussian2DFunction

    double background = params[Gaussian2DFunction.BACKGROUND];
    if (background == 0 && !zeroBackground) {
      // Get background
      background = getBackground(y, maxx, maxy, npeaks);
      params[Gaussian2DFunction.BACKGROUND] = background;

      // Input estimates are either signal or amplitude, both of which are above background.
      // The code below is appropriate if the input estimate is for absolute peak height inc.
      // background.

      // // For a single peak, check the height is above background
      // if (npeaks == 1 && amplitudeEstimate[0] && params[Gaussian2DFunction.SIGNAL] < background)
      // {
      // // Set the background to the min value in the data using the multiple peak option
      // background = getBackground(data, maxx, maxy, 2);
      //
      // // Check if still below background
      // if (params[Gaussian2DFunction.SIGNAL] < background)
      // {
      // // Set the height to the max value in the data
      // double yMax = y[0];
      // for (int i = 1; i < ySize; i++)
      // if (yMax < y[i])
      // yMax = y[i];
      // params[Gaussian2DFunction.SIGNAL] = (float) yMax;
      // }
      // }
      //
      // // Lower heights to get amplitude (if appropriate)
      // for (int j = Gaussian2DFunction.SIGNAL, i = 0; j < params.length; j += paramsPerPeak, i++)
      // {
      // if (amplitudeEstimate[i])
      // params[j] -= background;
      // }
    }

    double[] initialParams = Arrays.copyOf(params, params.length);

    // Check all the heights are valid first
    int zeroHeight = 0;
    for (int j = Gaussian2DFunction.SIGNAL; j < params.length; j += paramsPerPeak) {
      if (params[j] <= 0) {
        zeroHeight++;
      }
    }

    // Check there are peaks to fit
    if (zeroHeight == npeaks) {
      return new FitResult(FitStatus.BAD_PARAMETERS, 0, Double.NaN, initialParams, null, null,
          npeaks, 0, null, 0, 0);
    }

    // Set all zero height peaks to a fraction of the maximum to allow fitting
    if (zeroHeight > 0) {
      // Amplitude estimates
      double max = 0;
      for (int j = Gaussian2DFunction.SIGNAL, i = 0; j < params.length; j += paramsPerPeak, i++) {
        if (amplitudeEstimate[i] && max < params[j]) {
          max = params[j];
        }
      }
      if (max == 0) {
        max = y[0];
        for (int i = maxx * maxy; --i > 0;) {
          if (max < y[i]) {
            max = y[i];
          }
        }
        max -= getBackground(y, maxx, maxy, 2);
      }
      max *= 0.1; // Use fraction of the max peak
      for (int j = Gaussian2DFunction.SIGNAL, i = 0; j < params.length; j += paramsPerPeak, i++) {
        if (amplitudeEstimate[i] && params[j] <= 0) {
          params[j] = max;
        }
      }

      // Signal estimates
      max = 0;
      for (int j = Gaussian2DFunction.SIGNAL, i = 0; j < params.length; j += paramsPerPeak, i++) {
        if (!amplitudeEstimate[i] && max < params[j]) {
          max = params[j];
        }
      }
      if (max == 0) {
        for (int i = maxx * maxy; --i > 0;) {
          max += y[i];
        }
        max -= ySize * getBackground(y, maxx, maxy, 2);
      }
      max *= 0.1; // Use fraction of the max peak
      for (int j = Gaussian2DFunction.SIGNAL, i = 0; j < params.length; j += paramsPerPeak, i++) {
        if (!amplitudeEstimate[i] && params[j] <= 0) {
          params[j] = max;
        }
      }
    }

    if (!checkParameters(maxx, maxy, npeaks, params, amplitudeEstimate, ySize, y, paramsPerPeak,
        background, initialParams)) {
      return new FitResult(FitStatus.FAILED_TO_ESTIMATE_WIDTH, 0, Double.NaN, initialParams, null,
          null, npeaks, 0, null, 0, 0);
    }

    // Re-copy the parameters now they have all been set
    initialParams = params.clone();

    // -----------------------
    // Use alternative fitters
    // -----------------------

    fitConfiguration.initialise(npeaks, maxx, maxy, initialParams);
    solver = fitConfiguration.getFunctionSolver();

    // Bounds are more restrictive than constraints
    if (solver.isBounded()) {
      // Input configured bounds
      setBounds(maxx, maxy, npeaks, params, y, ySize, paramsPerPeak, this.lower, this.upper);
    } else if (solver.isConstrained()) {
      setConstraints(maxx, maxy, npeaks, params, y, ySize, paramsPerPeak);
    }

    final double[] paramsDev =
        (fitConfiguration.isComputeDeviations()) ? new double[params.length] : null;
    if (computeResiduals) {
      residuals = new double[ySize];
    }
    FitStatus result = solver.fit(y, residuals, params, paramsDev);

    // -----------------------

    if (result == FitStatus.OK) {
      // For debugging
      // double[] initialParams2 = initialParams.clone();
      // double[] params2 = params.clone();

      // Re-assemble all the parameters
      if (fitConfiguration.isXSDFitting()) {
        if (fitConfiguration.isYSDFitting()) {
          // Ensure widths are positive
          for (int i = 0, j = 0; i < npeaks; i++, j += paramsPerPeak) {
            params[j + Gaussian2DFunction.X_SD] = Math.abs(params[j + Gaussian2DFunction.X_SD]);
            params[j + Gaussian2DFunction.Y_SD] = Math.abs(params[j + Gaussian2DFunction.Y_SD]);
          }
        } else {
          // Ensure Y width is updated with the fitted X width
          for (int i = 0, j = 0; i < npeaks; i++, j += paramsPerPeak) {
            // Ensure width is positive
            params[j + Gaussian2DFunction.X_SD] = Math.abs(params[j + Gaussian2DFunction.X_SD]);
            params[j + Gaussian2DFunction.Y_SD] = params[j + Gaussian2DFunction.X_SD];
            if (paramsDev != null) {
              paramsDev[j + Gaussian2DFunction.Y_SD] = paramsDev[j + Gaussian2DFunction.X_SD];
            }
          }
        }
      }
      if (fitConfiguration.isAngleFitting()) {
        // Ensure the angle is within the correct bounds
        for (int i = Gaussian2DFunction.ANGLE; i < params.length; i += paramsPerPeak) {
          correctAngle(i, params, paramsDev);
        }
      }

      Object statusData = null;
      if (fitConfiguration.isFitValidation()) {
        result = fitConfiguration.validateFit(npeaks, initialParams, params, paramsDev);
        statusData = fitConfiguration.getValidationData();
      }

      if (residuals != null) {
        for (int i = 0; i < residuals.length; i++) {
          residuals[i] = y[i] - residuals[i];
        }
      }

      return new FitResult(result, FastMath.max(ySize - solver.getNumberOfFittedParameters(), 0),
          solver.getValue(), initialParams, params, paramsDev, npeaks,
          solver.getNumberOfFittedParameters(), statusData, solver.getIterations(),
          solver.getEvaluations());
    }
    residuals = null;
    return new FitResult(result, 0, Double.NaN, initialParams, null, null, npeaks,
        solver.getNumberOfFittedParameters(), null, solver.getIterations(),
        solver.getEvaluations());
  }

  /**
   * Check parameters.
   *
   * @param maxx the maxx
   * @param maxy the maxy
   * @param npeaks the npeaks
   * @param params the params
   * @param amplitudeEstimate the amplitude estimate
   * @param ySize the y size
   * @param y the y
   * @param paramsPerPeak the params per peak
   * @param background the background
   * @param initialParams the initial params
   * @return true, if successful
   */
  protected boolean checkParameters(final int maxx, final int maxy, final int npeaks,
      double[] params, final boolean[] amplitudeEstimate, final int ySize, final double[] y,
      final int paramsPerPeak, double background, double[] initialParams) {
    final int[] dim = new int[] {maxx, maxy};
    final int[] position = new int[2];
    final int[] cumul_region = new int[] {1, maxx, ySize};
    for (int i = 0, j = 0; i < npeaks; i++, j += paramsPerPeak) {
      // ----
      // Check all input parameters and estimate them if necessary
      // ----

      // Get the parameters
      double signal = params[j + Gaussian2DFunction.SIGNAL];
      double xpos = params[j + Gaussian2DFunction.X_POSITION];
      double ypos = params[j + Gaussian2DFunction.Y_POSITION];

      // Set-up for estimating peak width at half maximum
      position[0] = (int) Math.round(xpos);
      position[1] = (int) Math.round(ypos);
      final int index = position[1] * maxx + position[0];

      double sx;
      double sy;
      double angle;
      if (fitConfiguration.isZFitting()) {
        // Use the widths at z=0.
        // These are used to determine the centre-of-mass range search.
        // It does not matter if they are negative as we use max(1, sx+sy)
        // to search for the centre. It will also effect the conversion of
        // amplitudes to signal.
        sx = fitConfiguration.getInitialXSD();
        sy = fitConfiguration.getInitialYSD();
        angle = 0;
      } else {

        sx = params[j + Gaussian2DFunction.X_SD];
        sy = params[j + Gaussian2DFunction.Y_SD];
        angle = params[j + Gaussian2DFunction.ANGLE];

        if (sx == 0) {
          if (fitConfiguration.getInitialXSD() > 0) {
            sx = fitConfiguration.getInitialXSD();
          } else {
            // Fail if the width cannot be estimated due to out of bounds
            if (position[0] < 0 || position[0] > maxx || position[1] < 0 || position[1] > maxy) {
              return false;
            }

            sx = fwhm2sd(half_max_linewidth(y, index, position, dim, 0, cumul_region, background));
          }
        }

        if (sy == 0) {
          if (fitConfiguration.isYSDFitting()) {
            if (fitConfiguration.getInitialYSD() > 0) {
              sy = fitConfiguration.getInitialYSD();
            } else {
              // Fail if the width cannot be estimated
              if (position[0] < 0 || position[0] > maxx || position[1] < 0 || position[1] > maxy) {
                return false;
              }

              sy = fwhm2sd(
                  half_max_linewidth(y, index, position, dim, 1, cumul_region, background));
            }
          } else {
            sy = sx;
          }
        }

        // Guess the initial angle if input angle is out-of-bounds
        if (angle == 0) {
          if (fitConfiguration.isAngleFitting() && fitConfiguration.getInitialAngle() >= -Math.PI
              && fitConfiguration.getInitialAngle() <= -Math.PI) {
            if (sx != sy) {
              // There is no angle gradient information if the widths are equal. Zero and it will be
              // ignored
              angle = fitConfiguration.getInitialAngle();
            }
          }
        }

      }

      // If the position is on the integer grid then use a centre-of-mass approximation
      if (npeaks == 1 && xpos == position[0] && ypos == position[1]) {
        // Estimate using centre of mass around peak index
        // Use 2 * SD estimate to calculate the range around the index that should be considered.
        // SD = (sx+sy)/2 => Range = sx+sy
        final int range = Math.max(1, (int) Math.ceil(sx + sy));
        final double[] com = findCentreOfMass(y, dim, range, position);
        xpos = com[0];
        ypos = com[1];
      }

      // Convert amplitudes to signal
      if (amplitudeEstimate[i]) {
        signal *= 2 * Math.PI * sx * sy;
      }

      // Set all the parameters
      params[j + Gaussian2DFunction.SIGNAL] = signal;
      params[j + Gaussian2DFunction.X_POSITION] = xpos;
      params[j + Gaussian2DFunction.Y_POSITION] = ypos;
      params[j + Gaussian2DFunction.X_SD] = sx;
      params[j + Gaussian2DFunction.Y_SD] = sy;
      params[j + Gaussian2DFunction.ANGLE] = angle;
      // Leave the z-position (i.e. do not reset to zero)
      // params[j + Gaussian2DFunction.Z_POSITION] = 0;
    }

    return true;
  }

  /**
   * Finds the centre of the image using the centre of mass within the given range of the specified
   * centre-of-mass.
   *
   * @param subImage the sub image
   * @param dimensions the dimensions
   * @param range the range
   * @param centre the centre
   * @return the centre
   */
  static double[] findCentreOfMass(final double[] subImage, final int[] dimensions, final int range,
      final int[] centre) {
    final int[] min = new int[2];
    final int[] max = new int[2];
    for (int i = 2; i-- > 0;) {
      min[i] = centre[i] - range;
      max[i] = centre[i] + range;
      if (min[i] < 0) {
        min[i] = 0;
      }
      if (max[i] >= dimensions[i] - 1) {
        max[i] = dimensions[i] - 1;
      }
    }

    final double[] newCom = new double[2];
    double sum = 0;
    for (int y = min[1]; y <= max[1]; y++) {
      int index = dimensions[0] * y + min[0];
      for (int x = min[0]; x <= max[0]; x++, index++) {
        final double value = subImage[index];
        sum += value;
        newCom[0] += x * value;
        newCom[1] += y * value;
      }
    }

    for (int i = 2; i-- > 0;) {
      newCom[i] /= sum;
    }

    return newCom;
  }

  /**
   * Sets the bounds for the fitted parameters.
   *
   * @param maxx The x range of the data
   * @param maxy The y range of the data
   * @param npeaks The number of peaks
   * @param params The estimated parameters
   * @param y The data
   * @param ySize The size of the data
   * @param paramsPerPeak The number of parameters per peak
   * @param lower2 the input lower bounds
   * @param upper2 the input upper bounds
   */
  protected void setBounds(final int maxx, final int maxy, final int npeaks, final double[] params,
      final double[] y, final int ySize, final int paramsPerPeak, double[] lower2,
      double[] upper2) {
    // Create appropriate bounds for the parameters
    final double[] lower = new double[params.length];
    final double[] upper = new double[lower.length];
    double yMax = y[0];
    double yMin = y[0];
    for (int i = 1; i < ySize; i++) {
      if (yMax < y[i]) {
        yMax = y[i];
      } else if (yMin > y[i]) {
        yMin = y[i];
      }
    }
    if (fitConfiguration.isBackgroundFitting()) {
      if (yMax > params[Gaussian2DFunction.BACKGROUND]) {
        upper[Gaussian2DFunction.BACKGROUND] = yMax;
      } else {
        upper[Gaussian2DFunction.BACKGROUND] =
            params[Gaussian2DFunction.BACKGROUND] + (params[Gaussian2DFunction.BACKGROUND] - yMax);
      }

      if (yMin < params[Gaussian2DFunction.BACKGROUND]) {
        lower[Gaussian2DFunction.BACKGROUND] = yMin;
      } else {
        lower[Gaussian2DFunction.BACKGROUND] =
            params[Gaussian2DFunction.BACKGROUND] - (yMin - params[Gaussian2DFunction.BACKGROUND]);
      }

      if (lower[Gaussian2DFunction.BACKGROUND] < 0) {
        // This is a problem for MLE fitting
        if (solver.isStrictlyPositiveFunction()) {
          lower[Gaussian2DFunction.BACKGROUND] = 0;
        }
      }
    }

    // Configure the bounds for the width.
    // The factors are less strict than the fit configuration to allow some search space when
    // fitting close to the limits.
    final double min_wf;
    final double max_wf;
    final boolean isZFitting = fitConfiguration.isZFitting();
    if (isZFitting) {
      min_wf = 0;
      max_wf = Double.MAX_VALUE;
    } else {
      min_wf = getMinWidthFactor();
      max_wf = getMaxWidthFactor();
    }

    // Get the upper bounds for the width factor. This is just used to estimate the upper bounds for
    // the signal
    // So it does not matter if it is too wrong.
    final double wf = (max_wf < Double.MAX_VALUE) ? fitConfiguration.getMaxWidthFactor() : 3;

    if (npeaks == 1) {
      // Allow the signal to explain all the data. This assumes the data window entirely covers the
      // spot.
      double sum = 0;
      for (int i = 1; i < ySize; i++) {
        sum += y[i];
      }
      // Increase sum by 2 to allow for error
      upper[Gaussian2DFunction.SIGNAL] = 2 * sum - yMin * ySize;
    } else {
      final double height = yMax - yMin;
      // Signal = height * 2 * pi * sd0 * sd1
      // Allow a maximum using the width factor that defines the bounds on the width.
      // Increase the height by 2 to allow for error.
      final double factor = 2 * height * 2 * Math.PI * wf * wf;
      // System.out.printf("%f or %f\n", upper[Gaussian2DFunction.SIGNAL], factor *
      // params[Gaussian2DFunction.X_SD] *
      // params[Gaussian2DFunction.Y_SD]);
      for (int i = 0, j = 0; i < npeaks; i++, j += paramsPerPeak) {
        upper[j + Gaussian2DFunction.SIGNAL] =
            factor * params[j + Gaussian2DFunction.X_SD] * params[j + Gaussian2DFunction.Y_SD];
      }
    }

    for (int i = 0, j = 0; i < npeaks; i++, j += paramsPerPeak) {
      // Check if the signal bounds are appropriate
      if (params[j + Gaussian2DFunction.SIGNAL] < lower[j + Gaussian2DFunction.SIGNAL]) {
        lower[j + Gaussian2DFunction.SIGNAL] = params[j + Gaussian2DFunction.SIGNAL]
            - (lower[j + Gaussian2DFunction.SIGNAL] - params[j + Gaussian2DFunction.SIGNAL]);
      }
      if (lower[j + Gaussian2DFunction.SIGNAL] <= 0) {
        // This is a problem for MLE fitting
        if (solver.isStrictlyPositiveFunction()) {
          // If this is zero it causes problems when computing gradients since the
          // Gaussian function may not exist. So use a small value instead.
          lower[j + Gaussian2DFunction.SIGNAL] = 0; // .1;
        }
      }
      if (params[j + Gaussian2DFunction.SIGNAL] > upper[j + Gaussian2DFunction.SIGNAL]) {
        upper[j + Gaussian2DFunction.SIGNAL] = params[j + Gaussian2DFunction.SIGNAL]
            + (params[j + Gaussian2DFunction.SIGNAL] - upper[j + Gaussian2DFunction.SIGNAL]);
      }

      // All functions evaluate the x and y position.
      // Lower bounds on these will be zero when the array is initialised.
      // We may have an estimate outside the bounds (if including neighbours).
      upper[j + Gaussian2DFunction.X_POSITION] = Math.max(maxx,
          params[j + Gaussian2DFunction.X_POSITION] + params[j + Gaussian2DFunction.X_SD]);
      upper[j + Gaussian2DFunction.Y_POSITION] = Math.max(maxy,
          params[j + Gaussian2DFunction.Y_POSITION] + params[j + Gaussian2DFunction.Y_SD]);

      lower[j + Gaussian2DFunction.X_POSITION] = Math.min(0,
          params[j + Gaussian2DFunction.X_POSITION] - params[j + Gaussian2DFunction.X_SD]);
      lower[j + Gaussian2DFunction.Y_POSITION] = Math.min(0,
          params[j + Gaussian2DFunction.Y_POSITION] - params[j + Gaussian2DFunction.Y_SD]);

      if (fitConfiguration.isAngleFitting()) {
        lower[j + Gaussian2DFunction.ANGLE] = -Math.PI;
        upper[j + Gaussian2DFunction.ANGLE] = Math.PI;
      }
      // TODO - Add support for z-depth fitting. Currently this is unbounded.
      if (isZFitting) {
        lower[j + Gaussian2DFunction.Z_POSITION] = Double.NEGATIVE_INFINITY;
        upper[j + Gaussian2DFunction.Z_POSITION] = Double.POSITIVE_INFINITY;
        // The widths are not fit but set simple limits
        // to avoid errors when checking the limits
        upper[j + Gaussian2DFunction.X_SD] = params[j + Gaussian2DFunction.X_SD];
        upper[j + Gaussian2DFunction.Y_SD] = params[j + Gaussian2DFunction.Y_SD];
      } else {
        lower[j + Gaussian2DFunction.X_SD] = params[j + Gaussian2DFunction.X_SD] * min_wf;
        upper[j + Gaussian2DFunction.X_SD] = params[j + Gaussian2DFunction.X_SD] * max_wf;
        lower[j + Gaussian2DFunction.Y_SD] = params[j + Gaussian2DFunction.Y_SD] * min_wf;
        upper[j + Gaussian2DFunction.Y_SD] = params[j + Gaussian2DFunction.Y_SD] * max_wf;
      }
    }

    if (solver.isStrictlyPositiveFunction()) {
      // If the lower bounds are zero it causes problems when computing gradients since
      // the Gaussian function may not exist. So use a small value instead.
      for (int i = 0, j = 0; i < npeaks; i++, j += paramsPerPeak) {
        if (lower[j + Gaussian2DFunction.SIGNAL] <= 0) {
          lower[j + Gaussian2DFunction.SIGNAL] = 0.1;
        }
        if (!isZFitting) {
          if (lower[j + Gaussian2DFunction.X_SD] <= 0) {
            lower[j + Gaussian2DFunction.X_SD] = 0.01;
          }
          if (lower[j + Gaussian2DFunction.Y_SD] <= 0) {
            lower[j + Gaussian2DFunction.Y_SD] = 0.01;
          }
        }
      }
    }

    // Check against the configured bounds
    if (lower2 != null) {
      for (int i = Math.max(lower.length, lower2.length); i-- > 0;) {
        if (lower[i] < lower2[i]) {
          lower[i] = lower2[i];
        }
      }
    }
    if (upper2 != null) {
      for (int i = Math.max(upper.length, upper2.length); i-- > 0;) {
        if (upper[i] > upper2[i]) {
          upper[i] = upper2[i];
        }
      }
    }

    // Debug check
    for (int i = (fitConfiguration.isBackgroundFitting()) ? 0 : 1; i < params.length; i++) {
      if (params[i] < lower[i]) {
        System.out.printf("Param %d (%s) too low %f < %f\n", i, solver.getName(i), params[i],
            lower[i]);
        lower[i] = params[i] - (lower[i] - params[i]);
      }
      if (params[i] > upper[i]) {
        System.out.printf("Param %d (%s) too high %f > %f\n", i, solver.getName(i), params[i],
            upper[i]);
        upper[i] = params[i] + (params[i] - upper[i]);
      }
    }

    solver.setBounds(lower, upper);
  }

  /**
   * Gets the min width factor for the lower bounds.
   *
   * @return the min width factor
   */
  private double getMinWidthFactor() {
    if (fitConfiguration.getMinWidthFactor() < 1 && fitConfiguration.getMinWidthFactor() >= 0) {
      // Add some buffer to allow fitting to go past then come back to the limit
      // This also allows fitting to fail if the spot is definitely smaller than the configured
      // limits.
      return fitConfiguration.getMinWidthFactor() * 0.7;
    }
    return 0;
  }

  /**
   * Gets the max width factor for the upper bounds.
   *
   * @return the max width factor
   */
  private double getMaxWidthFactor() {
    if (fitConfiguration.getMaxWidthFactor() > 1) {
      // Add some buffer to allow fitting to go past then come back to the limit.
      // This also allows fitting to fail if the spot is definitely bigger than the configured
      // limits.
      return fitConfiguration.getMaxWidthFactor() * 1.5;
    }
    return Double.MAX_VALUE;
  }

  /**
   * Sets the constraints for the fitted parameters. This functions set the lower bounds of the
   * signal to zero and background to zero (or negative if the background estimate is {@code < 0}).
   *
   * @param maxx The x range of the data
   * @param maxy The y range of the data
   * @param npeaks The number of peaks
   * @param params The estimated parameters
   * @param y The data
   * @param ySize The size of the data
   * @param paramsPerPeak The number of parameters per peak
   */
  protected void setConstraints(final int maxx, final int maxy, final int npeaks,
      final double[] params, final double[] y, final int ySize, final int paramsPerPeak) {
    // Create appropriate bounds for the parameters
    final double[] lower = new double[params.length];
    final double[] upper = new double[lower.length];
    Arrays.fill(lower, Float.NEGATIVE_INFINITY);
    Arrays.fill(upper, Float.POSITIVE_INFINITY);

    lower[Gaussian2DFunction.BACKGROUND] = 0;
    // If the bias is subtracted then we may have negative data and a background estimate that is
    // negative
    if (params[Gaussian2DFunction.BACKGROUND] < 0) {
      double yMin = 0;
      for (int i = 0; i < ySize; i++) {
        if (yMin > y[i]) {
          yMin = y[i];
        }
      }
      if (yMin < params[Gaussian2DFunction.BACKGROUND]) {
        lower[Gaussian2DFunction.BACKGROUND] = yMin;
      } else {
        lower[Gaussian2DFunction.BACKGROUND] =
            params[Gaussian2DFunction.BACKGROUND] - (yMin - params[Gaussian2DFunction.BACKGROUND]);
      }
    }

    for (int j = Gaussian2DFunction.SIGNAL; j < lower.length; j += paramsPerPeak) {
      lower[j] = 0;
    }
    solver.setConstraints(lower, upper);
  }

  /**
   * Swap the axes so that the major axis is the X axis. Correct the fit angle to lie within the
   * 0-180 degree domain from the major-axis.
   *
   * @param i The angle position within the parameter array
   * @param params the params
   * @param paramsDev the params deveations
   */
  protected void correctAngle(final int i, final double[] params, final double[] paramsDev) {
    double angle = params[i];

    final double twicePI = 2 * Math.PI;
    double fixed = (angle + Math.PI) % twicePI;
    if (fixed < 0) {
      fixed += twicePI;
    }
    angle = fixed - Math.PI;

    // // Angle should now be in -180 - 180 degrees domain
    // if (angle < -Math.PI || angle > Math.PI)
    // {
    // System.out.printf("angle error %g != %g\n", angle, Math.asin(Math.sin(params[i])));
    // }

    // Commented out as this interferes with the PSF Estimator
    final int ix = i + Gaussian2DFunction.X_SD - Gaussian2DFunction.ANGLE;
    final int iy = i + Gaussian2DFunction.Y_SD - Gaussian2DFunction.ANGLE;
    final double xWidth = params[ix];
    final double yWidth = params[iy];
    // The fit will compute the angle from the major axis.
    // Standardise so it is always from the X-axis
    if (yWidth > xWidth) {
      swap(ix, iy, params);
      if (paramsDev != null) {
        swap(ix, iy, paramsDev);
      }

      // Rotate 90 degrees
      angle += Math.PI / 2.0;
      // Check domain
      if (angle > Math.PI) {
        angle -= 2 * Math.PI;
      }
    }

    // Return in 0 - 180 degrees domain since the Gaussian has 2-fold symmetry,
    // i.e. angle -10 == 170
    params[i] = (angle < 0) ? angle + Math.PI : angle;
  }

  private static void swap(final int i, final int j, final double[] params) {
    final double tmp = params[i];
    params[i] = params[j];
    params[j] = tmp;
  }

  /**
   * Convert the Full-Width at Half-Maximum to the Standard Deviation.
   *
   * @param fwhm the fwhm
   * @return sd
   */
  public static double fwhm2sd(double fwhm) {
    return fwhm / (2 * Math.sqrt(2 * Math.log(2)));
  }

  /**
   * Convert the Standard Deviation to the Full-Width at Half-Maximum.
   *
   * @param sd the sd
   * @return fwhm
   */
  public static double sd2fwhm(final double sd) {
    return sd * 2 * Math.sqrt(2 * Math.log(2));
  }

  /**
   * @return the residuals from the last successful fit. If fitting failed then this is null.
   */
  public double[] getResiduals() {
    return residuals;
  }

  /**
   * @return the computeResiduals.
   */
  public boolean isComputeResiduals() {
    return computeResiduals;
  }

  /**
   * @param computeResiduals Set to true to compute the residuals
   */
  public void setComputeResiduals(final boolean computeResiduals) {
    this.computeResiduals = computeResiduals;
  }

  /**
   * Return true if the last call to a fit(...) method created a function solver. This allows the
   * properties to be accessed for the last fit. Otherwise the properties will return zero.
   *
   * @return True if the last call to a fit(...) method created a function solver
   */
  public boolean solvedLastFit() {
    return (solver != null);
  }

  /**
   * Gets the function solver from the last call to a fit(...) method.
   *
   * @return the function solver
   */
  public FunctionSolver getFunctionSolver() {
    return solver;
  }

  /**
   * Set the maximum width factor to use to set the bounds for width parameter fitting. If the fit
   * configuration has a smaller width factor then that will be used instead.
   *
   * <p>The bounds are set using the initial estimate w in the range w*min to w*max.
   *
   * @param maximumWidthFactor the maximum width factor to use to set the bounds for width parameter
   *        fitting (must be above 1)
   * @throws IllegalArgumentException if the arguments are outside the valid ranges
   */
  public void setWidthFactor(
      // double minimumWidthFactor,
      double maximumWidthFactor) {
    // if (minimumWidthFactor < 0 || minimumWidthFactor > 1)
    // throw new IllegalArgumentException("MinWidth factor must be in the range 0-1: " +
    // minimumWidthFactor);
    if (maximumWidthFactor < 1) {
      throw new IllegalArgumentException(
          "MaxWidth factor must be 1 or above: " + maximumWidthFactor);
      // this.minimumWidthFactor = minimumWidthFactor;
      // this.maximumWidthFactor = maximumWidthFactor;
    }
  }

  /**
   * @return the optimised function value for the last fit.
   */
  public double getValue() {
    return (solver != null) ? solver.getValue() : 0;
  }

  /**
   * Sets the bounds on the parameter array in the next call to the fit() method.
   *
   * <p>Bounds must be the same length as the parameter array otherwise they are ignored. Bounds
   * should be cleared when finished by passing in null.
   *
   * @param lower the lower bounds
   * @param upper the upper bounds
   */
  public void setBounds(double[] lower, double[] upper) {
    this.lower = lower;
    this.upper = upper;
  }

  /**
   * Evaluate the function using the configured fit solver. No fit is attempted. The input
   * parameters are assumed to be the parameters of the Gaussian.
   *
   * @param data The data to fit
   * @param maxx The data size in the x dimension
   * @param maxy The data size in the y dimension
   * @param npeaks The number of peaks
   * @param params The parameters of the peaks, e.g. from a previous call to fit(...)
   * @return True if the function was evaluated
   */
  public boolean evaluate(double[] data, int maxx, int maxy, int npeaks, double[] params) {
    final int ySize = maxx * maxy;
    // The solvers require y to be the correct length
    final double[] y = (data.length == ySize) ? data : Arrays.copyOf(data, ySize);
    residuals = (computeResiduals) ? new double[ySize] : null;

    // -----------------------
    // Use alternative fitters
    // -----------------------

    fitConfiguration.initialise(npeaks, maxx, maxy, params);
    solver = fitConfiguration.getFunctionSolver();

    // Note: Do not apply bounds and constraints as it is assumed the input parameters are good

    final boolean result = solver.evaluate(y, residuals, params);

    if (result) {
      if (residuals != null) {
        for (int i = 0; i < residuals.length; i++) {
          residuals[i] = y[i] - residuals[i];
        }
      }
    }

    return result;
  }

  /**
   * Compute the deviations of the parameters of the function using the configured fit solver. No
   * fit is attempted. The input parameters are assumed to be the parameters of the Gaussian.
   *
   * @param data The data to fit
   * @param maxx The data size in the x dimension
   * @param maxy The data size in the y dimension
   * @param npeaks The number of peaks
   * @param params The parameters of the peaks, e.g. from a previous call to fit(...)
   * @param paramsDev The parameters deviations (output)
   * @return True if the function was evaluated
   */
  public boolean computeDeviations(double[] data, int maxx, int maxy, int npeaks, double[] params,
      double[] paramsDev) {
    final int ySize = maxx * maxy;
    // The solvers require y to be the correct length
    final double[] y = (data.length == ySize) ? data : Arrays.copyOf(data, ySize);

    // -----------------------
    // Use alternative fitters
    // -----------------------

    fitConfiguration.initialise(npeaks, maxx, maxy, params);
    solver = fitConfiguration.getFunctionSolver();

    // Note: Do not apply bounds and constraints as it is assumed the input parameters are good

    final boolean result = solver.computeDeviations(y, params, paramsDev);

    if (result) {
      // Ensure the deviations are copied for a symmetric Gaussian
      if (!fitConfiguration.isYSDFitting() && fitConfiguration.isXSDFitting()) {
        // Ensure Y deviation is updated with the X deviation
        for (int i = 0, j = 0; i < npeaks; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK) {
          paramsDev[j + Gaussian2DFunction.Y_SD] = paramsDev[j + Gaussian2DFunction.X_SD];
        }
      }
    }

    return result;
  }
}
