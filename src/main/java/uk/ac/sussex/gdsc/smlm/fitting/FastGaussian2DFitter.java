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

/**
 * Fits a 2-dimensional Gaussian function for the specified peak. Can optionally fit an elliptical
 * Gaussian function.
 *
 * <p>Performs fitting using the configured algorithm.
 *
 * <p>This is based on the Gaussian2DFitter class but does not support estimating the width of each
 * Gaussian. Widths must be provided in the fit configuration. Settings from the fit configuration
 * are cached and thus updates to the configuration after construction may be ignored.
 */
public class FastGaussian2DFitter extends Gaussian2DFitter {
  // Cache the fitting defaults
  private final boolean isZFitting;
  private final boolean isWidth1Fitting;
  private final boolean isAngleFitting;
  private final double angle;
  private final double sx;
  private final double sy;

  /**
   * Constructor.
   *
   * @param fitConfiguration the fit configuration
   * @throws IllegalArgumentException If the configuration is missing information, e.g. initial
   *         widths
   */
  public FastGaussian2DFitter(Gaussian2DFitConfiguration fitConfiguration) {
    super(fitConfiguration);

    // Note: Even if the Gaussian is z fitting there will be an initial estimate for
    // the width (i.e. at z=0).

    // Cache the estimate for the Gaussian
    if (fitConfiguration.getInitialXSD() > 0) {
      sx = fitConfiguration.getInitialXSD();
    } else {
      throw new IllegalArgumentException("No initial width0 estimate");
    }

    isWidth1Fitting = fitConfiguration.isYSDFitting();
    if (isWidth1Fitting) {
      if (fitConfiguration.getInitialYSD() > 0) {
        sy = fitConfiguration.getInitialYSD();
      } else {
        throw new IllegalArgumentException("No initial width1 estimate");
      }
    } else {
      sy = sx;
    }

    isZFitting = fitConfiguration.isZFitting();
    if (isZFitting) {
      // Z determines the width. Only support non-rotated 2D gaussian.
      angle = 0;
      isAngleFitting = false;
      // No cache of initial estimate for z. This is assumed to be zero.
    } else {
      isAngleFitting = fitConfiguration.isAngleFitting();
      if (isAngleFitting) {
        if (fitConfiguration.getInitialAngle() >= -Math.PI
            && fitConfiguration.getInitialAngle() <= Math.PI) {
          angle = fitConfiguration.getInitialAngle();
        } else {
          throw new IllegalArgumentException("No initial angle estimate");
        }
      } else {
        angle = 0;
      }
    }
  }

  @Override
  protected boolean checkParameters(final int maxx, final int maxy, final int npeaks,
      double[] params, final boolean[] amplitudeEstimate, final int ySize, final double[] y,
      final int paramsPerPeak, double background, double[] initialParams) {
    final int[] dim = new int[] {maxx, maxy};
    final int[] position = new int[2];
    for (int i = 0, j = 0; i < npeaks; i++, j += paramsPerPeak) {
      // ----
      // Check all input parameters and uses the default values if necessary
      // ----

      // Get the parameters
      double signal = params[j + Gaussian2DFunction.SIGNAL];
      double xpos = params[j + Gaussian2DFunction.X_POSITION];
      double ypos = params[j + Gaussian2DFunction.Y_POSITION];

      double sx;
      double sy;
      double angle;
      if (isZFitting) {
        // Use the widths at z=0.
        // These are used to determine the centre-of-mass range search.
        sx = this.sx;
        sy = this.sy;
        angle = 0;
      } else {
        sx = params[j + Gaussian2DFunction.X_SD];
        sy = params[j + Gaussian2DFunction.Y_SD];
        angle = params[j + Gaussian2DFunction.ANGLE];

        if (sx == 0) {
          sx = this.sx;
        }

        if (isWidth1Fitting) {
          if (sy == 0) {
            sy = this.sy;
          }
        } else {
          sy = sx;
        }

        // Guess the initial angle if input angle is out-of-bounds
        if (isAngleFitting) {
          if (angle == 0) {
            angle = this.angle;
          }
        }
      }

      // Set-up for estimating peak width at half maximum
      position[0] = (int) Math.round(xpos);
      position[1] = (int) Math.round(ypos);

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
        signal *= 6.283185307 * sx * sy; // 2 * Math.PI * sx * sy
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
}
