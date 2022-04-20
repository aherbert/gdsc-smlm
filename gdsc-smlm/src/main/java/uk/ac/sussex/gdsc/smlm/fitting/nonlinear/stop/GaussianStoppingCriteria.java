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

package uk.ac.sussex.gdsc.smlm.fitting.nonlinear.stop;

import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.lang3.ArrayUtils;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.StoppingCriteria;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;

/**
 * Defines the stopping criteria for the
 * {@link uk.ac.sussex.gdsc.smlm.fitting.nonlinear.NonLinearFit} class.
 *
 * <p>Stop when successive iterations with a reduced error move the fitted X,Y coordinates by less
 * than a specified distance (delta).
 *
 * <p>The criteria also ensure that signal, coordinates and peak-widths are held positive, otherwise
 * fitting is stopped.
 */
public class GaussianStoppingCriteria extends StoppingCriteria {
  private double delta = 0.01;

  /** The number of peaks. */
  protected int peaks;

  /** The function. */
  protected Gaussian2DFunction func;

  private double minimumSignal = Float.NEGATIVE_INFINITY;
  private double[] minimumPosition;
  private double[] maximumPosition;
  private double[] minimumSd;
  private double[] maximumSd;

  /**
   * Instantiates a new gaussian stopping criteria.
   *
   * @param func The Gaussian function
   */
  public GaussianStoppingCriteria(Gaussian2DFunction func) {
    this.func = func;
    this.peaks = func.getNPeaks();
  }

  @Override
  public void evaluate(double oldError, double newError, double[] a) {
    final Logger l = log;
    final StringBuilder sb;
    if (l != null && l.isLoggable(Level.INFO)) {
      sb = logParameters(oldError, newError, a);
    } else {
      sb = null;
    }

    if (newError > oldError) {
      // Fit is worse
      increment(a, false);
    } else {
      // Fit is improved - Check if the movement is negligible
      if (noCoordinateChange(a)) {
        // // Check if all params are within 2sf
        // FloatEquality eq = new FloatEquality(3, 1e-10f);
        // for (int i = 0; i < a.length; i++)
        // {
        // if (!eq.almostEqualComplement(bestA[i], a[i]))
        // System.out.printf("Stopping when still moving: %function => %function (%g)\n%s\n%s\n",
        // bestA[i], a[i], FloatEquality.relativeError(bestA[i], a[i]),
        // Arrays.toString(bestA), Arrays.toString(a));
        // }

        areAchieved = true;
        notSatisfied = false;
      }

      // Check the parameters are still valid
      if (invalidCoordinates(a)) {
        notSatisfied = false;
        if (sb != null) {
          sb.append(" Bad Coords: ").append(Arrays.toString(a));
        }
      }

      increment(a, true);
    }

    if (l != null && l.isLoggable(Level.INFO)) {
      sb.append(" Continue=").append(notSatisfied).append('\n');
      l.info(sb.toString());
    }
  }

  /**
   * Creates a string representation of the peak parameters for logging.
   *
   * @param oldError the old error
   * @param newError the new error
   * @param a The parameters
   * @return The string builder
   */
  protected StringBuilder logParameters(double oldError, double newError, double[] a) {
    final StringBuilder sb = new StringBuilder(128);
    sb.append("iter = ").append(getIteration() + 1).append(", error = ").append(oldError)
        .append(" -> ").append(newError);
    if (newError <= oldError) {

      for (int i = 0; i < peaks; i++) {
        sb.append(", Peak").append(i + 1).append("=[");
        int param = i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_POSITION;
        for (int j = 0; j < 2; j++, param++) {
          if (j > 0) {
            sb.append(',');
          }
          sb.append(a[param] - bestA[param]);
        }
        sb.append(']');
      }
    }
    return sb;
  }

  /**
   * Check if there was no coordinate change.
   *
   * @param a the parameters
   * @return true, if successful
   */
  protected boolean noCoordinateChange(double[] a) {
    for (int i = 0; i < peaks; i++) {
      int param = i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_POSITION;
      for (int j = 0; j < 2; j++, param++) {
        // Check if the coordinates have moved less than the delta limit
        if (Math.abs(bestA[param] - a[param]) > delta) {
          return false;
        }
      }
    }
    return true;
  }

  private boolean invalidCoordinates(double[] a) {
    for (int i = 0; i < peaks; i++) {
      if (a[i * Gaussian2DFunction.PARAMETERS_PER_PEAK
          + Gaussian2DFunction.SIGNAL] < minimumSignal) {
        return true;
      }

      if (isBelow(minimumPosition, a,
          i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_POSITION)) {
        return true;
      }
      if (isAbove(maximumPosition, a,
          i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_POSITION)) {
        return true;
      }

      if (func.evaluatesSD0()) {
        if (isBelow(minimumSd, a,
            i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_SD)) {
          return true;
        }
        if (isAbove(maximumSd, a,
            i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_SD)) {
          return true;
        }
      }

      if (func.evaluatesSD1()) {
        if (isBelow(minimumSd, a,
            i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.Y_SD)) {
          return true;
        }
        if (isAbove(maximumSd, a,
            i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.Y_SD)) {
          return true;
        }
      }
    }
    return false;
  }

  /**
   * Check if any dimension is below the threshold.
   *
   * @param threshold the threshold
   * @param params the params
   * @param paramIndex the param index
   * @return true, if is below
   */
  private static boolean isBelow(double[] threshold, double[] params, int paramIndex) {
    if (threshold != null) {
      for (int j = 0; j < 2; j++, paramIndex++) {
        if (params[paramIndex] < threshold[j]) {
          return true;
        }
      }
    }
    return false;
  }

  /**
   * Check if any dimension is above the threshold.
   *
   * @param threshold the threshold
   * @param params the params
   * @param paramIndex the param index
   * @return true, if is above
   */
  private static boolean isAbove(double[] threshold, double[] params, int paramIndex) {
    if (threshold != null) {
      for (int j = 0; j < 2; j++, paramIndex++) {
        if (params[paramIndex] > threshold[j]) {
          return true;
        }
      }
    }
    return false;
  }

  /**
   * Set the change in error that defines a negligible amount.
   *
   * @param delta the delta to set
   */
  public void setDelta(double delta) {
    this.delta = delta;
  }

  /**
   * Gets the delta.
   *
   * @return the delta.
   */
  public double getDelta() {
    return delta;
  }

  /**
   * Sets the minimum signal.
   *
   * @param minimumSignal the minimum signal
   */
  public void setMinimumSignal(double minimumSignal) {
    if (func.evaluatesSignal()) {
      this.minimumSignal = minimumSignal;
    }
  }

  /**
   * Gets the minimum signal.
   *
   * @return the minimum signal.
   */
  public double getMinimumSignal() {
    return minimumSignal;
  }

  /**
   * Sets the minimum position for each dimension.
   *
   * @param minimumPosition the minimum position for each dimension
   */
  public void setMinimumPosition(double[] minimumPosition) {
    if (func.evaluatesPosition()) {
      this.minimumPosition = checkArray(minimumPosition);
    }
  }

  /**
   * Gets the minimum position for each dimension.
   *
   * @return the minimum position for each dimension.
   */
  public double[] getMinimumPosition() {
    return minimumPosition;
  }

  /**
   * Sets the maximum position for each dimension.
   *
   * @param maximumPosition the maximum position for each dimension
   */
  public void setMaximumPosition(double[] maximumPosition) {
    if (func.evaluatesPosition()) {
      this.maximumPosition = checkArray(maximumPosition);
    }
  }

  /**
   * Gets the maximum position for each dimension.
   *
   * @return the maximum position for each dimension.
   */
  public double[] getMaximumPosition() {
    return maximumPosition;
  }

  /**
   * Sets the minimum standard deviation (SD) for each dimension.
   *
   * @param minimumSd the minimum SD for each dimension
   */
  public void setMinimumSd(double[] minimumSd) {
    if (func.evaluatesSD0()) {
      this.minimumSd = checkArray(minimumSd);
    }
  }

  /**
   * Gets the minimum standard deviation (SD) for each dimension.
   *
   * @return the minimum SD for each dimension.
   */
  public double[] getMinimumSd() {
    return minimumSd;
  }

  /**
   * Sets the maximum standard deviation (SD) for each dimension.
   *
   * @param maximumSd the maximum SD for each dimension
   */
  public void setMaximumSd(double[] maximumSd) {
    if (func.evaluatesSD0()) {
      this.maximumSd = checkArray(maximumSd);
    }
  }

  /**
   * Gets the maximum standard deviation (SD) for each dimension.
   *
   * @return the maximum SD for each dimension.
   */
  public double[] getMaximumSd() {
    return maximumSd;
  }

  private static double[] checkArray(double[] array) {
    return (ArrayUtils.getLength(array) == 2) ? array : null;
  }
}
