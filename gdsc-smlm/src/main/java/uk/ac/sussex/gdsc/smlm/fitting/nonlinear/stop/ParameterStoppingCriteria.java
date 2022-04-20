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

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;

/**
 * Defines the stopping criteria for the
 * {@link uk.ac.sussex.gdsc.smlm.fitting.nonlinear.NonLinearFit} class.
 *
 * <p>Stop when successive iterations with a reduced error move the fitted coordinates by less than
 * a specified distance.
 *
 * <p>The criteria also ensure that signal, coordinates and peak-widths are held positive, otherwise
 * fitting is stopped.
 */
public class ParameterStoppingCriteria extends GaussianStoppingCriteria {
  private int significantBits = 10;
  private double angleLimit = 1e-3f;

  private final DoubleEquality eq;

  /**
   * Instantiates a new parameter stopping criteria.
   *
   * @param func The Gaussian function
   */
  public ParameterStoppingCriteria(Gaussian2DFunction func) {
    super(func);
    eq = new DoubleEquality(DoubleEquality.getRelativeEpsilon(significantBits), 1e-16);
  }

  @Override
  protected StringBuilder logParameters(double oldError, double newError, double[] a) {
    final StringBuilder sb = new StringBuilder(158);
    sb.append("iter = ").append(getIteration() + 1).append(", error = ").append(oldError)
        .append(" -> ").append(newError);
    if (newError <= oldError) {
      if (func.evaluatesBackground()) {
        sb.append(", Back=[");
        sb.append(DoubleEquality.relativeError(bestA[0], a[0]));
        sb.append(']');
      }

      for (int i = 0; i < peaks; i++) {
        sb.append(", Peak").append(i + 1).append("=[");
        sb.append(DoubleEquality.relativeError(
            bestA[i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.SIGNAL],
            a[i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.SIGNAL]));
        sb.append(',');

        if (func.evaluatesAngle()) {
          final double x =
              bestA[i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.ANGLE];
          final double y = a[i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.ANGLE];
          sb.append(relativeAngle(x, y));
        } else {
          sb.append(0);
        }

        int param = i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_POSITION;
        for (int j = 0; j < 2 * 2; j++, param++) {
          sb.append(',');
          sb.append(DoubleEquality.relativeError(bestA[param], a[param]));
        }
        sb.append(']');
      }
    }
    return sb;
  }

  @Override
  protected boolean noCoordinateChange(double[] a) {
    // Old code does not correctly compute difference in angles. This is ignored for now.
    // return eq.almostEqualRelativeOrAbsolute(bestA, a);

    if (func.evaluatesBackground()
        && !eq.almostEqualRelativeOrAbsolute(bestA[Gaussian2DFunction.BACKGROUND],
            a[Gaussian2DFunction.BACKGROUND])) {
      return false;
    }

    for (int i = 0; i < peaks; i++) {
      if (!eq.almostEqualRelativeOrAbsolute(
          bestA[i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.SIGNAL],
          a[i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.SIGNAL])) {
        return false;
      }

      // Calculate the smallest angle between the two angles. This should be in the range 0 - 90
      // degrees.
      // Use this to compare if the angle has changed significantly relative to the maximum it could
      // change.
      if (func.evaluatesAngle()) {
        final double x =
            bestA[i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.ANGLE];
        final double y = a[i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.ANGLE];
        if (relativeAngle(x, y) > angleLimit) {
          return false;
        }
      }

      int param = i * Gaussian2DFunction.PARAMETERS_PER_PEAK + Gaussian2DFunction.X_POSITION;
      for (int j = 0; j < 2; j++, param++) {
        if (!eq.almostEqualRelativeOrAbsolute(bestA[param], a[param])) {
          return false;
        }
      }
    }

    return true;
  }

  private static double relativeAngle(double x, double y) {
    final double angle = Math.atan2(Math.sin(x - y), Math.cos(x - y));
    final double halfPi = Math.PI / 2;
    return Math.abs(angle / halfPi);
  }

  /**
   * Set the change in parameters that defines a negligible amount. Specified using the number of
   * binary significant digits (range [1, 52]).
   *
   * @param significantBits the significant bits to set
   */
  public void setSignificantBits(int significantBits) {
    this.significantBits = significantBits;
    eq.setMaxRelativeError(DoubleEquality.getRelativeEpsilon(significantBits));
    angleLimit = 1.0 / Math.pow(10, significantBits - 1);
  }

  /**
   * Gets the significant bits.
   *
   * @return the significant bits.
   */
  public int getSignificantBits() {
    return significantBits;
  }
}
