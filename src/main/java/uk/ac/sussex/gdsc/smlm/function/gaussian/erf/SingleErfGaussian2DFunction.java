/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.function.gaussian.erf;

import uk.ac.sussex.gdsc.smlm.function.ValueProcedure;

/**
 * Abstract base class for a 2-dimensional Gaussian function for a configured number of peaks.
 *
 * <p>The function will calculate the value of the Gaussian and evaluate the gradient of a set of
 * parameters. The class can specify which of the following parameters the function will
 * evaluate:<br> background, signal, z-depth, position0, position1, sd0, sd1
 *
 * <p>The class provides an index of the position in the parameter array where the parameter is
 * expected.
 */
public abstract class SingleErfGaussian2DFunction extends ErfGaussian2DFunction {
  // Required for the PSF

  /** The intensity. */
  // CHECKSTYLE.OFF: MemberName
  protected double tI;
  // CHECKSTYLE.ON: MemberName

  /**
   * Instantiates a new erf gaussian 2D function.
   *
   * @param maxx The maximum x value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   * @param maxy The maximum y value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   */
  public SingleErfGaussian2DFunction(int maxx, int maxy) {
    super(1, maxx, maxy);
  }

  @Override
  public int getNPeaks() {
    return 1;
  }

  /**
   * Evaluates a 2-dimensional Gaussian function for a single peak.
   *
   * @param x Input predictor
   * @return The Gaussian value
   */
  @Override
  public double eval(final int x) {
    // Unpack the predictor into the dimensions
    final int yy = x / maxx;
    final int xx = x % maxx;

    return tb + tI * deltaEx[xx] * deltaEy[yy];
  }

  /**
   * Evaluates a 2-dimensional Gaussian function for a single peak.
   *
   * @param x Input predictor
   * @param duda Partial gradient of function with respect to each coefficient
   * @return The predicted value
   */
  @Override
  public abstract double eval(final int x, final double[] duda);

  /**
   * Evaluates a 2-dimensional Gaussian function for a single peak.
   *
   * @param x Input predictor
   * @param duda Partial first gradient of function with respect to each coefficient
   * @param d2uda2 Partial second gradient of function with respect to each coefficient
   * @return The predicted value
   */
  @Override
  public abstract double eval2(final int x, final double[] duda, final double[] d2uda2);

  @Override
  public void forEach(ValueProcedure procedure) {
    if (tb == 0) {
      for (int y = 0; y < maxy; y++) {
        final double tI_deltaEy = tI * deltaEy[y];
        for (int x = 0; x < maxx; x++) {
          procedure.execute(tI_deltaEy * deltaEx[x]);
        }
      }
    } else {
      for (int y = 0; y < maxy; y++) {
        final double tI_deltaEy = tI * deltaEy[y];
        for (int x = 0; x < maxx; x++) {
          procedure.execute(tb + tI_deltaEy * deltaEx[x]);
        }
      }
    }
  }

  @Override
  public double[] computeValues(double[] variables) {
    initialise0(variables);
    final double[] values = new double[size()];
    if (tb == 0) {
      for (int y = 0, i = 0; y < maxy; y++) {
        final double tI_deltaEy = tI * deltaEy[y];
        for (int x = 0; x < maxx; x++) {
          values[i++] = tI_deltaEy * deltaEx[x];
        }
      }
    } else {
      for (int y = 0, i = 0; y < maxy; y++) {
        final double tI_deltaEy = tI * deltaEy[y];
        for (int x = 0; x < maxx; x++) {
          values[i++] = tb + tI_deltaEy * deltaEx[x];
        }
      }
    }
    return values;
  }

  // Force implementation
  @Override
  public abstract int getNumberOfGradients();

  // Force implementation
  @Override
  public abstract double integral(double[] a);
}
