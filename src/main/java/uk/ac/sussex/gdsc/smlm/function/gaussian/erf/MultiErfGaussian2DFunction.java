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

package uk.ac.sussex.gdsc.smlm.function.gaussian.erf;

import uk.ac.sussex.gdsc.smlm.function.ValueProcedure;

/**
 * Abstract base class for an 2-dimensional Gaussian function for a configured number of peaks.
 *
 * <p>The function will calculate the value of the Gaussian and evaluate the gradient of a set of
 * parameters. The class can specify which of the following parameters the function will
 * evaluate:<br> background, signal, z-depth, position0, position1, sd0, sd1
 *
 * <p>The class provides an index of the position in the parameter array where the parameter is
 * expected.
 */
public abstract class MultiErfGaussian2DFunction extends ErfGaussian2DFunction {
  /** The number of peaks. */
  protected final int numberOfPeaks;

  /** The gradient indices. */
  protected final int[] gradientIndices;

  /** The target intensity for each peak. */
  protected final double[] tI;

  /**
   * Instantiates a new multi-peak erf gaussian 2D function.
   *
   * @param numberOfPeaks The number of peaks
   * @param maxx The maximum x value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   * @param maxy The maximum y value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   */
  public MultiErfGaussian2DFunction(int numberOfPeaks, int maxx, int maxy) {
    super(numberOfPeaks, maxx, maxy);
    this.numberOfPeaks = numberOfPeaks;
    this.gradientIndices = createGradientIndices();
    tI = new double[numberOfPeaks];
  }

  @Override
  public int getNPeaks() {
    return numberOfPeaks;
  }

  /**
   * Creates the gradient indices.
   *
   * @return the the gradient indices
   */
  protected abstract int[] createGradientIndices();

  /**
   * Replicate the gradient indices from multiple peaks for the configured number of peaks.
   *
   * @param singleGradientIndices the single gradient indices
   * @return the multi gradient indices
   */
  protected int[] replicateGradientIndices(int[] singleGradientIndices) {
    final int start = (evaluatesBackground() ? 1 : 0);
    final int m = singleGradientIndices.length;
    final int[] indices = new int[start + numberOfPeaks * (m - start)];
    int p = 0;
    if (evaluatesBackground()) {
      indices[p++] = 0;
    }
    for (int n = 0, i = 0; n < numberOfPeaks; n++, i += PARAMETERS_PER_PEAK) {
      for (int j = start; j < m; j++) {
        indices[p++] = i + singleGradientIndices[j];
      }
    }
    return indices;
  }

  @Override
  public int[] gradientIndices() {
    return gradientIndices;
  }

  @Override
  public int getNumberOfGradients() {
    return gradientIndices.length;
  }

  /**
   * Evaluates an 2-dimensional Gaussian function for multiple peaks.
   *
   * @param i Input predictor
   * @return The Gaussian value
   */
  @Override
  public double eval(final int i) {
    // Unpack the predictor into the dimensions
    int yy = i / maxx;
    int xx = i % maxx;

    double I = tB;
    for (int n = 0; n < numberOfPeaks; n++, xx += maxx, yy += maxy) {
      I += tI[n] * deltaEx[xx] * deltaEy[yy];
    }
    return I;
  }

  /**
   * Evaluates an 2-dimensional Gaussian function for multiple peaks.
   *
   * @param i Input predictor
   * @param duda Partial gradient of function with respect to each coefficient
   * @return The predicted value
   */
  @Override
  public abstract double eval(final int i, final double[] duda);

  /**
   * Evaluates an 2-dimensional Gaussian function for multiple peaks.
   *
   * @param i Input predictor
   * @param duda Partial first gradient of function with respect to each coefficient
   * @param d2uda2 Partial second gradient of function with respect to each coefficient
   * @return The predicted value
   */
  @Override
  public abstract double eval(final int i, final double[] duda, final double[] d2uda2);

  @Override
  public void forEach(ValueProcedure procedure) {
    // Unroll for the number of peaks
    if (numberOfPeaks == 2) {
      if (tB == 0) {
        for (int y = 0; y < maxy; y++) {
          // Pre-compute
          final double tI_deltaEy0 = tI[0] * deltaEy[y];
          final double tI_deltaEy1 = tI[1] * deltaEy[y + maxy];

          for (int x = 0; x < maxx; x++) {
            procedure.execute(tI_deltaEy0 * deltaEx[x] + tI_deltaEy1 * deltaEx[x + maxx]);
          }
        }
      } else {
        for (int y = 0; y < maxy; y++) {
          // Pre-compute
          final double tI_deltaEy0 = tI[0] * deltaEy[y];
          final double tI_deltaEy1 = tI[1] * deltaEy[y + maxy];

          for (int x = 0; x < maxx; x++) {
            procedure.execute(tB + tI_deltaEy0 * deltaEx[x] + tI_deltaEy1 * deltaEx[x + maxx]);
          }
        }
      }
    } else {
      final double[] tI_deltaEy = new double[numberOfPeaks];
      for (int y = 0; y < maxy; y++) {
        // Pre-compute
        for (int n = 0, yy = y; n < numberOfPeaks; n++, yy += maxy) {
          tI_deltaEy[n] = tI[n] * deltaEy[yy];
        }

        for (int x = 0; x < maxx; x++) {
          double I = tB;
          for (int n = 0, xx = x; n < numberOfPeaks; n++, xx += maxx) {
            I += tI_deltaEy[n] * deltaEx[xx];
          }
          procedure.execute(I);
        }

        // No pre-compute
        // for (int x = 0; x < maxx; x++)
        // {
        // double I = tB;
        // for (int n = 0, xx = x, yy = y; n < numberOfPeaks; n++, xx += maxx, yy += maxy)
        // I += tI[n] * deltaEx[xx] * deltaEy[yy];
        // procedure.execute(I);
        // }
      }
    }
  }

  @Override
  public double[] computeValues(double[] variables) {
    initialise0(variables);
    final double[] values = new double[size()];
    // Unroll for the number of peaks
    if (numberOfPeaks == 2) {
      if (tB == 0) {
        for (int y = 0, i = 0; y < maxy; y++) {
          // Pre-compute
          final double tI_deltaEy0 = tI[0] * deltaEy[y];
          final double tI_deltaEy1 = tI[1] * deltaEy[y + maxy];

          for (int x = 0; x < maxx; x++) {
            values[i++] = (tI_deltaEy0 * deltaEx[x] + tI_deltaEy1 * deltaEx[x + maxx]);
          }
        }
      } else {
        for (int y = 0, i = 0; y < maxy; y++) {
          // Pre-compute
          final double tI_deltaEy0 = tI[0] * deltaEy[y];
          final double tI_deltaEy1 = tI[1] * deltaEy[y + maxy];

          for (int x = 0; x < maxx; x++) {
            values[i++] = (tB + tI_deltaEy0 * deltaEx[x] + tI_deltaEy1 * deltaEx[x + maxx]);
          }
        }
      }
    } else {
      final double[] tI_deltaEy = new double[numberOfPeaks];
      for (int y = 0, i = 0; y < maxy; y++) {
        // Pre-compute
        for (int n = 0, yy = y; n < numberOfPeaks; n++, yy += maxy) {
          tI_deltaEy[n] = tI[n] * deltaEy[yy];
        }

        for (int x = 0; x < maxx; x++) {
          double I = tB;
          for (int n = 0, xx = x; n < numberOfPeaks; n++, xx += maxx) {
            I += tI_deltaEy[n] * deltaEx[xx];
          }
          values[i++] = I;
        }
      }
    }
    return values;
  }

  // Force implementation
  @Override
  public abstract double integral(double[] a);
}
