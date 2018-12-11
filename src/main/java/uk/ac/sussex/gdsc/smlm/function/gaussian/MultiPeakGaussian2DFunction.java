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
 * Abstract base class for an N-dimensional Gaussian function for a configured number of peaks.
 *
 * <p>The function will calculate the value of the Gaussian and evaluate the gradient of a set of
 * parameters. The class can specify which of the following parameters the function will
 * evaluate:<br> background, amplitude, angle[N-1], position[N], sd[N]
 *
 * <p>The class provides the number of peaks and the gradient indices.
 */
public abstract class MultiPeakGaussian2DFunction extends Gaussian2DFunction {
  /** The number of peaks. */
  protected final int npeaks;

  /** The gradient indices. */
  protected final int[] gradientIndices;

  /**
   * Instantiates a new multi peak gaussian 2D function.
   *
   * @param npeaks The number of peaks
   * @param maxx The maximum x value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   * @param maxy The maximum y value of the 2-dimensional data (used to unpack a linear index into
   *        coordinates)
   */
  public MultiPeakGaussian2DFunction(int npeaks, int maxx, int maxy) {
    super(maxx, maxy);
    this.npeaks = npeaks;
    this.gradientIndices = createGradientIndices(npeaks);
  }

  /** {@inheritDoc} */
  @Override
  public int getNPeaks() {
    return npeaks;
  }

  /** {@inheritDoc} */
  @Override
  public int[] gradientIndices() {
    return gradientIndices;
  }
}
