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

package uk.ac.sussex.gdsc.smlm.filters;

/**
 * Computes the mean using a circular mask.
 *
 * <p>Adapted from ij.plugin.filter.RankFilters
 */
public class CircularMeanFilter extends CircularFilter {
  /** {@inheritDoc} */
  @Override
  protected Normaliser computeWeightedNormaliser(double radius) {
    final float[] npoints = weights.clone();
    final CircularSumFilter sum = new CircularSumFilter();
    sum.convolve(npoints, weightWidth, weightHeight, radius);
    return new PerPixelNormaliser(npoints);
  }

  /** {@inheritDoc} */
  @Override
  protected Normaliser computeNormaliser(int npoints) {
    return new FixedNormaliser(npoints);
  }

  /** {@inheritDoc} */
  @Override
  public CircularMeanFilter clone() {
    final CircularMeanFilter o = (CircularMeanFilter) super.clone();
    return o;
  }
}
