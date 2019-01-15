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

package uk.ac.sussex.gdsc.smlm.filters;

/**
 * Computes the sum using a circular mask.
 *
 * <p>Adapted from ij.plugin.filter.RankFilters
 */
public class CircularSumFilter extends CircularFilter {

  /**
   * Instantiates a new circular sum filter.
   */
  public CircularSumFilter() {
    // Do nothing
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected CircularSumFilter(CircularSumFilter source) {
    super(source);
  }

  /**
   * Create a copy.
   *
   * @return the copy
   */
  public CircularSumFilter copy() {
    return new CircularSumFilter(this);
  }

  @Override
  protected Normaliser computeWeightedNormaliser(double radius) {
    return NonNormaliser.INSTANCE;
  }

  @Override
  protected Normaliser computeNormaliser(int npoints) {
    return NonNormaliser.INSTANCE;
  }
}
