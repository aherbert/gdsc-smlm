/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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

@SuppressWarnings({"javadoc"})
public abstract class IntDataFilter {
  /** The name. */
  final String name;

  /** Flag to indicate that the filter can support non-integer box sizes. */
  final boolean isInterpolated;

  /** The min box size. */
  final int minBoxSize;

  /**
   * Instantiates a new data filter.
   *
   * @param name the name
   * @param isInterpolated Flag to indicate that the filter can support non-integer box sizes
   */
  public IntDataFilter(String name, boolean isInterpolated) {
    this(name, isInterpolated, 0);
  }

  /**
   * Instantiates a new data filter.
   *
   * @param name the name
   * @param isInterpolated Flag to indicate that the filter can support non-integer box sizes
   * @param minBoxSize the min box size
   */
  public IntDataFilter(String name, boolean isInterpolated, int minBoxSize) {
    this.name = name;
    this.isInterpolated = isInterpolated;
    this.minBoxSize = minBoxSize;
  }

  public abstract void filter(int[] data, int width, int height, int boxSize);

  public abstract void filterInternal(int[] data, int width, int height, int boxSize);

  public abstract void setWeights(float[] weights, int width, int height);
}
