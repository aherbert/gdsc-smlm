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

package uk.ac.sussex.gdsc.smlm.filters;

/**
 * Contains common functionality for weighted filters.
 */
public abstract class BaseWeightedFilter {
  /** The weights. */
  protected float[] weights;

  /** The width of the weights. */
  protected int weightWidth;

  /** The height of the weights. */
  protected int weightHeight;

  /**
   * Instantiates a new base weighted filter.
   */
  protected BaseWeightedFilter() {
    // Do nothing
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected BaseWeightedFilter(BaseWeightedFilter source) {
    // Should be thread safe
    this.weights = source.weights;
    // this.weights = ArrayUtils.clone(source.weights)
    this.weightWidth = source.weightWidth;
    this.weightHeight = source.weightHeight;
  }

  /**
   * Sets the weights of the data. This should be called before filtering data samples.
   *
   * @param weights the weights of the data (can be null)
   * @param width The width of the data
   * @param height The height of the data
   */
  public void setWeights(final float[] weights, final int width, final int height) {
    this.weights = weights;
    this.weightWidth = width;
    this.weightHeight = height;
    newWeights();
  }

  /**
   * Checks for weights.
   *
   * @return true, if successful
   */
  public boolean hasWeights() {
    return weights != null;
  }

  /**
   * Signal that new weight parameters have been set. Sub-classes can re-initialise for the new
   * weights.
   */
  protected abstract void newWeights();
}
