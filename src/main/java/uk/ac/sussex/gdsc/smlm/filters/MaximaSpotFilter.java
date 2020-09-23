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

import java.util.ArrayList;
import java.util.List;
import uk.ac.sussex.gdsc.core.filters.NonMaximumSuppression;

/**
 * Identifies candidate spots (local maxima) in an image using non-maximum suppression.
 */
public abstract class MaximaSpotFilter extends SpotFilter {
  private final int search;
  private final int border;
  private final NonMaximumSuppression nms;
  private float[] data2;

  /**
   * Create the spot filter.
   *
   * @param search The search width for non-maximum suppression
   * @param border The border to ignore for maxima
   * @throws IllegalArgumentException if search is below 1 or border is below zero
   */
  public MaximaSpotFilter(int search, int border) {
    if (search < 1) {
      throw new IllegalArgumentException("Search width must be 1 or above");
    }
    if (border < 0) {
      throw new IllegalArgumentException("Border must be 0 or above");
    }
    this.search = search;
    this.border = border;
    nms = new NonMaximumSuppression();
    // Do a neighbour check when using a low block size
    nms.setNeighbourCheck(search < 3);
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected MaximaSpotFilter(MaximaSpotFilter source) {
    search = source.search;
    border = source.border;
    nms = source.nms.copy();
  }

  @Override
  protected Spot[] find(final float[] data, final int width, final int height) {
    data2 = preprocessData(data, width, height);

    final int[] maxIndices = getMaxima(data2, width, height);
    if (maxIndices.length == 0) {
      return null;
    }

    final Spot[] spots = new Spot[maxIndices.length];
    for (int n = 0; n < maxIndices.length; n++) {
      final int y = maxIndices[n] / width;
      final int x = maxIndices[n] % width;
      final float intensity = data2[maxIndices[n]];
      spots[n] = new Spot(x, y, intensity);
    }
    return spots;
  }

  @Override
  public float[] getPreprocessedData() {
    return data2;
  }

  /**
   * Pre-process the data before finding local maxima.
   *
   * @param data the data
   * @param width the width
   * @param height the height
   * @return The pre-processed data
   */
  public abstract float[] preprocessData(final float[] data, final int width, final int height);

  /**
   * Find the indices of the maxima using the currently configured parameters.
   *
   * <p>Data must be arranged in yx block order, i.e. height rows of width.
   *
   * @param data the data
   * @param width the width
   * @param height the height
   * @return Indices of the maxima
   */
  protected int[] getMaxima(float[] data, int width, int height) {
    // Check upper limits are safe
    final int n = Math.min(search, Math.min(width, height));
    final int validBorder = Math.min(this.border, Math.min(width, height) / 2);
    return nms.blockFindInternal(data, width, height, n, validBorder);
  }

  /**
   * Gets the search width for maxima (maximum must be the highest point in a 2n+1 region).
   *
   * @return the search width for maxima
   */
  public int getSearch() {
    return search;
  }

  /**
   * Gets the border at the edge to ignore for maxima.
   *
   * @return the border.
   */
  public int getBorder() {
    return border;
  }

  @Override
  public List<String> getParameters() {
    final ArrayList<String> list = new ArrayList<>();
    list.add("search = " + search);
    list.add("border = " + border);
    return list;
  }
}
