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

package uk.ac.sussex.gdsc.smlm.search;

/**
 * Calculate the score of points within a search space.
 *
 * @param <T> the generic type
 */
public interface FullScoreFunction<T extends Comparable<T>> extends ScoreFunction<T> {
  /**
   * Return the score of the input points.
   *
   * @param points the points
   * @return the scores
   */
  SearchResult<T>[] score(double[][] points);

  /**
   * Cut the list of scores down to the given size by selecting only the best results. The input
   * list may not be sorted. The results should contain the best result at position 0 in the output
   * array.
   *
   * @param scores The scores
   * @param size The size
   * @return The reduced list
   */
  SearchResult<T>[] cut(SearchResult<T>[] scores, int size);
}
