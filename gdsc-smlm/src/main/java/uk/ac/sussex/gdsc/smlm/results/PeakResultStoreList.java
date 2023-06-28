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

package uk.ac.sussex.gdsc.smlm.results;

import java.util.Comparator;
import org.apache.commons.rng.UniformRandomProvider;
import uk.ac.sussex.gdsc.smlm.results.sort.FrameIdPeakResultComparator;

/**
 * Stores peak results with list access.
 */
public interface PeakResultStoreList extends PeakResultStore {
  /**
   * Gets the result.
   *
   * @param index the index
   * @return the peak result
   */
  PeakResult get(int index);

  /**
   * Removes the result.
   *
   * @param index the index
   * @return the peak result removed
   * @throws IndexOutOfBoundsException If the index is invalid
   */
  PeakResult remove(int index);

  /**
   * Removes a range of results.
   *
   * @param fromIndex the from index
   * @param toIndex the to index (inclusive)
   * @throws IndexOutOfBoundsException If the index is invalid
   */
  void remove(int fromIndex, int toIndex);

  /**
   * Sort the results.
   */
  default void sort() {
    sort(FrameIdPeakResultComparator.INSTANCE);
  }

  /**
   * Sort the results.
   *
   * @param comparator the comparator
   */
  void sort(Comparator<PeakResult> comparator);

  /**
   * Shuffle the results.
   *
   * @param randomSource the random source
   */
  void shuffle(UniformRandomProvider randomSource);

  /**
   * Returns the index of the first occurrence of the specified result in this store, or -1 if this
   * list does not contain the element. More formally, returns the lowest index {@code i} such
   * that {@code (result==null ? get(i)==null : result.equals(get(i)))}, or -1
   * if there is no such index.
   *
   * @param result the result
   * @return the index (or -1)
   */
  int indexOf(PeakResult result);

  /**
   * Returns the index of the last occurrence of the specified result in this store, or -1 if this
   * list does not contain the element. More formally, returns the lowest index {@code i} such
   * that {@code (result==null ? get(i)==null : result.equals(get(i)))}, or -1
   * if there is no such index.
   *
   * @param result the result
   * @return the index (or -1)
   */
  int lastIndexOf(PeakResult result);
}
