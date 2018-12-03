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
package uk.ac.sussex.gdsc.smlm.results;

/**
 * Stores peak results using a fixed size array.
 */
public class FixedPeakResultList {
  /** The results. */
  public PeakResult[] results;

  /** The size. */
  public int size;

  /**
   * Instantiates a new fixed peak result list.
   *
   * @param capacity the capacity
   */
  public FixedPeakResultList(int capacity) {
    this.results = new PeakResult[Math.max(capacity, 0)];
    this.size = 0;
  }

  /**
   * Adds the result.
   *
   * @param result the result
   */
  public void add(PeakResult result) {
    results[size++] = result;
  }

  /**
   * Get the result. <p> Note: This does not check against the current size so can return stale
   * data.
   *
   * @param index the index
   * @return the peak result
   */
  public PeakResult get(int index) {
    return results[index];
  }

  /**
   * Get the size.
   *
   * @return the size
   */
  public int size() {
    return size;
  }

  /**
   * Checks if is empty.
   *
   * @return true, if is empty
   */
  public boolean isEmpty() {
    return size == 0;
  }

  /**
   * Checks if is not empty.
   *
   * @return true, if is not empty
   */
  public boolean isNotEmpty() {
    return size != 0;
  }

  /**
   * Clear the list <p> Note: This does not remove the references to the underlying data or
   * reallocate storage thus {@link #get(int)} can return stale data.
   */
  public void clear() {
    size = 0;
  }

  /**
   * Convert to a new array using the current size.
   *
   * @return the peak result array
   */
  public PeakResult[] toArray() {
    final PeakResult[] array = new PeakResult[size];
    System.arraycopy(results, 0, array, 0, size);
    return array;
  }
}
