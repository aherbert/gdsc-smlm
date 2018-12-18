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

import uk.ac.sussex.gdsc.smlm.results.predicates.PeakResultPredicate;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;

import java.util.Collection;

/**
 * Stores peak results.
 */
public interface PeakResultStore {
  /**
   * Get the size.
   *
   * @return the size
   */
  int size();

  /**
   * Add a result. Not synchronized.
   *
   * @param result the result
   * @return true if the store is changed
   */
  boolean add(PeakResult result);

  /**
   * Add all results.
   *
   * @param results the results
   * @return true if the store is changed
   */
  boolean addCollection(Collection<PeakResult> results);

  /**
   * Add all results.
   *
   * @param results the results
   * @return true if the store is changed
   */
  boolean addArray(PeakResult[] results);

  /**
   * Adds the results.
   *
   * @param results the results
   * @return true if the store is changed
   */
  boolean addStore(PeakResultStore results);

  /**
   * Remove a result. Not synchronized.
   *
   * @param result the result
   * @return true if the store is changed
   */
  boolean remove(PeakResult result);

  /**
   * Remove all results.
   *
   * @param results the results
   * @return true if the store is changed
   */
  boolean removeCollection(Collection<PeakResult> results);

  /**
   * Remove all results.
   *
   * @param results the results
   * @return true if the store is changed
   */
  boolean removeArray(PeakResult[] results);

  /**
   * Removes the results.
   *
   * @param results the results
   * @return true if the store is changed
   */
  boolean removeStore(PeakResultStore results);

  /**
   * Retain all results.
   *
   * @param results the results
   * @return true if the store is changed
   */
  boolean retainCollection(Collection<PeakResult> results);

  /**
   * Retain all results.
   *
   * @param results the results
   * @return true if the store is changed
   */
  boolean retainArray(PeakResult[] results);

  /**
   * Retains the results.
   *
   * @param results the results
   * @return true if the store is changed
   */
  boolean retainStore(PeakResultStore results);

  /**
   * Clear the results.
   */
  void clear();

  /**
   * Trims the capacity of this instance to be the current size. An application can use this
   * operation to minimize the storage of an instance.
   */
  void trimToSize();

  /**
   * Convert to an array. This is a new allocation of storage space.
   *
   * @return the peak result array
   */
  PeakResult[] toArray();

  /**
   * Copy the results.
   *
   * @return the copy
   */
  PeakResultStore copy();

  /**
   * Copy the results.
   *
   * @param deepCopy Set to true to perform a deep copy
   * @return the copy
   */
  PeakResultStore copy(boolean deepCopy);

  /**
   * Removes the result if it matches the filter. If objects are removed then the order of elements
   * may change.
   *
   * @param filter the filter
   * @return true, if any were removed
   */
  boolean removeIf(PeakResultPredicate filter);

  /**
   * Execute the procedure on each result in the store.
   *
   * @param procedure the procedure
   */
  void forEach(PeakResultProcedure procedure);

  /**
   * Get a subset of the results if they match the filter.
   *
   * @param filter the filter
   * @return the results
   */
  PeakResult[] subset(PeakResultPredicate filter);

  /**
   * Returns <tt>true</tt> if this store contains the specified result. More formally, returns
   * <tt>true</tt> if and only if this store contains at least one element <tt>e</tt> such that
   * <tt>(result==null&nbsp;?&nbsp;e==null&nbsp;:&nbsp;result.equals(e))</tt>.
   *
   * @param result the result
   * @return <tt>true</tt> if this list contains the specified result true, if successful
   */
  boolean contains(PeakResult result);
}
