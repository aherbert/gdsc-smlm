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

package uk.ac.sussex.gdsc.smlm.results;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.function.Predicate;
import org.apache.commons.rng.UniformRandomProvider;
import uk.ac.sussex.gdsc.core.utils.rng.JdkRandomAdaptor;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;

/**
 * Stores peak results using an ArrayList.
 */
public class ArrayListPeakResultStore implements PeakResultStoreList, PeakResultStoreCollection {
  /** The results. */
  private final ArrayList<PeakResult> results;

  /**
   * Instantiates a new array list peak results store.
   *
   * @param capacity the capacity
   */
  public ArrayListPeakResultStore(int capacity) {
    this.results = new ArrayList<>(capacity);
  }

  /**
   * Instantiates a new array list peak result store.
   *
   * @param store the store to copy
   */
  public ArrayListPeakResultStore(ArrayListPeakResultStore store) {
    this.results = new ArrayList<>(store.results);
  }

  @Override
  public PeakResult get(int index) {
    return results.get(index);
  }

  @Override
  public int size() {
    return results.size();
  }

  @Override
  public boolean add(PeakResult result) {
    return results.add(result);
  }

  @Override
  public boolean addCollection(Collection<PeakResult> results) {
    return this.results.addAll(results);
  }

  @Override
  public boolean addArray(PeakResult[] results) {
    return this.results.addAll(Arrays.asList(results));
  }

  @Override
  public boolean addStore(PeakResultStore results) {
    if (results instanceof PeakResultStoreCollection) {
      return this.results.addAll(((PeakResultStoreCollection) results).getCollectionReference());
    }
    return addArray(results.toArray());
  }

  @Override
  public PeakResult remove(int index) {
    return results.remove(index);
  }

  @Override
  public void remove(int fromIndex, int toIndex) {
    if (fromIndex > toIndex) {
      throw new IllegalArgumentException("fromIndex must be <= toIndex");
    }
    for (int i = toIndex; i >= fromIndex; i--) {
      results.remove(i);
    }
  }

  @Override
  public boolean remove(PeakResult result) {
    return results.remove(result);
  }

  @Override
  public boolean removeCollection(Collection<PeakResult> results) {
    return this.results.removeAll(results);
  }

  @Override
  public boolean removeArray(PeakResult[] results) {
    return this.results.removeAll(Arrays.asList(results));
  }

  @Override
  public boolean removeStore(PeakResultStore results) {
    if (results instanceof PeakResultStoreCollection) {
      return this.results.removeAll(((PeakResultStoreCollection) results).getCollectionReference());
    }
    return removeArray(results.toArray());
  }

  @Override
  public boolean retainCollection(Collection<PeakResult> results) {
    return this.results.retainAll(results);
  }

  @Override
  public boolean retainArray(PeakResult[] results) {
    return this.results.retainAll(Arrays.asList(results));
  }

  @Override
  public boolean retainStore(PeakResultStore results) {
    if (results instanceof PeakResultStoreCollection) {
      return this.results.retainAll(((PeakResultStoreCollection) results).getCollectionReference());
    }
    return retainArray(results.toArray());
  }

  @Override
  public void clear() {
    results.clear();
  }

  @Override
  public void trimToSize() {
    results.trimToSize();
  }

  @Override
  public void sort(Comparator<PeakResult> comparator) {
    Collections.sort(results, comparator);
  }

  @Override
  public PeakResult[] toArray() {
    return results.toArray(new PeakResult[0]);
  }

  @Override
  public PeakResultStore copy() {
    return new ArrayListPeakResultStore(this);
  }

  @Override
  public PeakResultStore copy(boolean deepCopy) {
    if (deepCopy) {
      final ArrayListPeakResultStore copy = new ArrayListPeakResultStore(size());
      for (int i = 0, size = size(); i < size; i++) {
        copy.add(results.get(i).copy());
      }
      return copy;
    }
    return copy();
  }

  @Override
  public boolean removeIf(final Predicate<PeakResult> filter) {
    return results.removeIf(filter);
  }

  @Override
  public void forEach(PeakResultProcedure procedure) {
    for (int i = 0, size = size(); i < size; i++) {
      procedure.execute(results.get(i));
    }
  }

  @Override
  public PeakResult[] subset(Predicate<PeakResult> filter) {
    final ArrayPeakResultStore list = new ArrayPeakResultStore(10);
    for (int i = 0, size = size(); i < size; i++) {
      if (filter.test(results.get(i))) {
        list.add(results.get(i));
      }
    }
    return list.toArray();
  }

  @Override
  public void shuffle(UniformRandomProvider randomSource) {
    Collections.shuffle(results, new JdkRandomAdaptor(randomSource));
  }

  @Override
  public int indexOf(PeakResult result) {
    return results.indexOf(result);
  }

  @Override
  public int lastIndexOf(PeakResult result) {
    return results.lastIndexOf(result);
  }

  @Override
  public boolean contains(PeakResult result) {
    return results.contains(result);
  }

  @Override
  @SuppressWarnings("unchecked")
  public Collection<PeakResult> getCollection() {
    return (Collection<PeakResult>) results.clone();
  }

  @Override
  public Collection<PeakResult> getCollectionReference() {
    return results;
  }
}
