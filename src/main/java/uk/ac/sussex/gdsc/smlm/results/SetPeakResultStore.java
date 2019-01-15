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

package uk.ac.sussex.gdsc.smlm.results;

import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.function.Predicate;

/**
 * Stores peak results using a set.
 *
 * <p>Note that the {@link PeakResult} object does not implement {@link Object#hashCode()} or
 * {@link Object#equals(Object)} and so this stores all unique result references.
 */
public class SetPeakResultStore implements PeakResultStore, PeakResultStoreCollection {
  /** The results. */
  private final HashSet<PeakResult> results;

  /**
   * Instantiates a new set peak results store.
   *
   * @param capacity the capacity
   */
  public SetPeakResultStore(int capacity) {
    this.results = new HashSet<>(capacity);
  }

  /**
   * Instantiates a new array list peak result store.
   *
   * @param store the store to copy
   */
  public SetPeakResultStore(SetPeakResultStore store) {
    this.results = new HashSet<>(store.results);
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
  public void trimToSize() {}

  @Override
  public PeakResult[] toArray() {
    return results.toArray(new PeakResult[size()]);
  }

  @Override
  public PeakResultStore copy() {
    return new SetPeakResultStore(this);
  }

  @Override
  public PeakResultStore copy(boolean deepCopy) {
    if (deepCopy) {
      final SetPeakResultStore copy = new SetPeakResultStore(size());
      for (final PeakResult r : results) {
        copy.add(r.copy());
      }
      return copy;
    }
    return copy();
  }

  @Override
  public boolean removeIf(final Predicate<PeakResult> filter) {
    // Delegate to the list implementation
    final ArrayPeakResultStore list = new ArrayPeakResultStore(10);
    for (final PeakResult r : results) {
      if (filter.test(r)) {
        list.add(r);
      }
    }
    return this.results.removeAll(Arrays.asList(list.toArray()));
  }

  @Override
  public void forEach(PeakResultProcedure procedure) {
    for (final PeakResult r : results) {
      procedure.execute(r);
    }
  }

  @Override
  public PeakResult[] subset(Predicate<PeakResult> filter) {
    final ArrayPeakResultStore list = new ArrayPeakResultStore(10);
    for (final PeakResult r : results) {
      if (filter.test(r)) {
        list.add(r);
      }
    }
    return list.toArray();
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
