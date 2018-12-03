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

import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;

/**
 * Stores peak results using a set. This is similar to an HashSet but does not have concurrency
 * checking.
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

  /** {@inheritDoc} */
  @Override
  public int size() {
    return results.size();
  }

  /** {@inheritDoc} */
  @Override
  public boolean add(PeakResult result) {
    return results.add(result);
  }

  /** {@inheritDoc} */
  @Override
  public boolean addCollection(Collection<PeakResult> results) {
    return this.results.addAll(results);
  }

  /** {@inheritDoc} */
  @Override
  public boolean addArray(PeakResult[] results) {
    return this.results.addAll(Arrays.asList(results));
  }

  /** {@inheritDoc} */
  @Override
  public boolean addStore(PeakResultStore results) {
    if (results instanceof PeakResultStoreCollection) {
      return this.results.addAll(((PeakResultStoreCollection) results).getCollectionReference());
    }
    return addArray(results.toArray());
  }

  /** {@inheritDoc} */
  @Override
  public boolean remove(PeakResult result) {
    return results.remove(result);
  }

  /** {@inheritDoc} */
  @Override
  public boolean removeCollection(Collection<PeakResult> results) {
    return this.results.removeAll(results);
  }

  /** {@inheritDoc} */
  @Override
  public boolean removeArray(PeakResult[] results) {
    return this.results.removeAll(Arrays.asList(results));
  }

  /** {@inheritDoc} */
  @Override
  public boolean removeStore(PeakResultStore results) {
    if (results instanceof PeakResultStoreCollection) {
      return this.results.removeAll(((PeakResultStoreCollection) results).getCollectionReference());
    }
    return removeArray(results.toArray());
  }

  /** {@inheritDoc} */
  @Override
  public boolean retainCollection(Collection<PeakResult> results) {
    return this.results.retainAll(results);
  }

  /** {@inheritDoc} */
  @Override
  public boolean retainArray(PeakResult[] results) {
    return this.results.retainAll(Arrays.asList(results));
  }

  /** {@inheritDoc} */
  @Override
  public boolean retainStore(PeakResultStore results) {
    if (results instanceof PeakResultStoreCollection) {
      return this.results.retainAll(((PeakResultStoreCollection) results).getCollectionReference());
    }
    return retainArray(results.toArray());
  }

  /** {@inheritDoc} */
  @Override
  public void clear() {
    results.clear();
  }

  /** {@inheritDoc} */
  @Override
  public void trimToSize() {
    // results.trimToSize();
  }

  /** {@inheritDoc} */
  @Override
  public PeakResult[] toArray() {
    return results.toArray(new PeakResult[size()]);
  }

  /** {@inheritDoc} */
  @Override
  public PeakResultStore copy() {
    return new SetPeakResultStore(this);
  }

  /** {@inheritDoc} */
  @Override
  public PeakResultStore copy(boolean deepCopy) {
    if (deepCopy) {
      final SetPeakResultStore copy = new SetPeakResultStore(size());
      for (final PeakResult r : results) {
        copy.add(r.clone());
      }
      return copy;
    }
    return copy();
  }

  /** {@inheritDoc} */
  @Override
  public boolean removeIf(final PeakResultPredicate filter) {
    // Delegate to the list implementation
    final ArrayPeakResultStore list = new ArrayPeakResultStore(10);
    for (final PeakResult r : results) {
      if (filter.test(r)) {
        list.add(r);
      }
    }
    return this.results.removeAll(Arrays.asList(list.toArray()));
  }

  /** {@inheritDoc} */
  @Override
  public void forEach(PeakResultProcedure procedure) {
    for (final PeakResult r : results) {
      procedure.execute(r);
    }
  }

  /** {@inheritDoc} */
  @Override
  public PeakResult[] subset(PeakResultPredicate filter) {
    final ArrayPeakResultStore list = new ArrayPeakResultStore(10);
    for (final PeakResult r : results) {
      if (filter.test(r)) {
        list.add(r);
      }
    }
    return list.toArray();
  }

  /** {@inheritDoc} */
  @Override
  public boolean contains(PeakResult result) {
    return results.contains(result);
  }

  /** {@inheritDoc} */
  @Override
  @SuppressWarnings("unchecked")
  public Collection<PeakResult> getCollection() {
    return (Collection<PeakResult>) results.clone();
  }

  /** {@inheritDoc} */
  @Override
  public Collection<PeakResult> getCollectionReference() {
    return results;
  }
}
