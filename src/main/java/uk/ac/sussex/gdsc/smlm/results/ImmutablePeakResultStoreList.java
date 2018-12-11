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

import uk.ac.sussex.gdsc.core.data.DataException;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.rng.UniformRandomProvider;

import java.util.Comparator;

/**
 * Stores peak results and prevents modification.
 */
public class ImmutablePeakResultStoreList extends ImmutablePeakResultStore
    implements PeakResultStoreList {
  private final PeakResultStoreList store;

  /**
   * Instantiates a new immutable peak result store.
   *
   * @param store the store
   */
  public ImmutablePeakResultStoreList(PeakResultStoreList store) {
    super(store);
    this.store = store;
  }

  @Override
  public PeakResult get(int index) {
    return new ImmutablePeakResult(store.get(index));
  }

  @Override
  public PeakResult remove(int index) {
    throw new DataException("This result store is immutable");
  }

  @Override
  public void remove(int fromIndex, int toIndex) {
    throw new DataException("This result store is immutable");
  }

  @Override
  public void sort() {
    // store.sort();
    throw new DataException("This result store is immutable");
  }

  @Override
  public void sort(Comparator<PeakResult> comparator) {
    // store.sort(comparator);
    throw new DataException("This result store is immutable");
  }

  @Override
  public PeakResultStoreList copy() {
    return new ImmutablePeakResultStoreList((PeakResultStoreList) store.copy());
  }

  @Override
  public PeakResultStoreList copy(boolean deepCopy) {
    return new ImmutablePeakResultStoreList((PeakResultStoreList) store.copy(deepCopy));
  }

  @Override
  public void shuffle(final RandomGenerator randomSource) {
    // store.shuffle(randomGenerator);
    throw new DataException("This result store is immutable");
  }

  @Override
  public void shuffle(UniformRandomProvider randomSource) {
    // store.shuffle(randomSource);
    throw new DataException("This result store is immutable");
  }

  @Override
  public int indexOf(PeakResult result) {
    return store.indexOf(result);
  }

  @Override
  public int lastIndexOf(PeakResult result) {
    return store.lastIndexOf(result);
  }
}
