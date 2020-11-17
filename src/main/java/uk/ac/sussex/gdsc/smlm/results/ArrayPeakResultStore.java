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

package uk.ac.sussex.gdsc.smlm.results;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.Objects;
import java.util.function.Predicate;
import org.apache.commons.rng.UniformRandomProvider;
import uk.ac.sussex.gdsc.core.utils.MemoryUtils;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;

/**
 * Stores peak results using an array.
 */
public class ArrayPeakResultStore implements PeakResultStoreList, Serializable {
  private static final long serialVersionUID = 20190319L;

  /** The results. */
  private PeakResult[] results;

  /** The size. */
  private int size;

  /**
   * Instantiates a new array list peak results store.
   *
   * @param capacity the capacity
   */
  public ArrayPeakResultStore(int capacity) {
    this.results = new PeakResult[Math.max(capacity, 0)];
  }

  /**
   * Instantiates a new array peak result store.
   *
   * @param store the store to copy
   * @throws NullPointerException if the store is null
   */
  public ArrayPeakResultStore(ArrayPeakResultStore store) {
    this.results = store.toArray();
    this.size = store.size;
  }

  /**
   * Instantiates a new array peak result store.
   *
   * @param results the results
   * @throws NullPointerException if the results are null
   */
  public ArrayPeakResultStore(PeakResult[] results) {
    this.results = results;
    this.size = results.length;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Note: This does not check against the current size so can return stale data.
   */
  @Override
  public PeakResult get(int index) {
    return results[index];
  }

  @Override
  public int size() {
    return size;
  }

  /**
   * Increase the capacity to hold at least one more than the current capacity. This will reallocate
   * a new array and should only be called when the capacity has been checked and is known to be too
   * small.
   *
   * @return the new data array
   */
  private PeakResult[] increaseCapacity() {
    return results =
        Arrays.copyOf(results, MemoryUtils.createNewCapacity(results.length + 1, results.length));
  }

  /**
   * Increase the capacity to hold at least the minimum required capacity. This will reallocate a
   * new array and should only be called when the capacity has been checked and is known to be too
   * small.
   *
   * @param minCapacity the minimum required capacity
   * @return the new results array
   */
  private PeakResult[] increaseCapacity(final int minCapacity) {
    return results =
        Arrays.copyOf(results, MemoryUtils.createNewCapacity(minCapacity, results.length));
  }

  @Override
  public boolean add(PeakResult result) {
    final int s = size;
    PeakResult[] r = results;
    if (s == r.length) {
      r = increaseCapacity();
    }
    size = s + 1;
    r[s] = result;
    return true;
  }

  @Override
  public boolean addCollection(Collection<PeakResult> results) {
    return addArray(results.toArray(new PeakResult[0]));
  }

  @Override
  public boolean addArray(PeakResult[] results) {
    if (results == null) {
      return false;
    }
    return addArray(results, results.length);
  }

  private boolean addArray(PeakResult[] results, int length) {
    if (results == null || length == 0) {
      return false;
    }
    final int s = size;
    PeakResult[] r = this.results;
    // spare = r.length - size
    if (length > r.length - s) {
      r = increaseCapacity(s + length);
    }
    // Append
    System.arraycopy(results, 0, r, s, length);
    size = s + length;
    return true;
  }

  @Override
  public boolean addStore(PeakResultStore results) {
    if (results instanceof ArrayPeakResultStore) {
      final ArrayPeakResultStore store = (ArrayPeakResultStore) results;
      return addArray(store.results, store.size);
    }
    return addArray(results.toArray());
  }

  @Override
  public PeakResult remove(int index) {
    rangeCheck(index);
    final PeakResult oldValue = results[index];
    fastRemove(index);
    return oldValue;
  }

  @Override
  public void remove(int fromIndex, int toIndex) {
    if (fromIndex > toIndex) {
      throw new IllegalArgumentException("fromIndex must be <= toIndex");
    }
    rangeCheckWithLowerBounds(fromIndex);
    rangeCheck(toIndex); // This is above fromIndex so ignore lower bounds check
    toIndex++; // Make exclusive
    final int numMoved = size - toIndex;
    if (numMoved > 0) {
      System.arraycopy(results, toIndex, results, fromIndex, numMoved);
    }
    // Let gc do its work
    while (fromIndex++ < toIndex) {
      results[size--] = null;
    }
  }

  @Override
  public boolean remove(PeakResult result) {
    final int index = indexOf(result);
    if (index != -1) {
      fastRemove(index);
      return true;
    }
    return false;
  }

  /**
   * Checks if the given index is in range. If not, throws an appropriate runtime exception. This
   * method does *not* check if the index is negative: It is always used immediately prior to an
   * array access, which throws an ArrayIndexOutOfBoundsException if index is negative.
   */
  private void rangeCheck(int index) {
    if (index >= size) {
      throw new IndexOutOfBoundsException(outOfBoundsMsg(index));
    }
  }

  /**
   * A version of rangeCheck with lower bounds check.
   */
  private void rangeCheckWithLowerBounds(int index) {
    if (index > size || index < 0) {
      throw new IndexOutOfBoundsException(outOfBoundsMsg(index));
    }
  }

  /**
   * Constructs an IndexOutOfBoundsException detail message.
   */
  private String outOfBoundsMsg(int index) {
    return "Index: " + index + ", Size: " + size;
  }

  /*
   * Private remove method that skips bounds checking and does not return the value removed.
   */
  private void fastRemove(int index) {
    final int numMoved = size - index - 1;
    if (numMoved > 0) {
      System.arraycopy(results, index + 1, results, index, numMoved);
    }
    results[--size] = null; // Let gc do its work
  }

  @Override
  public boolean removeCollection(Collection<PeakResult> results) {
    return removeIf(results::contains);
  }

  @Override
  public boolean removeArray(PeakResult[] results) {
    if (results == null || results.length == 0) {
      return false;
    }
    return removeStore(new ArrayPeakResultStore(results));
  }

  @Override
  public boolean removeStore(PeakResultStore results) {
    return removeIf(results::contains);
  }

  @Override
  public boolean retainCollection(Collection<PeakResult> results) {
    return removeIf(e -> !results.contains(e));
  }

  @Override
  public boolean retainArray(PeakResult[] results) {
    if (results == null || results.length == 0) {
      final boolean result = size != 0;
      clear();
      return result;
    }
    return retainStore(new ArrayPeakResultStore(results));
  }

  @Override
  public boolean retainStore(PeakResultStore results) {
    return removeIf(e -> !results.contains(e));
  }

  /**
   * {@inheritDoc}
   *
   * <p>Note: This does not remove the references to the underlying data or reallocate storage thus
   * {@link #get(int)} can return stale data.
   */
  @Override
  public void clear() {
    size = 0;
  }

  @Override
  public void trimToSize() {
    if (size < results.length) {
      results = toArray();
    }
  }

  @Override
  public void sort(Comparator<PeakResult> comparator) {
    Arrays.sort(results, 0, size, comparator);
  }

  @Override
  public PeakResult[] toArray() {
    return Arrays.copyOf(results, size);
  }

  @Override
  public PeakResultStore copy() {
    return new ArrayPeakResultStore(this);
  }

  @Override
  public PeakResultStore copy(boolean deepCopy) {
    if (deepCopy) {
      final ArrayPeakResultStore copy = new ArrayPeakResultStore(size());
      for (int i = 0, max = size(); i < max; i++) {
        copy.add(results[i].copy());
      }
      return copy;
    }
    return copy();
  }

  /**
   * {@inheritDoc}
   *
   * <p>Functions as if testing all items:
   *
   * <pre>
   * int newSize = 0;
   * for (int i = 0; i < size; i++) {
   *   if (filter.test(data[i])) {
   *     // remove
   *     continue;
   *   }
   *   data[newSize++] = data[i];
   * }
   * size = newSize;
   * </pre>
   */
  @Override
  public boolean removeIf(Predicate<PeakResult> filter) {
    Objects.requireNonNull(filter);

    int index = 0;
    final PeakResult[] elements = results;
    final int length = size;

    // Find first item to filter.
    for (; index < length; index++) {
      if (filter.test(elements[index])) {
        break;
      }
    }
    if (index == length) {
      // Nothing to remove
      return false;
    }

    // The list has changed.
    // Do a single sweep across the remaining data copying elements if they are not filtered.
    // Note: The filter may throw and leave the list without a new size.
    // (E.g. is the filter is a collection that does not allow null).
    // So we set the new size in a finally block.
    int newSize = index;
    try {
      // We know the current index is identified by the filter so advance 1
      index++;

      // Scan the rest
      for (; index < length; index++) {
        final PeakResult e = elements[index];
        if (filter.test(e)) {
          continue;
        }
        elements[newSize++] = e;
      }
    } finally {
      // Ensure the length is correct
      if (index != length) {
        // Did not get to the end of the list (e.g. the filter may throw) so copy it verbatim.
        final int len = length - index;
        System.arraycopy(elements, index, elements, newSize, len);
        newSize += len;
      }
      // Clear old references
      for (int i = newSize; i < length; i++) {
        elements[i] = null;
      }
      size = newSize;
    }
    // The list was modified
    return true;
  }

  @Override
  public void forEach(PeakResultProcedure procedure) {
    for (int i = 0; i < size; i++) {
      procedure.execute(results[i]);
    }
  }

  @Override
  public PeakResult[] subset(Predicate<PeakResult> filter) {
    final ArrayPeakResultStore list = new ArrayPeakResultStore(10);
    final PeakResult[] r = results;
    final int length = size;
    for (int i = 0; i < length; i++) {
      if (filter.test(r[i])) {
        list.add(r[i]);
      }
    }
    return list.toArray();
  }

  @Override
  public void shuffle(UniformRandomProvider randomSource) {
    // Fisher-Yates shuffle
    final PeakResult[] r = results;
    for (int i = size; i-- > 1;) {
      final int j = randomSource.nextInt(i + 1);
      final PeakResult tmp = r[i];
      r[i] = r[j];
      r[j] = tmp;
    }
  }

  @Override
  public int indexOf(PeakResult result) {
    final PeakResult[] r = results;
    final int length = size;
    if (result == null) {
      for (int i = 0; i < length; i++) {
        if (r[i] == null) {
          return i;
        }
      }
    } else {
      for (int i = 0; i < length; i++) {
        if (result.equals(r[i])) {
          return i;
        }
      }
    }
    return -1;
  }

  @Override
  public int lastIndexOf(PeakResult result) {
    final PeakResult[] r = results;
    final int length = size;
    if (result == null) {
      for (int i = length; i-- > 0;) {
        if (r[i] == null) {
          return i;
        }
      }
    } else {
      for (int i = length; i-- > 0;) {
        if (result.equals(r[i])) {
          return i;
        }
      }
    }
    return -1;
  }

  @Override
  public boolean contains(PeakResult result) {
    return indexOf(result) != -1;
  }
}
