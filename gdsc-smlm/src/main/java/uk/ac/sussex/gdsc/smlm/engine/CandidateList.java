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

package uk.ac.sussex.gdsc.smlm.engine;

import java.util.Arrays;
import java.util.Comparator;

/**
 * Stores a list of candidates.
 */
class CandidateList {

  private int size;
  private Candidate[] list;

  /**
   * Simple interface for testing candidates.
   */
  public interface Predicate {
    /**
     * Test.
     *
     * @param candidate the candidate
     * @return true, if successful
     */
    boolean test(Candidate candidate);
  }

  /** Candidate Comparator. */
  private enum CandidateComparator implements Comparator<Candidate> {
    /** An instance of the comparator. */
    INSTANCE;

    @Override
    public int compare(Candidate o1, Candidate o2) {
      return o1.index - o2.index;
    }
  }

  /**
   * Instantiates a new candidate list.
   */
  CandidateList() {
    // Intentionally empty
  }

  /**
   * Instantiates a new candidate list.
   *
   * @param size the size
   * @param list the list
   */
  CandidateList(int size, Candidate[] list) {
    if (list != null) {
      if (list.length < size) {
        throw new IllegalArgumentException("List is smaller than size");
      }
    } else if (size != 0) {
      throw new IllegalArgumentException("List is null and size is not zero");
    }

    this.size = size;
    this.list = list;
  }

  /**
   * Instantiates a new candidate list.
   *
   * @param list the list
   */
  CandidateList(Candidate[] list) {
    if (list != null) {
      size = list.length;
    }
    this.list = list;
  }

  /**
   * Add a candidate.
   *
   * @param candidate the candidate
   */
  public void add(Candidate candidate) {
    if (list == null) {
      list = new Candidate[4];
    } else if (list.length == size) {
      // Allow creating a new list even if size is currently zero
      final Candidate[] list2 = new Candidate[1 + size * 2];
      System.arraycopy(list, 0, list2, 0, size);
      list = list2;
    }
    list[size++] = candidate;
  }

  /**
   * Sort in ascending order of Id.
   */
  public void sort() {
    if (size != 0) {
      Arrays.sort(list, 0, size, CandidateComparator.INSTANCE);
    }
  }

  /**
   * Gets the size.
   *
   * @return the size
   */
  public int getSize() {
    return size;
  }

  /**
   * Gets the length of the list. This may be larger than {@link #getSize()}. It is used when the
   * list of candidates is larger than the max candidate to process.
   *
   * @return the length of the list
   */
  int getLength() {
    return list == null ? 0 : list.length;
  }

  /**
   * Gets the candidate.
   *
   * @param index the index
   * @return the candidate
   */
  public Candidate get(int index) {
    return list[index];
  }

  /**
   * Copy this list.
   *
   * @return the new candidate list
   */
  public CandidateList copy() {
    if (size == 0) {
      return new CandidateList();
    }
    return new CandidateList(size, Arrays.copyOf(list, size));
  }

  /**
   * Removes candidates from the list if they fail the filter.
   *
   * @param filter the filter
   */
  public void removeIf(Predicate filter) {
    final int oldSize = size;
    // This in place resizes the list
    size = 0;
    for (int i = 0; i < oldSize; i++) {
      if (!filter.test(list[i])) {
        list[size++] = list[i];
      }
    }
  }

  /**
   * Copy to the destination list.
   *
   * @param list the list
   * @param position the position
   */
  public void copyTo(Candidate[] list, int position) {
    System.arraycopy(this.list, 0, list, position, size);
  }
}
