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

package uk.ac.sussex.gdsc.smlm.ga;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * Sorts chromosome using the fitness, highest fitness first.
 *
 * @param <T> the generic type
 */
public class ChromosomeComparator<T extends Comparable<T>>
    implements Comparator<Chromosome<T>>, Serializable {
  private static final long serialVersionUID = 1L;

  @Override
  public int compare(Chromosome<T> chromosome1, Chromosome<T> chromosome2) {
    return chromosome1.getFitness().compareTo(chromosome2.getFitness());
  }

  /**
   * Sort the list (highest fitness first).
   *
   * @param <T> the generic type
   * @param list the list
   */
  public static <T extends Comparable<T>> void sort(List<? extends Chromosome<T>> list) {
    Collections.sort(list, new ChromosomeComparator<T>());
  }

  /**
   * Sort the list (highest fitness first).
   *
   * @param <T> the generic type
   * @param list the list
   */
  public static <T extends Comparable<T>> void sort(Chromosome<T>[] list) {
    sort(list, 0, list.length);
  }

  /**
   * Sort the list (highest fitness first).
   *
   * @param <T> the generic type
   * @param list the list
   * @param fromIndex the from index
   * @param toIndex the to index
   */
  public static <T extends Comparable<T>> void sort(Chromosome<T>[] list, int fromIndex,
      int toIndex) {
    Arrays.sort(list, fromIndex, toIndex, new ChromosomeComparator<T>());
  }
}
