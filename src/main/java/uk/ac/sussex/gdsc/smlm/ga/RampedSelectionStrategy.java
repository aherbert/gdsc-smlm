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

package uk.ac.sussex.gdsc.smlm.ga;

import uk.ac.sussex.gdsc.core.data.ComputationException;

import org.apache.commons.math3.random.RandomDataGenerator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Selects the individuals using a linear ramp from the highest rank to the lowest to set the
 * probability of selection. Selection of breeding couples is uniform from the top n (n is
 * incremented from 1 each selection) crossed with a sample from the entire population using a
 * linear ramped weighting.
 *
 * @param <T> the generic type
 */
public class RampedSelectionStrategy<T extends Comparable<T>> extends SimpleSelectionStrategy<T>
    implements SelectionStrategy<T> {
  private List<? extends Chromosome<T>> sorted;
  private int numberSelected;
  private long[] sum;
  private long upper;

  /**
   * Instantiates a new ramped selection strategy.
   *
   * @param random the random
   * @param fraction The fraction of the individuals to select (set between 0 and 1)
   * @param max The maximum number of individuals to select
   */
  public RampedSelectionStrategy(RandomDataGenerator random, double fraction, int max) {
    super(random, fraction, max);
  }

  /**
   * Select the top individual and then the rest using a probability set by their rank. The
   * resulting subset will be at least size 2 (unless the input is smaller or there are not enough
   * valid individuals (fitness not null)).
   *
   * @param individuals the individuals
   * @return the subset
   */
  @Override
  public List<? extends Chromosome<T>> select(List<? extends Chromosome<T>> individuals) {
    if (individuals == null || individuals.size() < 3) {
      return individuals;
    }

    final ArrayList<Chromosome<T>> scored = new ArrayList<>(individuals.size());
    // Add only those with a fitness score
    for (final Chromosome<T> c : individuals) {
      if (c.getFitness() != null) {
        scored.add(c);
      }
    }
    if (scored.size() < 3) {
      return scored;
    }

    // Get the fraction relative to the input list size
    final int size = getSize(individuals.size());

    if (tracker != null) {
      tracker.progress(0);
    }

    // Sort the list
    ChromosomeComparator.sort(scored);

    // Create the output subset
    final ArrayList<Chromosome<T>> subset = new ArrayList<>(size);

    // Get the number of point available for selection:
    // n in this case is (size-1) since we include the top ranking individual.
    final int numberOfPoints = scored.size() - 1;

    // Add the top individual
    subset.add(scored.get(0));

    // Get the cumulative total of the rank: total = n(n+1)/2
    long cumulative = (numberOfPoints * (numberOfPoints + 1L)) / 2L;

    // Build an array of rank weighting. The highest ranked starts at n.
    // The first index is ignored since this has been included already.
    final int[] rank = new int[scored.size()];
    for (int i = 1; i < rank.length; i++) {
      rank[i] = numberOfPoints - i + 1;
    }

    // Now pick chromosomes using the cumulative as the upper limit
    while (subset.size() < size) {
      if (tracker != null) {
        tracker.progress(subset.size(), size);
      }

      // Used to check we pick something
      final long previous = cumulative;

      // Generate a random positive within the cumulative range - note the end points are inclusive
      final long next = random.nextLong(1L, cumulative);
      // Find the random position
      long total = 0;
      for (int i = 1; i < rank.length; i++) {
        total += rank[i];
        if (next <= total) {
          // Pick this point
          subset.add(scored.get(i));
          // Update the cumulative then eliminate from future selection
          cumulative -= rank[i];
          rank[i] = 0;
          break;
        }
      }

      // Check we chose something
      if (previous == cumulative) {
        throw new ComputationException(
            "Failed to select a candidate. Size = " + subset.size() + " / " + size);
      }
    }

    if (tracker != null) {
      tracker.progress(1);
    }

    return subset;
  }

  /** {@inheritDoc} */
  @Override
  public void initialiseBreeding(List<? extends Chromosome<T>> individuals) {
    sorted = null;
    if (individuals == null || individuals.size() < 2) {
      return;
    }

    // This method may be called when fitness is unknown.
    // Only sort those with a fitness score.
    ArrayList<Chromosome<T>> list = null;
    int count = 0;
    for (final Chromosome<T> c : individuals) {
      if (c.getFitness() == null) {
        count++;
      }
    }

    if (count != 0) {
      final int toSort = individuals.size() - count;
      if (toSort == 0) {
        // No sort possible. Revert to the simple selection strategy
        super.initialiseBreeding(individuals);
        return;
      }
      // A mixed population, some of which we can sort and some we cannot.
      list = new ArrayList<>(toSort);
      final ArrayList<Chromosome<T>> subset = new ArrayList<>(count);
      for (final Chromosome<T> c : individuals) {
        if (c.getFitness() == null) {
          subset.add(c);
        } else {
          list.add(c);
        }
      }

      ChromosomeComparator.sort(list);

      // Create a ramped sum for all those we can sort
      sum = createRampedSum(list.size());
      // Extend the sum linearly for those we cannot sort (i.e. they have the same selection chance)
      sum = extendSumLinearly(sum, subset.size());

      list.addAll(subset);
    } else {
      list = new ArrayList<>(individuals);
      ChromosomeComparator.sort(list);

      // Build a cumulative array of rank weighting. The highest ranked starts at n.
      sum = createRampedSum(list.size());
    }

    this.sorted = list;

    // Reset
    numberSelected = 0;

    // Bounds are inclusive so subtract 1
    upper = sum[sum.length - 1] - 1;
  }

  /**
   * Creates the sum.
   *
   * @param size the size
   * @return the sum
   */
  public static long[] createRampedSum(int size) {
    final long[] sum = new long[size];
    int rank = size;
    sum[0] = rank--;
    for (int i = 1; i < sum.length; i++) {
      sum[i] = rank-- + sum[i - 1];
    }
    return sum;
  }

  /**
   * Extend the sum.
   *
   * @param sum the sum
   * @param addition the addition
   * @return the new sum
   */
  public static long[] extendSumLinearly(long[] sum, int addition) {
    final long[] sum2 = Arrays.copyOf(sum, sum.length + addition);
    long rank = sum[sum.length - 1];
    for (int i = sum.length; i < sum2.length; i++) {
      sum2[i] = ++rank;
    }
    return sum2;
  }

  /**
   * Select pairs randomly from the population. The first is selected from the top n individuals
   * with n starting at 1 and incrementing for each call. The second is selected from the entire
   * population with the weighting equal to their ranking by fitness.
   */
  @Override
  public ChromosomePair<T> next() {
    if (sorted == null) {
      return super.next();
    }
    int first;
    int second;
    if (sorted.size() == 2) {
      first = 0;
      second = 1;
    } else {
      // Pick from the first n
      if (numberSelected == 0) {
        // This is the first call so select the top individual
        first = 0;
      } else {
        // Restrict the upper limit to the population range
        first = random.nextInt(0, Math.min(numberSelected, sorted.size() - 1));
      }
      numberSelected++;

      // Select the second from the ramped cumulative
      second = nextSample();
      // Avoid crossover with the same parent
      while (second == first) {
        second = nextSample();
      }
    }
    return new ChromosomePair<>(sorted.get(first), sorted.get(second));
  }

  private int nextSample() {
    final long next = random.nextLong(0, upper);
    // JUnit test shows that the simple scan through the array is faster for small arrays
    if (sum.length < 100) {
      return search(sum, next);
    }
    return binarySearch(sum, next);
  }

  /**
   * Find the index such that {@code sum[index-1] <= key < sum[index]}.
   *
   * @param sum the sum
   * @param key the key
   * @return the index (or -1)
   */
  public static int search(final long[] sum, final long key) {
    for (int i = 0; i < sum.length; i++) {
      if (key < sum[i]) {
        return i;
      }
    }
    return 0; // This should not happen
  }

  /**
   * Find the index such that {@code sum[index-1] <= key < sum[index]}.
   *
   * @param sum the sum
   * @param key the key
   * @return the index (or -1)
   */
  public static int binarySearch(final long[] sum, final long key) {
    if (key < sum[0]) {
      return 0;
    }

    // Adapted from Arrays.binarySearch which
    // finds the actual key value or returns a negative insertion point:

    // If we find the actual key return the next index above.
    // If we do not find the actual key return the insertion point

    int low = 0;
    int high = sum.length - 1;

    while (low <= high) {
      final int mid = (low + high) >>> 1;
      final long midVal = sum[mid];

      if (midVal < key) {
        low = mid + 1;
      } else if (midVal > key) {
        high = mid - 1;
      } else {
        // Changed from Arrays.binarySearch
        // return mid; // key found
        // If we find the actual key return the next index above.
        return mid + 1;
      }
    }

    // Changed from Arrays.binarySearch
    // return -(low + 1); // key not found.
    // If we do not find the actual key return the insertion point
    return low;
  }

  /** {@inheritDoc} */
  @Override
  public void finishBreeding() {
    // Free memory
    sorted = null;
    sum = null;
  }
}
