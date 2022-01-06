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

package uk.ac.sussex.gdsc.smlm.ga;

import java.util.Arrays;
import java.util.List;
import org.apache.commons.rng.UniformRandomProvider;
import uk.ac.sussex.gdsc.core.logging.TrackProgress;

/**
 * Selects the top individuals.
 *
 * @param <T> the generic type
 */
public class SimpleSelectionStrategy<T extends Comparable<T>> extends Randomiser
    implements SelectionStrategy<T> {
  /** The fraction of the individuals to select (set between 0 and 1). */
  final double fraction;
  /** The maximum number of individuals to select. */
  final int max;

  private List<? extends Chromosome<T>> individuals;

  /** The tracker. */
  TrackProgress tracker;

  /**
   * Instantiates a new simple selection strategy.
   *
   * @param random the random
   * @param fraction The fraction of the individuals to select (set between 0 and 1)
   * @param max The maximum number of individuals to select
   */
  public SimpleSelectionStrategy(UniformRandomProvider random, double fraction, int max) {
    super(random);
    if (fraction > 1) {
      fraction = 1;
    }
    this.fraction = fraction;
    this.max = max;
  }

  /**
   * Select the top individuals using the configured fraction. The resulting subset will be at least
   * size 2 (unless the input is smaller or there are not enough valid individuals (fitness not
   * null)).
   *
   * @param individuals the individuals
   * @return the subset
   */
  @Override
  public List<? extends Chromosome<T>> select(List<? extends Chromosome<T>> individuals) {
    if (individuals == null || individuals.size() < 2) {
      return individuals;
    }
    @SuppressWarnings("unchecked")
    final Chromosome<T>[] subset = new Chromosome[individuals.size()];
    int size = 0;
    // Add only those with a fitness score
    for (final Chromosome<T> c : individuals) {
      if (c.getFitness() != null) {
        subset[size++] = c;
      }
    }
    if (size < 3) {
      return Arrays.asList(Arrays.copyOf(subset, size));
    }
    if (tracker != null) {
      tracker.progress(0.5);
    }
    ChromosomeComparator.sort(subset, 0, size);

    // Get the fraction relative to the input list size
    // size = getSize(size);
    size = Math.min(size, getSize(individuals.size()));
    if (tracker != null) {
      tracker.progress(1);
    }
    return Arrays.asList(Arrays.copyOf(subset, size));
  }

  /**
   * Calculate the new size of the population after selection.
   *
   * @param size The current size of the population before selection
   * @return The new size of the population
   */
  protected int getSize(int size) {
    // Get the size using the fraction
    size = (int) Math.round(size * fraction);
    // Check against the max number to select
    if (max > 2 && size > max) {
      size = max;
    }
    // Check the size is at least 2
    if (size < 2) {
      size = 2;
    }
    return size;
  }

  @Override
  public void initialiseBreeding(List<? extends Chromosome<T>> individuals) {
    if (individuals != null && individuals.size() < 2) {
      individuals = null;
    }
    this.individuals = individuals;
  }

  /**
   * Select pairs randomly from the population.
   */
  @Override
  public ChromosomePair<T> next() {
    if (individuals == null) {
      return null;
    }
    int first;
    int second;
    if (individuals.size() == 2) {
      first = 0;
      second = 1;
    } else {
      final int upper = individuals.size();
      first = random.nextInt(upper);
      second = random.nextInt(upper);
      // Avoid crossover with the same parent
      while (second == first) {
        second = random.nextInt(upper);
      }
    }
    return new ChromosomePair<>(individuals.get(first), individuals.get(second));
  }

  @Override
  public void finishBreeding() {
    // Free memory
    individuals = null;
  }

  @Override
  public void setTracker(TrackProgress tracker) {
    this.tracker = tracker;
  }
}
