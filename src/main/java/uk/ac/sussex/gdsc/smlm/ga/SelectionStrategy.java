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

import java.util.List;
import uk.ac.sussex.gdsc.core.logging.TrackProgress;

/**
 * Defines a selection strategy of a population of individuals.
 *
 * @param <T> the generic type
 */
public interface SelectionStrategy<T extends Comparable<T>> {
  /**
   * Select a subset from the population using the fitness.
   *
   * @param individuals the individuals
   * @return a selection of individuals
   */
  List<? extends Chromosome<T>> select(List<? extends Chromosome<T>> individuals);

  /**
   * Initialise the selection of pairs for breeding using the fitness.
   *
   * @param individuals the population of individuals
   */
  void initialiseBreeding(List<? extends Chromosome<T>> individuals);

  /**
   * Get the next pair of individuals for breeding. Must be called after
   * {@link #initialiseBreeding(List)}.
   *
   * @return The next pair
   */
  ChromosomePair<T> next();

  /**
   * Finish selection of pairs for breeding.
   */
  void finishBreeding();

  /**
   * Set the tracker used to track progress. This should be used in the {@link #select(List)}
   * method. It is possible to know the end point when breeding and so it is not advised to use the
   * tracker in the {@link #next()} method.
   *
   * @param tracker the new tracker
   */
  void setTracker(TrackProgress tracker);
}
