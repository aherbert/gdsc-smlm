/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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

import uk.ac.sussex.gdsc.core.annotation.Nullable;

/**
 * Define the genetic sequence that can evolve.
 *
 * @param <T> the generic type
 */
public interface Chromosome<T extends Comparable<T>> {
  /**
   * Get the chromosome length.
   *
   * @return The chromosome length
   */
  int length();

  /**
   * Get the chromosome sequence.
   *
   * @return the chromosome sequence (must equal the length)
   */
  double[] sequence();

  /**
   * Create a new chromosome.
   *
   * @param sequence the chromosome sequence (must equal the current length)
   * @return A new chromosome with the given sequence
   */
  Chromosome<T> newChromosome(double[] sequence);

  /**
   * Get the range for mutation at each position in the sequence. This defines how far each position
   * in the sequence can mutate in a single step.
   *
   * @return The range for mutation at each position in the sequence (must equal length)
   */
  double[] mutationStepRange();

  /**
   * Get the lower limit at each position in the sequence. It is valid to return negative infinity
   * for any position or null for no limit.
   *
   * @return The lower limit for each position in the sequence (must equal length)
   */
  @Nullable
  double[] lowerLimit();

  /**
   * Get the upper limit at each position in the sequence. It is valid to return positive infinity
   * for any position or null for no limit.
   *
   * @return The upper limit for each position in the sequence (must equal length)
   */
  @Nullable
  double[] upperLimit();

  // Note: Default implementation of the getter/setter to store the double would require using Java
  // 8.

  /**
   * Set the fitness.
   *
   * @param fitness The fitness of the sequence
   */
  void setFitness(T fitness);

  /**
   * Get the fitness.
   *
   * <p>This should be null for an uninitialised score. The comparable should rank in ascending
   * order with the first item the fittest individual.
   *
   * @return The fitness of the sequence
   */
  T getFitness();

  /**
   * Calculate the distance to another chromosome.
   *
   * @param other the other chromosome
   * @return the distance (zero is a match)
   */
  double distance(Chromosome<T> other);

  /**
   * Calculate if equal to another chromosome.
   *
   * @param other the other chromosome
   * @return true if the same
   */
  boolean equalTo(Chromosome<T> other);
}
