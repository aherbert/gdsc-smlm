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

/**
 * Defines the fitness of a chromosome.
 *
 * <p>Note: this class has a natural ordering that is inconsistent with equals.
 *
 * @param <T> the generic type
 */
public class Fitness<T extends Comparable<T>> implements Comparable<Fitness<T>> {
  /** The comparable object. */
  final T object;

  /** The score. */
  final double score;

  /**
   * Instantiates a new fitness.
   *
   * @param object the comparable object
   * @param score the score
   */
  public Fitness(T object, double score) {
    this.object = object;
    this.score = score;
  }

  @Override
  public int compareTo(Fitness<T> other) {
    if (object == null) {
      return 1;
    }
    if (other.object == null) {
      return -1;
    }
    return object.compareTo(other.object);
  }
}
