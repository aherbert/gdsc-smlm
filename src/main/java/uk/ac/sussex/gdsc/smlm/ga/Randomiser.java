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

import org.apache.commons.rng.UniformRandomProvider;
import uk.ac.sussex.gdsc.core.utils.ValidationUtils;

/**
 * Base class for data generation using randomness.
 */
public abstract class Randomiser {
  /** The source of randomness. */
  final UniformRandomProvider random;

  /**
   * Instantiates a new randomiser.
   *
   * @param random the random data generator
   * @throws NullPointerException if the random generator is null
   */
  public Randomiser(UniformRandomProvider random) {
    this.random = ValidationUtils.checkNotNull(random, "Random generator cannot be null");
  }
}
