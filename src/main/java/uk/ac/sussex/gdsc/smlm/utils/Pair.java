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

package uk.ac.sussex.gdsc.smlm.utils;

/**
 * A generic pair.
 *
 * @param <A> the generic type
 * @param <B> the generic type
 */
public class Pair<A, B> {
  /** First item of the pair. */
  public final A item1;
  /** Second item of the pair. */
  public final B item2;

  /**
   * Instantiates a new pair.
   *
   * @param item1 the first item of the pair
   * @param item2 the second item of the pair
   */
  public Pair(A item1, B item2) {
    this.item1 = item1;
    this.item2 = item2;
  }
}
