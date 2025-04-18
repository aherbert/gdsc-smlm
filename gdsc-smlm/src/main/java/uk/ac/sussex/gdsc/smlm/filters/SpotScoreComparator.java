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

package uk.ac.sussex.gdsc.smlm.filters;

import java.util.Comparator;

/**
 * Sort {@link Spot} using the score in descending order.
 */
public enum SpotScoreComparator implements Comparator<Spot> {
  /** An instance of the comparator. */
  INSTANCE;

  /**
   * Gets an instance.
   *
   * @return an instance
   */
  public static SpotScoreComparator getInstance() {
    return INSTANCE;
  }

  @Override
  public int compare(Spot o1, Spot o2) {
    // Descending order
    return Double.compare(o2.getScore(), o1.getScore());
  }
}
