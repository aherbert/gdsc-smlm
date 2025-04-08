/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Package
 *
 * Software for single molecule localisation microscopy (SMLM) in ImageJ
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

package uk.ac.sussex.gdsc.smlm.ij.ij3d;

/**
 * Contains helper methods.
 */
final class ItemHelper {
  /** Threshold below which the transparency is zero. */
  static final float ZERO_TRANSPARENCY = 0.01f;

  /** No public construction. */
  private ItemHelper() {}

  /**
   * Check the actual size matches the expected size.
   *
   * @param actual the actual
   * @param expected the expected
   * @throws IllegalArgumentException if the size does not match
   */
  static void checkSize(final float actual, final int expected) {
    if (actual != expected) {
      throw new IllegalArgumentException("list of size " + expected + " expected");
    }
  }
}
