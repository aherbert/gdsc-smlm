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

package uk.ac.sussex.gdsc.smlm.results.filter;

/**
 * Contains Support direct filtering of PreprocessedPeakResult objects.
 *
 * <p>The decision to support for filtering as both a DirectFilter and Filter concurrently is left
 * to the implementing class. It is not a requirement.
 */
public final class FilterValidationOption {
  /**
   * Disable filtering using the width of the result.
   */
  public static final int NO_WIDTH = 0x000000001;

  /**
   * Disable filtering using the shift of the result.
   */
  public static final int NO_SHIFT = 0x000000002;

  /**
   * Enable filtering both X and Y widths.
   */
  public static final int XY_WIDTH = 0x000000004;

  /**
   * Disable Z filtering (use when not fitting in 3D).
   */
  public static final int NO_Z = 0x000000008;

  /** No public construction. */
  private FilterValidationOption() {}
}
