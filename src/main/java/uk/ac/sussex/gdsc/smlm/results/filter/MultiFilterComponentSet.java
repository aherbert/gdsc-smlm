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

package uk.ac.sussex.gdsc.smlm.results.filter;

/**
 * Contains a set of components of the multi filter.
 */
public interface MultiFilterComponentSet {
  /**
   * Gets the validation flags. These are possible return flags from the
   * {@link #validate(PreprocessedPeakResult)} method.
   *
   * @return the validation flags
   */
  int getValidationFlags();

  /**
   * Validate the peak.
   *
   * @param peak the peak
   * @return the result
   */
  int validate(final PreprocessedPeakResult peak);

  /**
   * Replace the first component.
   *
   * @param component the replacement component
   */
  abstract void replace0(MultiFilterComponent component);

  /**
   * Create a copy.
   *
   * @return the copy
   */
  MultiFilterComponentSet copy();
}
