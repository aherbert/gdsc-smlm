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

package uk.ac.sussex.gdsc.smlm.results;

/**
 * Gets a data value from a peak result.
 *
 * @param <E> the element type
 */
public interface PeakResultData<E> {
  /**
   * Gets the value of the result.
   *
   * @param result the result
   * @return the value
   */
  E getValue(PeakResult result);

  /**
   * Gets the name of the value.
   *
   * @return the name
   */
  String getValueName();

  /**
   * Gets the class type of the value.
   *
   * @return the name
   */
  Class<?> getValueClass();
}
