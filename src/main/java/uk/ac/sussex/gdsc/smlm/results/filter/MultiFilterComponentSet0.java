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
 * Contains a set of components of the multi filter.
 */
public class MultiFilterComponentSet0 implements MultiFilterComponentSet {
  /**
   * Instantiates a new multi filter component set for 0 components.
   *
   * @param components the components
   */
  MultiFilterComponentSet0(MultiFilterComponent[] components) {
    // Ignore the array
  }

  @Override
  public int getValidationFlags() {
    return 0;
  }

  @Override
  public int validate(final PreprocessedPeakResult peak) {
    return 0;
  }

  @Override
  public void replace0(MultiFilterComponent component) {
    // This set is empty so no replacement is possible. Don't throw an exception though!
  }

  @Override
  public MultiFilterComponentSet copy() {
    return new MultiFilterComponentSet0(null);
  }
}
