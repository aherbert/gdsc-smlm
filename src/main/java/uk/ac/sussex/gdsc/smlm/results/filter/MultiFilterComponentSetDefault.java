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
public class MultiFilterComponentSetDefault implements MultiFilterComponentSet {
  private final MultiFilterComponent[] components;

  /**
   * Instantiates a new multi filter component set default.
   *
   * @param components the components
   */
  public MultiFilterComponentSetDefault(MultiFilterComponent[] components) {
    this.components = components;
  }

  @Override
  public int getValidationFlags() {
    int flags = 0;
    for (int i = 0; i < components.length; i++) {
      flags |= components[i].getType();
    }
    return flags;
  }

  @Override
  public int validate(final PreprocessedPeakResult peak) {
    for (int i = 0; i < components.length; i++) {
      if (components[i].fail(peak)) {
        return components[i].getType();
      }
    }
    return 0;
  }

  @Override
  public void replace0(MultiFilterComponent component) {
    if (components.length > 0) {
      components[0] = component;
    }
  }

  @Override
  public MultiFilterComponentSet copy() {
    // Copy the array
    final MultiFilterComponent[] cmp = new MultiFilterComponent[components.length];
    if (cmp.length > 0) {
      System.arraycopy(components, 0, cmp, 0, cmp.length);
    }
    return new MultiFilterComponentSetDefault(cmp);
  }
}
