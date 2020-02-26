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

package uk.ac.sussex.gdsc.smlm.results.filter;

/**
 * Contains a set of components of the multi filter.
 */
public class MultiFilterComponentSet2 implements MultiFilterComponentSet {
  private MultiFilterComponent component0;
  private final MultiFilterComponent component1;

  /**
   * Instantiates a new multi filter component set for 2 components.
   *
   * @param components the components
   */
  MultiFilterComponentSet2(MultiFilterComponent[] components) {
    this.component0 = components[0];
    this.component1 = components[1];
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  MultiFilterComponentSet2(MultiFilterComponentSet2 source) {
    this.component0 = source.component0;
    this.component1 = source.component1;
  }

  @Override
  public int getValidationFlags() {
    return component0.getType() | component1.getType();
  }

  @Override
  public int validate(final PreprocessedPeakResult peak) {
    //@formatter:off
    if (component0.fail(peak)) {
      return component0.getType();
    }
    if (component1.fail(peak)) {
      return component1.getType();
    }
    //@formatter:on
    return 0;
  }

  @Override
  public void replace0(MultiFilterComponent component) {
    component0 = component;
  }

  @Override
  public MultiFilterComponentSet copy() {
    return new MultiFilterComponentSet2(this);
  }
}
