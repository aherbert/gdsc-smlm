/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
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
 * Filter results using Shift.
 */
public class MultiFilterShiftComponent extends MultiFilterComponent {
  private static final int type =
      IDirectFilter.V_X_RELATIVE_SHIFT | IDirectFilter.V_Y_RELATIVE_SHIFT;

  private final float offset;

  /**
   * Instantiates a new multi filter shift component.
   *
   * @param shift the shift
   */
  public MultiFilterShiftComponent(double shift) {
    this.offset = Filter.getUpperSquaredLimit(shift);
  }

  /** {@inheritDoc} */
  @Override
  public boolean fail(final PreprocessedPeakResult peak) {
    if (peak.getXRelativeShift2() > offset) {
      return true;
    }
    return (peak.getYRelativeShift2() > offset);
  }

  /** {@inheritDoc} */
  @Override
  public int getType() {
    return type;
  }
}
