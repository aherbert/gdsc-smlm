/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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
 * Filter results using Signal.
 */
public class MultiFilterSignalComponent implements MultiFilterComponent {
  private final float signal;

  /**
   * Instantiates a new multi filter signal component.
   *
   * @param signal the signal
   */
  public MultiFilterSignalComponent(double signal) {
    this.signal = (float) signal;
  }

  @Override
  public boolean fail(final PreprocessedPeakResult peak) {
    return peak.getSignal() < this.signal;
  }

  @Override
  public int getType() {
    return FilterValidationFlag.PHOTONS;
  }
}
