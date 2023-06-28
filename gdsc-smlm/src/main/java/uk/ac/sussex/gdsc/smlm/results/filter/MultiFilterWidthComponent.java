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
 * Filter results using Width. Assumes XY width are the same.
 */
public class MultiFilterWidthComponent implements MultiFilterComponent {
  private final float lowerSigmaThreshold;
  private final float upperSigmaThreshold;

  /**
   * Instantiates a new multi filter width component.
   *
   * @param minWidth the min width
   * @param maxWidth the max width
   */
  public MultiFilterWidthComponent(double minWidth, double maxWidth) {
    if (minWidth > 0 && minWidth < 1) {
      this.lowerSigmaThreshold = (float) minWidth;
    } else {
      lowerSigmaThreshold = 0;
    }
    this.upperSigmaThreshold = Filter.getUpperLimit(maxWidth);
  }

  @Override
  public boolean fail(final PreprocessedPeakResult peak) {
    final float xsdf = peak.getXSdFactor();
    return (xsdf > upperSigmaThreshold || xsdf < lowerSigmaThreshold);
  }

  @Override
  public int getType() {
    return FilterValidationFlag.X_SD_FACTOR;
  }
}
