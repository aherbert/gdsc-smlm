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
 * Filter results using Width. Assume XY width are different.
 */
public class MultiFilterXyWidthComponent implements MultiFilterComponent {
  private final float lowerSigmaThreshold;
  private final float upperSigmaThreshold;

  /**
   * Instantiates a new multi filter XY width component.
   *
   * @param minWidth the min width
   * @param maxWidth the max width
   */
  public MultiFilterXyWidthComponent(double minWidth, double maxWidth) {
    if (minWidth > 0 && minWidth < 1) {
      lowerSigmaThreshold = (float) (minWidth * minWidth);
    } else {
      lowerSigmaThreshold = 0;
    }
    upperSigmaThreshold = Filter.getUpperLimit(maxWidth * maxWidth);
  }

  @Override
  public boolean fail(final PreprocessedPeakResult peak) {
    final float s2 = peak.getXSdFactor() * peak.getYSdFactor();
    return (s2 > upperSigmaThreshold || s2 < lowerSigmaThreshold);
  }

  @Override
  public int getType() {
    return FilterValidationFlag.X_SD_FACTOR | FilterValidationFlag.Y_SD_FACTOR;
  }
}
