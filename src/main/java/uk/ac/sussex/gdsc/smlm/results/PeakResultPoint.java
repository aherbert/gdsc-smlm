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

package uk.ac.sussex.gdsc.smlm.results;

import java.util.Arrays;
import java.util.Objects;
import uk.ac.sussex.gdsc.core.match.BasePoint;

/**
 * A point that holds a reference to a PeakResult.
 */
public class PeakResultPoint extends BasePoint {
  /** The time. */
  private final int time;

  /** The peak result. */
  final PeakResult peakResult;

  /**
   * Instantiates a new peak result point.
   *
   * @param time the time
   * @param x the x
   * @param y the y
   * @param z the z
   * @param peakResult the peak result
   */
  public PeakResultPoint(int time, float x, float y, float z, PeakResult peakResult) {
    super(x, y, z);
    this.time = time;
    this.peakResult = peakResult;
  }

  /**
   * Gets the time.
   *
   * @return the time
   */
  public int getTime() {
    return time;
  }

  /**
   * Gets the peak result.
   *
   * @return the peak result
   */
  public PeakResult getPeakResult() {
    return peakResult;
  }

  @Override
  public boolean equals(Object object) {
    if (super.equals(object)) {
      // Super-class ensures the class is the same.
      // Compare new fields.
      final PeakResultPoint other = (PeakResultPoint) object;
      return time == other.time && peakResult == other.peakResult;
    }
    return false;
  }

  @Override
  public int hashCode() {
    return Arrays.hashCode(new int[] {super.hashCode(), time, Objects.hashCode(peakResult)});
  }
}
