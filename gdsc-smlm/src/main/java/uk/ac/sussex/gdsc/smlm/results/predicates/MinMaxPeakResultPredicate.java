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

package uk.ac.sussex.gdsc.smlm.results.predicates;

import java.util.Objects;
import java.util.function.Predicate;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResultValue;

/**
 * Accept results with a value within a min/max range.
 */
public class MinMaxPeakResultPredicate implements Predicate<PeakResult> {
  /** The min of the value range. */
  public final float min;

  /** The max of the value range. */
  public final float max;

  /** The value. */
  public final PeakResultValue value;

  /**
   * Instantiates a new min max peak result predicate.
   *
   * @param min the min of the value range
   * @param max the max of the value range
   * @param value the value
   * @throws IllegalArgumentException If the min/max are NaN, or if min is greater than max.
   */
  public MinMaxPeakResultPredicate(float min, float max, PeakResultValue value) {
    if (Float.isNaN(min) || Float.isNaN(max) || min > max) {
      throw new IllegalArgumentException("Min/Max range is invalid");
    }
    this.min = min;
    this.max = max;
    this.value = Objects.requireNonNull(value, "value");
  }

  @Override
  public boolean test(PeakResult peakResult) {
    final float v = value.getValue(peakResult);
    return v >= min && v <= max;
  }
}
