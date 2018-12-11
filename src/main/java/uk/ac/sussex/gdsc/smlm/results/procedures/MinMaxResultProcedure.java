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

package uk.ac.sussex.gdsc.smlm.results.procedures;

import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResultValue;

/**
 * Class to get min/max value in a set of results
 */
public class MinMaxResultProcedure implements PeakResultProcedure {
  private float min;
  private float max;
  private final PeakResultValue value;

  /**
   * Instantiates a new min max result procedure.
   *
   * @param results the results
   * @param value the value
   * @throws IllegalStateException If the results are empty
   */
  public MinMaxResultProcedure(MemoryPeakResults results, PeakResultValue value)
      throws IllegalStateException {
    this.value = value;
    min = max = value.getValue(results.getFirst());
    results.forEach(this);
  }

  /** {@inheritDoc} */
  @Override
  public void execute(PeakResult peakResult) {
    final float v = value.getValue(peakResult);
    if (min > v) {
      min = v;
    } else if (max < v) {
      max = v;
    }
  }

  /**
   * Gets the minimum.
   *
   * @return the minimum
   */
  public float getMinimum() {
    return min;
  }

  /**
   * Gets the maximum.
   *
   * @return the maximum
   */
  public float getMaximum() {
    return max;
  }
}
