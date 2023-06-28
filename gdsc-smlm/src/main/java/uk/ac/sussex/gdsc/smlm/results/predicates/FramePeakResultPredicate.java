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

package uk.ac.sussex.gdsc.smlm.results.predicates;

import java.util.function.Predicate;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Test a result using the frame.
 */
public class FramePeakResultPredicate implements Predicate<PeakResult> {
  /** The frame. */
  private final int frame;

  /**
   * Instantiates a new frame peak result predicate.
   *
   * @param frame the frame
   */
  public FramePeakResultPredicate(int frame) {
    this.frame = frame;
  }

  @Override
  public boolean test(PeakResult peakResult) {
    return peakResult.getFrame() == frame;
  }
}
