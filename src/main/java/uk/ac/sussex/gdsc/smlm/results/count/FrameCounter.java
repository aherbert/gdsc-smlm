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

package uk.ac.sussex.gdsc.smlm.results.count;

/**
 * Class to count per frame.
 */
public class FrameCounter extends Counter {
  private int current;
  private int previous;

  /**
   * Instantiates a new frame counter with the current value of -1.
   */
  public FrameCounter() {
    this(-1);
  }

  /**
   * Instantiates a new frame counter with the given frame.
   *
   * @param frame the current frame
   */
  public FrameCounter(int frame) {
    current = previous = frame;
  }

  /**
   * Advance the frame. If this is a new frame then the count is reset.
   *
   * @param frame the frame
   * @return true if the frame has changed
   */
  public boolean advanceAndReset(int frame) {
    if (current != frame) {
      previous = current;
      current = frame;
      reset();
      return true;
    }
    return false;
  }

  /**
   * Advance the frame.
   *
   * @param frame the frame
   * @return true if the frame has changed
   */
  public boolean advance(int frame) {
    if (current != frame) {
      previous = current;
      current = frame;
      return true;
    }
    return false;
  }

  /**
   * Get the current frame.
   *
   * @return the current frame
   */
  public int currentFrame() {
    return current;
  }

  /**
   * Gets the previous frame.
   *
   * <p>Note: If the counter has never been advanced to a new frame then the previous will be the
   * same as the current.
   *
   * @return the previous frame
   */
  public int previousFrame() {
    return previous;
  }
}
