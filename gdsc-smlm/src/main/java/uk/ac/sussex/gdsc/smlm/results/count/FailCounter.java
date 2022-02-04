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

package uk.ac.sussex.gdsc.smlm.results.count;

/**
 * Interface for the evaluation of whether to stop a sequential analysis.
 */
public interface FailCounter {
  /**
   * Gets the description of the fail counter.
   *
   * @return the description (including any parameter values)
   */
  String getDescription();

  /**
   * Called when the most recent event passed.
   */
  void pass();

  /**
   * Called when the n most recent events passed.
   *
   * <p>This method can be used when a series of events are known to pass.
   *
   * @param n the n
   */
  void pass(int n);

  /**
   * Called when the most recent event failed. It is expected that the result of {@link #isOk()} may
   * change after calling this method.
   */
  void fail();

  /**
   * Called when the n most recent event failed. It is expected that the result of {@link #isOk()}
   * may change after calling this method.
   *
   * <p>This method can be used when a series of events are known to fail.
   *
   * @param n the n
   */
  void fail(int n);

  /**
   * Checks if it is ok to continue the analysis. This is set to false when the analysis should
   * stop.
   *
   * @return true, if is ok to continue
   */
  boolean isOk();

  /**
   * Create a duplicate fail counter reset to the initialised state.
   *
   * @return the fail counter
   */
  FailCounter newCounter();

  /**
   * Reset the counter.
   */
  void reset();
}
