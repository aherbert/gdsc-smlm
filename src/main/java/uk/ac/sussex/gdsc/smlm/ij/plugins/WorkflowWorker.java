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

package uk.ac.sussex.gdsc.smlm.ij.plugins;

import uk.ac.sussex.gdsc.smlm.utils.Pair;

/**
 * A worker to update results based on new settings.
 *
 * @author Alex Herbert
 * @param <S> the generic type
 * @param <R> the generic type
 */
public abstract class WorkflowWorker<S, R> {
  /**
   * Compare the settings and return false if any settings that the work depends on have changed.
   *
   * <p>Both objects will not be null.
   *
   * @param current the current
   * @param previous the previous
   * @return true if settings have changed
   */
  public abstract boolean equalSettings(S current, S previous);

  /**
   * Compare the results and return false if any results that the work depends on have changed.
   *
   * <p>Either object could be null (if no results have yet been generated for this work).
   *
   * @param current the current
   * @param previous the previous
   * @return true if results have changed
   */
  public abstract boolean equalResults(R current, R previous);

  /**
   * Creates the results.
   *
   * @param work the work (the current settings and results)
   * @return the updated settings and results
   */
  public abstract Pair<S, R> doWork(Pair<S, R> work);

  /**
   * Called when there are new results in the current work. This can be used to reset the worker
   * before {@link #doWork(Pair)} is called.
   */
  protected void newResults() {
    // Do nothing
  }
}
