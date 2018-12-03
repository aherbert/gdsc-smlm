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
package uk.ac.sussex.gdsc.smlm.results.filter;

/**
 * Specifies the result of fitting a frame using different fitting methods. <p> The multi-path
 * results can be evaluated by the MultiPathFilter to determine which result from the different
 * paths should be accepted.
 */
public interface IMultiPathFitResults {
  /**
   * @return The frame containing the results.
   */
  public int getFrame();

  /**
   * Get the number of results. The {@link #getResult(int)} method should support being called with
   * any index up to the number of results (exclusive).
   *
   * @return The number of results
   */
  public int getNumberOfResults();

  /**
   * Gets the result.
   *
   * @param index the index
   * @return the result
   */
  public MultiPathFitResult getResult(int index);

  /**
   * Called when the results that would be returned by {@link #getResult(int)} are no longer
   * required
   *
   * @param index the index
   */
  public void complete(int index);

  /**
   * The total number of candidates. This may be greater than the size of the
   * {@link #getNumberOfResults()} if this is a subset of the results, i.e. has been prefiltered.
   *
   * @return the total candidates
   */
  public int getTotalCandidates();

  // Possible support for iteration
  // /**
  // * Begin. Called before a pass through the results using {@link #getResult(int)}.
  // *
  // * @return true, if a pass through the results is possible
  // */
  // public boolean begin();
  //
  // /**
  // * Called after a pass through the results. Returns a boolean indicating a repeat is possible.
  // *
  // * @return true, if another pass through the results is possible
  // */
  // public boolean end();
}
