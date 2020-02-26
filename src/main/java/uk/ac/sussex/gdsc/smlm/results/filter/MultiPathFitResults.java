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

package uk.ac.sussex.gdsc.smlm.results.filter;

/**
 * Specifies the result of fitting a frame using different fitting methods.
 *
 * <p>The multi-path results can be evaluated by the MultiPathFilter to determine which result from
 * the different paths should be accepted.
 *
 * <p>This class is used for benchmarking the fitting path options in the PeakFit algorithm.
 */
public class MultiPathFitResults implements IMultiPathFitResults {
  /** The frame containing the results. */
  private final int frame;

  /** The list of multi-path results. */
  private final MultiPathFitResult[] multiPathFitResultList;

  /**
   * The total number of candidates.
   */
  private final int totalCandidates;

  /**
   * The number of actual results in the frame.
   */
  private final int numberOfActualResults;

  /**
   * Instantiates a new multi path fit results.
   *
   * @param frame the frame
   * @param multiPathFitResults the multi path fit results
   */
  public MultiPathFitResults(int frame, MultiPathFitResult[] multiPathFitResults) {
    this(frame, multiPathFitResults, (multiPathFitResults == null) ? 0 : multiPathFitResults.length,
        0);
  }

  /**
   * Instantiates a new multi path fit results.
   *
   * @param frame the frame
   * @param multiPathFitResults the multi path fit results
   * @param totalCandidates the total candidates
   * @param numberOfActualResults the number of actual results in the frame
   */
  public MultiPathFitResults(int frame, MultiPathFitResult[] multiPathFitResults,
      int totalCandidates, int numberOfActualResults) {
    this.frame = frame;
    this.multiPathFitResultList = multiPathFitResults;
    this.totalCandidates = totalCandidates;
    this.numberOfActualResults = numberOfActualResults;
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected MultiPathFitResults(MultiPathFitResults source) {
    frame = source.frame;
    totalCandidates = source.totalCandidates;
    numberOfActualResults = source.numberOfActualResults;

    final MultiPathFitResult[] list = source.multiPathFitResultList;
    if (list != null) {
      multiPathFitResultList = new MultiPathFitResult[list.length];
      for (int i = 0; i < list.length; i++) {
        // Do a shallow copy. This matches the previous clone() functionality.
        multiPathFitResultList[i] = list[i].copy(false);
      }
    } else {
      multiPathFitResultList = null;
    }
  }

  /**
   * Create a copy.
   *
   * @return the copy
   */
  public MultiPathFitResults copy() {
    return new MultiPathFitResults(this);
  }

  @Override
  public int getFrame() {
    return frame;
  }

  @Override
  public int getNumberOfResults() {
    return multiPathFitResultList.length;
  }

  @Override
  public MultiPathFitResult getResult(int index) {
    return multiPathFitResultList[index];
  }

  @Override
  public void complete(int index) {
    // Do nothing
  }

  @Override
  public int getTotalCandidates() {
    return totalCandidates;
  }

  /**
   * Gets the multi path fit results.
   *
   * @return the multi path fit results
   */
  public MultiPathFitResult[] getMultiPathFitResults() {
    return multiPathFitResultList;
  }

  /**
   * Gets the number of actual results in the frame. Used during filter scoring.
   *
   * @return the number of actual results
   */
  public int getNumberOfActualResults() {
    return numberOfActualResults;
  }
}
