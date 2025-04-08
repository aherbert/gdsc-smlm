/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.engine;

import java.awt.Rectangle;
import java.util.List;
import uk.ac.sussex.gdsc.smlm.fitting.FitResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.filter.MultiPathFitResult;

/**
 * Specifies a job for peak fitting.
 */
public class ParameterisedFitJob extends FitJob {
  private final FitParameters parameters;
  private List<PeakResult> peakResults;
  private int[] indices = new int[0];
  private FitResult[] fitResults;
  private MultiPathFitResult[] multiPathResults;

  /**
   * Constructor with data. Exceptions are thrown if invalid bounds or data are passed.
   *
   * @param id the id
   * @param parameters the parameters
   * @param slice the slice
   * @param data the data
   * @param bounds the bounds
   */
  public ParameterisedFitJob(int id, FitParameters parameters, int slice, float[] data,
      Rectangle bounds) {
    super(id, slice, data, bounds);
    this.parameters = parameters;
  }

  /**
   * Constructor with data. Exceptions are thrown if invalid bounds or data are passed.
   *
   * @param parameters the parameters
   * @param slice the slice
   * @param data the data
   * @param bounds the bounds
   */
  public ParameterisedFitJob(FitParameters parameters, int slice, float[] data, Rectangle bounds) {
    super(slice, slice, data, bounds);
    this.parameters = parameters;
  }

  @Override
  public FitParameters getFitParameters() {
    return parameters;
  }

  @Override
  public void setResults(List<PeakResult> results) {
    this.peakResults = results;
  }

  /**
   * Gets the results.
   *
   * @return The results.
   */
  public List<PeakResult> getResults() {
    return peakResults;
  }

  @Override
  public void setIndices(int[] indices) {
    this.indices = indices;
  }

  @Override
  public void setFitResult(int n, FitResult fitResult) {
    if (fitResults == null) {
      fitResults = new FitResult[indices.length];
    }
    if (n < indices.length) {
      fitResults[n] = fitResult;
    }
  }

  @Override
  public void setMultiPathFitResult(int n, MultiPathFitResult fitResult) {
    if (multiPathResults == null) {
      multiPathResults = new MultiPathFitResult[indices.length];
    }
    if (n < indices.length) {
      multiPathResults[n] = fitResult;
    }
  }

  /**
   * Gets the indices of the data that were fitted.
   *
   * @return The indices of the data that were fitted.
   */
  public int[] getIndices() {
    return indices;
  }

  /**
   * The fit result of the specified index in the array of fitted indices.
   *
   * @param n the index
   * @return the fit result
   */
  public FitResult getFitResult(int n) {
    return fitResults[n];
  }

  /**
   * The fit result of the specified index in the array of fitted indices.
   *
   * @param n the index
   * @return the multi path fit result
   */
  public MultiPathFitResult getMultiPathFitResult(int n) {
    return multiPathResults[n];
  }
}
