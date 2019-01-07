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

package uk.ac.sussex.gdsc.smlm.results.procedures;

import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Contains core functionality to for result procedures.
 */
public abstract class AbstractResultProcedure {
  /** The results. */
  final MemoryPeakResults results;

  /** The counter for procedures. */
  protected int counter;

  /**
   * Instantiates a new abstract result procedure.
   *
   * @param results the results
   */
  public AbstractResultProcedure(MemoryPeakResults results) {
    if (results == null) {
      throw new IllegalArgumentException("results must not be null");
    }
    this.results = results;
  }

  /**
   * Get the size of the results.
   *
   * @return the size
   */
  public int size() {
    return results.size();
  }

  /**
   * Allocate the array to store the data based on the current results size.
   *
   * @param data the data
   * @return the (new) data
   */
  protected int[] allocate(int[] data) {
    return (data == null || data.length < size()) ? new int[size()] : data;
  }

  /**
   * Allocate the array to store the data based on the current results size.
   *
   * @param data the data
   * @return the (new) data
   */
  protected float[] allocate(float[] data) {
    return (data == null || data.length < size()) ? new float[size()] : data;
  }

  /**
   * Allocate the array to store the data based on the current results size.
   *
   * @param data the data
   * @return the (new) data
   */
  protected double[] allocate(double[] data) {
    return (data == null || data.length < size()) ? new double[size()] : data;
  }

  /**
   * Allocate the array to store the data based on the current results size.
   *
   * @param data the data
   * @return the (new) data
   */
  protected PeakResult[] allocate(PeakResult[] data) {
    return (data == null || data.length < size()) ? new PeakResult[size()] : data;
  }
}
