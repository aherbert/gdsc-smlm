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

import uk.ac.sussex.gdsc.smlm.filters.Spot;
import uk.ac.sussex.gdsc.smlm.results.filter.MultiPathFilter;

/**
 * Specifies additional parameters for the job.
 *
 * <p>Can be used to collect additional information during fitting.
 */
public class FitParameters {
  /**
   * The fit task.
   */
  public enum FitTask {
    //@formatter:off
    /** PSF fitting. */
    PSF_FITTING{ @Override
    public String getName() { return "PSF Fitting"; }},
    /** Maxima identification. */
    MAXIMA_IDENITIFICATION{ @Override
    public String getName() { return "Maxima Identification"; }},
    /** Benchmarking. */
    BENCHMARKING{ @Override
    public String getName() { return "Benchmarking"; }};
    //@formatter:on

    @Override
    public String toString() {
      return getName();
    }

    /**
     * Gets the name.
     *
     * @return the name
     */
    public abstract String getName();
  }

  /**
   * The noise for the image data.
   */
  public float noise = Float.NaN;

  /**
   * The spots to fit within the data.
   */
  public Spot[] spots;

  /**
   * The maximum candidate spot to fit. This should be equal to spots.length or less. It is used
   * when additional candidates have been added to the spots list that are neighbours of the primary
   * spot candidates.
   */
  public int maxCandidate;

  /**
   * The maxima to fit within the data.
   */
  public int[] maxIndices;

  /**
   * The task to perform.
   */
  public FitTask fitTask = FitTask.PSF_FITTING;

  /**
   * The frame acquisition end time for the input data. Used when data represents multiple frames.
   */
  public int endT = -1;

  /**
   * The filter used to pick the fitting path when benchmarking. This should be an instance for each
   * thread running fit jobs.
   *
   * <p>Note that during benchmarking all fitting paths will be computed. The current set of results
   * is then built by validating the results with this filter (in addition to the fit configuration
   * used to construct the FitWorker).
   */
  public MultiPathFilter benchmarkFilter;

  /**
   * The pass array.
   *
   * <p>If this is not null then it will be initialised to the length of the candidate list and used
   * to store the pass/fail flag for each candidate up to the point that fitting was halted.
   */
  public boolean[] pass;
}
