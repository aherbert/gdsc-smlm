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

package uk.ac.sussex.gdsc.smlm.fitting;

// TODO: Auto-generated Javadoc
/**
 * Define the status of a fit result.
 */
public enum FitStatus {
  /** OK. */
  OK("OK"),
  /** Singular non-linear model. */
  SINGULAR_NON_LINEAR_MODEL("Singular non-linear model"),
  /** Singular non-linear solution. */
  SINGULAR_NON_LINEAR_SOLUTION("Singular non-linear solution"),
  /** Invalid gradients. */
  INVALID_GRADIENTS("Invalid gradients"),
  /** Failed to converge. */
  FAILED_TO_CONVERGE("Failed to converge"),
  /** Too many iterations. */
  TOO_MANY_ITERATIONS("Too many iterations"),
  /** Too many evaluations. */
  TOO_MANY_EVALUATIONS("Too many evaluations"),
  /** Invalid likelihood. */
  INVALID_LIKELIHOOD("Invalid likelihood"),
  /** Bad parameters. */
  BAD_PARAMETERS("Bad parameters"),
  /** Failed to estimate width. */
  FAILED_TO_ESTIMATE_WIDTH("Failed to estimate width"),
  /** Coordinates moved. */
  COORDINATES_MOVED("Coordinates moved"),
  /** Outside fit region. */
  OUTSIDE_FIT_REGION("Outside fit region"),
  /** Insufficient signal. */
  INSUFFICIENT_SIGNAL("Insufficient signal"),
  /** Insufficient Signal-to-Noise Ratio (SNR). */
  INSUFFICIENT_SNR("Insufficient SNR"),
  /** Width diverged. */
  WIDTH_DIVERGED("Width diverged"),
  /** Z-coordinate moved. */
  Z_MOVED("Z-coordinate moved"),
  /** Insufficient precision. */
  INSUFFICIENT_PRECISION("Insufficient precision"),
  /** Neighbour overlap. */
  NEIGHBOUR_OVERLAP("Neighbour overlap"),
  /** Failed smart filter. */
  FAILED_SMART_FILTER("Failed smart filter"),
  /** Drift to another result. */
  DRIFT_TO_ANOTHER_RESULT("Drift to another result"),
  /** Failed validation. */
  FAILED_VALIDATION("Failed validation"),
  /** No model improvement. */
  NO_MODEL_IMPROVEMENT("No model improvement"),
  /** Line search error. */
  LINE_SEARCH_ERROR("Line search error"),
  /** Unknown. */
  UNKNOWN("Unknown");

  /** The nice name. */
  private String niceName;

  /**
   * Instantiates a new fit status.
   *
   * @param name the name
   */
  FitStatus(String name) {
    niceName = name;
  }

  @Override
  public String toString() {
    return getName();
  }

  /**
   * Gets the name.
   *
   * @return the name
   */
  public String getName() {
    return niceName;
  }
}
