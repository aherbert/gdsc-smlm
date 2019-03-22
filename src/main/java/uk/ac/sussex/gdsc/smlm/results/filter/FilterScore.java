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

package uk.ac.sussex.gdsc.smlm.results.filter;

/**
 * Store the filter score used in benchmarking.
 *
 * <p>Note: this class has a natural ordering that is inconsistent with equals.
 */
public class FilterScore implements Comparable<FilterScore> {
  /** The filter. */
  public final Filter filter;

  /** The score. */
  public final double score;

  /** The criteria. */
  public final double criteria;

  /** Flag to indicate if the criteria passed. */
  public final boolean criteriaPassed;

  /** Flag to indicate if the filters are all same type. */
  public final boolean allSameType;

  /**
   * Instantiates a new filter score.
   *
   * @param filter the filter
   * @param score the score
   * @param criteria the criteria
   * @param allSameType the all same type
   * @param criteriaPassed the criteria passed
   */
  public FilterScore(Filter filter, double score, double criteria, boolean allSameType,
      boolean criteriaPassed) {
    this.filter = filter;
    this.score = score;
    this.criteria = criteria;
    this.allSameType = allSameType;
    this.criteriaPassed = criteriaPassed;
  }

  @Override
  public int compareTo(FilterScore that) {
    if (that == null) {
      return -1;
    }

    if (this.criteriaPassed) {
      // Must pass criteria first
      if (!that.criteriaPassed) {
        return -1;
      }

      // Sort by the score then criteria, highest first
      int result = Double.compare(that.score, this.score);
      if (result != 0) {
        return result;
      }
      result = Double.compare(that.criteria, this.criteria);
      if (result != 0) {
        return result;
      }

      // If the same type then compare the parameters
      if (allSameType || this.filter.getType().equals(that.filter.getType())) {
        return compareParameters(that);
      }
      return 0;
    }

    // Must pass criteria first
    if (that.criteriaPassed) {
      return 1;
    }

    // Sort by how close we are to passing the criteria, highest first
    int result = Double.compare(that.criteria, this.criteria);
    if (result != 0) {
      return result;
    }
    result = Double.compare(that.score, this.score);
    if (result != 0) {
      return result;
    }

    // If the same type then compare the parameters
    if (allSameType || this.filter.getType().equals(that.filter.getType())) {
      return compareParameters(that);
    }
    return 0;
  }

  /**
   * Compare the parameters to the other score to count the number of strongest parameters.
   *
   * @param that the other filter score
   * @return the count difference
   */
  protected int compareParameters(FilterScore that) {
    // Get the filter with the strongest params
    return that.filter.weakestUnsafe(this.filter);
  }

  @Override
  public String toString() {
    // Add the score
    return String.format("%s : %.3f (%.3f)", filter.getName(), score, criteria);
  }
}
