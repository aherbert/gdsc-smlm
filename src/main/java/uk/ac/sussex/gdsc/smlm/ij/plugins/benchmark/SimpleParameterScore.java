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

package uk.ac.sussex.gdsc.smlm.ij.plugins.benchmark;

import uk.ac.sussex.gdsc.smlm.results.filter.DirectFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.FilterScore;

/**
 * Store the parameter score used in benchmarking.
 *
 * <p>Note: this class has a natural ordering that is inconsistent with equals.
 */
public class SimpleParameterScore extends FilterScore {
  /** The result. */
  final ParameterScoreResult result;

  /**
   * Instantiates a new simple parameter score.
   *
   * @param filter the filter
   * @param result the result
   * @param criteriaPassed the criteria passed
   */
  public SimpleParameterScore(DirectFilter filter, ParameterScoreResult result,
      boolean criteriaPassed) {
    super(filter, result.score, result.criteria, true, criteriaPassed);
    this.result = result;
  }

  @Override
  protected int compareParameters(FilterScore that) {
    // Compare the parameters and return the strongest, those most likely to restrict the output

    // 0 = failCount
    // 1 = residudalsThreshold
    // 2 = duplicateDistance
    final double[] p1 = this.result.parameters;
    final double[] p2 = ((SimpleParameterScore) that).result.parameters;

    // Lowest fail count
    if (p1[0] < p2[0]) {
      return -1;
    }
    if (p1[0] > p2[0]) {
      return 1;
    }

    // Lowest duplicate distance
    if (p1[2] < p2[2]) {
      return -1;
    }
    if (p1[2] > p2[2]) {
      return 1;
    }

    // Highest residuals threshold
    if (p1[2] > p2[2]) {
      return -1;
    }
    if (p1[2] < p2[2]) {
      return 1;
    }

    return 0;
  }
}
