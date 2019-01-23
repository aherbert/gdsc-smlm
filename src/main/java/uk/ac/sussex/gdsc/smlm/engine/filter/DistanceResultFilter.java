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

package uk.ac.sussex.gdsc.smlm.engine.filter;

import uk.ac.sussex.gdsc.smlm.fitting.FitResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Filter the results using the distance to a set of coordinates. Any fitted position within the
 * distance to the target coordinates is accepted.
 *
 * @deprecated Filtering of the results is no longer supported
 */
@Deprecated
public class DistanceResultFilter extends ResultFilter {
  /**
   * Instantiates a new distance result filter.
   *
   * @param filter the filter
   * @param distance the distance
   * @param numberOfMaxima the number of maxima
   */
  public DistanceResultFilter(List<float[]> filter, float distance, int numberOfMaxima) {
    super(filter, distance, numberOfMaxima);
    filteredFitResults = new FitResult[numberOfMaxima];
    filteredIndices = new int[numberOfMaxima];
    peakResults = new ArrayList<>(numberOfMaxima);
  }

  @Override
  public void filter(FitResult fitResult, int maxIndex, PeakResult... results) {
    boolean found = false;
    for (final PeakResult result : results) {
      if (result == null) {
        continue;
      }
      for (final float[] coord : filter) {
        final float dx = result.getXPosition() - coord[0];
        final float dy = result.getYPosition() - coord[1];
        if (dx * dx + dy * dy < d2) {
          found = true;
          peakResults.add(result);
          break;
        }
      }
    }
    if (found) {
      // Add the result and the fitted index to the filtered results
      filteredFitResults[filteredCount] = fitResult;
      filteredIndices[filteredCount] = maxIndex;
      filteredCount++;
    }
  }

  @Override
  public void filter(FitResult fitResult, int maxIndex, float x, float y) {
    boolean found = false;
    for (final float[] coord : filter) {
      final float dx = x - coord[0];
      final float dy = y - coord[1];
      if (dx * dx + dy * dy < d2) {
        found = true;
        break;
      }
    }
    if (found) {
      // Add the result and the fitted index to the filtered results
      filteredFitResults[filteredCount] = fitResult;
      filteredIndices[filteredCount] = maxIndex;
      filteredCount++;
    }
  }

  @Override
  public void finalise() {
    filteredFitResults = Arrays.copyOf(filteredFitResults, filteredCount);
    filteredIndices = Arrays.copyOf(filteredIndices, filteredCount);
  }
}
