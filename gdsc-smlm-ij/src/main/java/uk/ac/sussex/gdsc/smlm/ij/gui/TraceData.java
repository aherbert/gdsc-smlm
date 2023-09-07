/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Package
 *
 * Software for single molecule localisation microscopy (SMLM) in ImageJ
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.ij.gui;

import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import java.util.Arrays;
import java.util.stream.Collectors;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.Trace;

/**
 * Contains a trace and summary data.
 */
class TraceData {
  private final Trace trace;
  private final double cx;
  private final double cy;
  private final double cz;
  private final double sx;
  private final double sy;
  private final double sz;
  private final double msd;

  /**
   * Create an instance.
   *
   * @param trace the trace
   * @param useZ true to compute values for the z coordinate
   */
  TraceData(Trace trace, boolean useZ) {
    this.trace = trace;
    final double[] weights = getWeights(trace);
    cx = getCentroid(trace, weights, PeakResult.X);
    cy = getCentroid(trace, weights, PeakResult.Y);
    sx = getStandardDeviation(trace, weights, cx, PeakResult.X);
    sy = getStandardDeviation(trace, weights, cy, PeakResult.Y);
    if (useZ) {
      cz = getCentroid(trace, weights, PeakResult.Z);
      sz = getStandardDeviation(trace, weights, cz, PeakResult.Z);
    } else {
      cz = 0;
      sz = 0;
    }
    msd = getMsd(trace);
  }

  /**
   * Gets the weights. These will sum to the number of localisations.
   *
   * @param trace the trace
   * @return the weights
   */
  private static double[] getWeights(Trace trace) {
    final int n = trace.size();
    if (n == 1) {
      return new double[] {1};
    }
    final double[] w = new double[n];
    for (int i = 0; i < n; i++) {
      w[i] = trace.get(i).getIntensity();
    }
    // Make the weights sum to n
    final double norm = n / Arrays.stream(w).sum();
    for (int i = 0; i < n; i++) {
      w[i] *= norm;
    }
    return w;
  }

  private static double getCentroid(Trace trace, double[] weights, int index) {
    final int n = trace.size();
    if (n == 1) {
      return trace.getHead().getParameter(index);
    }
    // Signal weighted
    double c = 0;
    double sum = 0;
    for (int i = 0; i < n; i++) {
      final PeakResult result = trace.get(i);
      final double d = result.getParameter(index);
      final double w = weights[i];
      c += d * w;
      sum += w;
    }
    return c / sum;
  }

  private static double getStandardDeviation(Trace trace, double[] weights, double centre,
      int index) {
    final int n = trace.size();
    if (n == 1) {
      return 0;
    }
    // Signal weighted
    double sum = 0;
    for (int i = 0; i < n; i++) {
      final PeakResult result = trace.get(i);
      final double d = result.getParameter(index) - centre;
      sum += d * d * weights[i];
    }
    return Math.sqrt(sum / (n - 1.0));
  }

  /**
   * Gets the trace.
   *
   * @return the trace
   */
  Trace getTrace() {
    return trace;
  }

  /**
   * Gets the centroid in X.
   *
   * @return the cx
   */
  double getCx() {
    return cx;
  }

  /**
   * Gets the centroid in Y.
   *
   * @return the cy
   */
  double getCy() {
    return cy;
  }

  /**
   * Gets the centroid in Z.
   *
   * @return the cz
   */
  double getCz() {
    return cz;
  }

  /**
   * Gets the standard deviation from the centroid in X.
   *
   * @return the sx
   */
  double getSx() {
    return sx;
  }

  /**
   * Gets the standard deviation from the centroid in Y.
   *
   * @return the sy
   */
  double getSy() {
    return sy;
  }

  /**
   * Gets the standard deviation from the centroid in Z.
   *
   * @return the sz
   */
  double getSz() {
    return sz;
  }

  /**
   * Gets the mean squared deviation (2D).
   *
   * @return the msd
   */
  double getMsd() {
    return msd;
  }

  private static double getMsd(Trace trace) {
    final int n = trace.size();
    if (n == 1) {
      return 0;
    }
    trace.sort();
    double msd = 0;
    PeakResult p = trace.getHead();
    for (int i = 1; i < n; i++) {
      final PeakResult result = trace.get(i);
      msd += result.distance2(p);
      p = result;
    }
    return msd / (n - 1.0);
  }

  /**
   * Gets the category string. This is a sorted set of the category IDs in the trace.
   *
   * @return the category string
   */
  String getCategoryString() {
    final int n = trace.size();
    int cat1 = trace.getHead().getCategory();
    for (int i = 1; i < n; i++) {
      final PeakResult result = trace.get(i);
      if (result.getCategory() != cat1) {
        // Multiple categories
        IntOpenHashSet set = new IntOpenHashSet(Math.max(16, n / 2));
        set.add(cat1);
        set.add(result.getCategory());
        for (int j = i + 1; j < n; j++) {
          set.add(trace.get(j).getCategory());
        }
        return set.intStream().sorted().mapToObj(Integer::toString)
            .collect(Collectors.joining(", "));
      }
    }
    return Integer.toString(cat1);
  }
}
