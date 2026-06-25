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

package uk.ac.sussex.gdsc.smlm.fitting;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.numbers.gamma.Erfc;
import uk.ac.sussex.gdsc.core.data.VisibleForTesting;
import uk.ac.sussex.gdsc.smlm.math3.analysis.integration.CustomSimpsonIntegrator;

/**
 * Perform diffusion analysis.
 *
 * <p>Contains methods that compute the probability that a diffusing molecule remains within the
 * depth-of-field after a given time. This is based on the Spot-On model described in the methods
 * section of the paper:
 *
 * <p>Hansen, A.S., Woringer, M., Grimm, J.B., Lavis, L.D., Tjian, R., and Darzacq, X. (2018) Robust
 * model-based analysis of single-particle tracking experiments with Spot-On. eLife 7, e33125.
 * doi:10.7554/eLife.33125.
 */
public class DiffusionAnalysis {
  private DiffusionAnalysis() {}

  /**
   * Compute the probability that a molecule remains within bounds {@code [-dz/2, dz/2]} when
   * diffusing with a coefficient {@code d} for time {@code dt}.
   *
   * <p>Note: Parameters are not validated. It is assumed they are all strictly positive and finite.
   * If {@code dt} or {@code d} are zero the probability is 1. If {@code dz} is zero the computation
   * is invalid.
   *
   * @param dt the time delay
   * @param dz the depth of field
   * @param d the diffusion coefficient
   * @return the probability
   */
  public static double remaining(double dt, double dz, double d) {
    // The function is symmetric about z=0.
    // Integrate [0, z/2].
    // Each withinBound function is accurate to 1e-10 and the result is in [0, dz/2].
    // Use a maximum of n=200 evaluations.
    // When the molecule is diffusing fast out of bounds the evaluations is far less
    // using an integrator. If this fails we just use the last known result which is
    // more accurate than a sum of the same number of evenly space points.
    final int n = 200;
    final CustomSimpsonIntegrator in = new CustomSimpsonIntegrator(1e-8, 1e-8, 3, 63);
    final double denom = 1.0 / Math.sqrt(4 * d * dt);
    final UnivariateFunction fun = z -> DiffusionAnalysis.withinBound(z, dz, denom);
    try {
      in.integrate(n, fun, 0, dz * 0.5);
    } catch (TooManyEvaluationsException ignored) {
      // Ignore this.
      // Allow other exceptions to trickle up as they are actual errors in configuration.
    }
    // result in [0, 1]
    return 2 * in.getLastSum() / dz;
  }

  /**
   * Compute the probability that a molecule at position {@code z} remains within bounds
   * {@code [-dz/2, dz/2]} when diffusing with a coefficient {@code d} for time {@code dt}.
   *
   * @param z z position
   * @param dt the time delay (must be strictly positive)
   * @param dz the depth of field
   * @param d the diffusion coefficient
   * @return the probability
   */
  @VisibleForTesting
  static double withinBound(double z, double dt, double dz, double d) {
    final double denom = 1.0 / Math.sqrt(4 * d * dt);
    return withinBound(z, dz, denom);
  }

  /**
   * Compute the probability that a molecule at position {@code z} remains within bounds
   * {@code [-dz/2, dz/2]} when diffusing with a coefficient {@code d} for time {@code dt}.
   *
   * @param z z position (in range [-dz/2, dz/2])
   * @param dz the depth of field
   * @param denom 1 / sqrt(4*d*dt)
   * @return the probability
   */
  private static double withinBound(double z, double dz, double denom) {
    // As per Spot-On:
    // terms of the series evaluated until below 1e-10
    int n = 0;
    double sum = 0;
    double term;
    do {
      // erfc [ ((2n+1)dz / 2) - z ] / sqrt(4*d*dt)
      // erfc [ ((2n+1)dz / 2) + z ] / sqrt(4*d*dt)
      final double lo = (((n + 0.5) * dz) - z) * denom;
      final double hi = (((n + 0.5) * dz) + z) * denom;
      term = Erfc.value(lo) + Erfc.value(hi);
      // Alternating sum: (-1)^n * term
      sum += (n & 1) == 0 ? term : -term;
      n++;
    } while (term > 1e-10);
    // As |z| approaches dz/2 the sum approaches one.
    // Do not allow summation error of the alternating terms
    // of descending magnitude to exceed 1.
    return sum > 1 ? 0 : 1 - sum;
  }
}
