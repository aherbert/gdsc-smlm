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

package uk.ac.sussex.gdsc.smlm.math3.optim;

import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.util.FastMath;

/**
 * Check if the coordinates have converged.
 */
public abstract class CoordinateChecker
    implements OptimizationData, ConvergenceChecker<PointValuePair> {
  /** Relative tolerance threshold. */
  final double relative;
  /** Absolute tolerance threshold. */
  final double absolute;

  /**
   * Build an instance with specified thresholds.
   *
   * <p>In order to perform only relative checks, the absolute tolerance must be set to a negative
   * value. In order to perform only absolute checks, the relative tolerance must be set to a
   * negative value.
   *
   * @param relative relative tolerance threshold
   * @param absolute absolute tolerance threshold
   */
  public CoordinateChecker(double relative, double absolute) {
    this(relative, absolute, Integer.MAX_VALUE);
  }

  /**
   * Build an instance with specified thresholds.
   *
   * <p>In order to perform only relative checks, the absolute tolerance must be set to a negative
   * value. In order to perform only absolute checks, the relative tolerance must be set to a
   * negative value.
   *
   * @param relative the relative
   * @param absolute the absolute
   * @param fixedIterations the fixed number of iterations to signal convergence
   */
  public CoordinateChecker(double relative, double absolute, int fixedIterations) {
    this.relative = relative;
    this.absolute = absolute;
  }

  /**
   * Check if the coordinates have converged.
   *
   * @param p Previous
   * @param c Current
   * @return True if converged
   */
  public boolean converged(final double[] p, final double[] c) {
    for (int i = 0; i < p.length; ++i) {
      final double pi = p[i];
      final double ci = c[i];
      final double difference = Math.abs(pi - ci);
      final double size = FastMath.max(Math.abs(pi), Math.abs(ci));
      if (difference > size * relative && difference > absolute) {
        return false;
      }
    }
    return true;
  }

  /** {@inheritDoc} */
  @Override
  public boolean converged(int iteration, PointValuePair previous, PointValuePair current) {
    return converged(previous.getPointRef(), current.getPointRef());
  }
}
