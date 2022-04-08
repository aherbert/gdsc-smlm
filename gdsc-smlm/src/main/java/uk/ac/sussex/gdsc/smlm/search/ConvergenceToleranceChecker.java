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

package uk.ac.sussex.gdsc.smlm.search;

/**
 * Check if converged using a tolerance on the score and/or position change, and the number of
 * iterations.
 *
 * @param <T> the generic type
 */
public class ConvergenceToleranceChecker<T extends Comparable<T>> implements ConvergenceChecker<T> {
  /** The relative tolerance threshold. */
  public final double relative;
  /** The absolute tolerance threshold. */
  public final double absolute;
  /** The check score flag. */
  public final boolean checkScore;
  /** The check sequence flag. */
  public final boolean checkSequence;
  /** The max iterations. */
  public final int maxIterations;

  private int iterations;

  /**
   * Build an instance with specified thresholds. This only check convergence using the score.
   *
   * <p>In order to perform only relative checks, the absolute tolerance must be set to a negative
   * value. In order to perform only absolute checks, the relative tolerance must be set to a
   * negative value.
   *
   * <p>Warning: No checks are made to validate that convergence is possible.
   *
   * @param relative relative tolerance threshold
   * @param absolute absolute tolerance threshold
   */
  public ConvergenceToleranceChecker(double relative, double absolute) {
    this(relative, absolute, true, false, 0);
  }

  /**
   * Build an instance with specified thresholds.
   *
   * <p>In order to perform only relative checks, the absolute tolerance must be set to a negative
   * value. In order to perform only absolute checks, the relative tolerance must be set to a
   * negative value.
   *
   * <p>Warning: No checks are made to validate that convergence is possible.
   *
   * @param relative relative tolerance threshold
   * @param absolute absolute tolerance threshold
   * @param checkScore Set to true to check the score
   * @param checkSequence Set to true to check the position
   * @param maxIterations Set above zero to check the iterations (number of times
   *        {@link #converged(SearchResult, SearchResult)} is called)
   */
  public ConvergenceToleranceChecker(double relative, double absolute, boolean checkScore,
      boolean checkSequence, int maxIterations) {
    if (maxIterations < 0) {
      maxIterations = 0;
    }

    this.relative = relative;
    this.absolute = absolute;
    this.checkScore = checkScore;
    this.checkSequence = checkSequence;
    this.maxIterations = maxIterations;
  }

  /**
   * Check if the position has converged.
   *
   * @param previous Previous
   * @param current Current
   * @return True if converged
   */
  private boolean converged(final double[] previous, final double[] current) {
    for (int i = 0; i < previous.length; ++i) {
      if (!converged(previous[i], current[i])) {
        return false;
      }
    }
    return true;
  }

  /**
   * Check if the value has converged.
   *
   * @param previous Previous
   * @param current Current
   * @return True if converged
   */
  private boolean converged(final double previous, final double current) {
    final double difference = Math.abs(previous - current);
    if (difference <= absolute) {
      return true;
    }
    final double size = Math.max(Math.abs(previous), Math.abs(current));
    return (difference <= size * relative);
  }

  @Override
  public boolean converged(SearchResult<T> previous, SearchResult<T> current) {
    iterations++;
    if (maxIterations != 0 && iterations >= maxIterations) {
      return true;
    }
    if (checkScore && converged(previous.getScore(), current.getScore())) {
      return true;
    }
    return (checkSequence && converged(previous.getPoint(), current.getPoint()));
  }

  private boolean converged(T score, T score2) {
    return score.compareTo(score2) == 0;
  }

  /**
   * Gets the iterations.
   *
   * @return the iterations.
   */
  public int getIterations() {
    return iterations;
  }
}
