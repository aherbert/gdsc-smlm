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

package uk.ac.sussex.gdsc.smlm.results.count;

/**
 * Stop evaluating when a number of cumulative failures occurs. The failures count is reset to a
 * fraction of the current value for each pass.
 */
public class ResettingFailCounter extends BaseFailCounter {
  /** The fail count. */
  private int failCount;

  /** The number of allowed failures. */
  private final int allowedFailures;

  /** The fraction of the current failures count to reset to for a pass. */
  private final double resetFraction;

  /**
   * Instantiates a new resetting fail counter.
   *
   * @param allowedFailures the number of allowed failures
   * @param resetFraction the reset fraction
   */
  private ResettingFailCounter(int allowedFailures, double resetFraction) {
    this.allowedFailures = allowedFailures;
    this.resetFraction = resetFraction;
  }

  @Override
  protected String generateDescription() {
    return String.format("allowedFailures=%d;resetFraction=%f", allowedFailures, resetFraction);
  }

  /**
   * Instantiates a new resetting fail counter.
   *
   * @param allowedFailures the number of allowed failures
   * @param resetFraction The fraction of the current failures count to reset to for a pass.
   * @return the resetting fail counter
   */
  public static ResettingFailCounter create(int allowedFailures, double resetFraction) {
    if (!(resetFraction >= 0 && resetFraction <= 1)) {
      throw new IllegalArgumentException("Reset must be in the range 0-1");
    }
    return new ResettingFailCounter(Math.max(0, allowedFailures), resetFraction);
  }

  @Override
  public void pass() {
    failCount = (int) (failCount * resetFraction);
  }

  @Override
  public void pass(int n) {
    if (n < 0) {
      throw new IllegalArgumentException("Number of passes must be positive");
    }
    while (n-- > 0) {
      pass();
      if (failCount == 0) {
        break;
      }
    }
  }

  @Override
  public void fail() {
    if (failCount == Integer.MAX_VALUE) {
      throw new IllegalStateException("Unable to increment");
    }
    failCount++;
  }

  @Override
  public void fail(int n) {
    if (n < 0) {
      throw new IllegalArgumentException("Number of fails must be positive");
    }
    if (Integer.MAX_VALUE - n < failCount) {
      throw new IllegalStateException("Unable to increment");
    }
    failCount += n;
  }

  @Override
  public boolean isOk() {
    return failCount <= allowedFailures;
  }

  @Override
  public FailCounter newCounter() {
    return new ResettingFailCounter(allowedFailures, resetFraction);
  }

  @Override
  public void reset() {
    failCount = 0;
  }

  /**
   * Gets the fail count.
   *
   * @return the fail count
   */
  public long getFailCount() {
    return failCount;
  }

  /**
   * Gets the number of allowed failures.
   *
   * @return the number of allowed failures.
   */
  public int getAllowedFailures() {
    return allowedFailures;
  }

  /**
   * Gets the fraction of the current failures count to reset to for a pass.
   *
   * @return the reset fraction
   */
  public double getResetFraction() {
    return resetFraction;
  }
}
