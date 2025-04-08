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

package uk.ac.sussex.gdsc.smlm.results.count;

/**
 * Base class for fail counters.
 */
public abstract class BaseFailCounter implements FailCounter {
  private String description;

  @Override
  public String getDescription() {
    if (description == null) {
      description = generateDescription();
    }
    return description;
  }

  /**
   * Generate the description.
   *
   * @return the description
   */
  protected abstract String generateDescription();

  /**
   * Check the count is positive.
   *
   * @param count the count
   * @throws IllegalStateException if not positive
   */
  protected static void checkPositive(int count) {
    if (count < 0) {
      throw new IllegalStateException("Negative count: " + count);
    }
  }

  /**
   * Adds the result.
   *
   * @param pass Set to true if a pass
   */
  public void addResult(boolean pass) {
    if (pass) {
      pass();
    } else {
      fail();
    }
  }
}
