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

package uk.ac.sussex.gdsc.smlm.results.count;

import java.util.Objects;

/**
 * Combine the result of two fail counters.
 */
public abstract class CombinedFailCounter extends BaseFailCounter {
  /** The first fail counter . */
  protected final FailCounter c1;
  /** The second fail counter . */
  protected final FailCounter c2;

  /**
   * Instantiates a new combined fail counter.
   *
   * @param c1 the first counter
   * @param c2 the second counter
   */
  public CombinedFailCounter(FailCounter c1, FailCounter c2) {
    this.c1 = Objects.requireNonNull(c1, "Counter 1");
    this.c2 = Objects.requireNonNull(c2, "Counter 2");
  }

  @Override
  protected String generateDescription() {
    final StringBuilder sb = new StringBuilder();
    addDescription(sb, c1);
    sb.append(" ").append(getOperator()).append(" ");
    addDescription(sb, c2);
    return sb.toString();
  }

  private static void addDescription(StringBuilder sb, FailCounter counter) {
    if (counter instanceof CombinedFailCounter) {
      sb.append("(");
      sb.append(counter.getDescription());
      sb.append(")");
    } else {
      sb.append(counter.getDescription());
    }
  }

  /**
   * Get the string representation of the operator used to combine the two fail counters. This is
   * used in the filter name.
   *
   * @return The operator
   */
  protected abstract String getOperator();

  @Override
  public void pass() {
    c1.pass();
    c2.pass();
  }

  @Override
  public void pass(int n) {
    c1.pass(n);
    c2.pass(n);
  }

  @Override
  public void fail() {
    c1.fail();
    c2.fail();
  }

  @Override
  public void fail(int n) {
    c1.fail(n);
    c2.fail(n);
  }

  @Override
  public void reset() {
    c1.reset();
    c2.reset();
  }
}
