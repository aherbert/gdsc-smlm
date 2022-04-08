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

package uk.ac.sussex.gdsc.smlm.results;

import uk.ac.sussex.gdsc.core.annotation.NotNull;

/**
 * Contains the options to set before reading the results.
 */
public class ResultOption {
  /**
   * An empty immutable {@code ResultOption} array.
   */
  @NotNull
  public static final ResultOption[] EMPTY_ARRAY = new ResultOption[0];

  /** The id. */
  public final int id;

  /** The name. */
  public final String name;

  /** The set of valid values. This can be null. */
  public final Object[] values;

  /** The value. */
  private Object value;

  /**
   * Creates a new result option.
   *
   * @param id the id
   * @param name the name
   * @param value the value
   * @param values the values
   */
  ResultOption(int id, String name, Object value, Object[] values) {
    this.id = id;
    this.name = name;
    this.values = values;
    checkValue(value, values);
    this.value = value;
  }

  /**
   * Gets the value.
   *
   * @return the value
   */
  public Object getValue() {
    return value;
  }

  /**
   * Sets the value. If the list of valid values is not empty then an exception will be thrown if
   * the value is one of the valid values.
   *
   * @param value the new value
   * @throws IllegalArgumentException If the value is not in the list of valid values
   */
  public void setValue(Object value) {
    checkValue(value, values);
    this.value = value;
  }

  /**
   * Check the value is equal to a valid value.
   *
   * <p>If the set of valid values is null or empty then no exception is raised (any value is
   * valid).
   *
   * @param value the value
   * @param values the values
   * @throws IllegalArgumentException If the value is not in the list of valid values
   */
  static void checkValue(Object value, Object[] values) {
    if (values != null && values.length > 0) {
      for (int i = 0; i < values.length; i++) {
        if (values[i].equals(value)) {
          return;
        }
      }
      throw new IllegalArgumentException("Not a valid value: " + value);
    }
  }

  /**
   * Checks for valid values.
   *
   * @return true, if the list of valid values is not empty
   */
  public boolean hasValues() {
    return values != null && values.length > 0;
  }
}
