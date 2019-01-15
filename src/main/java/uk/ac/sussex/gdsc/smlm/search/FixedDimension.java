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

package uk.ac.sussex.gdsc.smlm.search;

import uk.ac.sussex.gdsc.core.utils.MathUtils;

/**
 * Specify the dimensions for a search.
 */
public class FixedDimension implements Cloneable, Dimension {
  /** The minimum of the range. */
  public final double min;

  /** The maximum of the range. */
  public final double max;

  /** The current lower bound of the range (will be clipped to min./max). */
  public final double lower;

  /** The current upper bound of the range (will be clipped to min./max). */
  public final double upper;

  /** The min increment to use around the centre. */
  public final double minIncrement;

  /**
   * Set to true if {@link #min} &lt; {@link #max}.
   */
  public final boolean active;

  /**
   * Instantiates a new inactive search dimension. The centre is set to zero.
   */
  public FixedDimension() {
    this(0);
  }

  /**
   * Instantiates a new inactive search dimension. The centre can be set to any value.
   *
   * @param centre the centre
   */
  public FixedDimension(double centre) {
    this(centre, centre, 0);
  }

  /**
   * Instantiates a new search dimension.
   *
   * @param min the minimum of the range
   * @param max the maximum of the range
   * @param minIncrement the min increment to use around the centre
   */
  public FixedDimension(double min, double max, double minIncrement) {
    this(min, max, minIncrement, min, max);
  }

  /**
   * Instantiates a new search dimension.
   *
   * @param min the minimum of the range
   * @param max the maximum of the range
   * @param minIncrement the min increment to use around the centre
   * @param lower the current lower bound of the range (will be clipped to min/max)
   * @param upper the current upper bound of the range (will be clipped to min/max)
   */
  public FixedDimension(double min, double max, double minIncrement, double lower, double upper) {
    if (isInvalid(min)) {
      throw new IllegalArgumentException("Min is not a valid number: " + min);
    }
    if (isInvalid(max)) {
      throw new IllegalArgumentException("Max is not a valid number: " + max);
    }
    if (isInvalid(lower)) {
      throw new IllegalArgumentException("Lower is not a valid number: " + lower);
    }
    if (isInvalid(upper)) {
      throw new IllegalArgumentException("Upper is not a valid number: " + upper);
    }
    if (isInvalid(minIncrement)) {
      throw new IllegalArgumentException("Min increment is not a valid number: " + minIncrement);
    }
    if (max < min) {
      throw new IllegalArgumentException("Max is less than min");
    }
    this.active = min < max;
    if (minIncrement < 0) {
      throw new IllegalArgumentException("Min increment is negative: " + minIncrement);
    }
    if (upper < lower) {
      throw new IllegalArgumentException("Upper is less than lower");
      // if (upper < min || upper > max)
      // throw new IllegalArgumentException("Upper is outside min/max range");
      // if (lower < min || lower > max)
      // throw new IllegalArgumentException("Lower is outside min/max range");
    }

    // Clip to range
    lower = Math.min(max, Math.max(lower, min));
    upper = Math.min(max, Math.max(upper, min));

    // We round to the min increment so that the values returned should be identical if the centre
    // is moved by a factor of the increment.
    this.minIncrement = minIncrement;
    this.min = round(min);
    this.max = round(max);
    this.lower = round(lower);
    this.upper = round(upper);
  }

  /**
   * Creates a new fixed dimension, respecting the current min/max and the increment settings. If
   * the current search dimension is not active then an inactive dimension is returned centred
   * between the lower and upper bounds.
   *
   * @param lower the lower
   * @param upper the upper
   * @return the fixed dimension
   */
  @Override
  public FixedDimension create(double lower, double upper) {
    if (!active) {
      return new FixedDimension((upper + lower) / 2);
    }
    if (lower < min) {
      lower = min;
    }
    if (upper > max) {
      upper = max;
    }
    return new FixedDimension(min, max, minIncrement, lower, upper);
  }

  /**
   * Creates a new search dimension, respecting the current settings.
   *
   * @param nIncrement the number of increments to use around the centre
   * @return the search dimension
   */
  public SearchDimension create(int nIncrement) {
    if (nIncrement <= 0) {
      // Compute the maximum number of increments to cover the range from the centre
      nIncrement = (int) Math.ceil(Math.ceil((max - min) / minIncrement) / 2);
    }
    return new SearchDimension(min, max, minIncrement, nIncrement, getLower(), getUpper());
  }

  private static boolean isInvalid(double d) {
    return Double.isNaN(d) || Double.isInfinite(d);
  }

  /**
   * Round the value to the nearest min increment. If min increment is zero no rounding is
   * performed.
   *
   * @param value the value
   * @return the rounded value
   */
  @Override
  public double round(double value) {
    if (canRound()) {
      return MathUtils.round(value, minIncrement);
    }
    return value;
  }

  /**
   * If the dimension is not active or min increment is zero no rounding is performed.
   *
   * @return true, if successful
   */
  @Override
  public boolean canRound() {
    return (active && minIncrement != 0);
  }

  /**
   * Gets the centre of the range in the dimension.
   *
   * @return the centre of the range in the dimension
   */
  @Override
  public double getCentre() {
    return round((upper + lower) / 2);
  }

  /**
   * Gets the current lower bound of the range.
   *
   * @return the current lower bound of the range
   */
  @Override
  public double getLower() {
    return lower;
  }

  /**
   * Gets the current upper bound of the range.
   *
   * @return the current upper bound of the range
   */
  @Override
  public double getUpper() {
    return upper;
  }

  @Override
  public boolean isAtBounds(double v) {
    return (v <= lower || v >= upper);
  }

  @Override
  public FixedDimension clone() {
    try {
      return (FixedDimension) super.clone();
    } catch (final CloneNotSupportedException ex) {
      return null;
    }
  }

  @Override
  public double getMin() {
    return min;
  }

  @Override
  public double getMax() {
    return max;
  }

  @Override
  public boolean isActive() {
    return active;
  }
}
