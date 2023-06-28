/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
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

package uk.ac.sussex.gdsc.smlm.search;

import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.ValidationUtils;

/**
 * Specify the dimensions for a search.
 */
public class FixedDimension implements Dimension {
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
   * @throws IllegalArgumentException if the centre is not finite
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
   * @throws IllegalArgumentException if the numbers are not finite or {@code max < min}
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
   * @throws IllegalArgumentException if the numbers are not finite or {@code max < min} or
   *         {@code upper < lower}
   */
  public FixedDimension(double min, double max, double minIncrement, double lower, double upper) {
    ValidationUtils.checkArgument(Double.isFinite(min), "Min is not finite: %s", min);
    ValidationUtils.checkArgument(Double.isFinite(max), "Max is not finite: %s", max);
    ValidationUtils.checkArgument(Double.isFinite(lower), "Lower is not finite: %s", lower);
    ValidationUtils.checkArgument(Double.isFinite(upper), "Upper is not finite: %s", upper);
    ValidationUtils.checkArgument(Double.isFinite(minIncrement), "Min increment is not finite: %s",
        minIncrement);
    ValidationUtils.checkArgument(min <= max, "Max (%s) is not greater than min (%s)", max, min);
    ValidationUtils.checkPositive(minIncrement, "Min increment");
    ValidationUtils.checkArgument(lower <= upper, "Upper (%s) is not greater than lower (%s)",
        upper, lower);

    this.active = min < max;

    // Clip to range
    lower = MathUtils.clip(min, max, lower);
    upper = MathUtils.clip(min, max, upper);

    // We round to the min increment so that the values returned should be identical if the centre
    // is moved by a factor of the increment.
    this.minIncrement = minIncrement;
    this.min = round(min);
    this.max = round(max);
    this.lower = round(lower);
    this.upper = round(upper);
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected FixedDimension(FixedDimension source) {
    this.min = source.min;
    this.max = source.max;
    this.lower = source.lower;
    this.upper = source.upper;
    this.minIncrement = source.minIncrement;
    this.active = source.active;
  }

  /**
   * Create a copy.
   *
   * @return the copy
   */
  public FixedDimension copy() {
    return new FixedDimension(this);
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
   * @param numberOfIncrements the number of increments to use around the centre
   * @return the search dimension
   */
  public SearchDimension create(int numberOfIncrements) {
    if (numberOfIncrements <= 0) {
      // Compute the maximum number of increments to cover the range from the centre
      numberOfIncrements = (int) Math.ceil(Math.ceil((max - min) / minIncrement) / 2);
    }
    return new SearchDimension(min, max, minIncrement, numberOfIncrements, getLower(), getUpper());
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
  public boolean isAtBounds(double value) {
    return (value <= lower || value >= upper);
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
