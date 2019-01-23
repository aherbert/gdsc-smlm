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

import uk.ac.sussex.gdsc.core.data.ComputationException;
import uk.ac.sussex.gdsc.core.utils.MathUtils;

import java.util.Arrays;

/**
 * Specify the dimensions for a search.
 */
public class SearchDimension implements Dimension {
  /** The minimum of the range. */
  public final double min;

  /** The maximum of the range. */
  public final double max;

  /** The min increment to use around the centre. */
  public final double minIncrement;

  /** The number of increments to use around the centre. */
  public final int increments;

  /**
   * Set to true if {@link #min} &lt; {@link #max}.
   */
  public final boolean active;

  private double centre;
  private double increment;
  private double reduceFactor = 0.5;
  private double[] values;
  private boolean pad;

  /**
   * Instantiates a new inactive search dimension. The centre can be set to any value, the default
   * is zero.
   */
  public SearchDimension() {
    this(0);
  }

  /**
   * Instantiates a new inactive search dimension. The centre can be set to any value.
   *
   * @param centre the centre
   */
  public SearchDimension(double centre) {
    this(0, 0, 0, 1);
    setCentre(centre);
  }

  /**
   * Instantiates a new search dimension.
   *
   * @param min the minimum of the range
   * @param max the maximum of the range
   * @param minIncrement the min increment to use around the centre
   * @param increments the number of increments to use around the centre
   */
  public SearchDimension(double min, double max, double minIncrement, int increments) {
    this(min, max, minIncrement, increments, min, max);
  }

  /**
   * Instantiates a new search dimension.
   *
   * @param min the minimum of the range
   * @param max the maximum of the range
   * @param minIncrement the min increment to use around the centre
   * @param increments the number of increments to use around the centre
   * @param lower the current lower bound of the range (will be clipped to min/max)
   * @param upper the current upper bound of the range (will be clipped to min/max)
   */
  public SearchDimension(double min, double max, double minIncrement, int increments, double lower,
      double upper) {
    if (!Double.isFinite(min)) {
      throw new IllegalArgumentException("Min is not a valid number: " + min);
    }
    if (!Double.isFinite(max)) {
      throw new IllegalArgumentException("Max is not a valid number: " + max);
    }
    if (!Double.isFinite(lower)) {
      throw new IllegalArgumentException("Lower is not a valid number: " + lower);
    }
    if (!Double.isFinite(upper)) {
      throw new IllegalArgumentException("Upper is not a valid number: " + upper);
    }
    if (!Double.isFinite(minIncrement)) {
      throw new IllegalArgumentException("Min increment is not a valid number: " + minIncrement);
    }
    if (max < min) {
      throw new IllegalArgumentException("Max is less than min");
    }
    this.active = min < max;
    if (active && increments < 1) {
      throw new IllegalArgumentException("Steps must be more than 0: " + increments);
    }
    if (minIncrement < 0) {
      throw new IllegalArgumentException("Min increment is negative: " + minIncrement);
    }
    if (upper < lower) {
      throw new IllegalArgumentException("Upper is less than lower");
    }

    // We round to the min increment so that the values returned should be identical if the centre
    // is moved by a factor of the increment.
    this.minIncrement = minIncrement;
    this.min = round(min);
    this.max = round(max);
    this.increments = increments;

    // Rounding changes the range so bring the upper and lower back within
    lower = MathUtils.clip(this.min, this.max, lower);
    upper = MathUtils.clip(this.min, this.max, upper);

    setCentre((upper + lower) / 2);
    setIncrement((upper - lower) / (2 * increments));
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected SearchDimension(SearchDimension source) {
    this.min = source.min;
    this.max = source.max;
    this.minIncrement = source.minIncrement;
    this.increments = source.increments;
    this.active = source.active;
    this.centre = source.centre;
    this.increment = source.increment;
    this.reduceFactor = source.reduceFactor;
    this.values = source.values;
    this.pad = source.pad;
  }

  /**
   * Create a copy.
   *
   * @return the copy
   */
  public SearchDimension copy() {
    return new SearchDimension(this);
  }

  /**
   * Creates a new search dimension, respecting the current min/max and the increment settings. If
   * the current search dimension is not active then an inactive dimension is returned centred
   * between the lower and upper bounds.
   *
   * @param lower the current lower bound of the range
   * @param upper the current upper bound of the range
   * @return the search dimension
   */
  @Override
  public SearchDimension create(double lower, double upper) {
    if (!active) {
      return new SearchDimension((upper + lower) / 2);
    }
    if (lower < min) {
      lower = min;
    }
    if (upper > max) {
      upper = max;
    }
    return new SearchDimension(min, max, minIncrement, increments, lower, upper);
  }

  /**
   * Creates a new search dimension, respecting the current settings and changing the number of
   * increments.
   *
   * @param increments the number of increments to use around the centre
   * @return the search dimension
   */
  public SearchDimension create(int increments) {
    return new SearchDimension(min, max, minIncrement, increments, getLower(), getUpper());
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
   */
  @Override
  public boolean canRound() {
    return (active && minIncrement != 0);
  }

  /**
   * Sets the centre.
   *
   * @param centre the new centre of the range in the dimension
   */
  public void setCentre(double centre) {
    if (active && (centre < min || centre > max)) {
      throw new IllegalArgumentException("Centre is outside min/max range");
    }
    this.centre = round(centre);
    values = null;
  }

  /**
   * Gets the centre of the range in the dimension.
   *
   * @return the centre of the range in the dimension
   */
  @Override
  public double getCentre() {
    return centre;
  }

  /**
   * Gets the increment.
   *
   * @return the increment
   */
  public double getIncrement() {
    return increment;
  }

  /**
   * Sets the increment.
   *
   * @param increment the new increment
   */
  public void setIncrement(double increment) {
    if (increment < minIncrement) {
      increment = minIncrement;
    }
    this.increment = round(increment);
    values = null;
  }

  /**
   * Gets the current lower bound of the range.
   *
   * @return the current lower bound of the range
   */
  @Override
  public double getLower() {
    return values()[0];
  }

  /**
   * Gets the current upper bound of the range.
   *
   * @return the current upper bound of the range
   */
  @Override
  public double getUpper() {
    values();
    return values[values.length - 1];
  }

  @Override
  public boolean isAtBounds(double value) {
    values();
    return (value <= values[0] || value >= values[values.length - 1]);
  }

  /**
   * Get the values of the dimension around the current centre using the configured increments.
   * Note: If the values are outside the min/max range then by default the number of values may be
   * reduced.
   *
   * <p>If the pad setting is enabled then the number of values should remain constant as the values
   * are padded in the opposite direction.
   *
   * @return the values
   */
  public double[] values() {
    if (values != null) {
      return values;
    }

    if (!active) {
      values = new double[] {centre};
      return values;
    }

    values = new double[getMaxLength()];
    int size = 0;

    double value = round(centre - increments * increment);
    if (value < min) {
      values[size++] = min;

      // Avoid further values below the min
      for (int i = increments - 1; i >= 1; i--) {
        value = round(centre - i * increment);
        if (value < min) {
          continue;
        }
        values[size++] = value;
      }

      if (centre != min) {
        values[size++] = centre;
      }
    } else {
      // Not at the limit
      for (int i = increments; i >= 1; i--) {
        values[size++] = round(centre - i * increment);
      }

      values[size++] = centre; // Already rounded and within range
    }

    for (int i = 1; i <= increments; i++) {
      value = round(centre + i * increment);
      if (value > max) {
        if (centre != max) {
          values[size++] = max;
        }
        // Avoid further values outside the range
        break;
      }
      values[size++] = value;
    }

    // Check for duplicates if at the limits
    if (size != values.length) {
      // Option to pad in the opposite direction
      if (pad) {
        if (values[0] == min) {
          if (values[size - 1] != max) {
            // Pad up
            for (int i = increments + 1; size < values.length; i++) {
              value = round(centre + i * increment);
              if (value > max) {
                values[size++] = max;
                break;
              }
              values[size++] = value;
            }
          }
        } else {
          // Pad down
          for (int i = increments + 1; size < values.length; i++) {
            value = round(centre - i * increment);
            if (value < min) {
              values[size++] = min;
              break;
            }
            values[size++] = value;
          }
          // Simple option is to sort.
          // A better option is to copy the values to the correct place.
          Arrays.sort(values, 0, size);
        }

        // In case we could not pad enough
        if (size != values.length) {
          values = Arrays.copyOf(values, size);
        }
      } else {
        // No padding so truncate
        values = Arrays.copyOf(values, size);
      }
    }

    return values;
  }

  /**
   * Gets the max length of the values array.
   *
   * @return the max length
   */
  public int getMaxLength() {
    return 2 * increments + 1;
  }

  /**
   * Gets the reduce factor.
   *
   * @return the reduce factor
   */
  public double getReduceFactor() {
    return reduceFactor;
  }

  /**
   * Set the reduce factor. A value of 1 will prevent the range being reduced by the
   * {@link #reduce()} method.
   *
   * @param reduceFactor the reduce factor (must be between 0 (exclusive) and 1 (inclusive))
   */
  public void setReduceFactor(double reduceFactor) {
    if (reduceFactor <= 0 || reduceFactor > 1) {
      throw new IllegalArgumentException("Reduce factor must be between 0 and 1 (inclusive)");
    }
    this.reduceFactor = reduceFactor;
  }

  /**
   * Reduce the size of the increment by multiplying by the reduce factor.
   */
  public void reduce() {
    setIncrement(increment * reduceFactor);
  }

  /**
   * Checks if the {@link #reduce()} method will result in a change to the range returned by
   * {@link #values()}.
   *
   * @return true, if the range can be reduced
   */
  public boolean canReduce() {
    return active && increment != minIncrement && reduceFactor < 1;
  }

  /**
   * Checks if padding the values in the opposite direction when the range overlaps the min/max.
   *
   * @return true, if padding the values
   */
  public boolean isPad() {
    return pad;
  }

  /**
   * Set to true if padding the values in the opposite direction when the range overlaps the
   * min/max.
   *
   * @param pad true, if padding the values in the opposite direction when the range overlaps the
   *        min/max
   */
  public void setPad(boolean pad) {
    this.pad = pad;
    values = null;
  }

  /**
   * Enumerate from the lower to the upper value using the configured increments (n) to define the
   * number of steps (2*n+1).
   *
   * <p>No range check is performed against the current min/max so the returned values can be
   * outside the allowed range. The min increment setting is respected so the number of actual steps
   * may be smaller.
   *
   * @param lower the lower
   * @param upper the upper
   * @return the double[]
   */
  public double[] enumerate(double lower, double upper) {
    return enumerate(lower, upper, 2 * increments + 1);
  }

  /**
   * Enumerate from the lower to the upper value using the number of steps.
   *
   * <p>No range check is performed against the current min/max so the returned values can be
   * outside the allowed range. The min increment setting is respected so the number of actual steps
   * may be smaller.
   *
   * @param lower the lower
   * @param upper the upper
   * @param steps the steps
   * @return the double[]
   */
  public double[] enumerate(double lower, double upper, int steps) {
    if (upper <= lower || steps < 2) {
      return new double[] {round(lower)};
    }

    lower = round(lower);
    upper = round(upper);

    double inc = (upper - lower) / (steps - 1);
    if (inc < minIncrement) {
      inc = minIncrement;
    }
    steps = 1 + (int) Math.ceil((upper - lower) / inc);

    final double[] values = new double[steps];
    int size = 0;
    for (int i = 0; i < values.length; i++) {
      final double v = round(lower + i * inc);
      if (v >= upper) {
        values[size++] = upper;
        break;
      }
      values[size++] = v;
    }

    // Check
    if (values[size - 1] != upper) {
      throw new ComputationException("enumeration is invalid");
    }

    return (size != values.length) ? Arrays.copyOf(values, size) : values;

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
