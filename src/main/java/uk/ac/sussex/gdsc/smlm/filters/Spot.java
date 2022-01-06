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

package uk.ac.sussex.gdsc.smlm.filters;

/**
 * Identify a candidate spot (local maximum).
 */
public class Spot {
  /** The x. */
  public final int x;

  /** The y. */
  public final int y;

  /** The intensity. */
  public final float intensity;

  /** The score. */
  private float score;

  /**
   * Constructor that sets the score (for sorting) equal to the intensity.
   *
   * @param x The x-coordinate
   * @param y The y-coordinate
   * @param intensity The intensity of the spot
   */
  public Spot(int x, int y, float intensity) {
    this.x = x;
    this.y = y;
    this.intensity = score = intensity;
  }

  /**
   * Instantiates a new spot.
   *
   * @param x The x-coordinate
   * @param y The y-coordinate
   * @param intensity The intensity of the spot
   * @param score The score used for sorting
   */
  public Spot(int x, int y, float intensity, float score) {
    this.x = x;
    this.y = y;
    this.intensity = intensity;
    this.score = score;
  }

  /**
   * Copy constructor.
   *
   * @param source the source
   */
  protected Spot(Spot source) {
    x = source.x;
    y = source.y;
    intensity = source.intensity;
    score = source.score;
  }

  /**
   * Get a copy.
   *
   * @return the copy
   */
  public Spot copy() {
    return new Spot(this);
  }

  /**
   * Get the distance between the two spots.
   *
   * @param other the other spot
   * @return The distance
   */
  public double distance(Spot other) {
    return Math.sqrt(distance2(other));
  }

  /**
   * Get the squared distance between the two spots.
   *
   * @param other the other spot
   * @return The squared distance
   */
  public double distance2(Spot other) {
    final double dx = (double) x - other.x;
    final double dy = (double) y - other.y;
    return dx * dx + dy * dy;
  }

  /**
   * Gets the score used for sorting the spots.
   *
   * @return the score.
   */
  public float getScore() {
    return score;
  }

  /**
   * Sets the score used for sorting the spots.
   *
   * @param score the score
   */
  public void setScore(float score) {
    this.score = score;
  }
}
