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

package uk.ac.sussex.gdsc.smlm.results.filter;

import uk.ac.sussex.gdsc.core.annotation.Nullable;

/**
 * Stores a set of results within a grid arrangement at a given resolution. Allows checking for
 * duplicates. The XY and Z resolution can be different.
 */
public interface CoordinateStore {
  /**
   * Gets the XY resolution of the store. If negative then nothing is stored.
   *
   * @return the XY resolution
   */
  double getXyResolution();

  /**
   * Gets the Z resolution of the store. If negative then this is ignored and the store behaves as
   * if processing 2D coordinates.
   *
   * @return the Z resolution
   */
  double getZResolution();

  /**
   * Queue a coordinate to the store.
   *
   * <p>It is not added to the store until flush is called.
   *
   * @param x the x
   * @param y the y
   * @param z the z
   */
  void addToQueue(double x, double y, double z);

  /**
   * Flush the queue to the store.
   */
  void flush();

  /**
   * Add a coordinate to the store. Assumes that the coordinates are within the size of the grid
   * otherwise they will be ignored.
   *
   * @param x the x
   * @param y the y
   * @param z the z
   */
  void add(double x, double y, double z);

  /**
   * Clear to the store.
   */
  void clear();

  /**
   * Check if the store contains the coordinates within the configured resolution.
   *
   * @param x the x
   * @param y the y
   * @param z the z
   * @return true, if the store contains another coordinate closer than the resolution
   */
  boolean contains(double x, double y, double z);

  /**
   * Find the closest coordinate within the configured resolution.
   *
   * @param x the x
   * @param y the y
   * @param z the z
   * @return the coordinate closer than the resolution (or null)
   */
  @Nullable
  double[] find(double x, double y, double z);

  /**
   * Create a new instance.
   *
   * @return the new coordinate store
   */
  CoordinateStore newInstance();

  /**
   * Resize to the given dimensions. If these match the existing dimensions the current store is
   * returned. Otherwise a new store is returned.
   *
   * @param minx the min x coordinate value
   * @param miny the min y coordinate value
   * @param width the width
   * @param height the height
   * @return the coordinate store
   */
  CoordinateStore resize(int minx, int miny, int width, int height);

  /**
   * Gets the the min x coordinate value.
   *
   * @return the min X
   */
  int getMinX();

  /**
   * Gets the the min y coordinate value.
   *
   * @return the min Y
   */
  int getMinY();

  /**
   * Gets the width.
   *
   * @return the width
   */
  int getWidth();

  /**
   * Gets the height.
   *
   * @return the height
   */
  int getHeight();
}
