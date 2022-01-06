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

package uk.ac.sussex.gdsc.smlm.results.filter;

/**
 * Stores a set of results within a grid arrangement at a given resolution. Allows checking for
 * duplicates.
 *
 * <p>Uses a block resolution of 1.
 */
public class GridCoordinateStore1 extends GridCoordinateStore {
  // Note: We have package level constructors so that the factory must be used to create an
  // instance.

  /**
   * Create a grid for coordinates.
   *
   * @param minx the min x coordinate value
   * @param miny the min y coordinate value
   * @param width the width
   * @param height the height
   * @param xyResolution the xy resolution
   * @param zResolution the z resolution
   */
  GridCoordinateStore1(int minx, int miny, int width, int height, double xyResolution,
      double zResolution) {
    super(minx, miny, width, height, xyResolution, zResolution);
    checkResolution(xyResolution);
  }

  private static void checkResolution(double xyResolution) {
    if (xyResolution > 1) {
      throw new IllegalArgumentException("XY Resolution must be <= 1");
    }
  }

  @Override
  public GridCoordinateStore newInstance(int minx, int miny, int width, int height) {
    return new GridCoordinateStore1(minx, miny, width, height, getXyResolution(), getZResolution());
  }

  @Override
  protected int getBlock(final double x) {
    // blockResolution is always 1
    return (int) x;
  }

  @Override
  public void changeXyResolution(double xyResolution) {
    checkResolution(xyResolution);
    super.changeXyResolution(xyResolution);
  }
}
