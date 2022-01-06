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

package uk.ac.sussex.gdsc.smlm.ij.ij3d;

import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Point3f;

/**
 * Interface for shape objects that represent a set of items.
 */
public interface ItemShape {
  /**
   * Gets the number of items.
   *
   * @return the size
   */
  int size();

  /**
   * Gets the coordinate of the specified item.
   *
   * @param index the index
   * @return the coordinate
   */
  Point3f getCoordinate(int index);

  /**
   * Sets the color for each item.
   *
   * @param color the new color
   */
  void setItemColor(final Color3f color);

  /**
   * Sets the color for each item.
   *
   * @param color the new color
   * @throws IllegalArgumentException if the number of colours is incorrect
   */
  void setItemColor(final Color3f[] color);
}
