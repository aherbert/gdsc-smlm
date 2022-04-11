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

import java.util.Arrays;
import java.util.Objects;
import org.scijava.java3d.Group;
import org.scijava.java3d.PointAttributes;
import org.scijava.java3d.PolygonAttributes;
import org.scijava.java3d.TransparencyAttributes;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Point3f;

/**
 * This class represents a list of repeated shapes in the universe that can be added to a scene. It
 * provides basic functionality for displayed Content.
 */
public abstract class ItemGroup extends Group implements ItemShape {
  /** The default colour. */
  protected static final Color3f DEFAULT_COLOUR = new Color3f(0, 1, 0);

  /** The list of points. */
  protected Point3f[] points;

  /**
   * Instantiates a new item geometry group.
   *
   * @param points the points
   */
  public ItemGroup(final Point3f[] points) {
    this.points = Objects.requireNonNull(points, "points");
  }

  /**
   * Gets the parent group to which all the shapes should be added.
   *
   * @return the parent group
   */
  protected Group getParentGroup() {
    return this;
  }

  /**
   * Gets the points.
   *
   * @return the points
   */
  public Point3f[] getPoints() {
    return points;
  }

  /**
   * Set the color of the points.
   *
   * @param color the new color
   */
  public void setColor(Color3f color) {
    setItemColor(color);
  }

  /**
   * Gets the color of the points.
   *
   * @return the color
   */
  public abstract Color3f getColor();

  /**
   * Set the transparency of the points. This is blended with a per-item alpha if present.
   *
   * @param transparency the new transparency
   */
  public void setTransparency(final float transparency) {
    final TransparencyAttributes ta = getTransparencyAttributes();
    if (ta == null) {
      return;
    }
    if (transparency == 0f) {
      ta.setTransparencyMode(TransparencyAttributes.NONE);
      ta.setTransparency(transparency);
    } else {
      ta.setTransparency(transparency);
      ta.setTransparencyMode(TransparencyAttributes.FASTEST);
    }
  }

  /**
   * Gets the global transparency attributes.
   *
   * @return the transparency attributes
   */
  protected abstract TransparencyAttributes getTransparencyAttributes();

  /**
   * Gets the transparency of the points. This is blended with a per-item alpha if present.
   *
   * @return the transparency
   */
  public float getTransparency() {
    final TransparencyAttributes ta = getTransparencyAttributes();
    if (ta == null) {
      return 0f;
    }
    return ta.getTransparency();
  }

  /**
   * Checks if any item is transparent.
   *
   * @return true, if is transparent
   */
  public boolean isTransparent() {
    final TransparencyAttributes ta = getTransparencyAttributes();
    if (ta == null) {
      return false;
    }
    return ta.getTransparency() != 0f && ta.getTransparencyMode() != TransparencyAttributes.NONE;
  }

  /**
   * Sets the shaded.
   *
   * @param shaded the new shaded
   */
  public void setShaded(boolean shaded) {
    final PolygonAttributes pa = getPolygonAttributes();
    if (pa == null) {
      return;
    }
    final int mode = (shaded) ? PolygonAttributes.POLYGON_FILL : PolygonAttributes.POLYGON_LINE;
    if (pa.getPolygonMode() != mode) {
      pa.setPolygonMode(mode);
    }
  }

  /**
   * Gets the global polygon attributes.
   *
   * @return the polygon attributes
   */
  protected abstract PolygonAttributes getPolygonAttributes();

  /**
   * Checks if is shaded.
   *
   * @return true, if is shaded
   */
  public boolean isShaded() {
    final PolygonAttributes pa = getPolygonAttributes();
    if (pa == null) {
      return false;
    }
    return pa.getPolygonMode() == PolygonAttributes.POLYGON_FILL;
  }

  /**
   * Sets the pointSize.
   *
   * @param pointSize the new pointSize
   */
  public void setPointSize(float pointSize) {
    final PointAttributes pa = getPointAttributes();
    if (pa == null) {
      return;
    }
    pa.setPointSize(pointSize);
  }

  /**
   * Gets the global point attributes.
   *
   * @return the point attributes
   */
  protected abstract PointAttributes getPointAttributes();

  /**
   * Gets the point size.
   *
   * @return the point size
   */
  public float getPointSize() {
    final PointAttributes pa = getPointAttributes();
    if (pa == null) {
      return 0f;
    }
    return pa.getPointSize();
  }

  /**
   * Returns a string representation of the underlying {@code List<Point3f>}.
   */
  @Override
  public String toString() {
    return Arrays.toString(points);
  }

  @Override
  public int size() {
    return points.length;
  }

  @Override
  public Point3f getCoordinate(int index) {
    return points[index];
  }
}
