/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Package
 *
 * Software for single molecule localisation microscopy (SMLM) in ImageJ
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

package uk.ac.sussex.gdsc.smlm.ij.ij3d;

import customnode.CustomPointMesh;
import java.util.Arrays;
import java.util.List;
import org.scijava.java3d.Appearance;
import org.scijava.java3d.Geometry;
import org.scijava.java3d.GeometryArray;
import org.scijava.java3d.GeometryUpdater;
import org.scijava.java3d.PointArray;
import org.scijava.java3d.PointAttributes;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Point3f;

/**
 * Create an object to represent a set of points.
 */
public class ItemPointMesh extends CustomPointMesh implements UpdateableItemShape {
  /**
   * Instantiates a new item point mesh.
   *
   * @param mesh the mesh
   */
  public ItemPointMesh(final List<Point3f> mesh) {
    super(mesh);
  }

  /**
   * Instantiates a new item point mesh.
   *
   * @param mesh the mesh
   * @param color the color
   * @param transparency the transparency
   */
  public ItemPointMesh(final List<Point3f> mesh, final Color3f color, final float transparency) {
    super(mesh, color, transparency);
  }

  @Override
  protected Appearance createAppearance() {
    final Appearance appearance = super.createAppearance();
    final PointAttributes pointAttributes = appearance.getPointAttributes();
    // This allows points to support transparency
    pointAttributes.setPointAntialiasingEnable(true);
    return appearance;
  }

  @Override
  protected GeometryArray createGeometry() {
    if (mesh == null || mesh.isEmpty()) {
      return null;
    }
    final int size = size();

    final Point3f[] coords = new Point3f[size];
    mesh.toArray(coords);

    final Color3f[] colors = new Color3f[size];
    Arrays.fill(colors, (color == null) ? DEFAULT_COLOR : color);

    final GeometryArray ta =
        new PointArray(size, GeometryArray.COORDINATES | GeometryArray.COLOR_3);

    ta.setValidVertexCount(size);

    ta.setCoordinates(0, coords);
    ta.setColors(0, colors);

    ta.setCapability(GeometryArray.ALLOW_COLOR_WRITE);
    ta.setCapability(GeometryArray.ALLOW_COORDINATE_WRITE);
    ta.setCapability(GeometryArray.ALLOW_COUNT_WRITE);
    ta.setCapability(GeometryArray.ALLOW_COUNT_READ);
    ta.setCapability(Geometry.ALLOW_INTERSECT);

    return ta;
  }

  @Override
  public void reorder(int[] indices) {
    checkIndices(indices, mesh.size());
    reorderFast(indices);
  }

  @Override
  public void reorderFast(int[] indices) {
    changed = true;

    final int oldSize = size();
    final int size = (indices == null) ? 0 : Math.min(oldSize, indices.length);

    if (size == 0 || indices == null) {
      mesh.clear();
      this.setGeometry(null);
      return;
    }

    // From here on we assume the current geometry will not be null
    // as this only happens when the original size is zero. Size has
    // been checked at this point to be the smaller of new and old.
    final GeometryArray ga = (GeometryArray) getGeometry();

    // Reorder all things in the geometry: coordinates and colour
    final Point3f[] oldCoords = mesh.toArray(new Point3f[0]);
    final float[] oldColors = new float[oldSize * 3];
    ga.getColors(0, oldColors);
    final Point3f[] coords = new Point3f[size];
    final float[] colors = new float[size * 3];
    for (int i = 0; i < size; i++) {
      final int j = indices[i];
      coords[i] = oldCoords[j];
      System.arraycopy(oldColors, j * 3, colors, i * 3, 3);
    }
    mesh = Arrays.asList(coords);

    ga.updateData(new GeometryUpdater() {
      @Override
      public void updateData(Geometry geometry) {
        final GeometryArray ga = (GeometryArray) geometry;
        // We re-use the geometry and just truncate the vertex count
        ga.setCoordinates(0, coords);
        ga.setColors(0, colors);
        ga.setValidVertexCount(coords.length);
      }
    });

    // this.setGeometry(ga);
  }

  /**
   * Check the indices contain a valid natural order of the specifed size.
   *
   * @param indices the indices
   * @param size the size
   */
  public static void checkIndices(int[] indices, int size) {
    if (indices == null || indices.length != size) {
      throw new IllegalArgumentException("Indices length do not match the size of the mesh");
    }

    // Check all indices are present.
    // Do a sort and then check it is a natural order
    final int[] check = indices.clone();
    Arrays.sort(check);
    for (int i = 0; i < check.length; i++) {
      if (check[i] != i) {
        throw new IllegalArgumentException("Indices do not contain a valid natural order");
      }
    }
  }

  @Override
  public int size() {
    return mesh.size();
  }

  @Override
  public Point3f getCoordinate(int index) {
    return mesh.get(index);
  }

  @Override
  public void setColor(Color3f color) {
    // Delegate this to the interface implementation.
    // Allows transparent version to only implement to the interface method.
    setItemColor(color);
  }

  @Override
  public void setItemColor(Color3f color) {
    super.setColor(color);
  }

  @Override
  public void setItemColor(Color3f[] color) {
    this.color = null;
    if (color.length != size()) {
      throw new IllegalArgumentException("list of size " + size() + " expected");
    }
    final GeometryArray ga = (GeometryArray) getGeometry();
    if (ga == null) {
      return;
    }
    ga.setColors(0, color);
    changed = true;
  }

  @Override
  public void calculateMinMaxCenterPoint(Point3f min, Point3f max, Point3f center) {
    final Point3f[] points = new Point3f[size()];
    mesh.toArray(points);
    CustomContentHelper.calculateMinMaxCenterPoint(min, max, center, points);
  }
}
