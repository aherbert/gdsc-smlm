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

import java.util.Arrays;
import java.util.List;
import org.jogamp.java3d.Geometry;
import org.jogamp.java3d.GeometryArray;
import org.jogamp.java3d.PointArray;
import org.jogamp.vecmath.Color3f;
import org.jogamp.vecmath.Color4f;
import org.jogamp.vecmath.Point3f;

/**
 * Create an object to represent a set of points.
 */
public class TransparentItemPointMesh extends ItemPointMesh implements TransparentItemShape {
  /**
   * Instantiates a new transparent item point mesh.
   *
   * @param mesh the mesh
   */
  public TransparentItemPointMesh(final List<Point3f> mesh) {
    super(mesh);
  }

  /**
   * Instantiates a new transparent item point mesh.
   *
   * @param mesh the mesh
   * @param color the color
   * @param transparency the transparency
   */
  public TransparentItemPointMesh(final List<Point3f> mesh, final Color3f color,
      final float transparency) {
    super(mesh, color, transparency);
  }

  @Override
  protected GeometryArray createGeometry() {
    if (mesh == null || mesh.isEmpty()) {
      return null;
    }
    final int size = size();

    final Point3f[] coords = new Point3f[size];
    mesh.toArray(coords);

    final Color4f[] colors = new Color4f[size];
    if (color == null) {
      color = DEFAULT_COLOR;
    }
    Arrays.fill(colors, new Color4f(color.x, color.y, color.z, 1));

    final GeometryArray ta =
        new PointArray(size, GeometryArray.COORDINATES | GeometryArray.COLOR_4);

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
  public void reorderFast(int[] indices) {
    changed = true;

    final int oldSize = size();
    final int size = (indices == null) ? 0 : Math.min(oldSize, indices.length);

    if (size == 0) {
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
    final float[] oldColors = new float[oldSize * 4];
    ga.getColors(0, oldColors);
    final Point3f[] coords = new Point3f[size];
    final float[] colors = new float[size * 4];
    for (int i = 0; i < size; i++) {
      final int j = indices[i];
      coords[i] = oldCoords[j];
      System.arraycopy(oldColors, j * 4, colors, i * 4, 4);
    }
    mesh = Arrays.asList(coords);

    ga.updateData(geometry -> {
      final GeometryArray geom = (GeometryArray) geometry;
      // We re-use the geometry and just truncate the vertex count
      geom.setCoordinates(0, coords);
      geom.setColors(0, colors);
      geom.setValidVertexCount(coords.length);
    });
  }

  @Override
  public void setItemColor(Color3f color) {
    if (color == null) {
      color = DEFAULT_COLOR;
    }
    this.color = color;
    final GeometryArray ga = (GeometryArray) getGeometry();
    if (ga == null) {
      return;
    }
    final int size = size();
    final float[] colors = new float[4 * size];
    ga.getColors(0, colors);
    int index = 0;
    while (index < colors.length) {
      colors[index++] = color.x;
      colors[index++] = color.y;
      colors[index++] = color.z;
      index++; // Skip over alpha
    }
    ga.setColors(0, colors);
    changed = true;
  }

  @Override
  public void setItemColor(Color3f[] color) {
    this.color = null;
    final int size = size();
    ItemHelper.checkSize(color.length, size);
    final GeometryArray ga = (GeometryArray) getGeometry();
    if (ga == null) {
      return;
    }
    final float[] colors = new float[4 * size];
    ga.getColors(0, colors);
    int index = 0;
    for (final Color3f c : color) {
      colors[index++] = c.x;
      colors[index++] = c.y;
      colors[index++] = c.z;
      index++; // Skip over alpha
    }
    ga.setColors(0, colors);
    changed = true;
  }

  @Override
  public void setItemColor4(Color4f[] color) {
    this.color = null;
    final int size = size();
    ItemHelper.checkSize(color.length, size);
    final GeometryArray ga = (GeometryArray) getGeometry();
    if (ga == null) {
      return;
    }
    ga.setColors(0, color);
    changed = true;
  }

  @Override
  public void setItemAlpha(float[] alpha) {
    final int size = size();
    ItemHelper.checkSize(alpha.length, size);
    final GeometryArray ga = (GeometryArray) getGeometry();
    if (ga == null) {
      return;
    }
    final float[] colors = new float[4 * size];
    ga.getColors(0, colors);
    for (int i = 0; i < size; i++) {
      // Set only alpha
      colors[i * 4 + 3] = alpha[i];
    }
    ga.setColors(0, colors);
    changed = true;
  }

  @Override
  public void setItemAlpha(float alpha) {
    final GeometryArray ga = (GeometryArray) getGeometry();
    if (ga == null) {
      return;
    }
    final int size = size();
    final float[] colors = new float[4 * size];
    ga.getColors(0, colors);
    for (int i = 0; i < size; i++) {
      // Set only alpha
      colors[i * 4 + 3] = alpha;
    }
    ga.setColors(0, colors);
    changed = true;
  }

  @Override
  public void getItemAlpha(float[] alpha) {
    final int size = size();
    ItemHelper.checkSize(alpha.length, size);
    final GeometryArray ga = (GeometryArray) getGeometry();
    if (ga == null) {
      return;
    }
    final float[] colors = new float[4 * size];
    ga.getColors(0, colors);
    for (int i = 0; i < size; i++) {
      // Get only alpha
      alpha[i] = colors[i * 4 + 3];
    }
  }
}
