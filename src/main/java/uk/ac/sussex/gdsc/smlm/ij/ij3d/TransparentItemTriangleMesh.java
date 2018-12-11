/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
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

import uk.ac.sussex.gdsc.core.logging.TrackProgress;

import org.scijava.java3d.Geometry;
import org.scijava.java3d.GeometryArray;
import org.scijava.java3d.GeometryUpdater;
import org.scijava.java3d.TriangleArray;
import org.scijava.java3d.utils.geometry.GeometryInfo;
import org.scijava.java3d.utils.geometry.NormalGenerator;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Color4f;
import org.scijava.vecmath.Point3f;
import org.scijava.vecmath.Vector3f;

import java.util.Arrays;

/**
 * Use a triangle mesh object to represent a set of points. The object is duplicated, scaled and
 * translated for each point.
 */
public class TransparentItemTriangleMesh extends ItemTriangleMesh implements TransparentItemShape {
  /**
   * Instantiates a new transparent item triangle mesh.
   *
   * <p>This will repeat the object for each input point. The object is assumed to be centred on the
   * origin. It will be scaled and translated for each input point.
   *
   * @param objectVertices the vertices of the object for a single point
   * @param points the points
   * @param sizes the size of each point
   * @param color the color
   * @param transp the transparency
   */
  public TransparentItemTriangleMesh(Point3f[] objectVertices, Point3f[] points, Point3f[] sizes,
      Color3f color, float transp) {
    super(objectVertices, points, sizes, color, transp);
  }

  /**
   * Instantiates a new transparent item triangle mesh.
   *
   * <p>This will repeat the object for each input point. The object is assumed to be centred on the
   * origin. It will be scaled and translated for each input point.
   *
   * <p>The crease angle is used to collapse facets normals at a vertex into a single normal for
   * smoothing shading. Set to 0 to draw the polygon with no shading.
   *
   * @param objectVertices the vertices of the object for a single point
   * @param points the points
   * @param sizes the size of each point
   * @param color the color
   * @param transp the transparency
   * @param creaseAngle the crease angle (in degrees). Set to negative to ignore. The default is 44.
   * @param progress the progress
   */
  public TransparentItemTriangleMesh(Point3f[] objectVertices, Point3f[] points, Point3f[] sizes,
      Color3f color, float transp, double creaseAngle, TrackProgress progress) {
    super(objectVertices, points, sizes, color, transp, creaseAngle, progress);
  }

  @Override
  protected GeometryArray createGeometry() {
    if (mesh == null || mesh.size() < 3) {
      return null;
    }
    final int vertexCount = mesh.size();

    final Point3f[] coords = new Point3f[vertexCount];
    mesh.toArray(coords);

    final Color4f colors[] = new Color4f[vertexCount];
    if (color == null) {
      color = DEFAULT_COLOR;
    }
    Arrays.fill(colors, new Color4f(color.x, color.y, color.z, 1));

    final GeometryArray ta = new TriangleArray(vertexCount,
        GeometryArray.COORDINATES | GeometryArray.COLOR_4 | GeometryArray.NORMALS);

    ta.setCoordinates(0, coords);
    ta.setColors(0, colors);

    // generate normals
    final GeometryArray result;

    if (dirty) {
      final GeometryInfo gi = new GeometryInfo(ta);
      final NormalGenerator ng = new NormalGenerator();
      ng.generateNormals(gi);
      result = gi.getGeometryArray();
    } else {
      // Use the same normals for each repeated object
      final Vector3f[] normals = new Vector3f[vertexCount];

      // Binary fill
      int fill = objectNormals.length;
      System.arraycopy(objectNormals, 0, normals, 0, fill);
      for (int i = 2; i < points.length; i *= 2) {
        System.arraycopy(normals, 0, normals, fill, fill);
        fill *= 2;
      }
      // Final fill
      System.arraycopy(normals, 0, normals, fill, normals.length - fill);

      ta.setNormals(0, normals);

      result = ta;
    }

    result.setCapability(GeometryArray.ALLOW_NORMAL_WRITE);
    result.setCapability(GeometryArray.ALLOW_COLOR_WRITE);
    result.setCapability(GeometryArray.ALLOW_COORDINATE_WRITE);
    result.setCapability(GeometryArray.ALLOW_COUNT_WRITE);
    result.setCapability(GeometryArray.ALLOW_COUNT_READ);
    result.setCapability(GeometryArray.ALLOW_FORMAT_READ);
    result.setCapability(Geometry.ALLOW_INTERSECT);
    result.setValidVertexCount(vertexCount);

    return result;
  }

  /** {@inheritDoc} */
  @Override
  public void reorderFast(int[] indices) throws IllegalArgumentException {
    if (dirty) {
      throw new IllegalArgumentException("Mesh has been modified");
    }

    changed = true;

    final int oldSize = size();
    final int size = (indices == null) ? 0 : Math.min(oldSize, indices.length);

    if (size == 0 || indices == null) {
      mesh.clear();
      points = new Point3f[0];
      sizes = new Point3f[0];
      this.setGeometry(null);
      return;
    }

    // From here on we assume the current geometry will not be null
    // as this only happens when the original size is zero. Size has
    // been checked at this point to be the smaller of new and old.
    final GeometryArray ga = (GeometryArray) getGeometry();

    points = reorder(points, indices);
    // Sizes could be null or a single size
    if (sizes != null && sizes.length == points.length) {
      sizes = reorder(sizes, indices);
    }

    // Reorder all things in the geometry: coordinates and colour
    // The normals can be copied as they are unchanged.
    // The mesh should contain the same coordinates as the geometry array.
    final int objectSize = objectVertices.length;
    final Point3f[] oldCoords = mesh.toArray(new Point3f[mesh.size()]);
    final float[] oldColors = new float[oldCoords.length * 4];
    ga.getColors(0, oldColors);
    final Point3f[] coords = new Point3f[size * objectSize];
    final float[] colors = new float[coords.length * 4];
    for (int i = 0; i < size; i++) {
      final int j = indices[i];

      final int ii = i * objectSize;
      final int jj = j * objectSize;
      System.arraycopy(oldCoords, jj, coords, ii, objectSize);
      System.arraycopy(oldColors, jj * 4, colors, ii * 4, objectSize * 4);
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

  /** {@inheritDoc} */
  @Override
  public void setItemColor(Color3f color) {
    if (color == null) {
      color = DEFAULT_COLOR;
    }
    this.color = color;
    final int size = size();
    final GeometryArray ga = (GeometryArray) getGeometry();
    if (ga == null) {
      return;
    }
    final int objectSize = objectVertices.length;
    final int N = objectSize * size;
    final float[] colors = new float[N * 4];
    ga.getColors(0, colors);
    int i = 0;
    while (i < colors.length) {
      colors[i++] = color.x;
      colors[i++] = color.y;
      colors[i++] = color.z;
      i++; // Skip over alpha
    }
    ga.setColors(0, colors);
    changed = true;
  }

  /** {@inheritDoc} */
  @Override
  public void setItemColor(Color3f[] color) throws IllegalArgumentException {
    this.color = null;
    final int size = size();
    if (color.length != size) {
      throw new IllegalArgumentException("list of size " + size + " expected");
    }
    final GeometryArray ga = (GeometryArray) getGeometry();
    if (ga == null) {
      return;
    }
    final int objectSize = objectVertices.length;
    final int N = objectSize * size;
    final float[] colors = new float[N * 4];
    ga.getColors(0, colors);
    int i = 0;
    for (final Color3f c : color) {
      for (int j = objectSize; j-- > 0;) {
        colors[i++] = c.x;
        colors[i++] = c.y;
        colors[i++] = c.z;
        i++; // Skip over alpha
      }
    }
    ga.setColors(0, colors);
    changed = true;
  }

  /** {@inheritDoc} */
  @Override
  public void setItemColor4(Color4f[] color) throws IllegalArgumentException {
    this.color = null;
    final int size = size();
    if (color.length != size) {
      throw new IllegalArgumentException("list of size " + size + " expected");
    }
    final GeometryArray ga = (GeometryArray) getGeometry();
    if (ga == null) {
      return;
    }
    final int objectSize = objectVertices.length;
    final int N = objectSize * size;
    final Color4f[] colors = new Color4f[N];
    int i = 0;
    for (final Color4f c : color) {
      for (int j = objectSize; j-- > 0;) {
        colors[i++] = c;
      }
    }
    ga.setColors(0, colors);
    changed = true;
  }

  /** {@inheritDoc} */
  @Override
  public void setItemAlpha(float[] alpha) throws IllegalArgumentException {
    final int size = size();
    if (alpha.length != size) {
      throw new IllegalArgumentException("list of size " + size + " expected");
    }
    final GeometryArray ga = (GeometryArray) getGeometry();
    if (ga == null) {
      return;
    }
    final int objectSize = objectVertices.length;
    final int N = objectSize * size;
    final float[] colors = new float[N * 4];
    ga.getColors(0, colors);
    int i = 3;
    for (final float f : alpha) {
      for (int j = objectSize; j-- > 0;) {
        colors[i] = f;
        i += 4;
      }
    }
    ga.setColors(0, colors);
    changed = true;
  }

  /** {@inheritDoc} */
  @Override
  public void setItemAlpha(float alpha) throws IllegalArgumentException {
    int size = size();
    final GeometryArray ga = (GeometryArray) getGeometry();
    if (ga == null) {
      return;
    }
    final int objectSize = objectVertices.length;
    final int N = objectSize * size;
    final float[] colors = new float[N * 4];
    ga.getColors(0, colors);
    int i = 3;
    while (size-- > 0) {
      for (int j = objectSize; j-- > 0;) {
        colors[i] = alpha;
        i += 4;
      }
    }
    ga.setColors(0, colors);
    changed = true;
  }

  /** {@inheritDoc} */
  @Override
  public void getItemAlpha(float[] alpha) throws IllegalArgumentException {
    final int size = size();
    if (alpha.length != size) {
      throw new IllegalArgumentException("list of size " + size + " expected");
    }
    final GeometryArray ga = (GeometryArray) getGeometry();
    if (ga == null) {
      return;
    }
    final int objectSize = objectVertices.length;
    final int N = objectSize * size;
    final float[] colors = new float[N * 4];
    ga.getColors(0, colors);
    for (int i = 0; i < size; i++) {
      // Get only alpha
      alpha[i] = colors[i * 4 * objectSize + 3];
    }
  }
}
