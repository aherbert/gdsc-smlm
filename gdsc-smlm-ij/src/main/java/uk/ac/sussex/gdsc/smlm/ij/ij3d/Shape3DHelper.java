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

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.lang3.tuple.Pair;
import org.scijava.java3d.Appearance;
import org.scijava.java3d.ColoringAttributes;
import org.scijava.java3d.GeometryArray;
import org.scijava.java3d.Material;
import org.scijava.java3d.PointArray;
import org.scijava.java3d.PointAttributes;
import org.scijava.java3d.PolygonAttributes;
import org.scijava.java3d.Shape3D;
import org.scijava.java3d.TriangleArray;
import org.scijava.java3d.utils.geometry.GeometryInfo;
import org.scijava.java3d.utils.geometry.NormalGenerator;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Point3f;
import org.scijava.vecmath.Vector3f;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.smlm.data.NamedObject;

/**
 * Create Shape3D objects.
 */
public final class Shape3DHelper {

  private static final float[][] TRI_VERTICES =
      {{sqrt(8d / 9), 0, 0}, {-sqrt(2d / 9), sqrt(2d / 3), 0}, {-sqrt(2d / 9), -sqrt(2d / 3), 0}};
  private static final int[][] TRI_FACES = {{0, 1, 2}};

  private static final float[][] SQUARE_VERTICES = {{1, 1, 0}, {-1, 1, 0}, {-1, -1, 0}, {1, -1, 0}};
  private static final int[][] SQUARE_FACES = {{0, 1, 3}, {3, 1, 2}};

  // https://en.m.wikipedia.org/wiki/Tetrahedron
  // based on alternated cube
  private static final float[][] TETRA_VERTICES =
      {{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}};
  private static final int[][] TETRA_FACES = {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}};

  // https://en.m.wikipedia.org/wiki/Octahedron
  private static final float[][] OCTA_VERTICES =
      {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1},};
  private static final int[][] OCTA_FACES =
      {{0, 3, 4}, {3, 1, 4}, {1, 2, 4}, {2, 0, 4}, {3, 0, 5}, {1, 3, 5}, {2, 1, 5}, {0, 2, 5},};

  private static final float[][] CUBE_VERTICES = {{1, 1, -1}, {-1, 1, -1}, {-1, -1, -1},
      {1, -1, -1}, {1, 1, 1}, {-1, 1, 1}, {-1, -1, 1}, {1, -1, 1},};
  private static final int[][] CUBE_FACES = {{0, 1, 3}, {3, 1, 2}, {0, 4, 7}, {0, 7, 3}, {1, 5, 6},
      {1, 6, 2}, {3, 7, 6}, {3, 6, 2}, {0, 4, 5}, {0, 5, 1}, {4, 5, 7}, {4, 6, 7},};
  private static final int[][] CUBE_FACES4 =
      {{0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}, {7, 4, 5, 6}, {2, 1, 0, 3},};

  //@formatter:off
  /**
   * The rendering.
   */
  public enum Rendering implements NamedObject {
    /** Point. */
    POINT {
      @Override public String getName() { return "Point"; }
      @Override public boolean is2D() { return true; }},
    /** Square. */
    SQUARE {
      @Override public String getName() { return "Square"; }
      @Override public boolean is2D() { return true; }},
    /** Hexagon. */
    HEXAGON{
      @Override public String getName() { return "Hexagon"; }
      @Override public boolean is2D() { return true; }},
    /** Low resolution circle. */
    LOW_RES_CIRCLE {
      @Override public String getName() { return "Low resolution circle"; }
      @Override public boolean is2D() { return true; }},
    /** High resolution circle. */
    HIGH_RES_CIRCLE {
      @Override public String getName() { return "High resolution circle"; }
      @Override public boolean is2D() { return true; }},
    /** Cube. */
    CUBE { @Override public String getName() { return "Cube"; }},
    /** Icosahedron. */
    ICOSAHEDRON  { @Override public String getName() { return "Icosahedron"; }},
    /** Low resolution sphere. */
    LOW_RES_SPHERE {
      @Override public String getName() { return "Low Resolution Sphere"; }
      @Override public boolean isHighResolution() { return true; }},
    /** High resolution sphere. */
    HIGH_RES_SPHERE  {
      @Override public String getName() { return "High Resolution Sphere"; }
      @Override public boolean isHighResolution() { return true; }},
    /** Super high resolution sphere. */
    SUPER_HIGH_RES_SPHERE {
      @Override public String getName() { return "Super-High Resolution Sphere"; }
      @Override public boolean isHighResolution() { return true; }},
    ;
    //@formatter:on

    /**
     * Checks if is 2d.
     *
     * @return true, if is 2d
     */
    public boolean is2D() {
      return false;
    }

    /**
     * Checks if is high resolution.
     *
     * @return true, if is high resolution
     */
    public boolean isHighResolution() {
      return false;
    }

    /**
     * For number.
     *
     * @param number the number
     * @return the rendering
     */
    public static Rendering forNumber(int number) {
      final Rendering[] values = Rendering.values();
      if (number < 0 || number >= values.length) {
        throw new IllegalArgumentException();
      }
      return values[number];
    }
  }

  private static int[] numberOfTriangles = new int[Rendering.values().length];

  /** No public construction. */
  private Shape3DHelper() {}

  /**
   * Creates the shape.
   *
   * @param rendering the rendering
   * @param colorDepth the color depth
   * @return the shape
   */
  public static Shape3D createShape(Rendering rendering, int colorDepth) {
    GeometryArray ga;
    final Appearance appearance = new Appearance();

    int vertexFormat = GeometryArray.COORDINATES;
    if (colorDepth == 3) {
      vertexFormat |= GeometryArray.COLOR_3;
    } else if (colorDepth == 4) {
      vertexFormat |= GeometryArray.COLOR_4;
    }

    // Support drawing as points ...
    if (rendering == Rendering.POINT) {
      ga = new PointArray(1, vertexFormat);

      final PointAttributes pa = new PointAttributes();
      pa.setPointAntialiasingEnable(true);
      appearance.setPointAttributes(pa);
    } else {
      ga = createGeometryArray(rendering, colorDepth);

      // // Test using the geometry from a sphere primitive
      // switch (r)
      // {
      // case HIGH_RES_SPHERE:
      // ga = ItemGeometryGroup.createSphere(50);
      // break;
      // case LOW_RES_SPHERE:
      // ga = ItemGeometryGroup.createSphere(16);
      // break;
      // }

      final PolygonAttributes pa = new PolygonAttributes();
      pa.setPolygonMode(PolygonAttributes.POLYGON_FILL);
      if (rendering.is2D()) {
        pa.setCullFace(PolygonAttributes.CULL_NONE);
        pa.setBackFaceNormalFlip(true);
      } else {
        pa.setCullFace(PolygonAttributes.CULL_BACK);
        pa.setBackFaceNormalFlip(false);
      }
      appearance.setPolygonAttributes(pa);

      final ColoringAttributes ca = new ColoringAttributes();
      // if (rendering.isHighResolution() || rendering.is2D())
      // // Smooth across vertices. Required to show 2D surfaces smoothly
      // ca.setShadeModel(ColoringAttributes.SHADE_GOURAUD);
      // else
      // // Faster polygon rendering with flat shading
      // ca.setShadeModel(ColoringAttributes.SHADE_FLAT);
      // For indexed models (with 1 normal per vertex) always use smooth shading
      ca.setShadeModel(ColoringAttributes.SHADE_GOURAUD);
      appearance.setColoringAttributes(ca);

      final Material m = new Material();
      m.setShininess(128f);
      m.setAmbientColor(0.1f, 0.1f, 0.1f);
      if (rendering.isHighResolution()) {
        // Allow shiny highlights on balls
        m.setSpecularColor(0.1f, 0.1f, 0.1f);
      } else {
        // For flat appearance
        m.setSpecularColor(0, 0, 0);
      }
      appearance.setMaterial(m);
    }

    return new Shape3D(ga, appearance);
  }

  /**
   * Gets the normals assuming triangle vertices.
   *
   * @param vertices the vertices
   * @param creaseAngle the crease angle (in degrees; 0=facet normals; 180=smooth shading)
   * @return the normals
   */
  public static Vector3f[] getNormals(Point3f[] vertices, double creaseAngle) {
    final int nVertices = vertices.length;
    final Vector3f[] normals = new Vector3f[nVertices];

    final GeometryArray ta = new TriangleArray(nVertices, GeometryArray.COORDINATES);
    ta.setCoordinates(0, vertices);
    final GeometryInfo gi = new GeometryInfo(ta);
    final NormalGenerator ng = new NormalGenerator();
    if (creaseAngle >= 0 && creaseAngle <= 180) {
      ng.setCreaseAngle(Math.toRadians(creaseAngle));
    }
    ng.generateNormals(gi);
    final Vector3f[] n = gi.getNormals();
    final int[] indices = gi.getNormalIndices();
    for (int i = 0; i < nVertices; i++) {
      normals[i] = n[indices[i]];
    }

    return normals;
  }

  /**
   * Gets the normals assuming triangle vertices.
   *
   * @param vertices the vertices
   * @param creaseAngle the crease angle (in degrees; 0=facet normals; 180=smooth shading)
   * @return the normals
   */
  public static Pair<Vector3f[], int[]> getIndexedNormals(Point3f[] vertices, double creaseAngle) {
    final int nVertices = vertices.length;
    final GeometryArray ta = new TriangleArray(nVertices, GeometryArray.COORDINATES);
    ta.setCoordinates(0, vertices);
    final GeometryInfo gi = new GeometryInfo(ta);
    final NormalGenerator ng = new NormalGenerator();
    if (creaseAngle >= 0 && creaseAngle <= 180) {
      ng.setCreaseAngle(Math.toRadians(creaseAngle));
    }
    ng.generateNormals(gi);
    return Pair.of(gi.getNormals(), gi.getNormalIndices());
  }

  /**
   * Creates the object used to draw a single localisation.
   *
   * @param rendering the rendering
   * @return the list of triangle vertices for the object
   */
  public static List<Point3f> createLocalisationObject(Rendering rendering) {
    switch (rendering) {
      case CUBE:
        return createCube();
      case SQUARE:
        return createSquare();
      case HEXAGON:
      case LOW_RES_CIRCLE:
      case HIGH_RES_CIRCLE:
        return createDisc(0, 0, 0, 0, 0, 1, 1, getDiscSubdivisions(rendering));

      // All handle the same way
      case ICOSAHEDRON:
      case LOW_RES_SPHERE:
      case HIGH_RES_SPHERE:
      case SUPER_HIGH_RES_SPHERE:
        return customnode.MeshMaker.createIcosahedron(getIcosahedronSubdivisions(rendering), 1f);

      case POINT:
      default:
        throw new IllegalStateException("Unknown rendering " + rendering);
    }
  }

  private static int getDiscSubdivisions(Rendering rendering) {
    switch (rendering) {
      case HIGH_RES_CIRCLE:
        return 20;
      case LOW_RES_CIRCLE:
        return 12;
      case HEXAGON:
        return 6;
      default:
        throw new IllegalStateException("Unsupporterd rendering: " + rendering);
    }
  }

  private static int getIcosahedronSubdivisions(Rendering rendering) {
    switch (rendering) {
      case ICOSAHEDRON:
        return 0;
      case LOW_RES_SPHERE:
        return 1;
      case HIGH_RES_SPHERE:
        return 2;
      case SUPER_HIGH_RES_SPHERE:
        return 3;
      default:
        throw new IllegalStateException("Unsupporterd rendering: " + rendering);
    }
  }

  /**
   * Creates the disc. This is copied from MeshMaker but the duplication of the vertices for both
   * sides on the disc is removed.
   *
   * @param x the x
   * @param y the y
   * @param z the z
   * @param nx the nx
   * @param ny the ny
   * @param nz the nz
   * @param radius the radius
   * @param edgePoints the edge points
   * @return the list
   */
  public static List<Point3f> createDisc(final double x, final double y, final double z,
      final double nx, final double ny, final double nz, final double radius,
      final int edgePoints) {
    double ax;
    double ay;
    double az;

    if (Math.abs(nx) >= Math.abs(ny)) {
      final double scale = 1 / Math.sqrt(nx * nx + nz * nz);
      ax = -nz * scale;
      ay = 0;
      az = nx * scale;
    } else {
      final double scale = 1 / Math.sqrt(ny * ny + nz * nz);
      ax = 0;
      ay = nz * scale;
      az = -ny * scale;
    }

    /*
     * Now to find the other vector in that plane, do the cross product of (ax,ay,az) with
     * (nx,ny,nz)
     */

    double bx = (ay * nz - az * ny);
    double by = (az * nx - ax * nz);
    double bz = (ax * ny - ay * nx);
    final double bScale = 1 / Math.sqrt(bx * bx + by * by + bz * bz);
    bx *= bScale;
    by *= bScale;
    bz *= bScale;

    final double[] circleX = new double[edgePoints + 1];
    final double[] circleY = new double[edgePoints + 1];
    final double[] circleZ = new double[edgePoints + 1];

    for (int i = 0; i < edgePoints + 1; ++i) {
      final double angle = (i * 2 * Math.PI) / edgePoints;
      final double c = Math.cos(angle);
      final double s = Math.sin(angle);
      circleX[i] = x + radius * c * ax + radius * s * bx;
      circleY[i] = y + radius * c * ay + radius * s * by;
      circleZ[i] = z + radius * c * az + radius * s * bz;
    }
    final LocalList<Point3f> list = new LocalList<>();
    final Point3f centre = new Point3f((float) x, (float) y, (float) z);
    for (int i = 0; i < edgePoints; ++i) {
      final Point3f t2 = new Point3f((float) circleX[i], (float) circleY[i], (float) circleZ[i]);
      final Point3f t3 =
          new Point3f((float) circleX[i + 1], (float) circleY[i + 1], (float) circleZ[i + 1]);
      list.add(centre);
      list.add(t2);
      list.add(t3);

      // We do not duplicate the triangle for both sides as we render the object as 2D
      // with setBackFaceNormalFlip(true)
      // list.add(centre);
      // list.add(t3);
      // list.add(t2);
    }
    return list;
  }

  // Note: The triangles are rendered using a right-hand coordinate system.
  // However for 2D shapes the handedness does matter as we set back-face cull off.
  // For polygons we check the handedness is
  // facing away from the centre so back-face cull can be on.

  private static float sqrt(double value) {
    return (float) Math.sqrt(value);
  }

  /**
   * Creates the triangle with vertices on a unit sphere.
   *
   * @return the list of vertices for the triangles
   */
  @SuppressWarnings("unused")
  private static List<Point3f> createTriangle() {
    return createSolid(TRI_VERTICES, TRI_FACES, true);
  }

  /**
   * Creates the square with vertices on a unit square.
   *
   * @return the list of vertices for the triangles
   */
  private static List<Point3f> createSquare() {
    return createSolid(SQUARE_VERTICES, SQUARE_FACES, true);
  }

  /**
   * Creates the tetrahedron with vertices on a unit sphere.
   *
   * @return the list of vertices for the triangles
   */
  @SuppressWarnings("unused")
  private static List<Point3f> createTetrahedron() {
    return createSolid(TETRA_VERTICES, TETRA_FACES, true);
  }

  /**
   * Creates the octahedron with vertices on a unit cube.
   *
   * @return the list of vertices for the triangles
   */
  @SuppressWarnings("unused")
  private static List<Point3f> createOctahedron() {
    // This is already normalised
    return createSolid(OCTA_VERTICES, OCTA_FACES, false);
  }

  /**
   * Creates the cube with vertices on a unit cube.
   *
   * @return the list of vertices for the triangles
   */
  private static List<Point3f> createCube() {
    return createSolid(CUBE_VERTICES, CUBE_FACES, true);
  }

  /**
   * Creates the solid with the defined faces and vertices on a unit sphere.
   *
   * @param vertices the vertices
   * @param faces the faces
   * @param normalise the normalise
   * @return the list of vertices for the triangles
   */
  private static List<Point3f> createSolid(float[][] vertices, int[][] faces, boolean normalise) {
    final List<Point3f> ps = new LocalList<>();
    for (final int[] face : faces) {
      for (int k = 0; k < 3; k++) {
        ps.add(new Point3f(vertices[face[k]]));
      }
    }
    // Project all vertices to the surface of a sphere of radius 1
    if (normalise) {
      normalise(ps);
    }

    return ps;
  }

  private static void normalise(List<Point3f> ps) {
    final Vector3f v = new Vector3f();
    for (final Point3f p : ps) {
      v.set(p);
      v.normalize();
      p.set(v);
    }
  }

  /**
   * Creates the object used to outline a single localisation.
   *
   * @param rendering the rendering
   * @return the list of triangle vertices for the object
   */
  public static List<Point3f> createLocalisationObjectOutline(Rendering rendering) {
    switch (rendering) {
      case SQUARE:
        return createSquareOutline();
      case HEXAGON:
      case LOW_RES_CIRCLE:
      case HIGH_RES_CIRCLE:
        return createDiscOutline(0, 0, 0, 0, 0, 1, 1, getDiscSubdivisions(rendering));

      default:
        return createLocalisationObject(rendering);
    }
  }

  /**
   * Creates the triangle with vertices on a unit sphere.
   *
   * @return the list of vertices for the triangles
   */
  @SuppressWarnings("unused")
  private static List<Point3f> createTriangleOutline() {
    return createSolidOutline(TRI_VERTICES, true);
  }

  /**
   * Creates the square with vertices on a unit sphere.
   *
   * @return the list of vertices for the triangles
   */
  private static List<Point3f> createSquareOutline() {
    return createSolidOutline(SQUARE_VERTICES, true);
  }

  /**
   * Creates the solid with the defined faces and vertices on a unit sphere.
   *
   * @param vertices the vertices
   * @param normalise the normalise
   * @return the list of vertices for the triangles
   */
  private static List<Point3f> createSolidOutline(float[][] vertices, boolean normalise) {
    final List<Point3f> ps = new LocalList<>();
    for (final float[] vertex : vertices) {
      ps.add(new Point3f(vertex));
    }
    // Make continuous
    ps.add(new Point3f(vertices[0]));
    if (normalise) {
      normalise(ps);
    }
    return ps;
  }

  /**
   * Creates the disc outline.
   *
   * @param x the x
   * @param y the y
   * @param z the z
   * @param nx the nx
   * @param ny the ny
   * @param nz the nz
   * @param radius the radius
   * @param edgePoints the edge points
   * @return the list
   */
  public static List<Point3f> createDiscOutline(final double x, final double y, final double z,
      final double nx, final double ny, final double nz, final double radius,
      final int edgePoints) {
    return createDiscEdge(false, x, y, z, nx, ny, nz, radius, edgePoints);
  }

  /**
   * Creates the disc outline with the first point at the centre, i.e. suitable for use as a
   * triangle fan array.
   *
   * @param x the x
   * @param y the y
   * @param z the z
   * @param nx the nx
   * @param ny the ny
   * @param nz the nz
   * @param radius the radius
   * @param edgePoints the edge points
   * @return the list
   */
  public static List<Point3f> createDiscFan(final double x, final double y, final double z,
      final double nx, final double ny, final double nz, final double radius,
      final int edgePoints) {
    return createDiscEdge(true, x, y, z, nx, ny, nz, radius, edgePoints);
  }

  /**
   * Creates the disc. This is copied from MeshMaker but the duplication of the vertices for both
   * sides on the disc is removed.
   *
   * @param x the x
   * @param y the y
   * @param z the z
   * @param nx the nx
   * @param ny the ny
   * @param nz the nz
   * @param radius the radius
   * @param edgePoints the edge points
   * @return the list
   */
  private static List<Point3f> createDiscEdge(boolean includeCenter, final double x, final double y,
      final double z, final double nx, final double ny, final double nz, final double radius,
      final int edgePoints) {
    double ax;
    double ay;
    double az;

    if (Math.abs(nx) >= Math.abs(ny)) {
      final double scale = 1 / Math.sqrt(nx * nx + nz * nz);
      ax = -nz * scale;
      ay = 0;
      az = nx * scale;
    } else {
      final double scale = 1 / Math.sqrt(ny * ny + nz * nz);
      ax = 0;
      ay = nz * scale;
      az = -ny * scale;
    }

    /*
     * Now to find the other vector in that plane, do the cross product of (ax,ay,az) with
     * (nx,ny,nz)
     */

    double bx = (ay * nz - az * ny);
    double by = (az * nx - ax * nz);
    double bz = (ax * ny - ay * nx);
    final double bScale = 1 / Math.sqrt(bx * bx + by * by + bz * bz);
    bx *= bScale;
    by *= bScale;
    bz *= bScale;

    final LocalList<Point3f> list = new LocalList<>();
    if (includeCenter) {
      list.add(new Point3f((float) x, (float) y, (float) z));
    }
    for (int i = edgePoints + 1; i-- > 0;) {
      // For consistency with the rendering of the icosahedron
      // we rotate by a quarter turn. The icosahedron projected
      // flat then looks like the hexagon.

      final double angle = Math.PI / 2 + (i * 2 * Math.PI) / edgePoints;
      final double c = Math.cos(angle);
      final double s = Math.sin(angle);
      final float px = (float) (x + radius * c * ax + radius * s * bx);
      final float py = (float) (y + radius * c * ay + radius * s * by);
      final float pz = (float) (z + radius * c * az + radius * s * bz);
      list.add(new Point3f(px, py, pz));
    }
    return list;
  }

  /**
   * Creates the object used to draw a single localisation.
   *
   * @param rendering the rendering
   * @param colorDepth the color depth
   * @return the geometry array
   */
  public static GeometryArray createGeometryArray(Rendering rendering, int colorDepth) {
    final GeometryInfo gi = createGeometryInfo(rendering, colorDepth);
    final boolean useCoordIndexOnly = gi.getUseCoordIndexOnly();
    return rendering.is2D() ? gi.getGeometryArray()
        : gi.getIndexedGeometryArray(false, false, false, useCoordIndexOnly, false);
  }

  /**
   * Creates the object used to draw a single localisation.
   *
   * @param rendering the rendering
   * @param colorDepth the color depth
   * @return the geometry info
   */
  public static GeometryInfo createGeometryInfo(Rendering rendering, int colorDepth) {
    List<Point3f> coords;
    int primitive = GeometryInfo.TRIANGLE_ARRAY;
    int[] stripsCounts = null;
    final Vector3f normal = new Vector3f();
    boolean normalise = false;
    switch (rendering) {
      // 2D
      case SQUARE:
        primitive = GeometryInfo.QUAD_ARRAY;
        coords = new LocalList<>();
        for (int i = 0; i < 4; i++) {
          coords.add(new Point3f(CUBE_VERTICES[i][0], CUBE_VERTICES[i][1], 0));
        }
        normalise = true;
        normal.set(0, 0, 1);
        break;

      case HEXAGON:
      case LOW_RES_CIRCLE:
      case HIGH_RES_CIRCLE:
        primitive = GeometryInfo.TRIANGLE_FAN_ARRAY;
        coords = createDiscFan(0, 0, 0, 0, 0, 1, 1, getDiscSubdivisions(rendering));
        stripsCounts = new int[] {coords.size()};
        normal.set(0, 0, 1);
        break;

      // 3D
      case CUBE:
        primitive = GeometryInfo.QUAD_ARRAY;
        coords = new LocalList<>();
        final Point3f[] vertices = new Point3f[8];
        for (int i = 0; i < 8; i++) {
          vertices[i] = new Point3f(CUBE_VERTICES[i][0], CUBE_VERTICES[i][1], CUBE_VERTICES[i][2]);
        }
        for (final int[] face : CUBE_FACES4) {
          for (int i = 0; i < 4; i++) {
            coords.add(vertices[face[i]]);
          }
        }
        normalise = true;
        break;

      // All spheres based on icosahedron for speed
      case ICOSAHEDRON:
      case LOW_RES_SPHERE:
      case HIGH_RES_SPHERE:
      case SUPER_HIGH_RES_SPHERE:
        coords = customnode.MeshMaker.createIcosahedron(getIcosahedronSubdivisions(rendering), 1f);
        break;

      case POINT:
      default:
        throw new IllegalStateException("Unknown rendering " + rendering);
    }

    final GeometryInfo gi = new GeometryInfo(primitive);
    gi.setUseCoordIndexOnly(true);

    if (normalise) {
      normalise(coords);
    }

    if (rendering.is2D()) {
      gi.setCoordinates(coords.toArray(new Point3f[0]));
      gi.setStripCounts(stripsCounts);

      if (colorDepth == 3) {
        gi.setColors3(new float[coords.size() * 3]);
      } else if (colorDepth == 4) {
        gi.setColors4(new float[coords.size() * 4]);
      }

      // Normals generated with the normal generator add extra normals for indexed arrays
      // so we do it manually
      final Vector3f[] normals = new Vector3f[coords.size()];
      Arrays.fill(normals, normal);
      gi.setNormals(normals);
    } else {
      // Generate indexes manually.
      // The GeometryInfo somehow does not do this correctly for all rendering modes.
      // E.g. the icosahedron gets extra indexes to normals that are 0,0,0.
      final Pair<Point3f[], int[]> p = createIndexedObject(coords);

      final Point3f[] iCoords = p.getKey();
      gi.setCoordinates(iCoords);
      gi.setCoordinateIndices(p.getValue());

      // Normals are just the vector from 0,0,0 to the vertex
      final Vector3f[] normals = new Vector3f[iCoords.length];
      for (int i = 0; i < normals.length; i++) {
        normal.set(iCoords[i]);
        normal.normalize();
        normals[i] = new Vector3f(normal);
      }

      gi.setNormals(normals);
      // gi.setNormalIndices(p.b);

      if (colorDepth == 3 || colorDepth == 4) {
        if (colorDepth == 3) {
          gi.setColors3(new float[iCoords.length * 3]);
        } else {
          gi.setColors4(new float[iCoords.length * 4]);
          // gi.setColorIndices(p.b);
        }
      }

      // final NormalGenerator ng = new NormalGenerator();
      // double creaseAngle = (rendering.isHighResolution()) ? Math.PI : 0;
      // ng.setCreaseAngle(creaseAngle);
      // //if (rendering== Rendering.ICOSAHEDRON)
      // // ng.setCreaseAngle(180);
      // ng.generateNormals(gi);
    }

    return gi;
  }

  /**
   * Creates an indexed object from a list of vertices.
   *
   * @param list the list of vertices
   * @return the vertices and indices of the object
   */
  public static Pair<Point3f[], int[]> createIndexedObject(List<Point3f> list) {
    // Compact the vertices to a set of vertices and faces
    final Object2IntOpenHashMap<Point3f> m = new Object2IntOpenHashMap<>(list.size());
    m.defaultReturnValue(-1);
    final LocalList<Point3f> vertices = new LocalList<>(list.size());
    final IntArrayList faces = new IntArrayList(list.size());
    int index = 0;
    // Process triangles
    for (final Point3f p : list) {
      int value = m.putIfAbsent(p, index);
      if (value == -1) {
        // Store the points in order
        vertices.add(p);
        value = index++;
      }
      faces.add(value);
    }

    return Pair.of(vertices.toArray(new Point3f[0]), faces.toIntArray());
  }

  /**
   * Gets the number of triangles.
   *
   * @param ga the geometry array
   * @return the number of triangles
   */
  public static int getNumberOfTriangles(GeometryArray ga) {
    final GeometryInfo gi = new GeometryInfo(ga);
    gi.convertToIndexedTriangles();
    return gi.getCoordinateIndices().length / 3;
  }

  /**
   * Gets the number of triangles.
   *
   * @param rendering the rendering
   * @return the number of triangles
   */
  public static int getNumberOfTriangles(Rendering rendering) {
    if (rendering == Rendering.POINT) {
      return 0;
    }
    final int index = rendering.ordinal();
    if (numberOfTriangles[index] == 0) {
      final GeometryInfo gi = createGeometryInfo(rendering, 0);
      // Remove so the conversion is faster
      gi.setColors((Color3f[]) null);
      gi.setColorIndices(null);
      gi.setNormals((Vector3f[]) null);
      gi.setNormalIndices(null);
      // this is required before conversion
      gi.setUseCoordIndexOnly(false);
      gi.convertToIndexedTriangles();
      numberOfTriangles[index] = gi.getCoordinateIndices().length / 3;
    }
    return numberOfTriangles[index];
  }
}
