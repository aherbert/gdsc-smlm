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

import customnode.CustomMesh;
import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;
import java.util.Objects;
import java.util.logging.Level;
import org.jogamp.java3d.Appearance;
import org.jogamp.java3d.ColoringAttributes;
import org.jogamp.java3d.Geometry;
import org.jogamp.java3d.GeometryArray;
import org.jogamp.java3d.GeometryStripArray;
import org.jogamp.java3d.GeometryUpdater;
import org.jogamp.java3d.IndexedGeometryArray;
import org.jogamp.java3d.IndexedGeometryStripArray;
import org.jogamp.java3d.IndexedPointArray;
import org.jogamp.java3d.Material;
import org.jogamp.java3d.PointArray;
import org.jogamp.java3d.PointAttributes;
import org.jogamp.java3d.PolygonAttributes;
import org.jogamp.java3d.TransparencyAttributes;
import org.jogamp.vecmath.Color3f;
import org.jogamp.vecmath.Color4f;
import org.jogamp.vecmath.Point3f;
import uk.ac.sussex.gdsc.core.data.NotImplementedException;
import uk.ac.sussex.gdsc.core.ij.ImageJPluginLoggerHelper;
import uk.ac.sussex.gdsc.core.utils.BitFlagUtils;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;

/**
 * Use a mesh object to represent a set of points. The object is duplicated, scaled and translated
 * for each point.
 *
 * <p>Note: TransparentItemShape is only supported if {@link #hasColor4()} is true, i.e. the input
 * geometry has per vertex colours with alpha.
 */
public class ItemMesh extends CustomMesh implements UpdateableItemShape, TransparentItemShape {
  private static int transparencyMode = TransparencyAttributes.FASTEST;

  /** The vertex count of the original geometry array. */
  protected final int vertexCount;

  /** The vertex format of the original geometry array. */
  protected final int vertexFormat;

  /** The index count of the original geometry array. */
  protected final int indexCount;

  /** The points. */
  protected Point3f[] points;

  /** The size of each point. */
  protected Point3f[] sizes;

  /** Set to true if this is a point array. */
  protected final boolean isPointArray;

  /** Set to true if this is a strip geometry array. */
  protected final boolean isStripGeometryArray;

  /**
   * If there are no per-vertex colours this is set to true to colour the object using the material,
   * false uses the ColoringAttributes.
   */
  protected boolean isColorByMaterial;

  /** The color updater used to update the colour array for a single item. */
  protected ArrayColorUpdater colorUpdater;

  /**
   * Instantiates a new item mesh.
   *
   * <p>This will repeat the object for each input point. The object is assumed to be centred on the
   * origin. It will be scaled and translated for each input point.
   *
   * <p>The input geometry array vertex format is checked and unsupported formats throw an
   * exception. Currently this only supports coordinates, normals and color. Strip and Fan arrays
   * are supported.
   *
   * @param point the point
   * @param ga the geometry array. If null then a default will be used. Assumed to be centred on the
   *        origin.
   * @param appearance the default appearance of the shape. PolygonAttributes, Material and
   *        TransparencyAttributes are used.
   * @param sizes the sizes of each point. Can be null (no scaling); length=1 (fixed scaling); or
   *        points.length.
   * @param color the color
   * @param transparency the transparency
   */
  public ItemMesh(Point3f point, GeometryArray ga, Appearance appearance, Point3f[] sizes,
      final Color3f color, final float transparency) {
    this(new Point3f[] {point}, ga, appearance, sizes, color, transparency);
  }

  /**
   * Instantiates a new item mesh.
   *
   * <p>This will repeat the object for each input point. The object is assumed to be centred on the
   * origin. It will be scaled and translated for each input point.
   *
   * <p>The input geometry array vertex format is checked and unsupported formats throw an
   * exception. Currently this only supports coordinates, normals and color. Strip and Fan arrays
   * are supported.
   *
   * @param points the points
   * @param ga the geometry array. If null then a default will be used. Assumed to be centred on the
   *        origin.
   * @param appearance the default appearance of the shape. PolygonAttributes, Material and
   *        TransparencyAttributes are used.
   * @param sizes the sizes of each point. Can be null (no scaling); length=1 (fixed scaling); or
   *        points.length.
   * @param color the color
   * @param transparency the transparency
   */
  public ItemMesh(Point3f[] points, GeometryArray ga, Appearance appearance, Point3f[] sizes,
      final Color3f color, final float transparency) {
    // Create empty
    super(null, color, transparency);

    if (sizes != null && points.length != sizes.length && sizes.length != 1) {
      throw new IllegalArgumentException("Points and sizes must be the same length");
    }

    if (ga.getVertexAttrCount() > 0) {
      throw new IllegalArgumentException("Vertex attributes are not supported");
    }
    if (ga.getTexCoordSetCount() > 0) {
      throw new IllegalArgumentException("Texture coordinates are not supported");
    }

    vertexCount = ga.getValidVertexCount();
    vertexFormat = ga.getVertexFormat();
    if (ga instanceof IndexedGeometryArray) {
      indexCount = ((IndexedGeometryArray) ga).getValidIndexCount();
    } else {
      indexCount = 0;
    }
    isPointArray = (ga instanceof PointArray || ga instanceof IndexedPointArray);
    isStripGeometryArray =
        (ga instanceof IndexedGeometryStripArray || ga instanceof GeometryStripArray);

    // Set the flags we support and error if there are others
    int flags = GeometryArray.COORDINATES;
    if ((vertexFormat & GeometryArray.NORMALS) != 0) {
      flags |= GeometryArray.NORMALS;
    }
    if (hasColor4()) {
      flags |= GeometryArray.COLOR_4;
    } else if (hasColor()) {
      flags |= GeometryArray.COLOR_3;
    }
    if ((vertexFormat & GeometryArray.USE_COORD_INDEX_ONLY) != 0) {
      flags |= GeometryArray.USE_COORD_INDEX_ONLY;
    }

    final int extra = BitFlagUtils.unset(vertexFormat, flags);
    if (extra != 0) {
      throw new IllegalArgumentException("Unsupported vertex format flags: " + extra);
    }

    this.points = points;
    this.sizes = sizes;

    // Now build the actual vertices by repeating the points.
    final float[] objectCoords = new float[vertexCount * 3];
    ga.getCoordinates(0, objectCoords);

    final float[] allCoords = new float[objectCoords.length * points.length];

    boolean sameSize = false;
    if (sizes == null || (sameSize = sameSize(sizes))) {
      if (sameSize) {
        // Scale the input object
        Objects.requireNonNull(sizes, "sizes should not be null here");
        final Point3f s = sizes[0];
        final float sx = s.x;
        final float sy = s.y;
        final float sz = s.z;
        for (int j = 0; j < objectCoords.length; j += 3) {
          objectCoords[j] *= sx;
          objectCoords[j + 1] *= sy;
          objectCoords[j + 2] *= sz;
        }
      }

      // Translate
      for (int i = 0, k = 0; i < points.length; i++) {
        final Point3f p = points[i];
        final float dx = p.x;
        final float dy = p.y;
        final float dz = p.z;
        for (int j = 0; j < objectCoords.length; j += 3) {
          allCoords[k++] = objectCoords[j] + dx;
          allCoords[k++] = objectCoords[j + 1] + dy;
          allCoords[k++] = objectCoords[j + 2] + dz;
        }
      }
    } else {
      // Translate and scale
      for (int i = 0, k = 0; i < points.length; i++) {
        final Point3f p = points[i];
        final float dx = p.x;
        final float dy = p.y;
        final float dz = p.z;
        final Point3f s = sizes[i];
        final float sx = s.x;
        final float sy = s.y;
        final float sz = s.z;
        for (int j = 0; j < objectCoords.length; j += 3) {
          allCoords[k++] = objectCoords[j] * sx + dx;
          allCoords[k++] = objectCoords[j + 1] * sy + dy;
          allCoords[k++] = objectCoords[j + 2] * sz + dz;
        }
      }
    }

    // Update the geometry.
    // Do this in a method to allow sub-classes to change the geometry.
    this.setGeometry(createGeometry(allCoords, ga));

    // Create a default appearance
    setAppearance(createAppearance(appearance, ga));

    // Initialise
    if (hasColor4()) {
      setItemAlpha(1f);
    }
    setColor(color);
    setTransparency(transparency);
  }

  @Override
  public void update() {
    // Ignore this
  }

  /**
   * Checks for normals.
   *
   * @return true, if successful
   */
  public boolean hasNormals() {
    return ((vertexFormat & GeometryArray.NORMALS) != 0);
  }

  /**
   * Checks for color.
   *
   * @return true, if successful
   */
  public boolean hasColor() {
    return ((vertexFormat & GeometryArray.COLOR_3) == GeometryArray.COLOR_3);
  }

  /**
   * Checks for color 3.
   *
   * @return true, if successful
   */
  public boolean hasColor3() {
    return hasColor() && !hasColor4();
  }

  /**
   * Checks for color 4.
   *
   * @return true, if successful
   */
  public boolean hasColor4() {
    return ((vertexFormat & GeometryArray.COLOR_4) == GeometryArray.COLOR_4);
  }

  /**
   * Gets the number of colours per vertex (0, 3, or 4).
   *
   * @return the number of colors per vertex
   */
  public int getNumberOfColorsPerVertex() {
    return (hasColor4()) ? 4 : (hasColor()) ? 3 : 0;
  }

  /**
   * Checks if is index geometry array.
   *
   * @return true, if is index geometry array
   */
  public boolean isIndexGeometryArray() {
    return (indexCount > 0);
  }

  /**
   * Gets the vertices per item.
   *
   * @return the vertices per item
   */
  public int getVerticesPerItem() {
    return (indexCount > 0) ? indexCount : vertexCount;
  }

  /**
   * Gets the vertices count per item. This may be lower than the vertices per item if using an
   * indexed array.
   *
   * @return the vertices count per item
   */
  public int getVerticesCountPerItem() {
    return vertexCount;
  }

  /**
   * Checks if is point array.
   *
   * @return true, if is point array
   */
  public boolean isPointArray() {
    return isPointArray;
  }

  /**
   * Checks if is strip geometry array.
   *
   * @return true, if is strip geometry array
   */
  public boolean isStripGeometryArray() {
    return isStripGeometryArray;
  }

  /**
   * Check if all the points are the same size.
   *
   * @param sizes the sizes (must be an array of at least 1)
   * @return true, if successful
   */
  public static boolean sameSize(Point3f[] sizes) {
    if (sizes.length == 1) {
      return true;
    }
    final Point3f s = sizes[0];
    for (int j = 1; j < sizes.length; j++) {
      if (!sizes[j].equals(s)) {
        return false;
      }
    }
    return true;
  }

  @Override
  public void setCoordinate(final int index, final Point3f point) {
    throw new NotImplementedException();
  }

  @Override
  public void setCoordinates(final int[] indices, final Point3f point) {
    throw new NotImplementedException();
  }

  @Override
  protected void addVertices(Point3f[] vertices) {
    throw new NotImplementedException();
  }

  @Override
  protected void addVerticesToGeometryArray(Point3f[] vertices) {
    throw new NotImplementedException();
  }

  @Override
  protected void addVerticesToGeometryStripArray(Point3f[] vertices) {
    throw new NotImplementedException();
  }

  @Override
  protected void removeVertices(int[] indices) {
    throw new NotImplementedException();
  }

  @Override
  protected GeometryArray createGeometry() {
    throw new NotImplementedException();
  }

  /**
   * Creates the geometry.
   *
   * @param coords the coords
   * @param sourceGa the source geometry array
   * @return the geometry array
   */
  protected GeometryArray createGeometry(float[] coords, GeometryArray sourceGa) {
    final GeometryArray ga = createGeometryArray(sourceGa, 0);

    ga.setCoordinates(0, coords);
    ga.setCapability(GeometryArray.ALLOW_COORDINATE_WRITE);

    // Handle normals
    boolean doNormals = hasNormals();

    // Handle colors
    if (hasColor()) {
      ga.setCapability(GeometryArray.ALLOW_COLOR_READ);
      ga.setCapability(GeometryArray.ALLOW_COLOR_WRITE);
      colorUpdater = ArrayColorUpdater.create(vertexCount, hasColor4());
    }

    // Fan are extensions of GeometryStripArray so do not need extra code.

    // Handle indexed array
    if (isIndexGeometryArray()) {
      final IndexedGeometryArray sourceIga = (IndexedGeometryArray) sourceGa;
      final IndexedGeometryArray iga = (IndexedGeometryArray) ga;
      final int objectIndexCount = sourceIga.getValidIndexCount();
      final int[] objectIndices = new int[objectIndexCount];
      final int[] allIndices = new int[objectIndices.length * points.length];
      sourceIga.getCoordinateIndices(0, objectIndices);
      duplicateIndices(objectIndices, allIndices);
      iga.setCoordinateIndices(0, allIndices);

      // Check if we need the color and normal indices
      if ((vertexFormat & GeometryArray.USE_COORD_INDEX_ONLY) == 0) {
        if (hasNormals()) {
          // Use the same index for all vertices for normals
          sourceIga.getNormalIndices(0, objectIndices);
          duplicate(objectIndices, 0, objectIndices.length, points.length, allIndices, 0);
          iga.setNormalIndices(0, allIndices);

          final float[] normals = new float[(MathUtils.max(objectIndices) + 1) * 3];
          sourceIga.getNormals(0, normals);
          iga.setNormals(0, normals);

          doNormals = false;
        }

        if (hasColor()) {
          // Use a single index per item for vertex colour
          for (int i = 0, k = 0; i < points.length; i++) {
            for (int j = 0; j < objectIndexCount; j++) {
              allIndices[k++] = i;
            }
          }
          iga.setColorIndices(0, allIndices);
          // Only have to update a single colour per item
          colorUpdater = ArrayColorUpdater.create(1, hasColor4());
        }
      }
    }

    if (doNormals) {
      final float[] objectNormals = new float[vertexCount * 3];
      sourceGa.getNormals(0, objectNormals);
      final float[] allNormals = new float[objectNormals.length * points.length];
      duplicate(objectNormals, 0, objectNormals.length, points.length, allNormals, 0);
      ga.setNormals(0, allNormals);
    }

    return ga;
  }

  /**
   * Creates the geometry array.
   *
   * @param sourceGa the source geometry array
   * @param format the format
   * @return the geometry array
   */
  protected GeometryArray createGeometryArray(GeometryArray sourceGa, int format) {
    // Create using reflection
    final GeometryArray ga;
    try {
      final Class<?> clazz = sourceGa.getClass();
      // clazz = clazz.asSubclass(clazz);

      final LocalList<Class<?>> paramTypes = new LocalList<>(4);
      final LocalList<Object> paramValues = new LocalList<>(4);

      paramTypes.add(int.class);
      paramTypes.add(int.class);
      paramValues.add(vertexCount * points.length);
      paramValues.add(vertexFormat | format);

      if (isIndexGeometryArray()) {
        paramTypes.add(int.class);
        paramValues.add(indexCount * points.length);
      }

      // Handle strips
      int numStrips = 0;
      int[] objectStripCounts = null;
      int[] allStripCounts;
      if (sourceGa instanceof IndexedGeometryStripArray) {
        final IndexedGeometryStripArray igsa = (IndexedGeometryStripArray) sourceGa;
        numStrips = igsa.getNumStrips();
        objectStripCounts = new int[numStrips];
        igsa.getStripIndexCounts(objectStripCounts);
      } else if (sourceGa instanceof GeometryStripArray) {
        final GeometryStripArray gsa = (GeometryStripArray) sourceGa;
        numStrips = gsa.getNumStrips();
        objectStripCounts = new int[numStrips];
        gsa.getStripVertexCounts(objectStripCounts);
      }

      if (objectStripCounts != null) {
        allStripCounts = new int[numStrips * points.length];
        duplicate(objectStripCounts, 0, numStrips, points.length, allStripCounts, 0);
        paramTypes.add(int[].class);
        paramValues.add(allStripCounts);
      }

      final Class<?>[] paramTypes2 = paramTypes.toArray(new Class<?>[0]);
      final Object[] paramValues2 = paramValues.toArray();
      ga = (GeometryArray) clazz.getConstructor(paramTypes2).newInstance(paramValues2);
    } catch (final NoSuchMethodException | InstantiationException | IllegalAccessException
        | IllegalArgumentException | InvocationTargetException | SecurityException ex) {
      ImageJPluginLoggerHelper.getDefaultLogger().log(Level.WARNING, "Failed to create geometry",
          ex);
      return null;
    }

    ga.setCapability(GeometryArray.ALLOW_COUNT_WRITE);
    ga.setCapability(GeometryArray.ALLOW_COUNT_READ);
    ga.setCapability(GeometryArray.ALLOW_FORMAT_READ);
    ga.setCapability(Geometry.ALLOW_INTERSECT);

    return ga;
  }

  /**
   * Duplicate the source into the destination n times.
   *
   * @param source the source
   * @param from the from
   * @param length the length
   * @param n the number of times
   * @param dest the dest
   * @param to the to
   */
  protected void duplicate(Object source, int from, int length, int n, Object dest, int to) {
    // Binary fill
    int fill = length;
    System.arraycopy(source, from, dest, to, fill);
    for (int i = 2; i < n; i *= 2) {
      System.arraycopy(dest, to, dest, to + fill, fill);
      fill *= 2;
    }
    // Final fill
    System.arraycopy(dest, to, dest, to + fill, n * length - fill);
  }

  /**
   * Duplicate the objects indices into the all indices array for each point.
   *
   * @param objectIndices the object indices
   * @param allIndices the all indices
   * @see #size()
   */
  protected void duplicateIndices(int[] objectIndices, int[] allIndices) {
    final int nIndices = MathUtils.max(objectIndices) + 1;
    for (int i = 0, k = 0; i < points.length; i++) {
      final int offset = i * nIndices;
      for (final int objectIndex : objectIndices) {
        allIndices[k++] = objectIndex + offset;
      }
    }
  }

  /**
   * Sets the transparency mode.
   *
   * @param mode the new transparency mode
   * @throws IllegalArgumentException If the mode is not valid
   * @see TransparencyAttributes#setTransparencyMode(int)
   */
  public static void setTransparencyMode(int mode) {
    if ((mode < TransparencyAttributes.FASTEST) || (mode > TransparencyAttributes.NONE)) {
      throw new IllegalArgumentException("Not a valid transparency mode");
    }
    transparencyMode = mode;
  }

  /**
   * Gets the transparency mode.
   *
   * @return the transparency mode
   * @see TransparencyAttributes#setTransparencyMode(int)
   */
  public static int getTransparencyMode() {
    return transparencyMode;
  }

  @Override
  public void setTransparency(final float transparency) {
    // We want to use a different transparency from the ij3d default which is FASTEST
    // so override this method.
    final Appearance appearance = getAppearance();
    final TransparencyAttributes ta = appearance.getTransparencyAttributes();
    if (transparency <= ItemHelper.ZERO_TRANSPARENCY) {
      this.transparency = 0;

      // For a strange reason if this is set to NONE before the mesh is added to
      // a scene then the transparency cannot be adjusted. So set to FASTEST.
      if (ta.isLive()) {
        ta.setTransparencyMode(TransparencyAttributes.NONE);
      } else {
        ta.setTransparencyMode(TransparencyAttributes.FASTEST);
      }
    } else {
      this.transparency = transparency;
      ta.setTransparencyMode(transparencyMode);
      // ta.setTransparencyMode(TransparencyAttributes.BLENDED);
    }
    ta.setTransparency(this.transparency);
  }

  @Override
  protected Appearance createAppearance() {
    throw new NotImplementedException();
  }

  /**
   * Creates the appearance.
   *
   * @param appearance the appearance
   * @param ga the geometry array
   * @return the appearance
   */
  protected Appearance createAppearance(Appearance appearance, GeometryArray ga) {
    // Create a suitable appearance for points or 3D shapes.
    if (appearance == null) {
      appearance = new Appearance();
    }

    appearance.setCapability(Appearance.ALLOW_TRANSPARENCY_ATTRIBUTES_READ);
    appearance.setCapability(Appearance.ALLOW_MATERIAL_READ);
    appearance.setCapability(Appearance.ALLOW_COLORING_ATTRIBUTES_READ);

    // Ensure we have the ability to colour the object
    if (!hasColor()) {
      isColorByMaterial = !isPointArray;
    }

    if (isPointArray) {
      shaded = false;
      appearance.setPolygonAttributes(null);
      appearance.setMaterial(null);

      PointAttributes pointAttributes = appearance.getPointAttributes();
      if (pointAttributes == null) {
        pointAttributes = new PointAttributes();
        pointAttributes.setPointAntialiasingEnable(true);
        appearance.setPointAttributes(pointAttributes);
      }
      pointAttributes.setCapability(PointAttributes.ALLOW_ANTIALIASING_WRITE);
      pointAttributes.setCapability(PointAttributes.ALLOW_SIZE_WRITE);

      if (hasColor()) {
        // We use the coordinates for the colour
        appearance.setColoringAttributes(null);
      } else {
        ColoringAttributes ca = appearance.getColoringAttributes();
        if (ca == null) {
          ca = new ColoringAttributes();
          ca.setShadeModel(ColoringAttributes.SHADE_FLAT);
          appearance.setColoringAttributes(ca);
        }
        ca.setCapability(ColoringAttributes.ALLOW_COLOR_WRITE);
      }
    } else {
      appearance.setPointAttributes(null);

      // These are the defaults. We may need them if we want to support mesh
      // display when the polygon mode is Line
      PolygonAttributes polygonAttributes = appearance.getPolygonAttributes();
      if (polygonAttributes == null) {
        polygonAttributes = new PolygonAttributes();
        polygonAttributes.setPolygonMode(PolygonAttributes.POLYGON_FILL);
        appearance.setPolygonAttributes(polygonAttributes);
        shaded = true;
      } else {
        shaded = polygonAttributes.getPolygonMode() == PolygonAttributes.POLYGON_FILL;
      }
      polygonAttributes.setCapability(PolygonAttributes.ALLOW_MODE_WRITE);

      ColoringAttributes ca = appearance.getColoringAttributes();
      if (ca == null) {
        ca = new ColoringAttributes();
        ca.setShadeModel(ColoringAttributes.SHADE_GOURAUD);
        appearance.setColoringAttributes(ca);
      }
      ca.setCapability(ColoringAttributes.ALLOW_SHADE_MODEL_WRITE);

      Material material = appearance.getMaterial();
      if (material == null) {
        material = new Material();
        material.setAmbientColor(0.1f, 0.1f, 0.1f);
        material.setSpecularColor(0.1f, 0.1f, 0.1f);
        material.setDiffuseColor(DEFAULT_COLOR);
        appearance.setMaterial(material);
      }
      // Ensure per vertex colours replace the diffuse colour
      material.setColorTarget(Material.DIFFUSE);
      if (isColorByMaterial) {
        material.setCapability(Material.ALLOW_COMPONENT_WRITE);
      }
    }

    // We require transparency attributes for global transparency
    TransparencyAttributes tr = appearance.getTransparencyAttributes();
    if (tr == null) {
      tr = new TransparencyAttributes();
      tr.setTransparencyMode(TransparencyAttributes.NONE);
      tr.setTransparency(0f);
      appearance.setTransparencyAttributes(tr);
    }
    tr.setCapability(TransparencyAttributes.ALLOW_VALUE_WRITE);
    tr.setCapability(TransparencyAttributes.ALLOW_MODE_WRITE);

    return appearance;
  }

  @Override
  public void reorder(int[] indices) {
    ItemPointMesh.checkIndices(indices, points.length);
    reorderFast(indices);
  }

  @Override
  public void reorderFast(int[] indices) {
    changed = true;

    final int oldSize = size();
    final int size = (indices == null) ? 0 : Math.min(oldSize, indices.length);

    if (size == 0 || indices == null) {
      points = new Point3f[0];
      sizes = new Point3f[0];
      this.setGeometry(null);
      return;
    }

    // From here on we assume the current geometry will not be null
    // as this only happens when the original size is zero. Size has
    // been checked at this point to be the smaller of new and old.
    final GeometryArray ga = (GeometryArray) getGeometry();

    points = reorderPoints(points, indices);
    // Sizes could be null or a single size
    if (sizes != null && sizes.length == points.length) {
      sizes = reorderPoints(sizes, indices);
    }

    // Reorder all things in the geometry: coordinates and colour.
    // The normals, indices, strip counts are are unchanged.
    // int objectSize = vertexCount;

    int countPerObject = vertexCount * 3;
    final float[] oldCoords = new float[oldSize * countPerObject];
    ga.getCoordinates(0, oldCoords);
    final float[] coords = new float[size * countPerObject];
    for (int i = 0; i < size; i++) {
      final int j = indices[i];
      final int ii = i * countPerObject;
      final int jj = j * countPerObject;
      System.arraycopy(oldCoords, jj, coords, ii, countPerObject);
    }

    final float[] colors;
    if (hasColor()) {
      countPerObject = colorUpdater.size();
      final int colorSize = oldSize * countPerObject;
      final float[] oldColors = (colorSize < oldCoords.length) ? oldCoords : new float[colorSize];
      ga.getColors(0, oldColors);
      colors = new float[size * countPerObject];
      for (int i = 0; i < size; i++) {
        final int j = indices[i];
        final int ii = i * countPerObject;
        final int jj = j * countPerObject;
        System.arraycopy(oldColors, jj, colors, ii, countPerObject);
      }
    } else {
      colors = null;
    }

    ga.updateData(new GeometryUpdater() {
      @Override
      public void updateData(Geometry geometry) {
        final GeometryArray ga = (GeometryArray) geometry;
        // We re-use the geometry and just truncate the vertex count
        ga.setCoordinates(0, coords);
        if (colors != null) {
          ga.setColors(0, colors);
        }

        if (size != oldSize) {
          if (isIndexGeometryArray()) {
            if (isStripGeometryArray()) {
              int[] indices = new int[indexCount * oldSize];
              ((IndexedGeometryStripArray) ga).getStripIndexCounts(indices);
              indices = Arrays.copyOf(indices, indexCount * size);
              ((IndexedGeometryStripArray) ga).setStripIndexCounts(indices);
            } else {
              ((IndexedGeometryArray) ga).setValidIndexCount(size * indexCount);
            }
          } else if (isStripGeometryArray()) {
            int[] indices = new int[vertexCount * oldSize];
            ((GeometryStripArray) ga).getStripVertexCounts(indices);
            indices = Arrays.copyOf(indices, vertexCount * size);
            ((GeometryStripArray) ga).setStripVertexCounts(indices);
          } else {
            ga.setValidVertexCount(size * vertexCount);
          }
        }
      }
    });
  }

  /**
   * Reorder the points using the indices.
   *
   * @param points the points
   * @param indices the indices
   * @return the new points
   */
  static Point3f[] reorderPoints(Point3f[] points, int[] indices) {
    final Point3f[] c = new Point3f[indices.length];
    for (int i = indices.length; i-- > 0;) {
      c[i] = points[indices[i]];
    }
    return c;
  }

  @Override
  public int size() {
    return points.length;
  }

  @Override
  public Point3f getCoordinate(int index) {
    return points[index];
  }

  @Override
  public void setColor(Color3f color) {
    // Delegate this to the interface implementation.
    // Allows transparent version to only implement to the interface method.
    setItemColor(color);
  }

  @Override
  public void setItemColor(Color3f color) {
    if (color == null) {
      color = DEFAULT_COLOR;
    }
    this.color = color;
    if (!hasColor()) {
      if (isColorByMaterial) {
        getAppearance().getMaterial().setDiffuseColor(color);
      } else {
        getAppearance().getColoringAttributes().setColor(color);
      }
      return;
    }
    final GeometryArray ga = (GeometryArray) getGeometry();
    if (ga == null) {
      return;
    }
    final int size = size();
    final int n = colorUpdater.size();
    final float[] colors = new float[size * n];
    if (hasColor3()) {
      final float[] tmp = new float[3];
      color.get(tmp);
      duplicate(tmp, 0, 3, colors.length / 3, colors, 0);
      ga.setColors(0, colors);
    } else {
      // Preserve alpha
      ga.getColors(0, colors);
      for (int i = 0; i < colors.length; i += 4) {
        colors[i] = color.x;
        colors[i + 1] = color.y;
        colors[i + 2] = color.z;
      }
      ga.setColors(0, colors);
    }
    changed = true;
  }

  @Override
  public void setItemColor(Color3f[] color) {
    if (!hasColor()) {
      setItemColor(color[0]);
      return;
    }
    this.color = null;
    final int size = size();
    ItemHelper.checkSize(color.length, size);
    final GeometryArray ga = (GeometryArray) getGeometry();
    if (ga == null) {
      return;
    }
    final int n = colorUpdater.size();
    final float[] colors = new float[size * n];
    if (hasColor3()) {
      for (int i = 0; i < color.length; i++) {
        System.arraycopy(colorUpdater.getColors(color[i]), 0, colors, i * n, n);
      }
      ga.setColors(0, colors);
    } else {
      // Preserve alpha
      ga.getColors(0, colors);
      for (int i = 0; i < color.length; i++) {
        final int offset = i * n;
        colorUpdater.getColors(color[i], colors[offset + 3]);
        System.arraycopy(colorUpdater.pointColor, 0, colors, offset, n);
      }
      ga.setColors(0, colors);
    }
    changed = true;
  }

  @Override
  public void calculateMinMaxCenterPoint(Point3f min, Point3f max, Point3f center) {
    CustomContentHelper.calculateMinMaxCenterPoint(min, max, center, points);
  }

  @Override
  public float getVolume() {
    return 0;
  }

  @Override
  public void setItemColor4(Color4f[] color) {
    checkPerItemAlpha();

    this.color = null;
    final int size = size();
    ItemHelper.checkSize(color.length, size);
    final GeometryArray ga = (GeometryArray) getGeometry();
    if (ga == null) {
      return;
    }
    final int n = colorUpdater.size();
    final float[] colors = new float[size * n];
    for (int i = 0; i < color.length; i++) {
      System.arraycopy(colorUpdater.getColors(color[i]), 0, colors, i * n, n);
    }
    ga.setColors(0, colors);
    changed = true;
  }

  @Override
  public void setItemAlpha(float[] alpha) {
    checkPerItemAlpha();

    final int size = size();
    ItemHelper.checkSize(alpha.length, size);
    final GeometryArray ga = (GeometryArray) getGeometry();
    if (ga == null) {
      return;
    }
    final int n = colorUpdater.size();
    final float[] colors = new float[size * n];
    // Preserve color
    ga.getColors(0, colors);
    for (int i = 0; i < size; i++) {
      final int offset = i * n;
      for (int j = 3; j < n; j += 4) {
        colors[j + offset] = alpha[i];
      }
    }
    ga.setColors(0, colors);
    changed = true;
  }

  @Override
  public void setItemAlpha(float alpha) {
    checkPerItemAlpha();

    final GeometryArray ga = (GeometryArray) getGeometry();
    if (ga == null) {
      return;
    }
    final int size = size();
    final int n = colorUpdater.size();
    final float[] colors = new float[size * n];
    // Preserve color
    ga.getColors(0, colors);
    for (int i = 0; i < size; i++) {
      final int offset = i * n;
      for (int j = 3; j < n; j += 4) {
        colors[j + offset] = alpha;
      }
    }
    ga.setColors(0, colors);
    changed = true;
  }

  @Override
  public void getItemAlpha(float[] alpha) {
    checkPerItemAlpha();

    final int size = size();
    ItemHelper.checkSize(alpha.length, size);
    final GeometryArray ga = (GeometryArray) getGeometry();
    if (ga == null) {
      return;
    }
    final int n = colorUpdater.size();
    final float[] colors = new float[size * n];
    ga.getColors(0, colors);
    for (int i = 0; i < size; i++) {
      // Get only alpha
      alpha[i] = colors[i * n + 3];
    }
  }

  @Override
  public void setShaded(boolean shaded) {
    if (!isPointArray) {
      super.setShaded(shaded);
    }
  }

  private void checkPerItemAlpha() {
    if (!hasColor4()) {
      throw new IllegalArgumentException("Per-item alpha not supported");
    }
  }
}
