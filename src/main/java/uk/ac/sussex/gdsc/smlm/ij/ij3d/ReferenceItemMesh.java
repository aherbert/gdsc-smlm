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

import java.util.Arrays;

import org.scijava.java3d.Appearance;
import org.scijava.java3d.Geometry;
import org.scijava.java3d.GeometryArray;
import org.scijava.java3d.GeometryStripArray;
import org.scijava.java3d.GeometryUpdater;
import org.scijava.java3d.IndexedGeometryArray;
import org.scijava.java3d.IndexedGeometryStripArray;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Color4f;
import org.scijava.vecmath.Point3f;

/**
 * Use a mesh object to represent a set of points. The object is duplicated, scaled and translated for
 * each point.
 * <p>
 * Note: TransparentItemShape is only supported if {@link #hasColor4()} is true, i.e. the input geometry has per vertex
 * colours with alpha.
 * <p>
 * The created geometry array is created by reference.
 */
public class ReferenceItemMesh extends ItemMesh
{
    /**
     * Instantiates a new item mesh.
     * <p>
     * This will repeat the object for each input point. The object
     * is assumed to be centred on the origin. It will be scaled and
     * translated for each input point.
     * <p>
     * The input geometry array vertex format is checked and unsupported formats throw an exception. Currently this only
     * supports coordinates, normals and color. Strip and Fan arrays are supported.
     *
     * @param point
     *            the point
     * @param ga
     *            the geometry array. If null then a default will be used. Assumed to be centred on the origin.
     * @param appearance
     *            the default appearance of the shape. PolygonAttributes, Material and TransparencyAttributes are used.
     * @param sizes
     *            the sizes of each point. Can be null (no scaling); length=1 (fixed scaling); or points.length.
     * @param color
     *            the color
     * @param transparency
     *            the transparency
     */
    public ReferenceItemMesh(Point3f point, GeometryArray ga, Appearance appearance, Point3f[] sizes,
            final Color3f color, final float transparency)
    {
        super(point, ga, appearance, sizes, color, transparency);
    }

    /**
     * Instantiates a new item mesh.
     * <p>
     * This will repeat the object for each input point. The object
     * is assumed to be centred on the origin. It will be scaled and
     * translated for each input point.
     * <p>
     * The input geometry array vertex format is checked and unsupported formats throw an exception. Currently this only
     * supports coordinates, normals and color. Strip and Fan arrays are supported.
     *
     * @param points
     *            the points
     * @param ga
     *            the geometry array. If null then a default will be used. Assumed to be centred on the origin.
     * @param appearance
     *            the default appearance of the shape. PolygonAttributes, Material and TransparencyAttributes are used.
     * @param sizes
     *            the sizes of each point. Can be null (no scaling); length=1 (fixed scaling); or points.length.
     * @param color
     *            the color
     * @param transparency
     *            the transparency
     */
    public ReferenceItemMesh(Point3f[] points, GeometryArray ga, Appearance appearance, Point3f[] sizes,
            final Color3f color, final float transparency)
    {
        super(points, ga, appearance, sizes, color, transparency);
    }

    @Override
    protected GeometryArray createGeometry(float[] coords, GeometryArray sourceGA)
    {
        final GeometryArray ga = createGeometryArray(sourceGA, GeometryArray.BY_REFERENCE);

        ga.setCoordRefFloat(coords);
        ga.setCapability(GeometryArray.ALLOW_REF_DATA_READ);
        ga.setCapability(GeometryArray.ALLOW_REF_DATA_WRITE);

        // Fan are extensions of GeometryStripArray so do not need extra code.

        // Handle indexed array
        if (isIndexGeometryArray())
        {
            final IndexedGeometryArray sourceIGA = (IndexedGeometryArray) sourceGA;
            final IndexedGeometryArray iga = (IndexedGeometryArray) ga;
            final int objectIndexCount = sourceIGA.getValidIndexCount();
            final int[] objectIndices = new int[objectIndexCount];
            final int[] allIndices = new int[objectIndices.length * points.length];
            sourceIGA.getCoordinateIndices(0, objectIndices);
            duplicateIndices(objectIndices, allIndices);
            iga.setCoordinateIndices(0, allIndices);

            // Check if we need the color and normal indices
            if ((vertexFormat & GeometryArray.USE_COORD_INDEX_ONLY) == 0)
            {
                // Duplicate the other indices
                if (hasNormals())
                {
                    sourceIGA.getNormalIndices(0, objectIndices);
                    duplicateIndices(objectIndices, allIndices);
                    iga.setNormalIndices(0, allIndices);
                }

                if (hasColor())
                {
                    sourceIGA.getColorIndices(0, objectIndices);
                    duplicateIndices(objectIndices, allIndices);
                    iga.setColorIndices(0, allIndices);
                }
            }
        }

        // Handle normals
        if (hasNormals())
        {
            final float[] objectNormals = new float[vertexCount * 3];
            sourceGA.getNormals(0, objectNormals);
            final float[] allNormals = new float[objectNormals.length * points.length];
            duplicate(objectNormals, 0, objectNormals.length, points.length, allNormals, 0);
            ga.setNormalRefFloat(allNormals);
        }

        // Handle colors
        if (hasColor())
        {
            colorUpdater = ArrayColorUpdater.create(vertexCount, hasColor4());
            ga.setColorRefFloat(new float[colorUpdater.size() * size()]);
        }

        return ga;
    }

    /** {@inheritDoc} */
    @Override
    public void reorderFast(int[] indices) throws IllegalArgumentException
    {
        changed = true;

        final int oldSize = size();
        final int size = (indices == null) ? 0 : Math.min(oldSize, indices.length);

        if (size == 0 || indices == null)
        {
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
        if (sizes != null && sizes.length == points.length)
            sizes = reorder(sizes, indices);

        // Reorder all things in the geometry: coordinates and colour.
        // The normals, indices, strip counts are are unchanged.
        //int objectSize = vertexCount;

        int n = vertexCount * 3;
        final float[] oldCoords = ga.getCoordRefFloat();
        final float[] coords = new float[size * n];
        for (int i = 0; i < size; i++)
        {
            final int j = indices[i];
            final int ii = i * n;
            final int jj = j * n;
            System.arraycopy(oldCoords, jj, coords, ii, n);
        }

        final float[] colors;
        if (hasColor())
        {
            n = colorUpdater.size();
            final float[] oldColors = ga.getColorRefFloat();
            colors = new float[size * n];
            for (int i = 0; i < size; i++)
            {
                final int j = indices[i];
                final int ii = i * n;
                final int jj = j * n;
                System.arraycopy(oldColors, jj, colors, ii, n);
            }
        }
        else
            colors = null;

        ga.updateData(new GeometryUpdater()
        {
            @Override
            public void updateData(Geometry geometry)
            {
                final GeometryArray ga = (GeometryArray) geometry;
                // We re-use the geometry and just truncate the vertex count
                ga.setCoordRefFloat(coords);
                if (colors != null)
                    ga.setColorRefFloat(colors);

                if (size != oldSize)
                    if (isIndexGeometryArray())
                    {
                        if (isStripGeometryArray())
                        {
                            int[] indices = new int[indexCount * oldSize];
                            ((IndexedGeometryStripArray) ga).getStripIndexCounts(indices);
                            indices = Arrays.copyOf(indices, indexCount * size);
                            ((IndexedGeometryStripArray) ga).setStripIndexCounts(indices);
                        }
                        else
                            ((IndexedGeometryArray) ga).setValidIndexCount(size * indexCount);
                    }
                    else if (isStripGeometryArray())
                    {
                        int[] indices = new int[vertexCount * oldSize];
                        ((GeometryStripArray) ga).getStripVertexCounts(indices);
                        indices = Arrays.copyOf(indices, vertexCount * size);
                        ((GeometryStripArray) ga).setStripVertexCounts(indices);
                    }
                    else
                        ga.setValidVertexCount(size * vertexCount);
            }
        });
    }

    /** {@inheritDoc} */
    @Override
    public void setItemColor(Color3f color)
    {
        if (color == null)
            color = DEFAULT_COLOR;
        this.color = color;
        if (!hasColor())
        {
            if (isColorByMaterial)
                getAppearance().getMaterial().setDiffuseColor(color);
            else
                getAppearance().getColoringAttributes().setColor(color);
            return;
        }
        final GeometryArray ga = (GeometryArray) getGeometry();
        if (ga == null)
            return;
        final float[] colors;
        if (hasColor3())
        {
            colors = new float[size() * colorUpdater.size()];
            final float[] tmp = new float[3];
            color.get(tmp);
            duplicate(tmp, 0, 3, colors.length / 3, colors, 0);
        }
        else
        {
            // Preserve alpha
            colors = ga.getColorRefFloat().clone();
            for (int i = 0; i < colors.length; i += 4)
            {
                colors[i] = color.x;
                colors[i + 1] = color.y;
                colors[i + 2] = color.z;
            }
        }
        ga.setColorRefFloat(colors);
        changed = true;
    }

    /** {@inheritDoc} */
    @Override
    public void setItemColor(Color3f[] color) throws IllegalArgumentException
    {
        if (!hasColor())
        {
            setItemColor(color[0]);
            return;
        }
        this.color = null;
        final int size = size();
        if (color.length != size)
            throw new IllegalArgumentException("list of size " + size + " expected");
        final GeometryArray ga = (GeometryArray) getGeometry();
        if (ga == null)
            return;
        final int n = colorUpdater.size();
        final float[] colors;
        if (hasColor3())
        {
            colors = new float[size() * n];
            for (int i = 0; i < color.length; i++)
                System.arraycopy(colorUpdater.getColors(color[i]), 0, colors, i * n, n);
        }
        else
        {
            // Preserve alpha
            colors = ga.getColorRefFloat().clone();
            for (int i = 0; i < color.length; i++)
            {
                final int offset = i * n;
                colorUpdater.getColors(color[i], colors[offset + 3]);
                System.arraycopy(colorUpdater.pointColor, 0, colors, offset, n);
            }
        }
        ga.setColorRefFloat(colors);
        changed = true;
    }

    /** {@inheritDoc} */
    @Override
    public void setItemColor4(Color4f[] color) throws IllegalArgumentException
    {
        if (!hasColor4())
            throw new IllegalArgumentException("Per-item alpha not supported");

        this.color = null;
        final int size = size();
        if (color.length != size)
            throw new IllegalArgumentException("list of size " + size + " expected");
        final GeometryArray ga = (GeometryArray) getGeometry();
        if (ga == null)
            return;
        final int n = colorUpdater.size();
        final float[] colors = new float[size() * n];
        for (int i = 0; i < color.length; i++)
            System.arraycopy(colorUpdater.getColors(color[i]), 0, colors, i * n, n);
        ga.setColorRefFloat(colors);
        changed = true;
    }

    /** {@inheritDoc} */
    @Override
    public void setItemAlpha(float[] alpha) throws IllegalArgumentException
    {
        if (!hasColor4())
            throw new IllegalArgumentException("Per-item alpha not supported");

        final int size = size();
        if (alpha.length != size)
            throw new IllegalArgumentException("list of size " + size + " expected");
        final GeometryArray ga = (GeometryArray) getGeometry();
        if (ga == null)
            return;
        final int n = colorUpdater.size();
        // Preserve color
        final float[] colors = ga.getColorRefFloat().clone();
        for (int i = 0; i < size; i++)
        {
            final int offset = i * n;
            for (int j = 3; j < n; j += 4)
                colors[j + offset] = alpha[i];
        }
        ga.setColorRefFloat(colors);
        changed = true;
    }

    /** {@inheritDoc} */
    @Override
    public void setItemAlpha(float alpha) throws IllegalArgumentException
    {
        if (!hasColor4())
            throw new IllegalArgumentException("Per-item alpha not supported");

        final int size = size();
        final GeometryArray ga = (GeometryArray) getGeometry();
        if (ga == null)
            return;
        final int n = colorUpdater.size();
        // Preserve color
        final float[] colors = ga.getColorRefFloat().clone();
        for (int i = 0; i < size; i++)
        {
            final int offset = i * n;
            for (int j = 3; j < n; j += 4)
                colors[j + offset] = alpha;
        }
        ga.setColorRefFloat(colors);
        changed = true;
    }

    /** {@inheritDoc} */
    @Override
    public void getItemAlpha(float[] alpha) throws IllegalArgumentException
    {
        if (!hasColor4())
            throw new IllegalArgumentException("Per-item alpha not supported");

        final int size = size();
        if (alpha.length != size)
            throw new IllegalArgumentException("list of size " + size + " expected");
        final GeometryArray ga = (GeometryArray) getGeometry();
        if (ga == null)
            return;
        final int n = colorUpdater.size();
        final float[] colors = ga.getColorRefFloat();
        for (int i = 0; i < size; i++)
            // Get only alpha
            alpha[i] = colors[i * n + 3];
    }
}
