package gdsc.smlm.ij.ij3d;

import org.scijava.java3d.Appearance;
import org.scijava.java3d.Geometry;
import org.scijava.java3d.GeometryArray;
import org.scijava.java3d.GeometryStripArray;
import org.scijava.java3d.GeometryUpdater;
import org.scijava.java3d.IndexedGeometryArray;
import org.scijava.java3d.IndexedGeometryStripArray;
import org.scijava.java3d.TransparencyAttributes;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Point3f;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import customnode.CustomMesh;
import gdsc.core.utils.BitFlags;
import gdsc.core.utils.NotImplementedException;

/**
 * Use a mesh object to represent a set of points. The object is duplicated, scaled and translated for
 * each point.
 */
public abstract class ItemMesh extends CustomMesh implements UpdateableItemShape
{
	/** The vertex count of the original geometry array. */
	final protected int vertexCount;

	/** The vertex format of the original geometry array. */
	final protected int vertexFormat;

	/** The points. */
	protected Point3f[] points;

	/** The size of each point. */
	protected Point3f[] sizes;

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
	 * @param colors
	 *            the per-item colors. Can be null. Must be the correct length for the number of vertices (including
	 *            alpha if using GeometryArray.COLOR_4)
	 */
	public ItemMesh(Point3f[] points, GeometryArray ga, Appearance appearance, Point3f[] sizes, float[] colors)
	{
		// Create empty
		super(null, null, 0f);

		if (sizes != null && points.length != sizes.length && sizes.length != 1)
			throw new IllegalArgumentException("Points and sizes must be the same length");

		if (ga.getVertexAttrCount() > 0)
			throw new IllegalArgumentException("Vertex attributes are not supported");
		if (ga.getTexCoordSetCount() > 0)
			throw new IllegalArgumentException("Texture coordinates are not supported");

		vertexCount = ga.getValidVertexCount();
		vertexFormat = ga.getVertexFormat();

		// Set the flags we support and error if there are others
		int flags = GeometryArray.COORDINATES;

		if (colors != null)
		{
			if (!hasColor())
				throw new IllegalArgumentException("Colours given but no colour format used in input geometry");
			flags |= GeometryArray.COLOR_3;
			if (hasColor4())
				flags |= GeometryArray.COLOR_4;
			int size = points.length * getNumberOfColorsPerVertex();
			if (colors.length != size)
				throw new IllegalArgumentException(
						"Incorrect number of colours used, input=" + colors.length + ", expected=" + size);
		}

		if ((vertexFormat & GeometryArray.NORMALS) != 0)
			flags |= GeometryArray.NORMALS;

		// Only support simple indexed arrays
		if (ga instanceof IndexedGeometryArray)
			flags |= GeometryArray.USE_COORD_INDEX_ONLY;

		int extra = BitFlags.unset(vertexFormat, flags);
		if (extra != 0)
			throw new IllegalArgumentException("Unsupported vertex format flags: " + extra);
		if (!BitFlags.areSet(vertexFormat, flags))
			throw new IllegalArgumentException("Unsupported vertex format: " + extra);

		this.points = points;
		this.sizes = sizes;

		// Now build the actual vertices by repeating the points.
		float[] objectCoords = new float[vertexCount * 3];
		ga.getCoordinates(0, objectCoords);

		float[] allCoords = new float[vertexCount * points.length];

		boolean sameSize = false;
		if (sizes == null || (sameSize = sameSize(sizes)))
		{
			if (sameSize)
			{
				// Scale the input object
				final Point3f s = sizes[0];
				final float sx = s.x;
				final float sy = s.y;
				final float sz = s.z;
				for (int j = 0; j < objectCoords.length; j += 3)
				{
					objectCoords[j] *= sx;
					objectCoords[j + 1] *= sy;
					objectCoords[j + 2] *= sz;
				}
			}

			// Translate
			for (int i = 0, k = 0; i < points.length; i++)
			{
				final Point3f p = points[i];
				final float dx = p.x;
				final float dy = p.y;
				final float dz = p.z;
				for (int j = 0; j < objectCoords.length; j += 3)
				{
					allCoords[k++] = objectCoords[j] + dx;
					allCoords[k++] = objectCoords[j + 1] + dy;
					allCoords[k++] = objectCoords[j + 2] + dz;
				}
			}
		}
		else
		{
			// Translate and scale
			for (int i = 0, k = 0; i < points.length; i++)
			{
				final Point3f p = points[i];
				final float dx = p.x;
				final float dy = p.y;
				final float dz = p.z;
				final Point3f s = sizes[i];
				final float sx = s.x;
				final float sy = s.y;
				final float sz = s.z;
				for (int j = 0; j < objectCoords.length; j += 3)
				{
					allCoords[k++] = objectCoords[j] * sx + dx;
					allCoords[k++] = objectCoords[j + 1] * sy + dy;
					allCoords[k++] = objectCoords[j + 2] * sz + dz;
				}
			}
		}

		// Update the geometry. 
		// Do this in a method to allow sub-classes to change the geometry.
		this.setGeometry(createGeometry(allCoords, colors, ga));

		// Create a default appearance
		setAppearance(appearance);
	}

	@Override
	public void update()
	{
		// Ignore this
	}

	/**
	 * Checks for normals.
	 *
	 * @return true, if successful
	 */
	public boolean hasNormals()
	{
		return ((vertexFormat & GeometryArray.NORMALS) != 0);
	}

	/**
	 * Checks for color.
	 *
	 * @return true, if successful
	 */
	public boolean hasColor()
	{
		return ((vertexFormat & GeometryArray.COLOR_3) == GeometryArray.COLOR_3);
	}

	/**
	 * Checks for color 3.
	 *
	 * @return true, if successful
	 */
	public boolean hasColor3()
	{
		return hasColor() && !hasColor4();
	}

	/**
	 * Checks for color 4.
	 *
	 * @return true, if successful
	 */
	public boolean hasColor4()
	{
		return ((vertexFormat & GeometryArray.COLOR_4) == GeometryArray.COLOR_4);
	}

	/**
	 * Gets the number of colours per vertex (0, 3, or 4).
	 *
	 * @return the number of colors per vertex
	 */
	public int getNumberOfColorsPerVertex()
	{
		return (hasColor4()) ? 4 : (hasColor()) ? 3 : 0;
	}

	/**
	 * Check if all the points are the same size.
	 *
	 * @param sizes
	 *            the sizes (must be an array of at least 1)
	 * @return true, if successful
	 */
	public static boolean sameSize(Point3f[] sizes)
	{
		if (sizes.length == 1)
			return true;
		final Point3f s = sizes[0];
		for (int j = 1; j < sizes.length; j++)
		{
			if (!sizes[j].equals(s))
				return false;
		}
		return true;
	}

	@Override
	public void setCoordinate(final int i, final Point3f p)
	{
		throw new NotImplementedException();
	}

	@Override
	public void setCoordinates(final int[] indices, final Point3f p)
	{
		throw new NotImplementedException();
	}

	@Override
	protected void addVertices(Point3f[] v)
	{
		throw new NotImplementedException();
	}

	@Override
	protected void addVerticesToGeometryArray(Point3f[] v)
	{
		throw new NotImplementedException();
	}

	@Override
	protected void addVerticesToGeometryStripArray(Point3f[] v)
	{
		throw new NotImplementedException();
	}

	@Override
	protected void removeVertices(int[] indices)
	{
		throw new NotImplementedException();
	}

	@Override
	protected GeometryArray createGeometry()
	{
		throw new NotImplementedException();
	}

	protected GeometryArray createGeometry(float[] coords, float[] colors, GeometryArray sourceGA)
	{
		// Create using reflection
		final GeometryArray ga;
		try
		{
			Class<?> clazz = sourceGA.getClass();
			//clazz = clazz.asSubclass(clazz);
			Class<?>[] paramTypes = { int.class, int.class };
			Object[] paramValues = { vertexCount * points.length, vertexFormat };

			ga = (GeometryArray) clazz.asSubclass(clazz).getConstructor(paramTypes).newInstance(paramValues);
		}
		catch (Exception e)
		{
			e.printStackTrace();
			return null;
		}

		ga.setCoordinates(0, coords);
		ga.setValidVertexCount(coords.length);
		ga.setCapability(GeometryArray.ALLOW_COORDINATE_WRITE);
		ga.setCapability(GeometryArray.ALLOW_COUNT_WRITE);
		ga.setCapability(GeometryArray.ALLOW_COUNT_READ);
		ga.setCapability(GeometryArray.ALLOW_FORMAT_READ);
		ga.setCapability(Geometry.ALLOW_INTERSECT);

		// Handle normals
		if (hasNormals())
		{
			float[] objectNormals = new float[vertexCount * 3];
			sourceGA.getNormals(0, objectNormals);
			float[] allNormals = new float[objectNormals.length * points.length];
			duplicate(objectNormals, 0, objectNormals.length, points.length, allNormals, 0);
			ga.setNormals(0, allNormals);
		}

		// Handle colors
		if (hasColor())
		{
			ga.setCapability(GeometryArray.ALLOW_COLOR_READ);
			ga.setCapability(GeometryArray.ALLOW_COLOR_WRITE);
			if (colors != null)
			{
				// Special case for a PointArray
				if (vertexCount == 1)
				{
					ga.setColors(0, colors);
				}
				else
				{
					int colorsPerVertex = getNumberOfColorsPerVertex();
					int coloursPerObject = colorsPerVertex * vertexCount;
					float[] all = new float[coloursPerObject * colors.length];
					for (int i = 0; i < colors.length; i += colorsPerVertex)
					{
						duplicate(colors, i, colorsPerVertex, vertexCount, all, i * vertexCount);
					}
					ga.setColors(0, all);
				}
			}
		}

		// Fan are extensions of GeometryStripArray so do not need extra code.

		// Handle indexed array
		if (sourceGA instanceof IndexedGeometryArray)
		{
			IndexedGeometryArray iga = (IndexedGeometryArray) sourceGA;
			int indexCount = iga.getValidIndexCount();
			int[] indices = new int[indexCount];
			iga.getCoordinateIndices(0, indices);
			int[] allIndices = new int[indices.length * points.length];
			for (int i = 0, k = 0; i < points.length; i++)
			{
				int offset = k;
				for (int j = 0; j < indices.length; j++)
					allIndices[k++] = indices[j] + offset;
			}
			((IndexedGeometryArray) ga).setCoordinateIndices(0, allIndices);

			// Handle strips
			if (sourceGA instanceof IndexedGeometryStripArray)
			{
				IndexedGeometryStripArray igsa = (IndexedGeometryStripArray) sourceGA;
				int numStrips = igsa.getNumStrips();
				int[] indexCounts = new int[numStrips];
				igsa.getStripIndexCounts(indexCounts);
				int[] allIndexCounts = new int[numStrips * points.length];
				duplicate(indexCounts, 0, numStrips, points.length, allIndexCounts, 0);
				((IndexedGeometryStripArray) ga).getStripIndexCounts(allIndexCounts);
			}
		}
		// Handle strips
		else if (sourceGA instanceof GeometryStripArray)
		{
			GeometryStripArray gsa = (GeometryStripArray) sourceGA;
			int numStrips = gsa.getNumStrips();
			int[] vertexCounts = new int[numStrips];
			gsa.getStripVertexCounts(vertexCounts);
			int[] allVertexCounts = new int[numStrips * points.length];
			duplicate(vertexCounts, 0, numStrips, points.length, allVertexCounts, 0);
			((GeometryStripArray) ga).getStripVertexCounts(allVertexCounts);
		}

		return ga;
	}

	private void duplicate(Object source, int from, int length, int n, Object dest, int to)
	{
		// Binary fill
		int fill = length;
		System.arraycopy(source, from, dest, to, fill);
		for (int i = 2; i < n; i *= 2)
		{
			System.arraycopy(dest, to, dest, to + fill, fill);
			fill *= 2;
		}
		// Final fill
		System.arraycopy(dest, to, dest, to + fill, n * length - fill);
	}

	private static int transparencyMode = TransparencyAttributes.FASTEST;

	/**
	 * Sets the transparency mode.
	 *
	 * @param mode
	 *            the new transparency mode
	 * @throws IllegalArgumentException
	 *             If the mode is not valid
	 * @see org.scijava.java3d.TransparencyAttributes.setTransparencyMode(int)
	 */
	public static void setTransparencyMode(int mode) throws IllegalArgumentException
	{
		if ((mode < TransparencyAttributes.FASTEST) || (mode > TransparencyAttributes.NONE))
		{
			throw new IllegalArgumentException("Not a valid transparency mode");
		}
		transparencyMode = mode;
	}

	/**
	 * Gets the transparency mode.
	 *
	 * @return the transparency mode
	 * @see org.scijava.java3d.TransparencyAttributes.setTransparencyMode(int)
	 */
	public static int getTransparencyMode()
	{
		return transparencyMode;
	}

	@Override
	public void setTransparency(final float transparency)
	{
		// We want to use a different transparency from the ij3d default which is FASTEST
		// so override this method.
		Appearance appearance = getAppearance();
		final TransparencyAttributes ta = appearance.getTransparencyAttributes();
		if (transparency <= .01f)
		{
			this.transparency = 0.0f;
			ta.setTransparencyMode(TransparencyAttributes.NONE);
		}
		else
		{
			this.transparency = transparency;
			ta.setTransparencyMode(transparencyMode);
		}
		ta.setTransparency(this.transparency);
	}

	@Override
	protected Appearance createAppearance()
	{
		throw new NotImplementedException();
	}

	protected Appearance createAppearance(GeometryArray ga)
	{
		// TOD0 - create a suitable appearance for points or 3D shapes.

		Appearance appearance = super.createAppearance();
		// Update the transparency to the default mode
		final TransparencyAttributes ta = appearance.getTransparencyAttributes();
		if (ta.getTransparencyMode() != TransparencyAttributes.NONE)
			ta.setTransparencyMode(transparencyMode);
		return appearance;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.UpdatedableItemMesh#reorder(int[])
	 */
	public void reorder(int[] indices) throws IllegalArgumentException
	{
		ItemPointMesh.checkIndices(indices, points.length);
		reorderFast(indices);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.UpdatedableItemMesh#reorderFast(int[])
	 */
	public void reorderFast(int[] indices) throws IllegalArgumentException
	{
		changed = true;

		int oldSize = size();
		int size = (indices == null) ? 0 : Math.min(oldSize, indices.length);

		if (size == 0)
		{
			points = new Point3f[0];
			sizes = new Point3f[0];
			this.setGeometry(null);
			return;
		}

		// From here on we assume the current geometry will not be null
		// as this only happens when the original size is zero. Size has 
		// been checked at this point to be the smaller of new and old. 
		GeometryArray ga = (GeometryArray) getGeometry();

		points = reorder(points, indices);
		// Sizes could be null or a single size
		if (sizes != null && sizes.length == points.length)
			sizes = reorder(sizes, indices);

		// Reorder all things in the geometry: coordinates and colour.
		// The normals, indices, strip counts are are unchanged.
		int objectSize = vertexCount;

		final float[] oldCoords = new float[oldSize * objectSize * 3];
		ga.getCoordinates(0, oldCoords);
		final float[] coords = new float[size * objectSize * 3];
		for (int i = 0; i < size; i++)
		{
			int j = indices[i];
			int ii = i * objectSize * 3;
			int jj = j * objectSize * 3;
			System.arraycopy(oldCoords, jj, coords, ii, objectSize);
		}

		final float[] colors;
		if (hasColor())
		{
			int n = getNumberOfColorsPerVertex();
			float[] oldColors = (n == 3) ? oldCoords : new float[oldSize * objectSize * n];
			ga.getColors(0, oldColors);
			colors = new float[size * objectSize * n];
			for (int i = 0; i < size; i++)
			{
				int j = indices[i];
				int ii = i * objectSize * n;
				int jj = j * objectSize * n;
				System.arraycopy(oldColors, jj, colors, ii, objectSize);
			}
		}
		else
		{
			colors = null;
		}

		ga.updateData(new GeometryUpdater()
		{
			public void updateData(Geometry geometry)
			{
				GeometryArray ga = (GeometryArray) geometry;
				// We re-use the geometry and just truncate the vertex count
				ga.setCoordinates(0, coords);
				if (colors != null)
					ga.setColors(0, colors);
				ga.setValidVertexCount(coords.length);
			}
		});
	}

	static Point3f[] reorder(Point3f[] p, int[] indices)
	{
		Point3f[] c = new Point3f[indices.length];
		for (int i = indices.length; i-- > 0;)
			c[i] = p[indices[i]];
		return c;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.ItemMesh#size()
	 */
	public int size()
	{
		return points.length;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.ItemShape#getCoordinate(int)
	 */
	public Point3f getCoordinate(int i)
	{
		return points[i];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see customnode.CustomMesh#setColor(org.scijava.vecmath.Color3f)
	 */
	@Override
	public void setColor(Color3f color)
	{
		// Delegate this to the interface implementation.
		// Allows transparent version to only implement to the interface method.
		setItemColor(color);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.ItemShape#setItemColor(org.scijava.vecmath.Color3f)
	 */
	public void setItemColor(Color3f color)
	{
		if (!hasColor())
		{
			// TODO - set colour in the Appearance, e.g. Material color or ColoringAttributes
			return;
		}
		if (color == null)
			color = DEFAULT_COLOR;
		this.color = color;
		int size = size();
		final GeometryArray ga = (GeometryArray) getGeometry();
		if (ga == null)
			return;
		int objectSize = vertexCount;
		int n = getNumberOfColorsPerVertex();
		final float[] colors = new float[objectSize * size * n];
		if (n == 3)
		{
			float[] tmp = new float[3];
			color.get(tmp);
			duplicate(tmp, 0, 3, size, colors, 0);
			ga.setColors(0, colors);
		}
		else
		{
			// Preserve alpha
			ga.getColors(0, colors);
			int i = 0;
			for (int j = objectSize; j-- > 0;)
			{
				colors[i++] = color.x;
				colors[i++] = color.y;
				colors[i++] = color.z;
				i++; // Skip over alpha
			}
			ga.setColors(0, colors);
		}
		changed = true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.ItemMesh#setItemColor(org.scijava.vecmath.Color3f[])
	 */
	public void setItemColor(Color3f[] color) throws IllegalArgumentException
	{
		if (!hasColor())
		{
			// TODO - set colour in the Appearance, e.g. Material color or ColoringAttributes
			return;
		}
		this.color = null;
		int size = size();
		if (color.length != size)
			throw new IllegalArgumentException("list of size " + size + " expected");
		final GeometryArray ga = (GeometryArray) getGeometry();
		if (ga == null)
			return;
		int objectSize = vertexCount;
		int n = getNumberOfColorsPerVertex();
		final float[] colors = new float[objectSize * size * n];
		if (n == 3)
		{
			int i = 0;
			for (Color3f c : color)
			{
				for (int j = objectSize; j-- > 0;)
				{
					colors[i++] = c.x;
					colors[i++] = c.y;
					colors[i++] = c.z;
				}
			}
			ga.setColors(0, colors);
		}
		else
		{
			// Preserve alpha
			ga.getColors(0, colors);
			int i = 0;
			for (Color3f c : color)
			{
				for (int j = objectSize; j-- > 0;)
				{
					colors[i++] = c.x;
					colors[i++] = c.y;
					colors[i++] = c.z;
					i++; // Skip over alpha
				}
			}
			ga.setColors(0, colors);
		}
		changed = true;
	}

	@Override
	public void calculateMinMaxCenterPoint(Point3f min, Point3f max, Point3f center)
	{
		CustomContentHelper.calculateMinMaxCenterPoint(min, max, center, points);
	}

	@Override
	public float getVolume()
	{
		return 0;
	}
}
