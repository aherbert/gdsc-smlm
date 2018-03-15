package gdsc.smlm.ij.ij3d;

import java.util.Arrays;

import org.scijava.java3d.Appearance;
import org.scijava.java3d.Geometry;
import org.scijava.java3d.GeometryArray;
import org.scijava.java3d.GeometryUpdater;
import org.scijava.java3d.TransparencyAttributes;
import org.scijava.java3d.TriangleArray;
import org.scijava.java3d.utils.geometry.GeometryInfo;
import org.scijava.java3d.utils.geometry.NormalGenerator;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Point3f;
import org.scijava.vecmath.Vector3f;

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

import customnode.CustomTriangleMesh;
import gdsc.core.logging.NullTrackProgress;
import gdsc.core.logging.Ticker;
import gdsc.core.logging.TrackProgress;

/**
 * Use a triangle mesh object to represent a set of points. The object is duplicated, scaled and translated for
 * each point.
 */
public class ItemTriangleMesh extends CustomTriangleMesh implements UpdateableItemMesh
{
	protected Point3f[] objectVertices;
	protected Vector3f[] objectNormals;
	protected Point3f[] points;
	protected Point3f[] sizes;
	protected boolean dirty = false;

	/**
	 * Instantiates a new item triangle mesh.
	 * <p>
	 * This will repeat the object for each input point. The object
	 * is assumed to be centred on the origin. It will be scaled and
	 * translated for each input point.
	 *
	 * @param objectVertices
	 *            the vertices of the object for a single point
	 * @param points
	 *            the points
	 * @param sizes
	 *            the size of each point
	 * @param color
	 *            the color
	 * @param transp
	 *            the transparency
	 */
	public ItemTriangleMesh(Point3f[] objectVertices, Point3f[] points, Point3f[] sizes, Color3f color, float transp)
	{
		this(objectVertices, points, sizes, color, transp, -1, NullTrackProgress.INSTANCE);
	}

	/**
	 * Instantiates a new item triangle mesh.
	 * <p>
	 * This will repeat the object for each input point. The object
	 * is assumed to be centred on the origin. It will be scaled and
	 * translated for each input point.
	 * <p>
	 * The crease angle is used to collapse facets normals at a vertex into a single normal for smoothing shading. Set
	 * to 0 to draw the polygon with no shading.
	 *
	 * @param objectVertices
	 *            the vertices of the object for a single point
	 * @param points
	 *            the points
	 * @param sizes
	 *            the size of each point
	 * @param color
	 *            the color
	 * @param transp
	 *            the transparency
	 * @param creaseAngle
	 *            the crease angle (in degrees). Set to negative to ignore. The default is 44.
	 * @param progress
	 *            the progress
	 */
	public ItemTriangleMesh(Point3f[] objectVertices, Point3f[] points, Point3f[] sizes, Color3f color, float transp,
			double creaseAngle, TrackProgress progress)
	{
		// Create empty 
		super(null, color, transp);

		if (sizes != null && points.length != sizes.length && sizes.length != 1)
			throw new IllegalArgumentException("Points and sizes must be the same length");

		progress = NullTrackProgress.createIfNull(progress);

		if (progress.isStatus())
			progress.status("Standardising object vertices");
		this.objectVertices = objectVertices;

		checkFacets(objectVertices);

		objectNormals = getNormals(objectVertices, creaseAngle);

		// Now build the actual vertices by repeating the points.
		this.points = points;
		this.sizes = sizes;

		Point3f[] vertices = new Point3f[objectVertices.length * points.length];

		final int n = objectVertices.length;
		if (progress.isStatus())
			progress.status("Computing vertices");
		Ticker ticker = Ticker.createStarted(progress, n, false);
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
				for (int j = 0; j < n; j++)
				{
					objectVertices[j].x *= sx;
					objectVertices[j].y *= sy;
					objectVertices[j].z *= sz;
				}
			}

			// Translate
			for (int i = 0, k = 0; i < points.length; i++)
			{
				final Point3f p = points[i];
				final float dx = p.x;
				final float dy = p.y;
				final float dz = p.z;
				for (int j = 0; j < n; j++)
				{
					vertices[k++] = new Point3f(objectVertices[j].x + dx, objectVertices[j].y + dy,
							objectVertices[j].z + dz);
				}
				ticker.tick();
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
				for (int j = 0; j < n; j++)
				{
					vertices[k++] = new Point3f(objectVertices[j].x * sx + dx, objectVertices[j].y * sy + dy,
							objectVertices[j].z * sz + dz);
				}
				ticker.tick();
			}
		}

		ticker.stop();

		this.mesh = Arrays.asList(vertices);

		// Update the geometry
		if (progress.isStatus())
			progress.status("Creating geometry");
		this.setGeometry(createGeometry());
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
		dirty = true;
		super.setCoordinate(i, p);
	}

	@Override
	public void setCoordinates(final int[] indices, final Point3f p)
	{
		dirty = true;
		super.setCoordinates(indices, p);
	}

	@Override
	public void addTriangle(Point3f p1, Point3f p2, Point3f p3)
	{
		dirty = true;
		super.addTriangle(p1, p2, p3);
	}

	@Override
	public void addTriangles(Point3f[] v)
	{
		dirty = true;
		super.addTriangles(v);
	}

	@Override
	public void removeTriangle(int index)
	{
		dirty = true;
		super.removeTriangle(index);
	}

	@Override
	public void removeTriangles(int[] indices)
	{
		dirty = true;
		super.removeTriangles(indices);
	}

	@Override
	protected void addVertices(Point3f[] v)
	{
		dirty = true;
		super.addVertices(v);
	}

	@Override
	protected void addVerticesToGeometryArray(Point3f[] v)
	{
		dirty = true;
		super.addVerticesToGeometryArray(v);
	}

	@Override
	protected void addVerticesToGeometryStripArray(Point3f[] v)
	{
		dirty = true;
		super.addVerticesToGeometryStripArray(v);
	}

	@Override
	protected void removeVertices(int[] indices)
	{
		dirty = true;
		super.removeVertices(indices);
	}

	@Override
	protected GeometryArray createGeometry()
	{
		if (mesh == null || mesh.size() < 3)
			return null;
		final int vertexCount = mesh.size();

		final Point3f[] coords = new Point3f[vertexCount];
		mesh.toArray(coords);

		// Do not try to get the colour back from the geometry as is done 
		// in the super-class. That will only work if the size is the same
		// and this method is likely to be called when the size changes.
		final Color3f colors[] = new Color3f[vertexCount];
		Arrays.fill(colors, (color == null) ? DEFAULT_COLOR : color);

		final GeometryArray ta = new TriangleArray(vertexCount,
				GeometryArray.COORDINATES | GeometryArray.COLOR_3 | GeometryArray.NORMALS);

		ta.setCoordinates(0, coords);
		ta.setColors(0, colors);

		// generate normals
		final GeometryArray result;

		if (dirty)
		{
			final GeometryInfo gi = new GeometryInfo(ta);
			final NormalGenerator ng = new NormalGenerator();
			ng.generateNormals(gi);
			result = gi.getGeometryArray();
		}
		else
		{
			// Use the same normals for each repeated object
			final Vector3f[] normals = new Vector3f[vertexCount];

			// Binary fill
			int fill = objectNormals.length;
			System.arraycopy(objectNormals, 0, normals, 0, fill);
			for (int i = 2; i < points.length; i *= 2)
			{
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

	private static int transparencyMode = TransparencyAttributes.BLENDED;

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
		Appearance appearance = super.createAppearance();
		// Update the transparency to the default mode
		final TransparencyAttributes ta = appearance.getTransparencyAttributes();
		if (ta.getTransparencyMode() != TransparencyAttributes.NONE)
			ta.setTransparencyMode(transparencyMode);
		return appearance;
	}

	/**
	 * Check the facet normals point out from the centre 0,0,0. If the normal points inwards then the vertices will be
	 * swapped.
	 *
	 * @param vertices
	 *            the vertices
	 * @return the count of the number swapped
	 */
	public static int checkFacets(Point3f[] vertices)
	{
		int count = 0;
		int nVertices = vertices.length;
		final Vector3f v1 = new Vector3f(), v2 = new Vector3f();
		for (int i = 0; i < nVertices; i += 3)
		{
			// Use the same order as that used to compute facet normals in 
			// org.scijava.java3d.utils.geometry.NormalGenerator
			v1.sub(vertices[i + 2], vertices[i + 1]);
			v2.sub(vertices[i], vertices[i + 1]);
			v1.cross(v1, v2);
			v1.normalize();

			// Project point (x,y,z) to plane with normal (a,b,c) and point (d,e,f)
			// t = (ad - ax + be - by + cd - cz) / (a^2 + b^2 + c^2)
			// projected point = (x+ta,y+tb,z+tc)

			// Project 0,0,0 to the facet
			double a = v1.x;
			double b = v1.y;
			double c = v1.z;
			double d = vertices[i].x;
			double e = vertices[i].y;
			double f = vertices[i].z;
			double t = a * d + b * e + c * f;
			if (t < 0)
			{
				count++;
				swap(vertices, i + 2, i);
			}
		}
		//System.out.printf("Swapped %d\n", count);
		return count;
	}

	private static void swap(Point3f[] vertices, int i, int j)
	{
		Point3f tmp = vertices[i];
		vertices[i] = vertices[j];
		vertices[j] = tmp;
	}

	/**
	 * Gets the normals assuming triangle vertices.
	 *
	 * @param vertices
	 *            the vertices
	 * @param creaseAngle
	 *            the crease angle (in degrees)
	 * @return the normals
	 */
	public static Vector3f[] getNormals(Point3f[] vertices, double creaseAngle)
	{
		int nVertices = vertices.length;
		Vector3f[] normals = new Vector3f[nVertices];

		final GeometryArray ta = new TriangleArray(nVertices, GeometryArray.COORDINATES | GeometryArray.NORMALS);
		ta.setCoordinates(0, vertices);
		final GeometryInfo gi = new GeometryInfo(ta);
		final NormalGenerator ng = new NormalGenerator();
		if (creaseAngle >= 0 && creaseAngle <= 180)
			ng.setCreaseAngle(creaseAngle * Math.PI / 180.0);
		ng.generateNormals(gi);
		Vector3f[] n = gi.getNormals();
		int[] indices = gi.getNormalIndices();
		for (int i = 0; i < nVertices; i++)
		{
			normals[i] = n[indices[i]];
		}

		return normals;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.UpdatedableItemMesh#reorder(int[])
	 */
	public void reorder(int[] indices) throws IllegalArgumentException
	{
		if (dirty)
			throw new IllegalArgumentException("Mesh has been modified");

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
		if (dirty)
			throw new IllegalArgumentException("Mesh has been modified");

		changed = true;

		int oldSize = size();
		int size = (indices == null) ? 0 : Math.min(oldSize, indices.length);

		if (size == 0)
		{
			mesh.clear();
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

		// Reorder all things in the geometry: coordinates and colour
		// The normals can be copied as they are unchanged.
		// The mesh should contain the same coordinates as the geometry array.
		int objectSize = objectVertices.length;
		Point3f[] oldCoords = mesh.toArray(new Point3f[mesh.size()]);
		float[] oldColors = new float[oldCoords.length * 3];
		ga.getColors(0, oldColors);
		final Point3f[] coords = new Point3f[size * objectSize];
		final float[] colors = new float[coords.length * 3];
		for (int i = 0; i < size; i++)
		{
			int j = indices[i];

			int ii = i * objectSize;
			int jj = j * objectSize;
			System.arraycopy(oldCoords, jj, coords, ii, objectSize);
			System.arraycopy(oldColors, jj * 3, colors, ii * 3, objectSize * 3);
		}
		mesh = Arrays.asList(coords);

		ga.updateData(new GeometryUpdater()
		{
			public void updateData(Geometry geometry)
			{
				GeometryArray ga = (GeometryArray) geometry;
				// We re-use the geometry and just truncate the vertex count
				ga.setCoordinates(0, coords);
				ga.setColors(0, colors);
				ga.setValidVertexCount(coords.length);
			}
		});

		//this.setGeometry(ga);
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
	 * @see gdsc.smlm.ij.ij3d.ItemMesh#setItemColor(org.scijava.vecmath.Color3f[])
	 */
	public void setItemColor(Color3f[] color) throws IllegalArgumentException
	{
		this.color = null;
		int size = size();
		if (color.length != size)
			throw new IllegalArgumentException("list of size " + size + " expected");
		final GeometryArray ga = (GeometryArray) getGeometry();
		if (ga == null)
			return;
		int objectSize = objectVertices.length;
		final int N = objectSize * size;
		final Color3f[] colors = new Color3f[N];
		int i = 0;
		for (Color3f c : color)
		{
			for (int j = objectSize; j-- > 0;)
				colors[i++] = c;
		}
		ga.setColors(0, colors);
		changed = true;
	}
}
