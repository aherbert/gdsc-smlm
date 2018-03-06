package gdsc.smlm.ij.ij3d;

import org.scijava.java3d.Appearance;
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

import customnode.CustomIndexedTriangleMesh;

/**
 * Use an indexed triangle mesh object to represent a set of points. The object is duplicated, scaled and translated for
 * each point.
 */
public class RepeatedIndexedTriangleMesh extends CustomIndexedTriangleMesh
{
	protected Point3f[] objectVertices;
	protected int[] objectFaces;
	protected Vector3f[] objectNormals;
	protected Point3f[] points;
	protected Point3f[] sizes;
	private boolean dirty = false;

	/**
	 * Instantiates a new repeated indexed triangle mesh.
	 * <p>
	 * This will repeat the object for each input point. The object
	 * is assumed to be centred on the origin. It will be scaled and
	 * translated for each input point.
	 *
	 * @param objectVertices
	 *            the vertices of the object for a single point
	 * @param objectFaces
	 *            the faces of the object for a single point
	 * @param points
	 *            the points
	 * @param sizes
	 *            the size of each point
	 * @param color
	 *            the color
	 * @param transp
	 *            the transparency
	 */
	public RepeatedIndexedTriangleMesh(Point3f[] objectVertices, int[] objectFaces, Point3f[] points, Point3f[] sizes,
			Color3f color, float transp)
	{
		// Create empty 
		super(new Point3f[0], new int[0], color, transp);

		if (sizes != null && points.length != sizes.length)
			throw new IllegalArgumentException("Points and sizes must be the same length");

		this.objectVertices = objectVertices;
		this.objectFaces = objectFaces;

		checkFacets(objectVertices, objectFaces);

		// Now build the actual vertices by repeating the points.
		this.points = points;
		this.sizes = sizes;

		vertices = new Point3f[objectVertices.length * points.length];
		faces = new int[objectFaces.length * points.length];
		this.nVertices = vertices.length;
		this.nFaces = faces.length;

		final int n = objectVertices.length;
		boolean sameSize = false;
		if (sizes == null || (sameSize = RepeatedTriangleMesh.sameSize(sizes)))
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
			}
		}

		final int m = objectFaces.length;
		for (int i = 0, k = 0; i < points.length; i++)
		{
			for (int j = 0, offset = i * n; j < m; j++)
			{
				// Offset the face index by the count of vertices beforehand
				faces[k++] = objectFaces[j] + offset;
			}
		}

		colors = new Color3f[this.vertices.length];
		for (int i = 0; i < nVertices; i++)
			// This may be a default colour set in the super constructor
			colors[i] = this.color;

		// Update the geometry
		this.setGeometry(createGeometry());
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
	protected Appearance createAppearance()
	{
		return super.createAppearance();

		//		final Appearance appearance = new Appearance();
		//		appearance.setCapability(Appearance.ALLOW_TRANSPARENCY_ATTRIBUTES_READ);
		//
		//		final PolygonAttributes polyAttrib = new PolygonAttributes();
		//		polyAttrib.setCapability(PolygonAttributes.ALLOW_MODE_WRITE);
		//		if (this.shaded)
		//			polyAttrib.setPolygonMode(PolygonAttributes.POLYGON_FILL);
		//		else
		//			polyAttrib.setPolygonMode(PolygonAttributes.POLYGON_LINE);
		//		polyAttrib.setCullFace(PolygonAttributes.CULL_NONE);
		//
		//		// This is what makes the polygons look the same on both sides!
		//		//polyAttrib.setBackFaceNormalFlip(true);
		//
		//		appearance.setPolygonAttributes(polyAttrib);
		//
		//		final ColoringAttributes colorAttrib = new ColoringAttributes();
		//		colorAttrib.setShadeModel(ColoringAttributes.SHADE_GOURAUD);
		//		if (null != color) // is null when colors are vertex-wise
		//			colorAttrib.setColor(color);
		//		appearance.setColoringAttributes(colorAttrib);
		//
		//		final TransparencyAttributes tr = new TransparencyAttributes();
		//		final int mode = TransparencyAttributes.FASTEST;
		//		tr.setCapability(TransparencyAttributes.ALLOW_VALUE_WRITE);
		//		tr.setCapability(TransparencyAttributes.ALLOW_MODE_WRITE);
		//		tr.setTransparencyMode(mode);
		//		tr.setTransparency(transparency);
		//		appearance.setTransparencyAttributes(tr);
		//
		//		final Material material = new Material();
		//		material.setCapability(Material.ALLOW_COMPONENT_WRITE);
		//		material.setAmbientColor(0.1f, 0.1f, 0.1f);
		//		material.setSpecularColor(0.1f, 0.1f, 0.1f);
		//		material.setDiffuseColor(0.1f, 0.1f, 0.1f);
		//		appearance.setMaterial(material);
		//		return appearance;
	}

	@Override
	public Vector3f[] getNormals()
	{
		// Take advantage of the fact that the normals will all be the 
		// same for the repeated object. This cannot be done if the 
		// vertices have been updated.
		if (dirty)
			super.getNormals();

		// Use the same functionality as the IndexedTriangleMesh
		if (objectNormals == null)
			objectNormals = getNormals(objectVertices, objectFaces);

		final Vector3f[] normals = new Vector3f[nVertices];

		for (int i = 0, k = 0; i < points.length; i++)
		{
			for (int j = 0; j < objectNormals.length; j++)
			{
				normals[k++] = objectNormals[j];
			}
		}

		return normals;
	}

	/**
	 * Check the facet normals point out from the centre 0,0,0. If the normal points inwards then the vertices will be
	 * swapped.
	 *
	 * @param vertices
	 *            the vertices
	 * @param faces
	 *            the faces
	 * @return the count of the number swapped
	 */
	public static int checkFacets(Point3f[] vertices, int[] faces)
	{
		int count = 0;
		int nFaces = faces.length;
		final Vector3f v1 = new Vector3f(), v2 = new Vector3f();
		for (int i = 0; i < nFaces; i += 3)
		{
			final int f1 = faces[i];
			final int f2 = faces[i + 1];
			final int f3 = faces[i + 2];

			// Use the same order as that used to compute facet normals in 
			// org.scijava.java3d.utils.geometry.NormalGenerator
			v1.sub(vertices[f3], vertices[f2]);
			v2.sub(vertices[f1], vertices[f2]);
			v1.cross(v1, v2);
			v1.normalize();

			// Project point (x,y,z) to plane with normal (a,b,c) and point (d,e,f)
			// t = (ad - ax + be - by + cd - cz) / (a^2 + b^2 + c^2)
			// projected point = (x+ta,y+tb,z+tc)

			// Project 0,0,0 to the facet
			double a = v1.x;
			double b = v1.y;
			double c = v1.z;
			double d = vertices[f1].x;
			double e = vertices[f1].y;
			double f = vertices[f1].z;
			double t = a * d + b * e + c * f;
			if (t < 0)
			{
				count++;
				swap(faces, i + 2, i);
			}
		}
		//System.out.printf("Swapped %d\n", count);
		return count;
	}

	private static void swap(int[] faces, int i, int j)
	{
		int tmp = faces[i];
		faces[i] = faces[j];
		faces[j] = tmp;
	}

	/**
	 * Gets the normals assuming triangle vertices on the given faces
	 *
	 * @param vertices
	 *            the vertices
	 * @param faces
	 *            the faces
	 * @return the normals
	 */
	public static Vector3f[] getNormals(Point3f[] vertices, int[] faces)
	{
		int nVertices = vertices.length;
		int nFaces = faces.length;

		final Vector3f[] normals = new Vector3f[nVertices];
		for (int i = 0; i < nVertices; i++)
			normals[i] = new Vector3f();

		//		final IndexedTriangleArray ta = new IndexedTriangleArray(vertices.length,
		//				GeometryArray.COORDINATES | GeometryArray.NORMALS, faces.length);
		//
		//		ta.setValidIndexCount(nFaces);
		//
		//		ta.setCoordinates(0, vertices);
		//		ta.setCoordinateIndices(0, faces);
		//		final GeometryInfo gi = new GeometryInfo(ta);
		//		final NormalGenerator ng = new NormalGenerator();
		//		ng.generateNormals(gi);
		//		final Vector3f[] n = gi.getNormals();
		//
		//		for (int i = 0; i < nFaces; i++)
		//		{
		//			normals[faces[i]].add(n[i]);
		//		}

		final Vector3f v1 = new Vector3f(), v2 = new Vector3f();
		for (int i = 0; i < nFaces; i += 3)
		{
			final int f1 = faces[i];
			final int f2 = faces[i + 1];
			final int f3 = faces[i + 2];

			// Use the same order as that used to compute facet normals in 
			// org.scijava.java3d.utils.geometry.NormalGenerator
			v1.sub(vertices[f3], vertices[f2]);
			v2.sub(vertices[f1], vertices[f2]);
			v1.cross(v1, v2);
			v1.normalize();

			normals[f1].add(v1);
			normals[f2].add(v1);
			normals[f3].add(v1);
		}

		for (int i = 0; i < nVertices; i++)
			normals[i].normalize();

		//		for (int i = 0; i < nVertices; i++)
		//		{
		//			System.out.printf("[%d] %s\n", i, normals[i]);
		//		}

		return normals;
	}
}
