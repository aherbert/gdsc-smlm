package gdsc.smlm.ij.ij3d;

import java.util.Arrays;
import java.util.List;

import org.scijava.java3d.Geometry;
import org.scijava.java3d.GeometryArray;
import org.scijava.java3d.GeometryUpdater;
import org.scijava.java3d.PointArray;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Color4f;
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

/**
 * Create an object to represent a set of points
 */
public class TransparentItemPointMesh extends ItemPointMesh implements TransparentItemShape
{
	/**
	 * Instantiates a new transparent item point mesh.
	 *
	 * @param mesh
	 *            the mesh
	 */
	public TransparentItemPointMesh(final List<Point3f> mesh)
	{
		super(mesh);
	}

	/**
	 * Instantiates a new transparent item point mesh.
	 *
	 * @param mesh
	 *            the mesh
	 * @param color
	 *            the color
	 * @param transparency
	 *            the transparency
	 */
	public TransparentItemPointMesh(final List<Point3f> mesh, final Color3f color, final float transparency)
	{
		super(mesh, color, transparency);
	}

	@Override
	protected GeometryArray createGeometry()
	{
		if (mesh == null || mesh.size() == 0)
			return null;
		final int size = size();

		final Point3f[] coords = new Point3f[size];
		mesh.toArray(coords);

		final Color4f[] colors = new Color4f[size];
		if (color == null)
			color = DEFAULT_COLOR;
		Arrays.fill(colors, new Color4f(color.x, color.y, color.z, 1));

		GeometryArray ta = new PointArray(size, GeometryArray.COORDINATES | GeometryArray.COLOR_4);

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
			mesh.clear();
			this.setGeometry(null);
			return;
		}

		// From here on we assume the current geometry will not be null
		// as this only happens when the original size is zero. Size has 
		// been checked at this point to be the smaller of new and old. 
		GeometryArray ga = (GeometryArray) getGeometry();

		// Reorder all things in the geometry: coordinates and colour
		Point3f[] oldCoords = mesh.toArray(new Point3f[oldSize]);
		float[] oldColors = new float[oldSize * 4];
		ga.getColors(0, oldColors);
		final Point3f[] coords = new Point3f[size];
		final float[] colors = new float[size * 4];
		for (int i = 0; i < size; i++)
		{
			int j = indices[i];
			coords[i] = oldCoords[j];
			System.arraycopy(oldColors, j * 4, colors, i * 4, 4);
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
		final float[] colors = new float[4 * size];
		ga.getColors(0, colors);
		int i = 0;
		for (Color3f c : color)
		{
			colors[i++] = c.x;
			colors[i++] = c.y;
			colors[i++] = c.z;
			i++; // Skip over alpha
		}
		ga.setColors(0, colors);
		changed = true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.TransparentItemMesh#setItemColor4(org.scijava.vecmath.Color4f[])
	 */
	public void setItemColor4(Color4f[] color) throws IllegalArgumentException
	{
		this.color = null;
		int size = size();
		if (color.length != size)
			throw new IllegalArgumentException("list of size " + size + " expected");
		final GeometryArray ga = (GeometryArray) getGeometry();
		if (ga == null)
			return;
		ga.setColors(0, color);
		changed = true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.TransparentItemMesh#setItemAlpha(float[])
	 */
	public void setItemAlpha(float[] alpha) throws IllegalArgumentException
	{
		final int size = size();
		if (alpha.length != size)
			throw new IllegalArgumentException("list of size " + size + " expected");
		final GeometryArray ga = (GeometryArray) getGeometry();
		if (ga == null)
			return;
		final float[] colors = new float[4 * size];
		ga.getColors(0, colors);
		for (int i = 0; i < size; i++)
		{
			// Set only alpha
			colors[i * 4 + 3] = alpha[i];
		}
		ga.setColors(0, colors);
		changed = true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.TransparentItemMesh#setItemAlpha(float)
	 */
	public void setItemAlpha(float alpha) throws IllegalArgumentException
	{
		final int size = size();
		final GeometryArray ga = (GeometryArray) getGeometry();
		if (ga == null)
			return;
		final float[] colors = new float[4 * size];
		ga.getColors(0, colors);
		for (int i = 0; i < size; i++)
		{
			// Set only alpha
			colors[i * 4 + 3] = alpha;
		}
		ga.setColors(0, colors);
		changed = true;
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.TransparentItemMesh#getItemAlpha(float[])
	 */
	public void getItemAlpha(float[] alpha) throws IllegalArgumentException
	{
		final int size = size();
		if (alpha.length != size)
			throw new IllegalArgumentException("list of size " + size + " expected");
		final GeometryArray ga = (GeometryArray) getGeometry();
		if (ga == null)
			return;
		final float[] colors = new float[4 * size];
		ga.getColors(0, colors);
		for (int i = 0; i < size; i++)
		{
			// Get only alpha
			alpha[i] = colors[i * 4 + 3];
		}
	}
}
