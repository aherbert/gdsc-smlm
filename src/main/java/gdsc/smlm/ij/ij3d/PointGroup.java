package gdsc.smlm.ij.ij3d;

import java.util.Arrays;

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

import org.scijava.java3d.Appearance;
import org.scijava.java3d.BoundingSphere;
import org.scijava.java3d.Bounds;
import org.scijava.java3d.BranchGroup;
import org.scijava.java3d.GeometryArray;
import org.scijava.java3d.Group;
import org.scijava.java3d.Material;
import org.scijava.java3d.PolygonAttributes;
import org.scijava.java3d.Shape3D;
import org.scijava.java3d.Transform3D;
import org.scijava.java3d.TransformGroup;
import org.scijava.java3d.TransparencyAttributes;
import org.scijava.java3d.utils.geometry.Primitive;
import org.scijava.java3d.utils.geometry.Sphere;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Point3d;
import org.scijava.vecmath.Point3f;
import org.scijava.vecmath.Vector3d;
import org.scijava.vecmath.Vector3f;

/**
 * This class represents a list as a number of objects in the universe.
 * <p>
 * This class is based on ideas in ij3d.pointlist.PointListShape.
 *
 * @author Alex Herbert
 */
public class PointGroup extends Group
{
	// Tips: https://webserver2.tecgraf.puc-rio.br/~ismael/Cursos/Cidade_CG/labs/Java3D/Java3D_onlinebook_selman/Htmls/3DJava_Ch04.htm#4

	// TODO -
	// Make the input require an object type, per item XYZ radius, colour and alpha.
	// Build transparency using a global transparency and per item alpha.

	// Create a list of the primitives appearance (Transparency/Material attributes) 
	// for direct access to them.
	// Use getAllChildren() to enumerate the actual children.
	// Add methods to change the transparency and colour
	// Transparency should be combined with per item alpha.

	// Basically extent TranparentItemMesh. Change the interface names of Mesh to Node, 
	// i.e. something that can be added to a scene graph. 

	// See if this can extend BranchGroup to allow picking the items.

	/** The list of points */
	private Point3f[] points;

	/** The default appearance when not using per-item colour/transparency */
	private Appearance defaultAppearance;

	/** The polygon attributes. These are shared. */
	private PolygonAttributes pa;

	/** The global transparency of the points. This may be combined with per item alpha. */
	private float transparency;

	/** The single colour of the material. */
	private Color3f color = new Color3f();

	/** The per-item colour. */
	private Color3f[] colors;

	/** The per-item alpha. */
	private float[] alphas;

	/** The per-item transparency attributes. */
	private TransparencyAttributes[] transparencyAttributes = null;

	/** The per-item material. */
	private Material[] material = null;

	public PointGroup(final Point3f[] points)
	{
		this(points, null, null, null, null, null);
	}

	public PointGroup(final Point3f[] points, GeometryArray ga, Appearance appearance, Point3f[] sizes,
			Color3f[] colors, float[] alphas)
	{
		if (points == null)
			throw new NullPointerException("Points must not be null");
		//setCapability(ALLOW_CHILDREN_EXTEND);
		//setCapability(ALLOW_CHILDREN_WRITE);
		this.points = points;
		this.defaultAppearance = createDefaultAppearance(appearance);
		this.transparency = defaultAppearance.getTransparencyAttributes().getTransparency();
		defaultAppearance.getMaterial().getDiffuseColor(this.color);

		// Initialise the geometry
		if (ga == null)
			ga = createSphere(6);
		Bounds bounds = new Shape3D(ga, null).getBounds();

		final boolean hasColor = colors != null && colors.length == points.length;
		final boolean hasAlpha = alphas != null && alphas.length == points.length;

		this.colors = colors;
		this.alphas = alphas;

		transparencyAttributes = new TransparencyAttributes[points.length];
		material = new Material[points.length];

		// Flag for creating per-item appearance
		final boolean perItem = hasColor || hasAlpha;

		float[] coordinates = new float[ga.getVertexCount() * 3];
		ga.getCoordinates(0, coordinates);
		float[] coordinates2 = new float[coordinates.length];

		// Handle a single scale
		final boolean hasSize;
		if (sizes == null || ItemTriangleMesh.sameSize(sizes))
		{
			hasSize = false;
			if (sizes != null)
			{
				// Scale the input object by the fixed scale
				final Point3f s = sizes[0];
				final float sx = s.x;
				final float sy = s.y;
				final float sz = s.z;
				for (int j = 0; j < coordinates.length; j += 3)
				{
					coordinates[j] *= sx;
					coordinates[j + 1] *= sy;
					coordinates[j + 2] *= sz;
				}
			}
		}
		else
			hasSize = true;

		// XXX - Write a new class that uses ordered group
		// to support custom sort of the displayed order.
		// OrderedPointGroup extends PointGroup
		Group parent = getParentGroup();

		Transform3D t3d = new Transform3D();
		Vector3f v3f = new Vector3f();
		final float alpha1 = 1 - transparency;

		for (int i = 0; i < points.length; i++)
		{
			v3f.set(points[i]);

			// Allow per-item appearance
			appearance = (Appearance) defaultAppearance.cloneNodeComponent(true);
			appearance.setPolygonAttributes(pa); // Shared
			if (perItem)
			{
				if (hasAlpha)
				{
					// Combine alphas to get the transparency
					float t = 1 - (alpha1 * alphas[i]);
					TransparencyAttributes tr = appearance.getTransparencyAttributes();
					final int mode = t == 0f ? TransparencyAttributes.NONE : TransparencyAttributes.FASTEST;
					//if (tr.getTransparencyMode() != mode)
					tr.setTransparencyMode(mode);
					//if (tr.getTransparency() != t)
					tr.setTransparency(t);
				}

				if (hasColor)
				{
					Material material = appearance.getMaterial();
					material.setDiffuseColor(colors[i]);
				}
			}

			// Store to allow fast update
			transparencyAttributes[i] = appearance.getTransparencyAttributes();
			material[i] = appearance.getMaterial();

			// Scale and translate
			GeometryArray ga2 = (GeometryArray) ga.cloneNodeComponent(true);
			if (hasSize)
			{
				final float sx = sizes[i].x;
				final float sy = sizes[i].y;
				final float sz = sizes[i].z;
				for (int j = 0; j < coordinates2.length; j += 3)
				{
					coordinates2[j] = coordinates[j] * sx + v3f.x;
					coordinates2[j + 1] = coordinates[j + 1] * sy + v3f.y;
					coordinates2[j + 2] = coordinates[j + 2] * sz + v3f.z;
				}
			}
			else
			{
				for (int j = 0; j < coordinates2.length; j += 3)
				{
					coordinates2[j] = coordinates[j] + v3f.x;
					coordinates2[j + 1] = coordinates[j + 1] + v3f.y;
					coordinates2[j + 2] = coordinates[j + 2] + v3f.z;
				}
			}
			ga2.setCoordinates(0, coordinates2);

			Shape3D shape = new Shape3D(ga2, appearance);

			// Transform the bounds
			t3d.set(v3f);
			Bounds bounds2 = (Bounds) bounds.clone();
			bounds2.transform(t3d);
			shape.setBounds(bounds2);
			shape.setBoundsAutoCompute(false);

			parent.addChild(shape);
		}
	}

	/**
	 * Gets the parent group to which all the shapes should be added.
	 *
	 * @return the parent group
	 */
	protected Group getParentGroup()
	{
		return this;
	}

	/**
	 * Gets the points.
	 *
	 * @return the points
	 */
	public Point3f[] getPoints()
	{
		return points;
	}

	/**
	 * Get the number of points.
	 *
	 * @return the size
	 */
	public int size()
	{
		return points.length;
	}

	/**
	 * Set the color of the points.
	 *
	 * @param color
	 *            the new color
	 */
	public void setColor(final Color3f color)
	{
		// Global colour
		for (int i = 0; i < material.length; i++)
			material[i].setDiffuseColor(color);
	}

	/**
	 * Set the transparency of the points. This is blended with a per-item alpha if present.
	 *
	 * @param transparency
	 *            the new transparency
	 */
	public void setTransparency(final float transparency)
	{
		// Store global transparency
		this.transparency = transparency;

		final boolean hasAlpha = alphas != null && alphas.length == points.length;
		if (hasAlpha)
		{
			// Combine with alpha
			final float alpha1 = 1 - transparency;
			for (int i = 0; i < transparencyAttributes.length; i++)
			{
				float t = 1 - (alpha1 * alphas[i]);
				TransparencyAttributes tr = transparencyAttributes[i];
				final int mode = t == 0f ? TransparencyAttributes.NONE : TransparencyAttributes.FASTEST;
				if (tr.getTransparencyMode() != mode)
					tr.setTransparencyMode(mode);
				if (tr.getTransparency() != t)
					tr.setTransparency(t);
			}
		}
		else
		{
			// Global transparency
			if (transparency == 0f)
			{
				for (int i = 0; i < transparencyAttributes.length; i++)
				{
					transparencyAttributes[i].setTransparencyMode(TransparencyAttributes.NONE);
				}
			}
			else
			{
				for (int i = 0; i < transparencyAttributes.length; i++)
				{
					transparencyAttributes[i].setTransparencyMode(TransparencyAttributes.FASTEST);
					transparencyAttributes[i].setTransparency(transparency);
				}
			}
		}
	}

	/**
	 * Sets the shaded.
	 *
	 * @param shaded
	 *            the new shaded
	 */
	public void setShaded(boolean shaded)
	{
		// TODO see if the PA can be shared across alll appearances 
		setPolygonMode((shaded) ? PolygonAttributes.POLYGON_FILL : PolygonAttributes.POLYGON_LINE);
	}

	private void setPolygonMode(int mode)
	{
		if (pa.getPolygonMode() != mode)
			pa.setPolygonMode(mode);
	}

	/**
	 * Create a default Appearance object. This will have the correct attributes and capability bits set to manipulate
	 * the material and transparency.
	 *
	 * @param appearance
	 *            the appearance
	 * @return the appearance
	 */
	private Appearance createDefaultAppearance(Appearance appearance)
	{
		if (appearance == null)
			appearance = new Appearance();
		appearance.setCapability(Appearance.ALLOW_TRANSPARENCY_ATTRIBUTES_READ);
		appearance.setCapability(Appearance.ALLOW_MATERIAL_READ);

		// These are the defaults. We may need them if we want to support mesh 
		// display when the polygon mode is Line
		pa = appearance.getPolygonAttributes();
		if (pa == null)
		{
			pa = new PolygonAttributes();
			pa.setCapability(PolygonAttributes.ALLOW_MODE_WRITE);
			pa.setPolygonMode(PolygonAttributes.POLYGON_FILL);
			//pa.setCullFace(PolygonAttributes.CULL_BACK);
			//pa.setBackFaceNormalFlip(false);
		}
		else
		{
			// Remove these to speed up cloning
			appearance.setPolygonAttributes(null);
		}

		// We may need this to choose a different shade model
		//		final ColoringAttributes colorAttrib = new ColoringAttributes();
		//		colorAttrib.setShadeModel(ColoringAttributes.SHADE_GOURAUD);
		//		// This is ignored. The material colour is used.
		//		colorAttrib.setColor(color);
		//		colorAttrib.setCapability(ColoringAttributes.ALLOW_COLOR_WRITE);
		//		appearance.setColoringAttributes(colorAttrib);

		TransparencyAttributes tr = appearance.getTransparencyAttributes();
		if (tr == null)
		{
			tr = new TransparencyAttributes();
			tr.setTransparencyMode(TransparencyAttributes.NONE);
			tr.setTransparency(0f);
			appearance.setTransparencyAttributes(tr);
		}
		tr.setCapability(TransparencyAttributes.ALLOW_VALUE_WRITE);
		tr.setCapability(TransparencyAttributes.ALLOW_MODE_WRITE);

		Material material = appearance.getMaterial();
		if (material == null)
		{
			material = new Material();
			material.setDiffuseColor(1, 1, 0);
			appearance.setMaterial(material);
		}
		material.setCapability(Material.ALLOW_COMPONENT_WRITE);

		return appearance;
	}

	/**
	 * Returns a string representation of the underlying List<Point3f>.
	 */
	@Override
	public String toString()
	{
		return Arrays.toString(points);
	}

	/**
	 * Creates the geometry for a sphere using the given number of divisions.
	 *
	 * @param divisions
	 *            the divisions
	 * @return the geometry array
	 */
	public static GeometryArray createSphere(int divisions)
	{
		int flags = Primitive.GENERATE_NORMALS | Primitive.ENABLE_APPEARANCE_MODIFY | Primitive.ENABLE_GEOMETRY_PICKING;
		final Sphere sphere = new Sphere(1, flags, divisions, null);
		return (GeometryArray) sphere.getShape(0).getGeometry();
	}
}
