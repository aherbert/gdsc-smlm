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
import org.scijava.java3d.BranchGroup;
import org.scijava.java3d.GeometryArray;
import org.scijava.java3d.Group;
import org.scijava.java3d.Material;
import org.scijava.java3d.Shape3D;
import org.scijava.java3d.Transform3D;
import org.scijava.java3d.TransformGroup;
import org.scijava.java3d.TransparencyAttributes;
import org.scijava.java3d.utils.geometry.Primitive;
import org.scijava.java3d.utils.geometry.Sphere;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Point3d;
import org.scijava.vecmath.Point3f;
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

	/** The global transparency of the points. This may be combined with per item alpha. */
	private float transparency;

	/** The default appearance when not using per-item colour/transparency */
	private Appearance defaultAppearance;

	/** The per-item transparency attributes. */
	private TransparencyAttributes[] transparencyAttributes = null;

	/** The per-item material. */
	private Material[] material = null;

	public PointGroup(final Point3f[] points)
	{
		this(points, null, null, 1, null, 0f);
	}

	public PointGroup(final Point3f[] points, GeometryArray ga, Appearance appearance, float scale, Color3f color,
			float transparency)
	{
		if (points == null)
			throw new NullPointerException("Points must not be null");
		//setCapability(ALLOW_CHILDREN_EXTEND);
		//setCapability(ALLOW_CHILDREN_WRITE);
		this.points = points;
		this.transparency = transparency;
		this.defaultAppearance = createDefaultAppearance(appearance);

		defaultAppearance = createAppearance(defaultAppearance, color, transparency);

		if (ga == null)
			ga = createSphere(6);

		Transform3D t3d = new Transform3D();
		Transform3D s3d = new Transform3D();
		Vector3f v3f = new Vector3f();

		float[] coordinates = new float[ga.getVertexCount() * 3];
		ga.getCoordinates(0, coordinates);
		float[] coordinates2 = new float[coordinates.length];

		Point3d centre = new Point3d();
		
		// XXX - Write a new class that uses ordered group
		// to support custom sort of the displayed order.
		// OrderedPointGroup extends PointGroup
		Group parent = this;
		
		for (int i = 0; i < points.length; i++)
		{
			v3f.set(points[i]);

			// TODO - Allow per-item appearance
			appearance = defaultAppearance;

			if (false)
			{
				// Use transforms
				t3d.set(v3f);
				final TransformGroup tg = new TransformGroup(t3d);

				// Note: Using a non-uniform scale will results in
				// reduced performanec. See the notes in TransformGroup javadoc.
				s3d.setScale(scale);
				final TransformGroup sg = new TransformGroup(s3d);
				tg.addChild(sg);

				Shape3D shape = new Shape3D(ga, appearance);
				shape.setBounds(new BoundingSphere(centre, scale));
				shape.setBoundsAutoCompute(false);

				sg.addChild(shape);
				
				parent.addChild(tg);
			}
			else
			{
				// Scale and translate
				GeometryArray ga2 = (GeometryArray) ga.cloneNodeComponent(true);
				for (int j = 0; j < coordinates2.length; j += 3)
				{
					coordinates2[j] = coordinates[j] * scale + v3f.x;
					coordinates2[j + 1] = coordinates[j + 1] * scale + v3f.y;
					coordinates2[j + 2] = coordinates[j + 2] * scale + v3f.z;
				}
				ga2.setCoordinates(0, coordinates2);

				// XXX - fix this to be generic
				Shape3D shape = new Shape3D(ga2, appearance);
				shape.setBounds(new BoundingSphere(new Point3d(points[i]), scale));
				shape.setBoundsAutoCompute(false);

				parent.addChild(shape);
			}
		}
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
	 * @param c
	 */
	public void setColor(final Color3f c)
	{
		// Shared colour
		defaultAppearance.getMaterial().setDiffuseColor(c);
	}

	/**
	 * Set the transparency of the points.
	 * 
	 * @param t
	 */
	public void setTransparency(final float transparency)
	{
		this.transparency = transparency;

		// Shared appearance
		final TransparencyAttributes tr = defaultAppearance.getTransparencyAttributes();
		final int mode = transparency == 0f ? TransparencyAttributes.NONE : TransparencyAttributes.FASTEST;

		tr.setTransparencyMode(mode);
		tr.setTransparency(transparency);

		//		for (int i = 0; i < numChildren(); i++)
		//		{
		//			final BranchGroup bg = (BranchGroup) getChild(i);
		//			final TransformGroup tg = (TransformGroup) ((Group) bg.getChild(0)).getChild(0);
		//			final Sphere s = (Sphere) tg.getChild(0);
		//			s.setAppearance(appearance);
		//		}
	}

	/**
	 * Create an Appearance object.
	 *
	 * @param appearance
	 *            the appearance
	 * @return the appearance
	 */
	private static Appearance createDefaultAppearance(Appearance appearance)
	{
		if (appearance == null)
			appearance = new Appearance();
		appearance.setCapability(Appearance.ALLOW_TRANSPARENCY_ATTRIBUTES_READ);
		appearance.setCapability(Appearance.ALLOW_MATERIAL_READ);

		// These are the defaults. We may need them if we want to support mesh 
		// display when the polygon mode is Line
		//final PolygonAttributes polyAttrib = new PolygonAttributes();
		////polyAttrib.setCapability(PolygonAttributes.ALLOW_MODE_WRITE);
		//polyAttrib.setPolygonMode(PolygonAttributes.POLYGON_FILL);
		//polyAttrib.setCullFace(PolygonAttributes.CULL_BACK);
		//polyAttrib.setBackFaceNormalFlip(false);

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
	 * Create an Appearance object.
	 *
	 * @param color
	 *            the color
	 * @param transparency
	 *            the transparency
	 * @return the appearance
	 */
	private static Appearance createAppearance(Appearance appearance, Color3f color, float transparency)
	{
		appearance = (Appearance) appearance.cloneNodeComponent(true);

		TransparencyAttributes tr = appearance.getTransparencyAttributes();
		final int mode = transparency == 0f ? TransparencyAttributes.NONE : TransparencyAttributes.FASTEST;
		if (tr.getTransparencyMode() != mode)
			tr.setTransparencyMode(mode);
		if (tr.getTransparency() != transparency)
			tr.setTransparency(transparency);

		if (color != null)
		{
			Material material = appearance.getMaterial();
			material.setDiffuseColor(color);
		}

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
