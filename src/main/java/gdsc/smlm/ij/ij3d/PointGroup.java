package gdsc.smlm.ij.ij3d;

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
import org.scijava.java3d.BranchGroup;
import org.scijava.java3d.Group;
import org.scijava.java3d.Material;
import org.scijava.java3d.Transform3D;
import org.scijava.java3d.TransformGroup;
import org.scijava.java3d.TransparencyAttributes;
import org.scijava.java3d.utils.geometry.Primitive;
import org.scijava.java3d.utils.geometry.Sphere;
import org.scijava.vecmath.Color3f;
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

	/** The color of the points */
	private Color3f color;

	/** The radius of the points */
	private float radius;

	/** The default appearance which is used for new points */
	private Appearance appearance;

	/**
	 * Constructor. Creates a new BranchGroup, and adds all points from the
	 * specified List<Point3f> as points to it.
	 * 
	 * @param points
	 */
	public PointGroup(final Point3f[] points)
	{
		this(points, 1, null);
	}

	/**
	 * Constructor. Creates a new BranchGroup, and adds all points from the
	 * specified List<Point3f> as points to it.
	 * 
	 * @param points
	 */
	public PointGroup(final Point3f[] points, float radius, Color3f color)
	{
		if (points == null)
			throw new NullPointerException("Points must not be null");
		//setCapability(ALLOW_CHILDREN_EXTEND);
		//setCapability(ALLOW_CHILDREN_WRITE);
		this.points = points;
		this.radius = radius;
		this.color = color;
		initAppearance();
		for (int i = 0; i < points.length; i++)
		{
			addPointToGeometry(points[i]);
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
	 * Gets the color.
	 *
	 * @return the color
	 */
	public Color3f getColor()
	{
		return color;
	}

	/**
	 * Returns the radius of the points.
	 *
	 * @return the radius
	 */
	public float getRadius()
	{
		return radius;
	}

	/**
	 * Set the radius of the points.
	 * 
	 * @param r
	 */
	public void setRadius(final float r)
	{
		this.radius = r;
		final Transform3D t3d = new Transform3D();
		t3d.setScale(r);
		for (int i = 0; i < numChildren(); i++)
		{
			final BranchGroup bg = (BranchGroup) getChild(i);
			final TransformGroup tg = (TransformGroup) bg.getChild(0);
			final TransformGroup sg = (TransformGroup) tg.getChild(0);
			//sg.getTransform(t3d);
			//t3d.setScale(radius);
			sg.setTransform(t3d);
		}
	}

	/**
	 * Set the color of the points.
	 * 
	 * @param c
	 */
	public void setColor(final Color3f c)
	{
		color = c;
		initAppearance();
		for (int i = 0; i < numChildren(); i++)
		{
			final BranchGroup bg = (BranchGroup) getChild(i);
			final TransformGroup tg = (TransformGroup) ((Group) bg.getChild(0)).getChild(0);
			final Sphere s = (Sphere) tg.getChild(0);
			s.setAppearance(appearance);
		}
	}

	/**
	 * Set the transparency of the points.
	 * 
	 * @param t
	 */
	public void setTransparency(final float transparency)
	{
		// Shared appearance
		final TransparencyAttributes tr = appearance.getTransparencyAttributes();
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

	/*
	 * *************************************************************
	 * Private methods for updating the scenegraph
	 *************************************************************/
	private final Transform3D t3d = new Transform3D();
	private final Transform3D s3d = new Transform3D();
	private final Vector3f v3f = new Vector3f();

	/**
	 * Add a new Point to the scenegraph
	 * 
	 * @param p
	 */
	private void addPointToGeometry(final Point3f p)
	{
		v3f.x = p.x;
		v3f.y = p.y;
		v3f.z = p.z;

		final BranchGroup bg = new BranchGroup();

		t3d.set(v3f);
		final TransformGroup tg = new TransformGroup(t3d);
		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		bg.addChild(tg);

		s3d.setScale(radius);
		final TransformGroup sg = new TransformGroup(s3d);
		sg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		tg.addChild(sg);

		// Allow different primitives
		int flags = Primitive.GENERATE_NORMALS | Primitive.ENABLE_APPEARANCE_MODIFY | Primitive.ENABLE_GEOMETRY_PICKING;
		//final Primitive primative = new Box(1, 1, 1, flags, appearance);
		int divisions = 6;
		final Primitive primative = new Sphere(1, flags, divisions, appearance);

		//primative.setPickable(false);

		sg.addChild(primative);

		addChild(bg);
	}

	/**
	 * Create a default Appearance object using the specified color.
	 * 
	 * @param color
	 */
	private void initAppearance()
	{
		Color3f color = this.color == null ? new Color3f(1, 1, 0) : this.color;

		appearance = new Appearance();
		appearance.setCapability(Appearance.ALLOW_TRANSPARENCY_ATTRIBUTES_READ);
		appearance.setCapability(Appearance.ALLOW_MATERIAL_READ);

		//final PolygonAttributes polyAttrib = new PolygonAttributes();
		////polyAttrib.setCapability(PolygonAttributes.ALLOW_MODE_WRITE);
		//polyAttrib.setPolygonMode(PolygonAttributes.POLYGON_FILL);
		//polyAttrib.setCullFace(PolygonAttributes.CULL_BACK);
		//polyAttrib.setBackFaceNormalFlip(false);

		//		// This is ignored. The material colour is used.
		//		final ColoringAttributes colorAttrib = new ColoringAttributes();
		//		colorAttrib.setShadeModel(ColoringAttributes.SHADE_GOURAUD);
		//		//colorAttrib.setColor(color);
		//		colorAttrib.setCapability(ColoringAttributes.ALLOW_COLOR_WRITE);
		//		appearance.setColoringAttributes(colorAttrib);

		final TransparencyAttributes tr = new TransparencyAttributes();
		float transparency = 0.6f;
		//float transparency = 0f;
		final int mode = transparency == 0f ? TransparencyAttributes.NONE : TransparencyAttributes.FASTEST;
		tr.setCapability(TransparencyAttributes.ALLOW_VALUE_WRITE);
		tr.setCapability(TransparencyAttributes.ALLOW_MODE_WRITE);
		tr.setTransparencyMode(mode);
		tr.setTransparency(transparency);
		appearance.setTransparencyAttributes(tr);

		final Material material = new Material();
		material.setCapability(Material.ALLOW_COMPONENT_WRITE);
		//material.setAmbientColor(color);
		//material.setSpecularColor(color);
		material.setDiffuseColor(color);
		material.setShininess(128f);
		appearance.setMaterial(material);
	}

	/**
	 * Returns a string representation of the underlying List<Point3f>.
	 */
	@Override
	public String toString()
	{
		return points.toString();
	}
}
