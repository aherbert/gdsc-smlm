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
import org.scijava.java3d.Bounds;
import org.scijava.java3d.GeometryArray;
import org.scijava.java3d.Group;
import org.scijava.java3d.Material;
import org.scijava.java3d.PointArray;
import org.scijava.java3d.PointAttributes;
import org.scijava.java3d.PolygonAttributes;
import org.scijava.java3d.Shape3D;
import org.scijava.java3d.Transform3D;
import org.scijava.java3d.TransparencyAttributes;
import org.scijava.java3d.utils.geometry.Primitive;
import org.scijava.java3d.utils.geometry.Sphere;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Color4f;
import org.scijava.vecmath.Point3f;
import org.scijava.vecmath.Vector3f;

/**
 * This class represents a list as a number of repeated shapes in the universe. The shape is defined using a geometry
 * array. Colouring is assumed to be done using the material diffuse colour. If the geometry has per vertex colours then
 * this class will not work.
 * <p>
 * A special exception is made for a PointArray as that has no surface to colour. In this case it must be created using
 * the flag
 *
 * @author Alex Herbert
 */
public class ItemGeometryGroup extends Group implements TransparentItemShape
{
	private final Color3f DEFAULT_COLOUR = new Color3f(1, 1, 0);

	// Tips: https://webserver2.tecgraf.puc-rio.br/~ismael/Cursos/Cidade_CG/labs/Java3D/Java3D_onlinebook_selman/Htmls/3DJava_Ch04.htm#4

	// See if this can extend BranchGroup to allow picking the items.

	/** The list of points */
	protected Point3f[] points;

	/**
	 * The default appearance when not using per-item colour/transparency. This is also the reference for shared
	 * attributes.
	 */
	protected Appearance defaultAppearance;

	/** The global transparency of the points. This may be combined with per item alpha. */
	protected float transparency;

	/** The default colour of the item. */
	protected Color3f color = new Color3f();

	/** The per-item alpha. */
	protected float[] alphas;

	/** Flag indicating if this is a point array. */
	protected final boolean isPointArray;

	/** The per-item transparency attributes. */
	protected TransparencyAttributes[] transparencyAttributes = null;

	/** The per-item material. Used for polygons. PointArrays use color4 coordinates. */
	protected Material[] material = null;

	/** The per-item material. Used for polygons. PointArrays use color4 coordinates. */
	protected GeometryArray[] geometryArray = null;

	private static abstract class PointArrayColorUpdater
	{
		/** The per-item colour. Used for PointArrays. PointArrays use color4 coordinates. */
		final float[] pointColor;

		PointArrayColorUpdater(int n)
		{
			pointColor = new float[4 * n];
		}

		/**
		 * Set the colour and alpha value and gets the current color array.
		 *
		 * @param alpha
		 *            the alpha
		 * @return the color
		 */
		public abstract float[] getColors(Color4f color4f);

		/**
		 * Set the colour and alpha value and gets the current color array.
		 *
		 * @param alpha
		 *            the alpha
		 * @return the color
		 */
		public abstract float[] getColors(Color3f color, float alpha);

		/**
		 * Set the colour value and gets the current color array.
		 *
		 * @param alpha
		 *            the alpha
		 * @return the color
		 */
		public abstract float[] getColors(Color3f color);

		/**
		 * Set the alpha value and gets the current color array.
		 *
		 * @param alpha
		 *            the alpha
		 * @return the color
		 */
		public abstract float[] getColors(float alpha);
	}

	private static class SinglePointArrayColorUpdater extends PointArrayColorUpdater
	{
		SinglePointArrayColorUpdater()
		{
			super(1);
		}

		public float[] getColors(Color4f color)
		{
			pointColor[0] = color.x;
			pointColor[1] = color.y;
			pointColor[2] = color.z;
			pointColor[3] = color.w;
			return pointColor;
		}

		public float[] getColors(Color3f color, float alpha)
		{
			pointColor[0] = color.x;
			pointColor[1] = color.y;
			pointColor[2] = color.z;
			pointColor[3] = alpha;
			return pointColor;
		}

		public float[] getColors(Color3f color)
		{
			pointColor[0] = color.x;
			pointColor[1] = color.y;
			pointColor[2] = color.z;
			return pointColor;
		}

		public float[] getColors(float alpha)
		{
			pointColor[3] = alpha;
			return pointColor;
		}
	}

	private static class MultiPointArrayColorUpdater extends PointArrayColorUpdater
	{
		MultiPointArrayColorUpdater(int n)
		{
			super(n);
		}

		public float[] getColors(Color4f color)
		{
			for (int i = 0; i < pointColor.length;)
			{
				pointColor[i++] = color.x;
				pointColor[i++] = color.y;
				pointColor[i++] = color.z;
				pointColor[i++] = color.w;
			}
			return pointColor;
		}

		public float[] getColors(Color3f color, float alpha)
		{
			for (int i = 0; i < pointColor.length;)
			{
				pointColor[i++] = color.x;
				pointColor[i++] = color.y;
				pointColor[i++] = color.z;
				pointColor[i++] = alpha;
			}
			return pointColor;
		}

		public float[] getColors(Color3f color)
		{
			for (int i = 0; i < pointColor.length;)
			{
				pointColor[i++] = color.x;
				pointColor[i++] = color.y;
				pointColor[i++] = color.z;
				i++;
			}
			return pointColor;
		}

		public float[] getColors(float alpha)
		{
			for (int i = 3; i < pointColor.length; i += 4)
			{
				pointColor[i] = alpha;
			}
			return pointColor;
		}
	}

	protected final PointArrayColorUpdater pointArrayColorUpdater;

	/**
	 * Instantiates a new item geometry group.
	 *
	 * @param points
	 *            the points
	 */
	public ItemGeometryGroup(final Point3f[] points)
	{
		this(points, null, null, null, null, null);
	}

	/**
	 * Instantiates a new item geometry group. The geometry is scaled and translated for each point. The default
	 * appearance is cloned per item and optionally can be updated with per-item material colour and transparency. The
	 * transparency in the appearance is blended with per-item alphas to create a per-item transparency.
	 *
	 * @param points
	 *            the points.
	 * @param ga
	 *            the geometry array. If null then a default will be used. Assumed to be centred on the origin.
	 * @param appearance
	 *            the default appearance of the shape. PolygonAttributes, Material and TransparencyAttributes are used.
	 * @param sizes
	 *            the sizes of each point. Can be null (no scaling); length=1 (fixed scaling); or points.length.
	 * @param colors
	 *            the per-item colors. Can be null.
	 * @param alphas
	 *            the per-item alphas. Can be null.
	 */
	public ItemGeometryGroup(final Point3f[] points, GeometryArray ga, Appearance appearance, Point3f[] sizes,
			Color3f[] colors, float[] alphas)
	{
		if (points == null)
			throw new NullPointerException("Points must not be null");

		if (ga == null)
		{
			// Default geometry
			ga = createSphere(6);
			isPointArray = false;
		}
		else
		{
			isPointArray = (ga instanceof PointArray);
			// Check input geometry
			int format = ga.getVertexFormat();
			if (isPointArray)
			{
				if ((format & GeometryArray.COLOR_4) == 0)
					throw new NullPointerException("PointArray must have COLOR_4 vertex type");
			}
			else if ((format & GeometryArray.COLOR_3) != 0)
			{
				throw new NullPointerException("GeometryArray must not have COLOR vertex type");
			}
		}

		this.points = points;
		this.defaultAppearance = createDefaultAppearance(appearance, ga);
		this.transparency = defaultAppearance.getTransparencyAttributes().getTransparency();
		if (isPointArray)
		{
			pointArrayColorUpdater = (ga.getValidVertexCount() == 1) ? new SinglePointArrayColorUpdater()
					: new MultiPointArrayColorUpdater(ga.getValidVertexCount());

			// Get the first color as the default
			ga.getColor(0, pointArrayColorUpdater.pointColor);
			// Only uses index [0,1,2] so ignores the transparency
			color.set(pointArrayColorUpdater.pointColor);
			if (color.x == 0 && color.y == 0 && color.z == 0)
			{
				this.color.set(DEFAULT_COLOUR);
				// Update the input to the default
				ga = (GeometryArray) ga.cloneNodeComponent(true);
				ga.setColors(0, pointArrayColorUpdater.getColors(color, 1));
			}
		}
		else
		{
			pointArrayColorUpdater = null;
			defaultAppearance.getMaterial().getDiffuseColor(this.color);
		}

		// Get the bounds so we can set the centroid and bounds for each object 
		Bounds bounds = new Shape3D(ga, null).getBounds();

		final boolean hasColor = colors != null && colors.length == points.length;
		final boolean hasAlpha = alphas != null && alphas.length == points.length;

		this.alphas = alphas;

		if (isPointArray)
		{
			geometryArray = new GeometryArray[points.length];
		}
		else
		{
			transparencyAttributes = new TransparencyAttributes[points.length];
			material = new Material[points.length];
		}

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

		// Allow use of an ordered group (in sub-class)
		// to support custom sort of the displayed order.
		Group parent = getParentGroup();

		Transform3D t3d = new Transform3D();
		Vector3f v3f = new Vector3f();
		final float alpha1 = 1 - transparency;

		TransparencyAttributes ta = defaultAppearance.getTransparencyAttributes();
		Material m = defaultAppearance.getMaterial();

		for (int i = 0; i < points.length; i++)
		{
			v3f.set(points[i]);

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

			// Store the point index in the geometry for intersection analysis
			ga2.setUserData(i);

			// Allow per-item appearance with shared attributes. 
			// This allows a global transparency for PointArray.
			appearance = (Appearance) defaultAppearance.cloneNodeComponent(false);
			if (!isPointArray)
			{
				// Not shared attributes
				appearance.setTransparencyAttributes((TransparencyAttributes) ta.cloneNodeComponent(true));
				appearance.setMaterial((Material) m.cloneNodeComponent(true));
			}

			if (perItem)
			{
				if (isPointArray)
				{
					Color3f c = (hasColor) ? colors[i] : color;
					float alpha = (hasAlpha) ? alphas[i] : 1f;
					ga2.setColors(0, pointArrayColorUpdater.getColors(c, alpha));
				}
				else
				{
					// Note that this entire class is based on an assumption that setting
					// the colour using attributes is faster/easier than if the input GA
					// has per-vertex colours. It is definitely easier to support more GA
					// formats and types (e.g. indexed/stripped) if appearance is used per 
					// item. No testing has been done on the speed the image is rendered. 

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
			}

			// Store to allow fast update
			if (isPointArray)
			{
				geometryArray[i] = ga2;
			}
			else
			{
				transparencyAttributes[i] = appearance.getTransparencyAttributes();
				material[i] = appearance.getMaterial();
			}

			Shape3D shape = new Shape3D(ga2, appearance);
			// Each object can be picked. Is this needed?
			//shape.setCapability(Shape3D.ENABLE_PICK_REPORTING);
			shape.setPickable(true);

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
	 * Set the color of the points.
	 *
	 * @param color
	 *            the new color
	 */
	public void setColor(Color3f color)
	{
		setItemColor(color);
	}

	/**
	 * Gets the default color of the points.
	 *
	 * @return the color
	 */
	public Color3f getColor()
	{
		return new Color3f(color);
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

		if (isPointArray)
		{
			// Global transparency
			defaultAppearance.getTransparencyAttributes().setTransparency(transparency);
		}
		else
		{

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
				// All items the same transparency
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
	}

	/**
	 * Gets the transparency of the points. This is blended with a per-item alpha if present.
	 *
	 * @return the transparency
	 */
	public float getTransparency()
	{
		return transparency;
	}

	/**
	 * Checks if any item is transparent.
	 *
	 * @return true, if is transparent
	 */
	public boolean isTransparent()
	{
		if (transparency != 0)
			return true;
		if (alphas != null && alphas.length == points.length)
		{
			for (int i = 0; i < alphas.length; i++)
				if (alphas[i] != 1)
					return true;
		}
		return false;
	}

	/**
	 * Sets the shaded.
	 *
	 * @param shaded
	 *            the new shaded
	 */
	public void setShaded(boolean shaded)
	{
		PolygonAttributes pa = defaultAppearance.getPolygonAttributes();
		if (pa == null)
			return;
		int mode = (shaded) ? PolygonAttributes.POLYGON_FILL : PolygonAttributes.POLYGON_LINE;
		if (pa.getPolygonMode() != mode)
			pa.setPolygonMode(mode);
	}

	/**
	 * Checks if is shaded.
	 *
	 * @return true, if is shaded
	 */
	public boolean isShaded()
	{
		PolygonAttributes pa = defaultAppearance.getPolygonAttributes();
		if (pa == null)
			return false;
		return pa.getPolygonMode() == PolygonAttributes.POLYGON_FILL;
	}

	/**
	 * Sets the pointSize.
	 *
	 * @param pointSize
	 *            the new pointSize
	 */
	public void setPointSize(float pointSize)
	{
		PointAttributes pa = defaultAppearance.getPointAttributes();
		if (pa == null)
			return;
		pa.setPointSize(pointSize);
	}

	/**
	 * Gets the point size.
	 *
	 * @return the point size
	 */
	public float getPointSize()
	{
		PointAttributes pa = defaultAppearance.getPointAttributes();
		if (pa == null)
			return 0f;
		return pa.getPointSize();
	}

	/**
	 * Create a default Appearance object. This will have the correct attributes and capability bits set to manipulate
	 * the material and transparency.
	 *
	 * @param appearance
	 *            the appearance
	 * @param ga
	 * @return the appearance
	 */
	private Appearance createDefaultAppearance(Appearance appearance, GeometryArray ga)
	{
		if (appearance == null)
			appearance = new Appearance();
		appearance.setCapability(Appearance.ALLOW_TRANSPARENCY_ATTRIBUTES_READ);
		appearance.setCapability(Appearance.ALLOW_MATERIAL_READ);

		if (ga instanceof PointArray)
		{
			appearance.setPolygonAttributes(null);
			appearance.setMaterial(null);

			PointAttributes pointAttributes = appearance.getPointAttributes();
			if (pointAttributes == null)
			{
				pointAttributes = new PointAttributes();
				pointAttributes.setCapability(PointAttributes.ALLOW_ANTIALIASING_WRITE);
				pointAttributes.setCapability(PointAttributes.ALLOW_SIZE_WRITE);
				pointAttributes.setPointAntialiasingEnable(true);
				appearance.setPointAttributes(pointAttributes);
			}

			// We use the coordinates for the colour
		}
		else
		{
			appearance.setPointAttributes(null);

			// These are the defaults. We may need them if we want to support mesh 
			// display when the polygon mode is Line
			PolygonAttributes polygonAttributes = appearance.getPolygonAttributes();
			if (polygonAttributes == null)
			{
				polygonAttributes = new PolygonAttributes();
				polygonAttributes.setCapability(PolygonAttributes.ALLOW_MODE_WRITE);
				polygonAttributes.setPolygonMode(PolygonAttributes.POLYGON_FILL);
				appearance.setPolygonAttributes(polygonAttributes);
			}

			// We require material attributes for colour
			Material material = appearance.getMaterial();
			if (material == null)
			{
				material = new Material();
				material.setDiffuseColor(DEFAULT_COLOUR);
				appearance.setMaterial(material);
			}
			material.setCapability(Material.ALLOW_COMPONENT_WRITE);
		}

		// We require transparency attributes for global transparency
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.ItemShape#size()
	 */
	public int size()
	{
		return points.length;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.ItemShape#setItemColor(org.scijava.vecmath.Color3f)
	 */
	public void setItemColor(Color3f color)
	{
		// Global colour
		if (color == null)
			color = this.color;
		if (isPointArray)
		{
			pointArrayColorUpdater.getColors(color, 1f);
			if (alphas != null)
			{
				for (int i = 0; i < geometryArray.length; i++)
				{
					geometryArray[i].setColors(i, pointArrayColorUpdater.getColors(alphas[i]));
				}
			}
			else
			{
				for (int i = 0; i < geometryArray.length; i++)
				{
					geometryArray[i].setColors(i, pointArrayColorUpdater.pointColor);
				}
			}
		}
		else
		{
			for (int i = 0; i < material.length; i++)
				material[i].setDiffuseColor(color);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.ItemShape#setItemColor(org.scijava.vecmath.Color3f[])
	 */
	public void setItemColor(Color3f[] color) throws IllegalArgumentException
	{
		if (color == null)
		{
			// Default
			setColor(null);
			return;
		}
		int size = size();
		if (color.length != size)
			throw new IllegalArgumentException("list of size " + size + " expected");
		if (isPointArray)
		{
			if (alphas != null)
			{
				for (int i = 0; i < geometryArray.length; i++)
				{
					geometryArray[i].setColors(i, pointArrayColorUpdater.getColors(color[i], alphas[i]));
				}
			}
			else
			{
				pointArrayColorUpdater.getColors(1f);
				for (int i = 0; i < geometryArray.length; i++)
				{
					geometryArray[i].setColors(i, pointArrayColorUpdater.getColors(color[i]));
				}
			}
		}
		else
		{
			for (int i = 0; i < material.length; i++)
				material[i].setDiffuseColor(color[i]);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.TransparentItemShape#setItemColor4(org.scijava.vecmath.Color4f[])
	 */
	public void setItemColor4(Color4f[] color) throws IllegalArgumentException
	{
		if (color == null)
		{
			this.alphas = null;

			// Default
			setColor(null);
			setTransparency(this.transparency);
			return;
		}
		int size = size();
		if (color.length != size)
			throw new IllegalArgumentException("list of size " + size + " expected");
		if (alphas == null)
			alphas = new float[size];
		if (isPointArray)
		{
			for (int i = 0; i < geometryArray.length; i++)
			{
				geometryArray[i].setColors(i, pointArrayColorUpdater.getColors(color[i]));
				alphas[i] = color[i].w;
			}
		}
		else
		{
			for (int i = 0; i < material.length; i++)
			{
				material[i].setDiffuseColor(color[i].x, color[i].y, color[i].z);
				alphas[i] = color[i].w;
			}
			setTransparency(this.transparency);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.TransparentItemShape#setItemAlpha(float[])
	 */
	public void setItemAlpha(float[] alpha) throws IllegalArgumentException
	{
		if (isPointArray)
		{
			final boolean hasAlpha;
			if (alpha != null)
			{
				int size = size();
				if (alpha.length != size)
					throw new IllegalArgumentException("list of size " + size + " expected");
				hasAlpha = true;
			}
			else
				hasAlpha = false;
			this.alphas = alpha;
			for (int i = 0; i < geometryArray.length; i++)
			{
				GeometryArray ga = geometryArray[i];
				ga.getColors(0, pointArrayColorUpdater.pointColor);
				pointArrayColorUpdater.getColors((hasAlpha) ? alpha[i] : 1f);
				ga.setColors(i, pointArrayColorUpdater.pointColor);
			}
		}
		else
		{
			setItemAlpha(alpha, this.transparency);
		}
	}

	/**
	 * Sets the item alpha and the global transparency in one operation.
	 *
	 * @param alpha
	 *            the alpha
	 * @param transparency
	 *            the transparency
	 * @throws IllegalArgumentException
	 *             the illegal argument exception
	 * @see #setItemAlpha(float[])
	 */
	public void setItemAlpha(float[] alpha, float transparency) throws IllegalArgumentException
	{
		if (alpha != null)
		{
			int size = size();
			if (alpha.length != size)
				throw new IllegalArgumentException("list of size " + size + " expected");
		}
		this.alphas = alpha;
		if (isPointArray)
		{
			// PointArray alpha must be updated  
			setItemAlpha(alpha);
		}
		setTransparency(this.transparency);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.TransparentItemShape#setItemAlpha(float)
	 */
	public void setItemAlpha(float alpha) throws IllegalArgumentException
	{
		if (isPointArray)
		{
			if (alphas == null)
				alphas = new float[size()];
			Arrays.fill(alphas, alpha);

			for (int i = 0; i < geometryArray.length; i++)
			{
				GeometryArray ga = geometryArray[i];
				ga.getColors(0, pointArrayColorUpdater.pointColor);
				pointArrayColorUpdater.getColors(alpha);
				ga.setColors(i, pointArrayColorUpdater.pointColor);
			}
		}
		else
		{
			setItemAlpha(alpha, this.transparency);
		}
	}

	/**
	 * Sets the item alpha and the global transparency in one operation.
	 *
	 * @param alpha
	 *            the alpha
	 * @param transparency
	 *            the transparency
	 * @throws IllegalArgumentException
	 *             the illegal argument exception
	 * @see #setItemAlpha(float)
	 */
	public void setItemAlpha(float alpha, float transparency) throws IllegalArgumentException
	{
		// Reuse current alpha storage
		if (alphas == null)
			alphas = new float[size()];
		Arrays.fill(alphas, alpha);
		if (isPointArray)
		{
			// PointArray alpha must be updated  
			for (int i = 0; i < geometryArray.length; i++)
			{
				GeometryArray ga = geometryArray[i];
				ga.getColors(0, pointArrayColorUpdater.pointColor);
				pointArrayColorUpdater.getColors(alpha);
				ga.setColors(i, pointArrayColorUpdater.pointColor);
			}
		}
		setTransparency(transparency);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.TransparentItemShape#getItemAlpha(float[])
	 */
	public void getItemAlpha(float[] alpha) throws IllegalArgumentException
	{
		int size = size();
		if (alpha.length != size)
			throw new IllegalArgumentException("list of size " + size + " expected");
		if (this.alphas == null)
			// No alpha
			Arrays.fill(alpha, 1f);
		else
			System.arraycopy(this.alphas, 0, alpha, 0, size);
	}
}
