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
import org.scijava.vecmath.Vector3d;
import org.scijava.vecmath.Vector3f;

/**
 * This class represents a list as a number of repeated shapes in the universe. The shape is defined using a geometry
 * array. Colouring is assumed to be done using the material diffuse colour. If the geometry has per vertex colours then
 * this class will not work.
 * <p>
 * A special exception is made for a PointArray as that has no surface to colour. In this case it must be created using
 * the GeometryArray.COLOR_4 flag.
 *
 * @author Alex Herbert
 */
public class ItemGeometryGroup extends ItemGroup implements TransparentItemShape
{
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
		super(points);

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

		// Get the bounds so we can set the centroid and bounds for each object 
		Bounds bounds = new Shape3D(ga, null).getBounds();

		Transform3D t3d = new Transform3D();
		Vector3f translate = new Vector3f();
		Vector3d scale = new Vector3d();

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

				// Scale the bounds
				scale.set(s);
				t3d.setScale(scale);
				bounds.transform(t3d);
			}
		}
		else
			hasSize = true;

		// Allow use of an ordered group (in sub-class)
		// to support custom sort of the displayed order.
		Group parent = getParentGroup();

		final float alpha1 = 1 - transparency;

		TransparencyAttributes ta = defaultAppearance.getTransparencyAttributes();
		Material m = defaultAppearance.getMaterial();

		for (int i = 0; i < points.length; i++)
		{
			translate.set(points[i]);

			// Scale and translate
			GeometryArray ga2 = (GeometryArray) ga.cloneNodeComponent(true);
			if (hasSize)
			{
				scale.set(sizes[i]); // For scaling the bounds
				final float sx = sizes[i].x;
				final float sy = sizes[i].y;
				final float sz = sizes[i].z;
				for (int j = 0; j < coordinates2.length; j += 3)
				{
					coordinates2[j] = coordinates[j] * sx + translate.x;
					coordinates2[j + 1] = coordinates[j + 1] * sy + translate.y;
					coordinates2[j + 2] = coordinates[j + 2] * sz + translate.z;
				}
			}
			else
			{
				for (int j = 0; j < coordinates2.length; j += 3)
				{
					coordinates2[j] = coordinates[j] + translate.x;
					coordinates2[j + 1] = coordinates[j + 1] + translate.y;
					coordinates2[j + 2] = coordinates[j + 2] + translate.z;
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
			t3d.set(translate);
			if (hasSize)
				t3d.setScale(scale);
			Bounds bounds2 = (Bounds) bounds.clone();
			bounds2.transform(t3d);
			shape.setBounds(bounds2);
			shape.setBoundsAutoCompute(false);

			parent.addChild(shape);
		}
	}

	@Override
	public Color3f getColor()
	{
		return new Color3f(color);
	}

	@Override
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

	@Override
	public float getTransparency()
	{
		return transparency;
	}

	@Override
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

	@Override
	protected TransparencyAttributes getTransparencyAttributes()
	{
		// Ignore this as transparency is handled by this class alone 
		return null;
	}
	
	@Override
	protected PolygonAttributes getPolygonAttributes()
	{
		return defaultAppearance.getPolygonAttributes();
	}

	@Override
	protected PointAttributes getPointAttributes()
	{
		return defaultAppearance.getPointAttributes();
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
					geometryArray[i].setColors(0, pointArrayColorUpdater.getColors(alphas[i]));
				}
			}
			else
			{
				for (int i = 0; i < geometryArray.length; i++)
				{
					geometryArray[i].setColors(0, pointArrayColorUpdater.pointColor);
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
					geometryArray[i].setColors(0, pointArrayColorUpdater.getColors(color[i], alphas[i]));
				}
			}
			else
			{
				pointArrayColorUpdater.getColors(1f);
				for (int i = 0; i < geometryArray.length; i++)
				{
					geometryArray[i].setColors(0, pointArrayColorUpdater.getColors(color[i]));
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
				geometryArray[i].setColors(0, pointArrayColorUpdater.getColors(color[i]));
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
		if (alpha != null)
		{
			int size = size();
			if (alpha.length != size)
				throw new IllegalArgumentException("list of size " + size + " expected");
			this.alphas = alpha.clone();
		}
		else
		{
			this.alphas = null;
		}
		if (isPointArray)
		{
			// PointArray alpha must be updated  
			final boolean hasAlpha = alpha != null;
			for (int i = 0; i < geometryArray.length; i++)
			{
				GeometryArray ga = geometryArray[i];
				ga.getColors(0, pointArrayColorUpdater.pointColor);
				pointArrayColorUpdater.getColors((hasAlpha) ? alpha[i] : 1f);
				ga.setColors(0, pointArrayColorUpdater.pointColor);
			}
		}
		else
		{
			// Global transparency is merged with alpha
			setTransparency(this.transparency);
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
			this.alphas = alpha.clone();
		}
		else
		{
			this.alphas = null;
		}
		if (isPointArray)
		{
			// PointArray alpha must be updated  
			final boolean hasAlpha = alpha != null;
			for (int i = 0; i < geometryArray.length; i++)
			{
				GeometryArray ga = geometryArray[i];
				ga.getColors(0, pointArrayColorUpdater.pointColor);
				pointArrayColorUpdater.getColors((hasAlpha) ? alpha[i] : 1f);
				ga.setColors(0, pointArrayColorUpdater.pointColor);
			}
		}
		setTransparency(transparency);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.TransparentItemShape#setItemAlpha(float)
	 */
	public void setItemAlpha(float alpha) throws IllegalArgumentException
	{
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
				ga.setColors(0, pointArrayColorUpdater.pointColor);
			}
		}
		else
		{
			// Global transparency is merged with alpha
			setTransparency(this.transparency);
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
				ga.setColors(0, pointArrayColorUpdater.pointColor);
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
