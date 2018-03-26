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

import org.scijava.java3d.Group;
import org.scijava.java3d.PointAttributes;
import org.scijava.java3d.PolygonAttributes;
import org.scijava.java3d.TransparencyAttributes;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Color4f;
import org.scijava.vecmath.Point3f;

/**
 * This class represents a list of repeated shapes in the universe that can be added to a scene. It provides basic
 * functionality for displayed Content.
 *
 * @author Alex Herbert
 */
public abstract class ItemGroup extends Group implements ItemShape
{
	protected final Color3f DEFAULT_COLOUR = new Color3f(0, 1, 0);

	public static abstract class PointArrayColorUpdater
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

	public static class SinglePointArrayColorUpdater extends PointArrayColorUpdater
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

	public static class MultiPointArrayColorUpdater extends PointArrayColorUpdater
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
	
	/** The list of points */
	protected Point3f[] points;

	/**
	 * Instantiates a new item geometry group.
	 *
	 * @param points
	 *            the points
	 */
	public ItemGroup(final Point3f[] points)
	{
		if (points == null)
			throw new NullPointerException("Points must not be null");
		this.points = points;
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
	 * Gets the color of the points.
	 *
	 * @return the color
	 */
	public abstract Color3f getColor();

	/**
	 * Set the transparency of the points. This is blended with a per-item alpha if present.
	 *
	 * @param transparency
	 *            the new transparency
	 */
	public void setTransparency(final float transparency)
	{
		TransparencyAttributes ta = getTransparencyAttributes();
		if (ta == null)
			return;
		if (transparency == 0f)
		{
			ta.setTransparencyMode(TransparencyAttributes.NONE);
			ta.setTransparency(transparency);
		}
		else
		{
			ta.setTransparency(transparency);
			ta.setTransparencyMode(TransparencyAttributes.FASTEST);
		}
	}

	/**
	 * Gets the global transparency attributes.
	 *
	 * @return the transparency attributes
	 */
	protected abstract TransparencyAttributes getTransparencyAttributes();

	/**
	 * Gets the transparency of the points. This is blended with a per-item alpha if present.
	 *
	 * @return the transparency
	 */
	public float getTransparency()
	{
		TransparencyAttributes ta = getTransparencyAttributes();
		if (ta == null)
			return 0f;
		return ta.getTransparency();
	}

	/**
	 * Checks if any item is transparent.
	 *
	 * @return true, if is transparent
	 */
	public boolean isTransparent()
	{
		TransparencyAttributes ta = getTransparencyAttributes();
		if (ta == null)
			return false;
		return ta.getTransparency() != 0f && ta.getTransparencyMode() != TransparencyAttributes.NONE;
	}

	/**
	 * Sets the shaded.
	 *
	 * @param shaded
	 *            the new shaded
	 */
	public void setShaded(boolean shaded)
	{
		PolygonAttributes pa = getPolygonAttributes();
		if (pa == null)
			return;
		int mode = (shaded) ? PolygonAttributes.POLYGON_FILL : PolygonAttributes.POLYGON_LINE;
		if (pa.getPolygonMode() != mode)
			pa.setPolygonMode(mode);
	}

	/**
	 * Gets the global polygon attributes.
	 *
	 * @return the polygon attributes
	 */
	protected abstract PolygonAttributes getPolygonAttributes();

	/**
	 * Checks if is shaded.
	 *
	 * @return true, if is shaded
	 */
	public boolean isShaded()
	{
		PolygonAttributes pa = getPolygonAttributes();
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
		PointAttributes pa = getPointAttributes();
		if (pa == null)
			return;
		pa.setPointSize(pointSize);
	}

	/**
	 * Gets the global point attributes.
	 *
	 * @return the point attributes
	 */
	protected abstract PointAttributes getPointAttributes();

	/**
	 * Gets the point size.
	 *
	 * @return the point size
	 */
	public float getPointSize()
	{
		PointAttributes pa = getPointAttributes();
		if (pa == null)
			return 0f;
		return pa.getPointSize();
	}

	/**
	 * Returns a string representation of the underlying List<Point3f>.
	 */
	@Override
	public String toString()
	{
		return Arrays.toString(points);
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
	 * @see gdsc.smlm.ij.ij3d.ItemShape#getCoordinate(int)
	 */
	public Point3f getCoordinate(int i)
	{
		return points[i];
	}
}
