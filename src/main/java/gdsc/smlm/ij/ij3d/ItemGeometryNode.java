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

import org.scijava.java3d.View;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Point3f;
import org.scijava.vecmath.Tuple3d;

import ij3d.ContentNode;

/**
 * Node to allow display of ItemGeometryGroup content.
 *
 * @author Alex Herbert
 */
public class ItemGeometryNode extends ContentNode
{
	private ItemGeometryGroup pointGroup;

	protected Point3f min, max, center;

	public ItemGeometryNode(final ItemGeometryGroup pointGroup)
	{
		this.pointGroup = pointGroup;
		calculateMinMaxCenterPoint();
		addChild(pointGroup);
	}

	public ItemGeometryGroup getItemGeometry()
	{
		return pointGroup;
	}

	@Override
	public void getMin(final Tuple3d min)
	{
		min.set(this.min);
	}

	@Override
	public void getMax(final Tuple3d max)
	{
		max.set(this.max);
	}

	@Override
	public void getCenter(final Tuple3d center)
	{
		center.set(this.center);
	}

	@Override
	public void channelsUpdated(final boolean[] channels)
	{
		// do nothing
	}

	@Override
	public void lutUpdated(final int[] r, final int[] g, final int[] b, final int[] a)
	{
		// do nothing
	}

	@Override
	public void colorUpdated(final Color3f color)
	{
		pointGroup.setColor(color);
	}

	@Override
	public void eyePtChanged(final View view)
	{
		// do nothing
	}

	@Override
	public float getVolume()
	{
		return 0;
	}

	@Override
	public void shadeUpdated(final boolean shaded)
	{
		pointGroup.setShaded(shaded);
	}

	@Override
	public void thresholdUpdated(final int threshold)
	{
		// do nothing
	}

	@Override
	public void transparencyUpdated(final float transparency)
	{
		pointGroup.setTransparency(transparency);
	}

	private void calculateMinMaxCenterPoint()
	{
		min = new Point3f();
		max = new Point3f();
		center = new Point3f();
		if (pointGroup.size() == 0)
			return;

		Point3f[] points = pointGroup.getPoints();
		CustomContentHelper.calculateMinMaxCenterPoint(min, max, center, points);
	}

	@Override
	public void restoreDisplayedData(final String path, final String name)
	{
		// do nothing
	}

	@Override
	public void swapDisplayedData(final String path, final String name)
	{
		// do nothing
	}

	@Override
	public void clearDisplayedData()
	{
		// do nothing
	}
}
