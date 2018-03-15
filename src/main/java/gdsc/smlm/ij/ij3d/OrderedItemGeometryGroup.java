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
import org.scijava.java3d.GeometryArray;
import org.scijava.java3d.Group;
import org.scijava.java3d.OrderedGroup;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Point3f;

/**
 * This class represents a list as a number of repeated shapes in the universe. The shape is defined using a geometry
 * array.
 * <p>
 * The order is fixed but can be user defined.
 *
 * @author Alex Herbert
 */
public class OrderedItemGeometryGroup extends ItemGeometryGroup implements UpdateableItemShape
{
	protected OrderedGroup orderedGroup;

	/**
	 * Instantiates a new ordered item geometry group.
	 *
	 * @param points
	 *            the points
	 */
	public OrderedItemGeometryGroup(final Point3f[] points)
	{
		super(points, null, null, null, null, null);
	}

	/**
	 * Instantiates a new ordered item geometry group. The geometry is scaled and translated for each point. The default
	 * appearance is cloned per item and optionally can be updated with per-item material colour and transparency. The
	 * transparency in the appearance is blended with per-item alphas to create a per-item transparency.
	 *
	 * @param points
	 *            the points.
	 * @param ga
	 *            the geometry array. If null then a default will be used.
	 * @param appearance
	 *            the default appearance of the shape. PolygonAttributes, Material and TransparencyAttributes are used.
	 * @param sizes
	 *            the sizes of each point. Can be null (no scaling); length=1 (fixed scaling); or points.length.
	 * @param colors
	 *            the per-item colors. Can be null.
	 * @param alphas
	 *            the per-item alphas. Can be null.
	 */
	public OrderedItemGeometryGroup(final Point3f[] points, GeometryArray ga, Appearance appearance, Point3f[] sizes,
			Color3f[] colors, float[] alphas)
	{
		super(points, ga, appearance, sizes, colors, alphas);
	}

	/**
	 * Gets the parent group to which all the shapes should be added.
	 *
	 * @return the parent group
	 */
	protected Group getParentGroup()
	{
		// Force the results to be ordered
		orderedGroup = new OrderedGroup();
		//orderedGroup.setCapability(OrderedGroup.ALLOW_CHILD_INDEX_ORDER_READ); // On by default
		orderedGroup.setCapability(OrderedGroup.ALLOW_CHILD_INDEX_ORDER_WRITE);
		addChild(orderedGroup);
		return orderedGroup;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.UpdateableItemShape#reorder(int[])
	 */
	public void reorder(int[] indices) throws IllegalArgumentException
	{
		reorderFast(indices);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.UpdateableItemShape#reorderFast(int[])
	 */
	public void reorderFast(int[] indices) throws IllegalArgumentException
	{
		orderedGroup.setChildIndexOrder(indices);
	}
}
