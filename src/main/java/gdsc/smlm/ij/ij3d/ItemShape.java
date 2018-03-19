package gdsc.smlm.ij.ij3d;

import org.scijava.vecmath.Color3f;
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
 * Interface for shape objects that represent a set of items.
 */
public interface ItemShape
{
	/**
	 * Gets the number of items.
	 *
	 * @return the size
	 */
	public int size();

	/**
	 * Gets the coordinate of the specified item.
	 *
	 * @param i
	 *            the index
	 * @return the coordinate
	 */
	public Point3f getCoordinate(int i);

	/**
	 * Sets the color for each item.
	 *
	 * @param color
	 *            the new color
	 */
	public void setItemColor(final Color3f color);

	/**
	 * Sets the color for each item.
	 *
	 * @param color
	 *            the new color
	 * @throws IllegalArgumentException
	 *             if the number of colours is incorrect
	 */
	public void setItemColor(final Color3f[] color) throws IllegalArgumentException;
}
