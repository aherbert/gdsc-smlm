package gdsc.smlm.ij.ij3d;

import java.util.List;

import org.scijava.vecmath.Color4f;

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
 * Interface for mesh objects that represent a set of items to have transparent colours.
 */
public interface TransparentItemMesh extends ItemMesh
{
	/**
	 * Sets the color for each item.
	 *
	 * @param color
	 *            the new color
	 * @throws IllegalArgumentException
	 *             if the number of colours is incorrect
	 */
	public void setItemColor4(final List<Color4f> color) throws IllegalArgumentException;

	/**
	 * Sets the alpha for each item.
	 *
	 * @param alpha
	 *            the new alpha
	 * @throws IllegalArgumentException
	 *             if the number of colours is incorrect
	 */
	public void setItemAlpha(final float[] alpha) throws IllegalArgumentException;
}
