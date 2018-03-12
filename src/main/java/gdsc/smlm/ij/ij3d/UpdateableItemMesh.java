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

/**
 * Interface allowing mesh objects that represent a set of items to be updated
 */
public interface UpdateableItemMesh
{
	/**
	 * Reorder the mesh items using the given indices. The number of indices must match the number of items in the mesh
	 * and contain all indices from 0 to size-1.
	 *
	 * @param indices
	 *            the indices
	 * @throws IllegalArgumentException
	 *             if the indices are not valid
	 */
	public void reorder(int[] indices) throws IllegalArgumentException;

	/**
	 * Reorder the mesh items using the given indices. The output number of items will be the minimum of indices.length
	 * or the current mesh size. There is no checking whether the indices are valid.
	 *
	 * @param indices
	 *            the indices
	 * @throws IllegalArgumentException
	 *             if the indices are not valid
	 */
	public void reorderFast(int[] indices) throws IllegalArgumentException;
}
