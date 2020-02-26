/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */

package uk.ac.sussex.gdsc.smlm.ij.ij3d;

/**
 * Interface allowing shape objects that represent a set of items to be updated.
 */
public interface UpdateableItemShape extends ItemShape {
  /**
   * Reorder the mesh items using the given indices. The number of indices must match the number of
   * items in the mesh and contain all indices from 0 to size-1.
   *
   * @param indices the indices
   * @throws IllegalArgumentException if the indices are not valid
   */
  void reorder(int[] indices);

  /**
   * Reorder the mesh items using the given indices. The output number of items will be the minimum
   * of indices.length or the current mesh size. There is no checking whether the indices are valid.
   *
   * @param indices the indices
   * @throws IllegalArgumentException if the indices are not valid
   */
  void reorderFast(int[] indices);
}
