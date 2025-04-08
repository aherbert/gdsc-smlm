/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Package
 *
 * Software for single molecule localisation microscopy (SMLM) in ImageJ
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.ij.gui;

import java.util.Arrays;
import java.util.BitSet;
import javax.swing.ListSelectionModel;
import org.apache.commons.lang3.ArrayUtils;

/**
 * A helper class for the ListSelectionModel.
 */
public final class ListSelectionModelHelper {
  /**
   * No public constructor.
   */
  private ListSelectionModelHelper() {}

  /**
   * Gets the selected indices from the selection model.
   *
   * <p>Adapted from javax.swing.JList.
   *
   * @param sm the selection model
   * @return the selected indices
   */
  public static int[] getSelectedIndices(ListSelectionModel sm) {
    final int iMin = sm.getMinSelectionIndex();
    final int iMax = sm.getMaxSelectionIndex();

    // Any negative
    if ((iMin | iMax) < 0) {
      return new int[0];
    }

    final int[] rvTmp = new int[1 + (iMax - iMin)];
    int count = 0;
    for (int i = iMin; i <= iMax; i++) {
      if (sm.isSelectedIndex(i)) {
        rvTmp[count++] = i;
      }
    }
    return Arrays.copyOf(rvTmp, count);
  }

  /**
   * Sets the selected indices.
   *
   * @param sm the selection model
   * @param indices the indices
   */
  public static void setSelectedIndices(ListSelectionModel sm, int[] indices) {
    setSelectedIndices(sm, indices, false);
  }

  /**
   * Sets the selected indices.
   *
   * <p>To continue with further adjustments set the {@code adjusting} flag to true. To propagate
   * the selection change event the value must be set to false using
   * {@link ListSelectionModel#setValueIsAdjusting(boolean)}.
   *
   * @param sm the selection model
   * @param indices the indices
   * @param adjusting the value for {@link ListSelectionModel#setValueIsAdjusting(boolean)}
   */
  static void setSelectedIndices(ListSelectionModel sm, int[] indices, boolean adjusting) {
    if (ArrayUtils.getLength(indices) == 0) {
      return;
    }

    sm.setValueIsAdjusting(true);
    sm.setSelectionInterval(indices[0], indices[0]);
    for (int i = 1; i < indices.length; i++) {
      sm.addSelectionInterval(indices[i], indices[i]);
    }
    sm.setValueIsAdjusting(adjusting);
  }

  /**
   * Invert the selection.
   *
   * @param size the size of the model data
   * @param sm the selection model
   */
  static void invertSelection(int size, ListSelectionModel sm) {
    if (sm == null) {
      return;
    }
    final int iMin = sm.getMinSelectionIndex();
    final int iMax = sm.getMaxSelectionIndex();
    // Any negative
    if ((iMin | iMax) < 0) {
      sm.setSelectionInterval(0, size - 1);
    } else {
      final BitSet selected = new BitSet(size);
      for (int i = iMin; i <= iMax; i++) {
        if (sm.isSelectedIndex(i)) {
          selected.set(i, true);
        }
      }
      // Select the inverted range
      sm.setValueIsAdjusting(true);
      sm.clearSelection();
      int i = selected.nextClearBit(0);
      while (i < size) {
        final int j = selected.nextSetBit(i);
        if (j < 0) {
          sm.addSelectionInterval(i, size - 1);
          break;
        }
        sm.addSelectionInterval(i, j - 1);
        i = selected.nextClearBit(j);
      }
      sm.setValueIsAdjusting(false);
    }
  }
}
