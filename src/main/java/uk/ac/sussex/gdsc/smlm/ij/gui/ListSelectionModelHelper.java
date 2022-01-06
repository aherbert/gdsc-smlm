/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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
    if (ArrayUtils.getLength(indices) == 0) {
      return;
    }

    sm.setValueIsAdjusting(true);
    sm.setSelectionInterval(indices[0], indices[0]);
    for (int i = 1; i < indices.length; i++) {
      sm.addSelectionInterval(indices[i], indices[i]);
    }
    sm.setValueIsAdjusting(false);
  }
}
