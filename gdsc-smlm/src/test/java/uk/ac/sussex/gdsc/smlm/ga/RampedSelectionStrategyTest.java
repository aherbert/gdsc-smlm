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

package uk.ac.sussex.gdsc.smlm.ga;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

@SuppressWarnings({"javadoc"})
class RampedSelectionStrategyTest {

  @Test
  void canSearchUsingActualKey() {
    final long[] sum = RampedSelectionStrategy.createRampedSum(10);

    for (int i = 0; i < sum.length - 1; i++) {
      final long key = sum[i];
      final int j = RampedSelectionStrategy.search(sum, key);
      Assertions.assertEquals(i + 1, j);
    }
  }

  @Test
  void canBinarySearchUsingActualKey() {
    final long[] sum = RampedSelectionStrategy.createRampedSum(10);

    for (int i = 0; i < sum.length - 1; i++) {
      final long key = sum[i];
      final int j = RampedSelectionStrategy.binarySearch(sum, key);
      Assertions.assertEquals(i + 1, j);
    }
  }

  @Test
  void canSearchUsingNotActualKey() {
    final long[] sum = RampedSelectionStrategy.createRampedSum(10);

    for (int i = 0; i < sum.length; i++) {
      final long key = sum[i] - 1;
      final int j = RampedSelectionStrategy.search(sum, key);
      Assertions.assertEquals(i, j);
    }
  }

  @Test
  void canBinarySearchUsingNotActualKey() {
    final long[] sum = RampedSelectionStrategy.createRampedSum(10);

    for (int i = 0; i < sum.length; i++) {
      final long key = sum[i] - 1;
      final int j = RampedSelectionStrategy.binarySearch(sum, key);
      Assertions.assertEquals(i, j);
    }
  }

  @Test
  void binarySearchEqualsSearch() {
    final long[] sum = RampedSelectionStrategy.createRampedSum(100);
    for (int key = (int) sum[sum.length - 1]; key-- > 0;) {
      final int i = RampedSelectionStrategy.search(sum, key);
      final int j = RampedSelectionStrategy.binarySearch(sum, key);
      Assertions.assertEquals(i, j);
    }
  }
}
