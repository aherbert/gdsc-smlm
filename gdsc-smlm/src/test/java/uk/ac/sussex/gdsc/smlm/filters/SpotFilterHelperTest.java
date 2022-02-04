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

package uk.ac.sussex.gdsc.smlm.filters;

import java.awt.Rectangle;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import uk.ac.sussex.gdsc.core.utils.rng.RandomUtils;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;

@SuppressWarnings({"javadoc"})
class SpotFilterHelperTest {
  private static Spot[] createData(UniformRandomProvider rg, int width, int height, int n) {
    if (n == 0) {
      return new Spot[0];
    }

    final int[] data = new int[width * height];
    for (int i = n; i-- > 0;) {
      data[i] = 1;
    }

    RandomUtils.shuffle(data, rg);

    final Spot[] spots = new Spot[n];
    for (int i = 0, j = 0; i < data.length; i++) {
      if (data[i] == 1) {
        spots[j++] = new Spot(i % width, i / width, 1);
        if (j == n) {
          break;
        }
      }
    }

    return spots;
  }

  @SeededTest
  void canCountNeighbours(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeed());

    final int width = 64;
    final int height = 64;
    final int size = width * height;
    int boxSize = 1; // Don't test simple case of no neighbours
    final SpotFilterHelper h = new SpotFilterHelper();
    while (boxSize < size / 4) {
      boxSize *= 2;
      for (int loop = 0; loop < 5; loop++) {
        final Spot[] spots = createData(rg, width, height, boxSize);
        for (final int box : new int[] {1, 2, 3, 4, 5}) {
          final int[] e = countNeighbours(spots, width, height, box);
          final int[] count = h.countNeighbours(spots, box);
          Assertions.assertArrayEquals(e, count);
          final int[] count2 = h.countNeighbours(spots, width, height, box);
          Assertions.assertArrayEquals(e, count2);
        }
      }
    }
  }

  private static int[] countNeighbours(Spot[] spots, int width, int height, int box) {
    final short[] data = new short[width * height];
    for (int i = 0; i < spots.length; i++) {
      data[spots[i].x + width * spots[i].y] = 1;
    }
    final Rectangle r = new Rectangle(width, height);
    final Rectangle bounds = new Rectangle(2 * box + 1, 2 * box + 1);
    final int[] count = new int[spots.length];
    for (int i = 0; i < spots.length; i++) {
      bounds.x = spots[i].x - box;
      bounds.y = spots[i].y - box;
      final Rectangle limits = r.intersection(bounds);
      int sum = -1;
      for (int y = limits.y; y < limits.y + limits.height; y++) {
        for (int x = limits.x, j = y * width + x; x < limits.x + limits.width; x++) {
          sum += data[j++];
        }
      }
      count[i] = sum;
    }
    return count;
  }
}
