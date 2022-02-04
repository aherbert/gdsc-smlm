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

package uk.ac.sussex.gdsc.smlm.results;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import java.awt.Rectangle;
import java.util.BitSet;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.function.IntDoubleConsumer;

@SuppressWarnings({"javadoc"})
class LocalDensityTest {
  @Test
  void estimateThrowsWithDimensionMismatch() {
    final int[] x = new int[2];
    final int[] y = new int[3];
    final int border = 1;
    Assertions.assertThrows(IllegalArgumentException.class,
        () -> LocalDensity.estimate(x, y, border));
  }

  @Test
  void estimateThrowsWithBorderOverflow() {
    final int[] x = new int[2];
    Assertions.assertThrows(IllegalArgumentException.class,
        () -> LocalDensity.estimate(x, x, Integer.MAX_VALUE));
    Assertions.assertThrows(IllegalArgumentException.class,
        () -> LocalDensity.estimate(x, x, 1 << 30));
    // OK
    final int border = (1 << 30) - 1;
    Assertions.assertEquals(Integer.MAX_VALUE, 2 * border + 1);
    LocalDensity.estimate(x, x, border);
  }

  @Test
  void estimateReturnsZeroForZeroPoints() {
    final int[] x = {};
    final int border = 1;
    final RegionConsumer consumer = new RegionConsumer();
    Assertions.assertEquals(0, LocalDensity.estimate(x, x, border));
    Assertions.assertEquals(0, LocalDensity.estimate(x, x, border, consumer::accept));
    Assertions.assertEquals(0, consumer.size());
  }

  @Test
  void estimateReturnsMinDensityForOnePoint() {
    final int[] x = {0};
    final RegionConsumer consumer = new RegionConsumer();
    final int[] counts = {1};
    for (int border = 0; border <= 3; border++) {
      final int size = 2 * border + 1;
      final double expected = 1.0 / (size * size);
      final double[] areas = {size * size};
      Assertions.assertEquals(expected, LocalDensity.estimate(x, x, border));
      Assertions.assertEquals(expected, LocalDensity.estimate(x, x, border, consumer::accept));
      Assertions.assertTrue(consumer.equals(counts, areas));
      consumer.clear();
    }
  }

  @Test
  void estimateReturnsMinDensityForNoPossibleInteraction() {
    final RegionConsumer consumer = new RegionConsumer();
    final int[] counts = {1, 1};
    for (int border = 0; border <= 2; border++) {
      final int size = 2 * border + 1;
      final double density = 1.0 / (size * size);
      final double[] areas = {size * size, size * size};
      Assertions.assertEquals(density,
          LocalDensity.estimate(new int[] {0, 5}, new int[] {0, 5}, border));
      Assertions.assertEquals(density,
          LocalDensity.estimate(new int[] {0, 5}, new int[] {0, 0}, border));
      Assertions.assertEquals(density,
          LocalDensity.estimate(new int[] {0, 0}, new int[] {0, 5}, border));
      Assertions.assertEquals(density,
          LocalDensity.estimate(new int[] {0, 5}, new int[] {0, 5}, border, consumer::accept));
      Assertions.assertTrue(consumer.equals(counts, areas));
      consumer.clear();
      Assertions.assertEquals(density,
          LocalDensity.estimate(new int[] {0, 5}, new int[] {0, 0}, border, consumer::accept));
      Assertions.assertTrue(consumer.equals(counts, areas));
      consumer.clear();
      Assertions.assertEquals(density,
          LocalDensity.estimate(new int[] {0, 0}, new int[] {0, 5}, border, consumer::accept));
      Assertions.assertTrue(consumer.equals(counts, areas));
      consumer.clear();
    }
  }

  @Test
  void estimateDensityForTwoPoints() {
    final int[] xs = {-1, 1};
    final int[] borders = {-1, 0, 1, 2};
    for (final int x : xs) {
      for (final int y : xs) {
        for (final int border : borders) {
          assertDensityForTwoPoints(0, 0, x, y, border);
        }
      }
    }

    // Overflow in difference between the coordinates
    assertDensityForTwoPoints(-1, 0, Integer.MAX_VALUE, 0, 5);
    assertDensityForTwoPoints(0, -1, 0, Integer.MAX_VALUE, 5);

    // Maximum square size
    assertDensityForTwoPoints(0, 0, 0, 0, (1 << 30) - 1);
  }

  private static void assertDensityForTwoPoints(int x0, int y0, int x1, int y1, int border) {
    final int size = 2 * Math.max(0, border) + 1;
    final Rectangle r0 = new Rectangle(x0, y0, size, size);
    final Rectangle r1 = new Rectangle(x1, y1, size, size);
    final Rectangle intersection = r0.intersection(r1);
    final double minArea = (double) size * size;
    int[] counts;
    double[] areas;
    if (!intersection.isEmpty()) {
      final double area = 2.0 * size * size - (double) intersection.width * intersection.height;
      counts = new int[] {2};
      areas = new double[] {area};
    } else {
      counts = new int[] {1, 1};
      areas = new double[] {minArea, minArea};
    }
    final double expected = 2 / MathUtils.sum(areas);
    Assertions.assertEquals(expected,
        LocalDensity.estimate(new int[] {x0, x1}, new int[] {y0, y1}, border));
    final RegionConsumer consumer = new RegionConsumer();
    Assertions.assertEquals(expected,
        LocalDensity.estimate(new int[] {x0, x1}, new int[] {y0, y1}, border, consumer::accept));
    Assertions.assertTrue(consumer.equals(counts, areas));
  }

  @Test
  void estimateDensityForThreePoints() {
    // No overlap
    assertDensityForThreePoints(0, 0, 5, 5, 10, 10, 3);
    // r0+r1;r2
    assertDensityForThreePoints(0, 0, 1, 2, 10, 10, 3);
    // r0+r1+r2
    assertDensityForThreePoints(0, 0, 1, 2, 3, 4, 3);
    // r0+r2+r1
    assertDensityForThreePoints(0, 0, 3, 4, 1, 2, 2);
    // r0+r2;r1
    assertDensityForThreePoints(11, 12, 0, 0, 13, 14, 3);
    // r1+r2;r0
    assertDensityForThreePoints(0, 0, 11, 12, 13, 14, 3);

    // Edge case: r2 intersects with the bounding rectangle of r0+r1, but not either rectangle
    // r0+r1;r2
    assertDensityForThreePoints(0, 0, 2, 2, -1, 4, 1);
  }

  private static void assertDensityForThreePoints(int x0, int y0, int x1, int y1, int x2, int y2,
      int border) {
    final int size = 2 * Math.max(0, border) + 1;
    // Determine the bounding box
    final Rectangle r0 = new Rectangle(x0, y0, size, size);
    final Rectangle r1 = new Rectangle(x1, y1, size, size);
    final Rectangle r2 = new Rectangle(x2, y2, size, size);
    final boolean r0r1 = r0.intersects(r1);
    final boolean r0r2 = r0.intersects(r2);
    final boolean r1r2 = r1.intersects(r2);
    // Add overlapping areas
    final double minArea = (double) size * size;
    int[] counts;
    double[] areas;
    if (r0r1) {
      if (r0r2 || r1r2) {
        counts = new int[] {3};
        areas = new double[] {area(r0, r1, r2)};
      } else {
        counts = new int[] {2, 1};
        areas = new double[] {area(r0, r1), minArea};
      }
    } else if (r0r2) {
      if (r1r2) {
        counts = new int[] {3};
        areas = new double[] {area(r0, r1, r2)};
      } else {
        counts = new int[] {2, 1};
        areas = new double[] {area(r0, r2), minArea};
      }
    } else if (r1r2) {
      counts = new int[] {2, 1};
      areas = new double[] {area(r1, r2), minArea};
    } else {
      counts = new int[] {1, 1, 1};
      areas = new double[] {minArea, minArea, minArea};
    }
    final double expected = 3 / MathUtils.sum(areas);
    Assertions.assertEquals(expected,
        LocalDensity.estimate(new int[] {x0, x1, x2}, new int[] {y0, y1, y2}, border));
    final RegionConsumer consumer = new RegionConsumer();
    Assertions.assertEquals(expected, LocalDensity.estimate(new int[] {x0, x1, x2},
        new int[] {y0, y1, y2}, border, consumer::accept));
    Assertions.assertTrue(consumer.equals(counts, areas));
  }

  private static double area(Rectangle... rectangles) {
    Rectangle union = rectangles[0].union(rectangles[1]);
    for (int i = 2; i < rectangles.length; i++) {
      union = union.union(rectangles[i]);
    }
    // Fill the area
    final int w = union.width;
    final int h = union.height;
    final int[] area = new int[w * h];
    for (final Rectangle r : rectangles) {
      final int ox = r.x - union.x;
      final int oy = r.y - union.y;
      final int maxx = ox + r.width;
      final int maxy = oy + r.height;
      for (int y = oy; y < maxy; y++) {
        for (int x = ox, index = y * w + ox; x < maxx; x++, index++) {
          area[index] = 1;
        }
      }
    }
    return MathUtils.sum(area);
  }

  @Test
  void estimateDensityForFivePoints() {
    // @formatter:off
    // ### ###     ### ###
    // #a###b#     #a# #b#      ###
    // ###e###  =  ### ###  +   #e#   +       +   ###
    //   #g#                    ###               #g#
    // #######     ### ###                        ###
    // #c###d#     #c# #d#               ###
    // ###f###     ### ###               #f#
    //   ###                             ###
    Assertions.assertEquals(7.0 / (7 * 7 - 1 - 2 - 2 + 3),
                                      // a  b  c  d  e  f  g
        LocalDensity.estimate(new int[] {0, 4, 0, 4, 2, 2, 2},
                              new int[] {0, 0, 4, 4, 1, 5, 2}, 1));
    // @formatter:on
  }

  @Test
  void estimateDensityForTwentyPoints() {
    // This tests increasing the capacity of the multi-square and adding
    // squares below the minimum bounds
    final int[] x = SimpleArrayUtils.newArray(20, 0, -2);
    Assertions.assertEquals(20.0 / (9 + 8 * 19), LocalDensity.estimate(x, x, 1));
  }

  /**
   * Simple consumer of the regions.
   */
  private static class RegionConsumer implements IntDoubleConsumer {
    TIntArrayList list1 = new TIntArrayList();
    TDoubleArrayList list2 = new TDoubleArrayList();

    @Override
    public void accept(int t, double value) {
      list1.add(t);
      list2.add(value);
    }

    int size() {
      return list1.size();
    }

    boolean equals(int[] ts, double[] values) {
      if (ts.length == list1.size()) {
        if (ts.length == 0) {
          return true;
        }
        final BitSet set = new BitSet(ts.length);
        for (int n = 0; n < ts.length; n++) {
          final int t = ts[n];
          final double value = values[n];
          for (int i = list1.size() - 1; i >= 0; i--) {
            if (!set.get(i) && list1.getQuick(i) == t && list2.getQuick(i) == value) {
              // Found, mark as used
              set.set(i);
              break;
            }
          }
        }
        return set.cardinality() == ts.length;
      }
      return false;
    }

    void clear() {
      list1.resetQuick();
      list2.resetQuick();
    }
  }
}
