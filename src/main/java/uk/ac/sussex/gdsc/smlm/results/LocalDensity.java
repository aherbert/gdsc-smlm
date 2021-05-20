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

package uk.ac.sussex.gdsc.smlm.results;

import java.awt.Rectangle;
import java.awt.geom.AffineTransform;
import java.awt.geom.Area;
import java.awt.geom.PathIterator;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.MemoryUtils;
import uk.ac.sussex.gdsc.core.utils.ValidationUtils;
import uk.ac.sussex.gdsc.core.utils.function.IntDoubleConsumer;

/**
 * Contains functions to compute the local density of 2D coordinates.
 */
public final class LocalDensity {
  /** The border limit. */
  private static final int BORDER_LIMIT = (1 << 30) - 1;

  /**
   * A construct of 1 or more overlapping squares of the same size.
   *
   * <p>This stores squares using the centre. It has methods to check if another square overlaps any
   * square within the multi-square; to add squares; to add multi-squares; and compute the area.
   */
  private static class MultiSquare {
    private static final int[] EMPTY = {};

    /** The size of the square. */
    private final int size;
    /** The list of x points. */
    private int[] xp = EMPTY;
    /** The list of y points. */
    private int[] yp = EMPTY;
    /**
     * The number of points.
     *
     * <p>Warning: This is zero when the multi square has only one point. In this case the centre is
     * recovered from the bounds and the size.
     */
    private int np;

    // Intersection testing

    /** The minimum x position of the multi-square minus the size. */
    private long minx;
    /** The minimum y position of the multi-square minus the size. */
    private long miny;
    /** The maximum x position of the multi-square plus the size. */
    private long maxx;
    /** The maximum y position of the multi-square plus the size. */
    private long maxy;

    /**
     * Create an instance.
     *
     * @param x the x
     * @param y the y
     * @param size the size
     */
    MultiSquare(int x, int y, int size) {
      this.size = size;
      minx = (long) x - size;
      miny = (long) y - size;
      maxx = (long) x + size;
      maxy = (long) y + size;
    }

    /**
     * Adds the square centred at the given point.
     *
     * @param x the x
     * @param y the y
     */
    void add(int x, int y) {
      createList();
      xp[np] = x;
      yp[np] = y;
      np++;
      // Expand bounds
      final long s = size;
      if (x - s < minx) {
        minx = x - s;
      } else if (x + s > maxx) {
        maxx = x + s;
      }
      if (y - s < miny) {
        miny = y - s;
      } else if (y + s > maxy) {
        maxy = y + s;
      }
    }

    /**
     * Adds the multi-square. Assumes the border is the same.
     *
     * @param ms the multi-square
     */
    void add(MultiSquare ms) {
      createList();
      if (ms.np != 0) {
        // Add all the squares
        increaseCapacity(np + ms.np);
        System.arraycopy(ms.xp, 0, xp, np, ms.np);
        System.arraycopy(ms.yp, 0, yp, np, ms.np);
        np += ms.np;
        // Expand bounds
        minx = Math.min(minx, ms.minx);
        maxx = Math.max(maxx, ms.maxx);
        miny = Math.min(miny, ms.miny);
        maxy = Math.max(maxy, ms.maxy);
      } else {
        // Only one square
        add((int) (ms.minx + size), (int) (ms.miny + size));
      }
    }

    /**
     * Creates the list of squares using the initial bounds. This should be called before adding to
     * the multi-square.
     */
    private void createList() {
      if (xp.length == np) {
        // Special handling of empty points as we must add the starting point
        if (np == 0) {
          xp = new int[10];
          yp = new int[10];
          xp[0] = (int) (minx + size);
          yp[0] = (int) (miny + size);
          np = 1;
        } else {
          // Increase by 50%
          increaseCapacity((np * 3) >>> 1);
        }
      }
    }

    /**
     * Increase capacity.
     *
     * @param minCapacity the min capacity
     */
    private void increaseCapacity(int minCapacity) {
      final int size = MemoryUtils.createNewCapacity(minCapacity, np);
      xp = Arrays.copyOf(xp, size);
      yp = Arrays.copyOf(yp, size);
    }

    /**
     * Checks whether or not a square centred at the given point intersects with any square in the
     * multi-square.
     *
     * @param x the x
     * @param y the y
     * @return true if the square intersects any square in the multi-square
     */
    boolean intersects(final long x, final long y) {
      // Bounding region check
      if (x > minx && x < maxx && y > miny && y < maxy) {
        if (np == 0) {
          // Single square
          return true;
        }
        // Check each square
        final long distance = size;
        for (int i = np - 1; i >= 0; i--) {
          if (LocalDensity.withinReach(xp[i], yp[i], x, y, distance)) {
            return true;
          }
        }
      }
      return false;
    }

    /**
     * Compute the area of the multi-square.
     *
     * @return the area
     */
    double area() {
      if (np == 0) {
        return (double) size * size;
      }
      if (np == 2) {
        return LocalDensity.getArea(xp[0], yp[0], xp[1], yp[1], size);
      }

      // Combine all squares
      final Area area = new Area(new Rectangle(xp[0], yp[0], size, size));
      for (int i = 1; i < np; i++) {
        area.add(new Area(new Rectangle(xp[i], yp[i], size, size)));
      }

      // Compute the area
      double sum = 0;
      // Trust that the geometry of the path iterator will not return more vertices
      // than the input squares.
      // Note: This allocation has the potential to overflow. Clip it and
      // leave the code to array index out of bounds if the path is that large.
      // In this case we cannot compute the area anyway.
      final double[] x = new double[MemoryUtils.createNewCapacity(np * 4, 0)];
      final double[] y = new double[x.length];
      int n = 0;
      final double[] coords = new double[6];
      final PathIterator pIter = area.getPathIterator(new AffineTransform());
      while (!pIter.isDone()) {
        final int segType = pIter.currentSegment(coords);
        // We are only interested in move-to and line-to operations
        switch (segType) {
          case PathIterator.SEG_MOVETO:
            // The start of a new path
            sum += getArea(x, y, n);
            n = 0;
            // CHECKSTYLE.OFF: FallThroughCheck
          case PathIterator.SEG_LINETO:
            // Add the point
            x[n] = coords[0];
            y[n] = coords[1];
            n++;
            // CHECKSTYLE.ON: FallThroughCheck
            break;
          default:
            // assume it is PathIterator.SEG_CLOSE
            break;
        }
        pIter.next();
      }
      return sum + getArea(x, y, n);
    }

    /**
     * Gets the area of a polygon using the Shoelace formula
     * (https://en.wikipedia.org/wiki/Shoelace_formula).
     *
     * <p>The area formula is valid for any non-self-intersecting (simple) polygon, which can be
     * convex or concave.
     *
     * @param x the x
     * @param y the y
     * @param n the count of points
     * @return the area
     */
    private static double getArea(double[] x, double[] y, int n) {
      double sum1 = 0;
      double sum2 = 0;
      for (int i = n, j = 0; i-- > 0; j = i) {
        sum1 += x[i] * y[j];
        sum2 += x[j] * y[i];
      }
      return Math.abs(sum1 - sum2) / 2;
    }

    /**
     * Get the size. This is the number of squares in the multi-square.
     *
     * @return the size
     */
    int size() {
      // np is zero when there is only one square
      return np == 0 ? 1 : np;
    }
  }

  /**
   * No public constructor.
   */
  private LocalDensity() {}

  /**
   * Estimate the local density of points. The border should be the extent around the point
   * considered to be close to the point, for example the extent of the point spread function (PSF).
   *
   * <p>This function calls {@link #estimate(int[], int[], int, IntDoubleConsumer) estimate} with a
   * {@code null} consumer for the regions.
   *
   * @param x the x positions
   * @param y the y positions
   * @param border the border
   * @return the local density
   * @see #estimate(int[], int[], int, IntDoubleConsumer)
   * @throw {@link IllegalArgumentException} if x and y lengths do not match; or if
   *        {@code 2 * border + 1} is greater than the maximum integer size.
   */
  public static double estimate(int[] x, int[] y, int border) {
    return estimate(x, y, border, null);
  }

  /**
   * Estimate the local density of points. The border should be the extent around the point
   * considered to be close to the point, for example the extent of the point spread function (PSF).
   *
   * <p>A bounding box is placed around each XY point defining the local region of interest. Any
   * other point's bounding box that overlap the bounding box is considered to be local. The
   * overlapping bounding box regions are combined to an area and the density computed as:
   *
   * <pre>
   * density = points / (sum local areas)
   * </pre>
   *
   * <p>The count of points and area of each local region can be obtained using the {@code regions}
   * consumer.
   *
   * <p>Any points that do not interact with other points have a minimum density of
   * {@code Math.pow(2 * border + 1, -2)}.
   *
   * <p>Note: This estimate requires that each square can be represented using the {@link Rectangle}
   * class. Thus the width and height of each square is limited to the maximum integer size.
   *
   * @param x the x positions
   * @param y the y positions
   * @param border the border (negative values are set to zero)
   * @param regions the regions (can be {@code null})
   * @return the local density
   * @throw {@link IllegalArgumentException} if x and y lengths do not match; or if
   *        {@code 2 * border + 1} is greater than the maximum integer size.
   */
  public static double estimate(int[] x, int[] y, int border, IntDoubleConsumer regions) {
    final int n = x.length;
    ValidationUtils.checkArgument(n == y.length, "xy length mismatch: %d != %d", n, y.length);
    // Check the border is positive and within the limit
    if (border < 0) {
      border = 0;
    } else {
      ValidationUtils.checkArgument(border <= BORDER_LIMIT, "border too large: %d", border);
    }
    final int size = 2 * border + 1;
    if (n == 0) {
      return 0;
    }
    if (n == 1) {
      final double area = (double) size * size;
      if (regions != null) {
        regions.accept(1, area);
      }
      return 1.0 / area;
    }
    // Edge case for 2 points
    if (n == 2) {
      // The points interact if they are within the size for x and y.
      if (withinReach(x[0], y[0], x[1], y[1], size)) {
        final double area = getArea(x[0], y[0], x[1], y[1], size);
        if (regions != null) {
          regions.accept(2, area);
        }
        return 2.0 / area;
      }
      // Minimum density
      final double area = (double) size * size;
      if (regions != null) {
        regions.accept(1, area);
        regions.accept(1, area);
      }
      return 1.0 / area;
    }

    // Construct bounding boxes and test for intersections.
    // Join intersecting squares to shapes.
    final LinkedList<MultiSquare> areas = new LinkedList<>();
    areas.add(new MultiSquare(x[0], y[0], size));
    for (int i = 1; i < n; i++) {
      // Use long as the intersects method uses long
      final long xx = x[i];
      final long yy = y[i];
      MultiSquare parent = null;
      final Iterator<MultiSquare> iter = areas.iterator();
      while (iter.hasNext()) {
        final MultiSquare mr = iter.next();
        if (mr.intersects(xx, yy)) {
          if (parent == null) {
            // First intersection, add the square to the multi-square
            parent = mr;
            parent.add(x[i], y[i]);
          } else {
            // Another intersection. This square bridges between to other multi-squares.
            // Combine the multi-squares and remove one from the list.
            parent.add(mr);
            iter.remove();
          }
        }
      }
      // Not found so start a new multi-square.
      if (parent == null) {
        areas.add(new MultiSquare(x[i], y[i], size));
      }
    }

    // Compute the density of interacting regions
    int sumCount = 0;
    double sumArea = 0;
    for (final MultiSquare mr : areas) {
      final int count = mr.size();
      final double area = mr.area();
      sumCount += count;
      sumArea += area;
      if (regions != null) {
        regions.accept(count, area);
      }
    }

    return MathUtils.div0(sumCount, sumArea);
  }

  /**
   * Return true if the points are within the reach distance.
   *
   * <pre>
   * {@code |x1-x2| < distance && |y1-y2| < distance}
   * </pre>
   *
   * @param x1 the first point x
   * @param y1 the first point y
   * @param x2 the second point x
   * @param y2 the second point y
   * @param distance the reach distance
   * @return true if within reach
   */
  static boolean withinReach(long x1, long y1, long x2, long y2, long distance) {
    return diff(x1, x2) < distance && diff(y1, y2) < distance;
  }

  /**
   * Compute the absolute difference between the two values.
   *
   * @param x value x
   * @param y value y
   * @return the difference |x-y|
   */
  private static long diff(long x, long y) {
    return x > y ? x - y : y - x;
  }

  /**
   * Gets the area of the overlapping squares.
   *
   * @param x1 the first point x
   * @param y1 the first point y
   * @param x2 the second point x
   * @param y2 the second point y
   * @param size the size
   * @return the area
   */
  static long getArea(int x1, int y1, int x2, int y2, final int size) {
    // Compute the union area.
    final Rectangle r0 = new Rectangle(x1, y1, size, size);
    final Rectangle r1 = new Rectangle(x2, y2, size, size);
    final Rectangle intersection = r0.intersection(r1);
    return 2L * size * size - (long) intersection.width * intersection.height;
  }
}
