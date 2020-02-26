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

package uk.ac.sussex.gdsc.smlm.ij.utils;

import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import java.util.Arrays;
import uk.ac.sussex.gdsc.core.annotation.Nullable;

/**
 * Find objects defined by contiguous pixels of the same value.
 */
public class ObjectAnalyzer {
  private static final int[] DIR_X_OFFSET = {0, 1, 0, -1, 1, 1, -1, -1};
  private static final int[] DIR_Y_OFFSET = {-1, 0, 1, 0, -1, 1, 1, -1};

  private final ImageProcessor ip;
  private boolean eightConnected;
  private int[] objectMask;
  private int maxObject;
  private int minObjectSize;

  /**
   * Instantiates a new object analyzer using 4N connected neighbours.
   *
   * @param ip the image
   */
  public ObjectAnalyzer(ImageProcessor ip) {
    this(ip, false);
  }

  /**
   * Instantiates a new object analyzer.
   *
   * @param ip the image
   * @param eightConnected Set to true to use 8N connected neighbours (the default is 4N)
   */
  public ObjectAnalyzer(ImageProcessor ip, boolean eightConnected) {
    this.ip = ip;
    this.eightConnected = eightConnected;
  }

  /**
   * Gets the object mask.
   *
   * @return A pixel array containing the object number for each pixel in the input image.
   */
  public int[] getObjectMask() {
    analyseObjects();
    return objectMask;
  }

  /**
   * Gets the object mask.
   *
   * @param objectNumber The object number
   * @return A byte mask of the object (objects pixels set to (byte)255)
   */
  public @Nullable byte[] getObjectMask(final int objectNumber) {
    if (objectNumber > getMaxObject()) {
      return null;
    }
    final byte[] pixels = new byte[objectMask.length];
    for (int i = 0; i < pixels.length; i++) {
      if (objectMask[i] == objectNumber) {
        pixels[i] = (byte) 255;
      }
    }
    return pixels;
  }

  /**
   * Gets the object processor.
   *
   * @param objectNumber The object number
   * @return A byte mask of the object
   */
  public ByteProcessor getObjectProcessor(final int objectNumber) {
    if (objectNumber > getMaxObject()) {
      return null;
    }
    return new ByteProcessor(ip.getWidth(), ip.getHeight(), getObjectMask(objectNumber));
  }

  /**
   * Gets the maximum object number.
   *
   * @return The maximum object number.
   */
  public int getMaxObject() {
    analyseObjects();
    return maxObject;
  }

  private void analyseObjects() {
    if (objectMask != null) {
      return;
    }

    final int[] maskImage = new int[ip.getPixelCount()];
    for (int i = 0; i < maskImage.length; i++) {
      maskImage[i] = ip.get(i);
    }

    // Perform a search for objects.
    // Expand any non-zero pixel value into all 8-connected pixels of the same value.
    objectMask = new int[maskImage.length];
    maxObject = 0;

    final int[][] ppList = new int[1][];
    ppList[0] = new int[100];
    initialise(ip);

    int[] sizes = new int[100];

    for (int i = 0; i < maskImage.length; i++) {
      // Look for non-zero values that are not already in an object
      if (maskImage[i] != 0 && objectMask[i] == 0) {
        maxObject++;
        final int size = expandObjectXy(maskImage, objectMask, i, maxObject, ppList);
        if (sizes.length == maxObject) {
          sizes = Arrays.copyOf(sizes, (int) (maxObject * 1.5));
        }
        sizes[maxObject] = size;
      }
    }

    // Remove objects that are too small
    if (minObjectSize > 0) {
      final int[] map = new int[maxObject + 1];
      maxObject = 0;
      for (int i = 1; i < map.length; i++) {
        if (sizes[i] >= minObjectSize) {
          map[i] = ++maxObject;
        }
      }

      for (int i = 0; i < objectMask.length; i++) {
        if (objectMask[i] != 0) {
          objectMask[i] = map[objectMask[i]];
        }
      }
    }
  }

  /**
   * Searches from the specified point to find all coordinates of the same value and assigns them to
   * given maximum ID.
   */
  private int expandObjectXy(final int[] image, final int[] objectMask, final int index0,
      final int id, int[][] ppList) {
    objectMask[index0] = id; // mark first point
    int listI = 0; // index of current search element in the list
    int listLen = 1; // number of elements in the list
    final int neighbours = (eightConnected) ? 8 : 4;

    // we create a list of connected points and start the list at the current point
    int[] pointList = ppList[0];
    pointList[listI] = index0;

    final int v0 = image[index0];

    do {
      final int index1 = pointList[listI];
      final int x1 = index1 % maxx;
      final int y1 = index1 / maxx;

      final boolean isInnerXy = (y1 != 0 && y1 != ylimit) && (x1 != 0 && x1 != xlimit);

      for (int d = neighbours; d-- > 0;) {
        if (isInnerXy || isWithinXy(x1, y1, d)) {
          final int index2 = index1 + offset[d];
          if (objectMask[index2] != 0) {
            // This has been done already, ignore this point
            continue;
          }

          final int v2 = image[index2];

          if (v2 == v0) {
            // Add this to the search
            pointList[listLen++] = index2;
            objectMask[index2] = id;
            if (pointList.length == listLen) {
              pointList = Arrays.copyOf(pointList, (int) (listLen * 1.5));
            }
          }
        }
      }

      listI++;

    } while (listI < listLen);

    ppList[0] = pointList;

    return listLen;
  }

  private int maxx;
  private int maxy;
  private int xlimit;
  private int ylimit;
  private int[] offset;

  /**
   * Creates the direction offset tables.
   */
  private void initialise(ImageProcessor ip) {
    maxx = ip.getWidth();
    maxy = ip.getHeight();

    xlimit = maxx - 1;
    ylimit = maxy - 1;

    // Create the offset table (for single array 3D neighbour comparisons)
    offset = new int[DIR_X_OFFSET.length];
    for (int d = offset.length; d-- > 0;) {
      offset[d] = maxx * DIR_Y_OFFSET[d] + DIR_X_OFFSET[d];
    }
  }

  /**
   * returns whether the neighbour in a given direction is within the image. NOTE: it is assumed
   * that the pixel x,y itself is within the image! Uses class variables xlimit, ylimit: (dimensions
   * of the image)-1
   *
   * @param x x-coordinate of the pixel that has a neighbour in the given direction
   * @param y y-coordinate of the pixel that has a neighbour in the given direction
   * @param direction the direction from the pixel towards the neighbour
   * @return true if the neighbour is within the image (provided that x, y is within)
   */
  private boolean isWithinXy(int x, int y, int direction) {
    switch (direction) {
      // 4-connected directions
      case 0:
        return (y > 0);
      case 1:
        return (x < xlimit);
      case 2:
        return (y < ylimit);
      case 3:
        return (x > 0);
      // Then remaining 8-connected directions
      case 4:
        return (y > 0 && x < xlimit);
      case 5:
        return (y < ylimit && x < xlimit);
      case 6:
        return (y < ylimit && x > 0);
      case 7:
        return (y > 0 && x > 0);
      default:
        return false;
    }
  }

  /**
   * Gets the image width.
   *
   * @return The image width.
   */
  public int getWidth() {
    return ip.getWidth();
  }

  /**
   * Gets the image height.
   *
   * @return The image height.
   */
  public int getHeight() {
    return ip.getHeight();
  }

  /**
   * Get the centre-of-mass and pixel count of each object. Data is stored indexed by the object
   * value so processing of results should start from 1.
   *
   * @return The centre-of-mass of each object (plus the pixel count) [object][cx,cy,n]
   */
  public double[][] getObjectCentres() {
    final int[] count = new int[maxObject + 1];
    final double[] sumx = new double[count.length];
    final double[] sumy = new double[count.length];
    final int maxy = getHeight();
    final int maxx = getWidth();
    for (int y = 0, i = 0; y < maxy; y++) {
      for (int x = 0; x < maxx; x++, i++) {
        final int value = objectMask[i];
        if (value != 0) {
          sumx[value] += x;
          sumy[value] += y;
          count[value]++;
        }
      }
    }
    final double[][] data = new double[count.length][3];
    for (int i = 1; i < count.length; i++) {
      data[i][0] = sumx[i] / count[i];
      data[i][1] = sumy[i] / count[i];
      data[i][2] = count[i];
    }
    return data;
  }

  /**
   * Gets the minimum object size. Objects below this are removed.
   *
   * @return The minimum object size
   */
  public int getMinObjectSize() {
    return minObjectSize;
  }

  /**
   * Sets the minimum object size. Objects below this are removed.
   *
   * @param minObjectSize The minimum object size
   */
  public void setMinObjectSize(int minObjectSize) {
    this.minObjectSize = minObjectSize;
  }

  /**
   * Checks if objects should use 8-connected pixels. The default is 4-connected.
   *
   * @return True if objects should use 8-connected pixels
   */
  public boolean isEightConnected() {
    return eightConnected;
  }

  /**
   * Sets to true if objects should use 8-connected pixels. The default is 4-connected.
   *
   * @param eightConnected True if objects should use 8-connected pixels
   */
  public void setEightConnected(boolean eightConnected) {
    this.eightConnected = eightConnected;
  }
}
