/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.model;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

/**
 * Samples uniformly from the specified mask. All non-zero pixels are sampled. The centre of the
 * mask corresponds to XY=0. The z coordinate is randomly sampled from the image depth and centred
 * on zero.
 *
 * <p>X coordinates are returned in the interval -width/2 to width/2. These can be converted to
 * different values using the scale parameter. Likewise for the Y coordinates. E.g. a mask of
 * 100x100 (range of -50:50) can be used to generate coordinates in the range -100:100 using a scale
 * of 2.
 *
 * <p>Sub pixel locations and z-depth are sampled from a uniform distribution. A Halton sequence is
 * used by default but this can be changed by setting a custom uniform distribution.
 */
public class MaskDistribution implements SpatialDistribution {
  // Used for a particle search
  private static final int[] DIR_X_OFFSET = {0, 1, 1, 1, 0, -1, -1, -1};
  private static final int[] DIR_Y_OFFSET = {-1, -1, 0, 1, 1, 1, 0, -1};

  private final RandomGenerator randomGenerator;
  private UniformDistribution uniformDistribution;
  private final int[] mask;
  private final int[] indices;
  private final int maxx;
  private final int maxy;
  private final double halfWidth;
  private final double halfHeight;
  private final double minDepth;
  private final double depth;
  private int particle;
  private final double scaleX;
  private final double scaleY;

  // Used for a particle search
  private int xlimit = -1;
  private int ylimit;
  private int[] offset;

  /**
   * Create a distribution from the mask image (packed in YX order).
   *
   * @param mask the mask
   * @param width The width of the mask in pixels
   * @param height the height of the mask in pixels
   * @param depth The mask depth
   * @param scaleX Used to scale the mask X-coordinate to a new value
   * @param scaleY Used to scale the mask Y-coordinate to a new value
   */
  public MaskDistribution(byte[] mask, int width, int height, double depth, double scaleX,
      double scaleY) {
    this(mask, width, height, depth, scaleX, scaleY, null);
  }

  /**
   * Create a distribution from the mask image (packed in YX order).
   *
   * @param mask the mask
   * @param width The width of the mask in pixels
   * @param height the height of the mask in pixels
   * @param depth The mask depth
   * @param scaleX Used to scale the mask X-coordinate to a new value
   * @param scaleY Used to scale the mask Y-coordinate to a new value
   */
  public MaskDistribution(int[] mask, int width, int height, double depth, double scaleX,
      double scaleY) {
    this(mask, width, height, depth, scaleX, scaleY, null);
  }

  /**
   * Create a distribution from the mask image (packed in YX order).
   *
   * @param mask the mask
   * @param width The width of the mask in pixels
   * @param height the height of the mask in pixels
   * @param depth The mask depth
   * @param scaleX Used to scale the mask X-coordinate to a new value
   * @param scaleY Used to scale the mask Y-coordinate to a new value
   * @param randomGenerator Used to pick random pixels in the mask
   */
  public MaskDistribution(byte[] mask, int width, int height, double depth, double scaleX,
      double scaleY, RandomGenerator randomGenerator) {
    this(convert(mask), width, height, depth, scaleX, scaleY, randomGenerator);
  }

  private static int[] convert(byte[] mask) {
    final int[] newMask = new int[mask.length];
    for (int i = 0; i < mask.length; i++) {
      if (mask[i] != 0) {
        newMask[i] = 1;
      }
    }
    return newMask;
  }

  /**
   * Create a distribution from the mask image (packed in YX order).
   *
   * @param mask the mask
   * @param width The width of the mask in pixels
   * @param height the height of the mask in pixels
   * @param depth The mask depth
   * @param scaleX Used to scale the mask X-coordinate to a new value
   * @param scaleY Used to scale the mask Y-coordinate to a new value
   * @param randomGenerator Used to pick random pixels in the mask
   */
  public MaskDistribution(int[] mask, int width, int height, double depth, double scaleX,
      double scaleY, RandomGenerator randomGenerator) {
    this(mask, width, height, depth, scaleX, scaleY, randomGenerator, null);
  }

  /**
   * Create a distribution from the mask image (packed in YX order).
   *
   * @param mask the mask
   * @param width The width of the mask in pixels
   * @param height the height of the mask in pixels
   * @param depth The mask depth
   * @param scaleX Used to scale the mask X-coordinate to a new value
   * @param scaleY Used to scale the mask Y-coordinate to a new value
   * @param randomGenerator Used to pick random pixels in the mask
   * @param uniformDistribution Used for sub-pixel location and z-depth
   */
  public MaskDistribution(int[] mask, int width, int height, double depth, double scaleX,
      double scaleY, RandomGenerator randomGenerator, UniformDistribution uniformDistribution) {
    this(mask, width, height, depth, scaleX, scaleY, randomGenerator, uniformDistribution, false);
  }

  /**
   * Create a distribution from the mask image (packed in YX order)
   *
   * <p>This is a package scope constructor allowing the mask to be created with all zero pixels.
   *
   * @param mask the mask
   * @param width The width of the mask in pixels
   * @param height the height of the mask in pixels
   * @param depth The mask depth
   * @param scaleX Used to scale the mask X-coordinate to a new value
   * @param scaleY Used to scale the mask Y-coordinate to a new value
   * @param randomGenerator Used to pick random pixels in the mask
   * @param uniformDistribution Used for sub-pixel location and z-depth
   * @param allowZeroMask Set to true to allow a mask with no non-zero pixels
   */
  MaskDistribution(int[] mask, int width, int height, double depth, double scaleX, double scaleY,
      RandomGenerator randomGenerator, UniformDistribution uniformDistribution,
      boolean allowZeroMask) {
    if (width < 1 || height < 1) {
      throw new IllegalArgumentException("Dimensions must be above zero");
    }
    if (scaleX < 0 || scaleY < 0) {
      throw new IllegalArgumentException("Scale must be above zero");
    }
    if (mask == null || mask.length < width * height) {
      throw new IllegalArgumentException(
          "Mask must not be null and must at least (width * height) in size");
    }
    if (randomGenerator == null) {
      randomGenerator = new Well19937c(System.currentTimeMillis() + System.identityHashCode(this));
    }

    this.randomGenerator = randomGenerator;
    setUniformDistribution(uniformDistribution);
    this.mask = mask;
    maxx = width;
    maxy = height;
    this.scaleX = scaleX;
    this.scaleY = scaleY;
    halfWidth = width * 0.5;
    halfHeight = height * 0.5;
    minDepth = depth * -0.5;
    this.depth = depth;

    final int size = width * height;
    int count = 0;
    for (int i = 0; i < size; i++) {
      if (mask[i] != 0) {
        count++;
      }
    }

    if (count == 0 && !allowZeroMask) {
      throw new IllegalArgumentException("Mask must have non-zero pixels");
    }

    indices = new int[count];
    count = 0;
    for (int i = 0; i < size; i++) {
      if (mask[i] != 0) {
        indices[count++] = i;
      }
    }

    // Fischer-Yates shuffle the indices to scramble the mask positions
    for (int i = indices.length; i-- > 1;) {
      final int j = randomGenerator.nextInt(i + 1);
      final int tmp = indices[i];
      indices[i] = indices[j];
      indices[j] = tmp;
    }
  }

  @Override
  public double[] next() {
    final int randomPosition = randomGenerator.nextInt(indices.length);
    final int y = indices[randomPosition] / maxx;
    final int x = indices[randomPosition] % maxx;
    final double[] d = uniformDistribution.nextUnit();
    // Ensure XY = 0 is the centre of the image by subtracting half the width/height
    d[0] = (x + d[0] - halfWidth) * scaleX;
    d[1] = (y + d[1] - halfHeight) * scaleY;
    d[2] = minDepth + d[2] * depth;
    return d;
  }

  @Override
  public boolean isWithin(double[] xyz) {
    if (!isWithinXy(xyz)) {
      return false;
    }
    return (xyz[2] >= minDepth && xyz[2] <= minDepth + depth);
  }

  @Override
  public boolean isWithinXy(double[] xyz) {
    // Ensure XY = 0 is the centre of the image
    final int index = getIndex(xyz);
    if (index < 0 || index >= mask.length || mask[index] == 0) {
      return false;
    }
    // Check if the search was initialised in a particle
    if (particle == 0) {
      // No starting particle so just accept the position
      return true;
    }
    // Must be in the same particle as the initial position
    return mask[index] == particle;
  }

  private int getIndex(double[] xyz) {
    final int x = (int) (halfWidth + xyz[0] / scaleX);
    final int y = (int) (halfHeight + xyz[1] / scaleY);
    if (x < 0 || x >= maxx || y < 0 || y >= maxy) {
      return -1;
    }
    return y * maxx + x;
  }

  @Override
  public void initialise(double[] xyz) {
    findParticles();

    // Now store the particle that contains the position
    final int index = getIndex(xyz);
    particle = (index < 0 || index >= mask.length) ? 0 : mask[index];
  }

  /**
   * Convert the mask to connected particles, each with a unique number. This allows the within
   * function to restrict movement to the particle of origin
   */
  private void findParticles() {
    // Check if already initialised
    if (xlimit != -1) {
      return;
    }

    xlimit = maxx - 1;
    ylimit = maxy - 1;

    // Create the offset table (for single array 2D neighbour comparisons)
    offset = new int[DIR_X_OFFSET.length];
    for (int d = offset.length; d-- > 0;) {
      offset[d] = maxx * DIR_Y_OFFSET[d] + DIR_X_OFFSET[d];
    }

    final int[] pointList = new int[mask.length];

    // Store all the non-zero positions
    final boolean[] binaryMask = new boolean[mask.length];
    for (int i = 0; i < mask.length; i++) {
      binaryMask[i] = (mask[i] != 0);
    }

    // Find particles
    int particles = 0;
    for (int i = 0; i < binaryMask.length; i++) {
      if (binaryMask[i]) {
        expandParticle(binaryMask, mask, pointList, i, ++particles);
      }
    }

    // Free memory
    offset = null;
  }

  /**
   * Searches from the specified point to find all connected points and assigns them to given
   * particle.
   */
  private void expandParticle(boolean[] binaryMask, int[] mask, int[] pointList, int index0,
      final int particle) {
    binaryMask[index0] = false; // mark as processed
    int listI = 0; // index of current search element in the list
    int listLen = 1; // number of elements in the list

    // we create a list of connected points and start the list at the particle start position
    pointList[listI] = index0;

    do {
      final int index1 = pointList[listI];
      // Mark this position as part of the particle
      mask[index1] = particle;

      // Search the 8-connected neighbours
      final int x1 = index1 % maxx;
      final int y1 = index1 / maxx;

      final boolean isInnerXy = (x1 != 0 && x1 != xlimit) && (y1 != 0 && y1 != ylimit);

      if (isInnerXy) {
        for (int d = 8; d-- > 0;) {
          final int index2 = index1 + offset[d];
          if (binaryMask[index2]) {
            binaryMask[index2] = false; // mark as processed
            // Add this to the search
            pointList[listLen++] = index2;
          }
        }
      } else {
        for (int d = 8; d-- > 0;) {
          if (isWithinDirection(x1, y1, d)) {
            final int index2 = index1 + offset[d];
            if (binaryMask[index2]) {
              binaryMask[index2] = false; // mark as processed
              // Add this to the search
              pointList[listLen++] = index2;
            }
          }
        }
      }

      listI++;

    }
    while (listI < listLen);
  }

  /**
   * Returns whether the neighbour in a given direction is within the image. NOTE: it is assumed
   * that the pixel x,y itself is within the image! Uses class variables xlimit, ylimit: (dimensions
   * of the image)-1
   *
   * @param x x-coordinate of the pixel that has a neighbour in the given direction
   * @param y y-coordinate of the pixel that has a neighbour in the given direction
   * @param direction the direction from the pixel towards the neighbour
   * @return true if the neighbour is within the image (provided that x, y is within)
   */
  private boolean isWithinDirection(int x, int y, int direction) {
    switch (direction) {
      case 0:
        return (y > 0);
      case 1:
        return (y > 0 && x < xlimit);
      case 2:
        return (x < xlimit);
      case 3:
        return (y < ylimit && x < xlimit);
      case 4:
        return (y < ylimit);
      case 5:
        return (y < ylimit && x > 0);
      case 6:
        return (x > 0);
      case 7:
        return (y > 0 && x > 0);
      default:
        return false;
    }
  }

  /**
   * Gets the number of non-zero pixels in the mask.
   *
   * @return The number of non-zero pixels in the mask.
   */
  public int getSize() {
    return indices.length;
  }

  /**
   * Gets the width of the mask in pixels.
   *
   * @return The width
   */
  public int getWidth() {
    return maxx;
  }

  /**
   * Gets the height of the mask in pixels.
   *
   * @return The height
   */
  public int getHeight() {
    return maxy;
  }

  /**
   * Gets the X-scale.
   *
   * @return The X-scale.
   */
  public double getScaleX() {
    return scaleX;
  }

  /**
   * Gets the Y-scale.
   *
   * @return The Y-scale.
   */
  public double getScaleY() {
    return scaleY;
  }

  /**
   * Gets the mask.
   *
   * @return The mask (packed in YX order).
   */
  protected int[] getMask() {
    return mask;
  }

  /**
   * The UniformDistribution to pick the sub pixel x,y coordinates and z-depth.
   *
   * @param uniformDistribution the uniformDistribution to set
   */
  public void setUniformDistribution(UniformDistribution uniformDistribution) {
    if (uniformDistribution == null) {
      uniformDistribution =
          new UniformDistribution(null, new double[] {1, 1, 1}, randomGenerator.nextInt());
    }
    this.uniformDistribution = uniformDistribution;
  }
}
