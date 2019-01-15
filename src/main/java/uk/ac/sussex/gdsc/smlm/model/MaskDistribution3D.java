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

import java.util.Arrays;
import java.util.List;

/**
 * Samples uniformly from the specified masks. All non-zero pixels are sampled. The centre of the
 * mask stack corresponds to XY=0. The z coordinate is randomly sampled from the slice depth offset
 * by the slice position in the stack. The distribution of Z is centred on zero.
 *
 * <p>X coordinates are returned in the interval -width/2 to width/2. These can be converted to
 * different values using the scale parameter. Likewise for the Y coordinates. E.g. a mask of
 * 100x100 (range of -50:50) can be used to generate coordinates in the range -100:100 using a scale
 * of 2.
 *
 * <p>Sub pixel locations and z-depth are sampled from a uniform distribution. A Halton sequence is
 * used by default but this can be changed by setting a custom uniform distribution.
 */
public class MaskDistribution3D implements SpatialDistribution {
  // Used for a particle search
  private static final int[] DIR_X_OFFSET =
      {0, 1, 1, 1, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1, -1, -1, 0, 0, 1, 1, 1, 0, -1, -1, -1, 0};
  private static final int[] DIR_Y_OFFSET =
      {-1, -1, 0, 1, 1, 1, 0, -1, -1, -1, 0, 1, 1, 1, 0, -1, 0, -1, -1, 0, 1, 1, 1, 0, -1, 0};
  private static final int[] DIR_Z_OFFSET =
      {0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  private final RandomGenerator randomGenerator;
  private UniformDistribution uniformDistribution;
  private final int[] mask;
  private int[] indices;
  private final int maxx;
  private final int maxy;
  private final int maxz;
  private final int maxxMaxy;
  private final double halfWidth;
  private final double halfHeight;
  private final double minDepth;
  private final double depth;
  private int particle;
  private final double scaleX;
  private final double scaleY;

  private final double sliceDepth;
  private MaskDistribution projection;

  // Used for a particle search
  private int xlimit = -1;
  private int ylimit;
  private int zlimit;
  private int[] offset;

  /**
   * Create a distribution from the stack of mask images (packed in YX order).
   *
   * @param masks the masks
   * @param width The width of the mask in pixels
   * @param height the height of the mask in pixels
   * @param sliceDepth The depth of each slice
   * @param scaleX Used to scale the mask X-coordinate to a new value
   * @param scaleY Used to scale the mask Y-coordinate to a new value
   */
  public MaskDistribution3D(List<int[]> masks, int width, int height, double sliceDepth,
      double scaleX, double scaleY) {
    this(masks, width, height, sliceDepth, scaleX, scaleY, null);
  }

  /**
   * Create a distribution from the stack of mask images (packed in YX order).
   *
   * @param masks the masks
   * @param width The width of the mask in pixels
   * @param height the height of the mask in pixels
   * @param sliceDepth The depth of each slice
   * @param scaleX Used to scale the mask X-coordinate to a new value
   * @param scaleY Used to scale the mask Y-coordinate to a new value
   * @param randomGenerator Used to pick random pixels in the mask
   */
  public MaskDistribution3D(List<int[]> masks, int width, int height, double sliceDepth,
      double scaleX, double scaleY, RandomGenerator randomGenerator) {
    this(masks, width, height, sliceDepth, scaleX, scaleY, randomGenerator, null);
  }

  /**
   * Create a distribution from the stack of mask images (packed in YX order).
   *
   * @param masks the masks
   * @param width The width of the mask in pixels
   * @param height the height of the mask in pixels
   * @param sliceDepth The depth of each slice
   * @param scaleX Used to scale the mask X-coordinate to a new value
   * @param scaleY Used to scale the mask Y-coordinate to a new value
   * @param randomGenerator Used to pick random pixels in the mask
   * @param uniformDistribution Used for sub-pixel location and slice z-depth
   */
  public MaskDistribution3D(List<int[]> masks, int width, int height, double sliceDepth,
      double scaleX, double scaleY, RandomGenerator randomGenerator,
      UniformDistribution uniformDistribution) {
    if (width < 1 || height < 1) {
      throw new IllegalArgumentException("Dimensions must be above zero");
    }
    if (sliceDepth <= 0) {
      throw new IllegalArgumentException("Slice depth must be above zero");
    }
    if (masks == null || masks.isEmpty()) {
      throw new IllegalArgumentException("Mask must not be null or empty");
    }
    if (randomGenerator == null) {
      randomGenerator = new Well19937c(System.currentTimeMillis() + System.identityHashCode(this));
    }

    this.randomGenerator = randomGenerator;
    setUniformDistribution(uniformDistribution);
    maxx = width;
    maxy = height;
    maxz = masks.size();
    this.scaleX = scaleX;
    this.scaleY = scaleY;
    halfWidth = width * 0.5;
    halfHeight = height * 0.5;
    this.sliceDepth = sliceDepth;
    depth = sliceDepth * maxz;
    minDepth = depth * -0.5;

    maxxMaxy = maxx * maxy;
    mask = new int[maxz * maxxMaxy];
    indices = new int[mask.length];

    int count = 0;
    int index = 0;
    for (final int[] mask : masks) {
      if (mask.length < maxxMaxy) {
        throw new IllegalArgumentException("Masks must be the same size");
      }
      for (int i = 0; i < maxxMaxy; i++) {
        this.mask[index] = mask[i];
        if (mask[i] != 0) {
          indices[count++] = index;
        }
        index++;
      }
    }

    if (count == 0) {
      throw new IllegalArgumentException("Mask must have non-zero pixels");
    }

    indices = Arrays.copyOf(indices, count);

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
    final int[] xyz = new int[3];
    getXyz(indices[randomPosition], xyz);
    final double[] d = uniformDistribution.nextUnit();

    // Ensure XY = 0 is the centre of the image by subtracting half the width/height
    d[0] = (xyz[0] + d[0] - halfWidth) * scaleX;
    d[1] = (xyz[1] + d[1] - halfHeight) * scaleY;
    d[2] = (xyz[2] + d[2]) * sliceDepth + minDepth;

    return d;
  }

  private int getIndex(double[] xyz) {
    final int x = (int) (halfWidth + xyz[0] / scaleX);
    final int y = (int) (halfHeight + xyz[1] / scaleY);
    final int z = (int) ((xyz[2] - minDepth) / sliceDepth);
    if (x < 0 || x >= maxx || y < 0 || y >= maxy || z < 0 || z >= maxz) {
      return -1;
    }
    return getIndex(x, y, z);
  }

  /**
   * Return the single index associated with the x,y,z coordinates.
   *
   * @param x the x
   * @param y the y
   * @param z the z
   * @return The index
   */
  private int getIndex(int x, int y, int z) {
    return (maxxMaxy) * z + maxx * y + x;
  }

  @Override
  public boolean isWithin(double[] xyz) {
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

  @Override
  public boolean isWithinXy(double[] xyz) {
    createProjection();
    return projection.isWithinXy(xyz);
  }

  private void createProjection() {
    // Create a projection of the masks
    if (projection == null) {
      final int[] mask2 = new int[maxxMaxy];
      for (int z = 0, index = 0; z < maxz; z++) {
        for (int j = 0; j < maxxMaxy; j++) {
          if (mask[index++] != 0) {
            mask2[j] = 1;
          }
        }
      }
      projection =
          new MaskDistribution(mask2, maxx, maxy, sliceDepth, scaleX, scaleY, randomGenerator);
    }
  }

  @Override
  public void initialise(double[] xyz) {
    findParticles();

    // Now store the particle that contains the position
    final int index = getIndex(xyz);
    particle = (index < 0 || index >= mask.length) ? 0 : mask[index];

    // Also initialise for isWithinXy()
    createProjection();
    projection.initialise(xyz);
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
    zlimit = maxz - 1;

    // Create the offset table (for single array 3D neighbour comparisons)
    offset = new int[DIR_X_OFFSET.length];
    for (int d = offset.length; d-- > 0;) {
      offset[d] = getIndex(DIR_X_OFFSET[d], DIR_Y_OFFSET[d], DIR_Z_OFFSET[d]);
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
   * Convert the single index into x,y,z coords, Input array must be length >= 3.
   *
   * @param index the index
   * @param xyz the xyz
   * @return The xyz array
   */
  private int[] getXyz(int index, int[] xyz) {
    xyz[2] = index / (maxxMaxy);
    final int mod = index % (maxxMaxy);
    xyz[1] = mod / maxx;
    xyz[0] = mod % maxx;
    return xyz;
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

    final int[] xyz = new int[3];

    // we create a list of connected points and start the list at the particle start position
    pointList[listI] = index0;

    do {
      final int index1 = pointList[listI];
      // Mark this position as part of the particle
      mask[index1] = particle;

      getXyz(index1, xyz);

      final int x1 = xyz[0];
      final int y1 = xyz[1];
      final int z1 = xyz[2];

      final boolean isInnerXy = (y1 != 0 && y1 != ylimit) && (x1 != 0 && x1 != xlimit);
      final boolean isInnerXyz = (zlimit == 0) ? isInnerXy : isInnerXy && (z1 != 0 && z1 != zlimit);

      // Search the neighbours
      for (int d = 26; d-- > 0;) {
        if (isInnerXyz || (isInnerXy && isWithinZ(z1, d)) || isWithinXyz(x1, y1, z1, d)) {
          final int index2 = index1 + offset[d];
          if (binaryMask[index2]) {
            binaryMask[index2] = false; // mark as processed
            // Add this to the search
            pointList[listLen++] = index2;
          }
        }
      }

      listI++;

    }
    while (listI < listLen);
  }

  /**
   * returns whether the neighbour in a given direction is within the image. NOTE: it is assumed
   * that the pixel x,y,z itself is within the image! Uses class variables xlimit, ylimit, zlimit:
   * (dimensions of the image)-1
   *
   * @param x x-coordinate of the pixel that has a neighbour in the given direction
   * @param y y-coordinate of the pixel that has a neighbour in the given direction
   * @param z z-coordinate of the pixel that has a neighbour in the given direction
   * @param direction the direction from the pixel towards the neighbour
   * @return true if the neighbour is within the image (provided that x, y, z is within)
   */
  private boolean isWithinXyz(int x, int y, int z, int direction) {
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
      case 8:
        return (z > 0 && y > 0);
      case 9:
        return (z > 0 && y > 0 && x < xlimit);
      case 10:
        return (z > 0 && x < xlimit);
      case 11:
        return (z > 0 && y < ylimit && x < xlimit);
      case 12:
        return (z > 0 && y < ylimit);
      case 13:
        return (z > 0 && y < ylimit && x > 0);
      case 14:
        return (z > 0 && x > 0);
      case 15:
        return (z > 0 && y > 0 && x > 0);
      case 16:
        return (z > 0);
      case 17:
        return (z < zlimit && y > 0);
      case 18:
        return (z < zlimit && y > 0 && x < xlimit);
      case 19:
        return (z < zlimit && x < xlimit);
      case 20:
        return (z < zlimit && y < ylimit && x < xlimit);
      case 21:
        return (z < zlimit && y < ylimit);
      case 22:
        return (z < zlimit && y < ylimit && x > 0);
      case 23:
        return (z < zlimit && x > 0);
      case 24:
        return (z < zlimit && y > 0 && x > 0);
      case 25:
        return (z < zlimit);
      default:
        return false;
    }
  }

  /**
   * returns whether the neighbour in a given direction is within the image. NOTE: it is assumed
   * that the pixel z itself is within the image! Uses class variables zlimit: (dimensions of the
   * image)-1
   *
   * @param z z-coordinate of the pixel that has a neighbour in the given direction
   * @param direction the direction from the pixel towards the neighbour
   * @return true if the neighbour is within the image (provided that z is within)
   */
  private boolean isWithinZ(int z, int direction) {
    // z = 0
    if (direction < 8) {
      return true;
    }
    // z = -1
    if (direction < 17) {
      return (z > 0);
    }
    // z = 1
    return z < zlimit;
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
   * @return The mask (packed in ZYX order).
   */
  protected int[] getMask() {
    return mask;
  }

  /**
   * The UniformDistribution to pick the sub pixel x,y coordinates and slice z-depth.
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
