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

import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import org.apache.commons.rng.simple.RandomSource;

import java.util.Arrays;

/**
 * Populates an image with well spaced unary or binary localisations.
 *
 * <p>Creates a grid layout using cells of the specified size within the image area. Each cell can
 * have one or two localisations. The first localisation is placed within the central 50% of the
 * cell. The second localisation (if present) is placed randomly up to a maximum distance away. When
 * all cells have been sampled then no more localisations are generated.
 */
public class GridDistribution implements SpatialDistribution {
  private final UniformRandomProvider rng;
  private final NormalizedGaussianSampler gauss;
  private final int size;
  private final int cellSize;
  private final double probBinary;
  private final double minBinaryDistance;
  private final double maxBinaryDistance;
  private final double min;
  private final double depth;

  private int cell = -1;
  private final int cellsPerRow;
  private final int totalCells;
  private double[] previous;

  /**
   * Create a distribution with the binary spots placed from 0 - distance.
   *
   * @param size the size
   * @param depth the depth
   * @param cellSize the cell size
   * @param probBinary the probability of a binary spot
   * @param binaryDistance the probability of a binary spot distance
   */
  public GridDistribution(int size, double depth, int cellSize, double probBinary,
      double binaryDistance) {
    this(size, depth, cellSize, probBinary, 0, binaryDistance, null);
  }

  /**
   * Create a distribution with the binary spots placed from min - max distance.
   *
   * @param size the size
   * @param depth the depth
   * @param cellSize the cell size
   * @param probBinary the probability of a binary spot
   * @param minBinaryDistance the min binary distance
   * @param maxBinaryDistance the max binary distance
   */
  public GridDistribution(int size, double depth, int cellSize, double probBinary,
      double minBinaryDistance, double maxBinaryDistance) {
    this(size, depth, cellSize, probBinary, minBinaryDistance, maxBinaryDistance, null);
  }

  /**
   * Create a distribution with the binary spots placed from min - max distance.
   *
   * @param size the size
   * @param depth the depth
   * @param cellSize the cell size
   * @param probBinary the probability of a binary spot
   * @param minBinaryDistance the min binary distance
   * @param maxBinaryDistance the max binary distance
   * @param randomGenerator the random generator
   */
  public GridDistribution(int size, double depth, int cellSize, double probBinary,
      double minBinaryDistance, double maxBinaryDistance, UniformRandomProvider randomGenerator) {
    if (size < 1) {
      throw new IllegalArgumentException("Size must be above zero");
    }
    if (size < cellSize) {
      throw new IllegalArgumentException("Size must be >= cell size");
    }
    if (probBinary < 0 || probBinary > 1) {
      throw new IllegalArgumentException("Probability must be between 0 and 1");
    }
    if (maxBinaryDistance < 0) {
      throw new IllegalArgumentException("Max distance must be positive");
    }
    if (minBinaryDistance > maxBinaryDistance) {
      throw new IllegalArgumentException("Min distance must be below max distance");
    }
    if (randomGenerator == null) {
      this.rng = RandomSource.create(RandomSource.XOR_SHIFT_1024_S);
    } else {
      this.rng = randomGenerator;
    }
    gauss = SamplerUtils.createNormalizedGaussianSampler(rng);
    this.size = size;
    this.min = -depth / 2;
    this.depth = depth;
    this.cellSize = cellSize;
    this.probBinary = probBinary;
    this.minBinaryDistance = minBinaryDistance;
    this.maxBinaryDistance = maxBinaryDistance;

    cellsPerRow = size / cellSize;
    totalCells = cellsPerRow * cellsPerRow;
  }

  @Override
  public double[] next() {
    // See if a binary localisation should be created near the previous spot
    if (previous != null && rng.nextDouble() < probBinary) {
      final double[] xyz = Arrays.copyOf(previous, 3);

      // Create a random unit vector

      double x = gauss.sample();
      double y = gauss.sample();
      double z = gauss.sample();
      final double length = Math.sqrt(x * x + y * y + z * z);
      if (length != 0) {
        // Shift by a random distance
        final double distance = (maxBinaryDistance == minBinaryDistance) ? maxBinaryDistance
            : nextUniform(minBinaryDistance, maxBinaryDistance);
        final double d = distance / length;
        x *= d;
        y *= d;
        z *= d;
      }
      xyz[0] += x;
      xyz[1] += y;
      xyz[2] += z;
      previous = null;
      return xyz;
    }
    previous = null;
    // See if any more localisations will fit in the grid
    if (++cell < totalCells) {
      final int cellx = cell % cellsPerRow;
      final int celly = cell / cellsPerRow;

      previous = new double[3];
      // Ensure the centre of the distribution is [0,0,0]
      previous[0] = cellx * cellSize - size / 2 + cellSize * nextUniform(0.25, 0.75);
      previous[1] = celly * cellSize - size / 2 + cellSize * nextUniform(0.25, 0.75);
      previous[2] = min + rng.nextDouble() * depth;
    }
    return previous;
  }

  /**
   * Create a uniform deviate between low and high inclusive.
   *
   * @param lo the low
   * @param hi the high
   * @return the double
   */
  private double nextUniform(double lo, double hi) {
    final double u = rng.nextDouble();
    return u * hi + (1 - u) * lo;
  }

  @Override
  public boolean isWithin(double[] xyz) {
    for (int i = 0; i < xyz.length; i++) {
      if (xyz[i] < 0 || xyz[i] > size) {
        return false;
      }
    }
    return true;
  }

  @Override
  public boolean isWithinXy(double[] xyz) {
    for (int i = 0; i < 2; i++) {
      if (xyz[i] < 0 || xyz[i] > size) {
        return false;
      }
    }
    return true;
  }

  @Override
  public void initialise(double[] xyz) {
    // Ignore
  }
}
