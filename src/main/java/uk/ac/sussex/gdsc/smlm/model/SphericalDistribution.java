/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
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

import org.apache.commons.math3.random.JDKRandomGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.FastMath;

/**
 * Samples uniformly from the specified spherical volume.
 */
public class SphericalDistribution implements SpatialDistribution {
  private final double radius, r2, range;
  private final RandomGenerator randomGenerator;
  private boolean useRejectionMethod = true;
  private final double[] origin = new double[3];

  /**
   * Instantiates a new spherical distribution.
   *
   * @param radius the radius
   */
  public SphericalDistribution(double radius) {
    this(radius, null);
  }

  /**
   * Instantiates a new spherical distribution.
   *
   * @param radius the radius
   * @param randomGenerator the random generator
   */
  public SphericalDistribution(double radius, RandomGenerator randomGenerator) {
    if (randomGenerator == null) {
      randomGenerator = new JDKRandomGenerator();
    }
    if (radius < 0) {
      throw new IllegalArgumentException("Radius must be positive: {0}");
    }
    this.radius = radius;
    this.r2 = radius * radius;
    this.range = 2 * radius;
    this.randomGenerator = randomGenerator;
  }

  /** {@inheritDoc} */
  @Override
  public double[] next() {
    final double[] xyz = new double[3];
    if (radius > 0) {
      if (useRejectionMethod) {
        // -=-=-=-
        // Rejection method:
        // Sample from a cube and then check if within a sphere
        // -=-=-=-
        double d2 = 0;
        do {
          for (int i = 0; i < 3; i++) {
            // xyz[i] = randomGenerator.nextDouble() * ((randomGenerator.nextBoolean()) ? -radius :
            // radius);
            // Avoid extra call to the random generator
            xyz[i] = randomGenerator.nextDouble() * range - radius;
          }
          d2 = xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2];
        }
        while (d2 > r2);
      } else {
        // -=-=-=-
        // Transformation method:
        // Generate a random point on the surface of the sphere and then sample within.
        // -=-=-=-

        // Generate a random unit vector: X1, X2, X3 sampled with mean 0 and variance 1
        for (int i = 0; i < 3; i++) {
          xyz[i] = randomGenerator.nextGaussian();
        }

        // Calculate the distance: RsU^1/3 / length
        final double d = (radius * FastMath.cbrt(randomGenerator.nextDouble()))
            / Math.sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
        for (int i = 0; i < 3; i++) {
          xyz[i] *= d;
        }
      }
    }
    return xyz;
  }

  /** {@inheritDoc} */
  @Override
  public boolean isWithin(double[] xyz) {
    final double[] delta = {xyz[0] - origin[0], xyz[1] - origin[1], xyz[2] - origin[2]};
    return (delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]) < r2;
  }

  /** {@inheritDoc} */
  @Override
  public boolean isWithinXY(double[] xyz) {
    final double[] delta = {xyz[0] - origin[0], xyz[1] - origin[1]};
    return (delta[0] * delta[0] + delta[1] * delta[1]) < r2;
  }

  /**
   * @return If true then sample from the distribution using the rejection method. The alternative
   *         is a transformation method.
   */
  public boolean isUseRejectionMethod() {
    return useRejectionMethod;
  }

  /**
   * @param useRejectionMethod If true then sample from the distribution using the rejection method.
   *        The alternative is a transformation method.
   */
  public void setUseRejectionMethod(boolean useRejectionMethod) {
    this.useRejectionMethod = useRejectionMethod;
  }

  /** {@inheritDoc} */
  @Override
  public void initialise(double[] xyz) {
    if (xyz != null && xyz.length > 2) {
      for (int i = 0; i < 3; i++) {
        origin[i] = xyz[i];
      }
    }
  }
}
