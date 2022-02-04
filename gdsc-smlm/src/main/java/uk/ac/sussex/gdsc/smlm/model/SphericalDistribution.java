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

package uk.ac.sussex.gdsc.smlm.model;

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.ObjectSampler;
import org.apache.commons.rng.sampling.shape.UnitBallSampler;
import uk.ac.sussex.gdsc.core.utils.ValidationUtils;
import uk.ac.sussex.gdsc.core.utils.rng.UniformRandomProviders;

/**
 * Samples uniformly from the specified spherical volume.
 */
public class SphericalDistribution implements SpatialDistribution {
  private final double r2;
  private final double[] origin = new double[3];
  private ObjectSampler<double[]> sampler;

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
  public SphericalDistribution(double radius, UniformRandomProvider randomGenerator) {
    ValidationUtils.checkPositive(radius, "Radius");
    if (randomGenerator == null) {
      randomGenerator = UniformRandomProviders.create();
    }
    if (radius > 0) {
      r2 = radius * radius;
      final ObjectSampler<double[]> ball = UnitBallSampler.of(randomGenerator, 3);
      final double r = radius;
      sampler = () -> {
        final double[] next = ball.sample();
        next[0] *= r;
        next[1] *= r;
        next[2] *= r;
        return next;
      };
    } else {
      r2 = 0;
      sampler = () -> new double[3];
    }
  }

  @Override
  public double[] next() {
    return sampler.sample();
  }

  @Override
  public boolean isWithin(double[] xyz) {
    final double[] delta = {xyz[0] - origin[0], xyz[1] - origin[1], xyz[2] - origin[2]};
    return (delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]) < r2;
  }

  @Override
  public boolean isWithinXy(double[] xyz) {
    final double[] delta = {xyz[0] - origin[0], xyz[1] - origin[1]};
    return (delta[0] * delta[0] + delta[1] * delta[1]) < r2;
  }

  @Override
  public void initialise(double[] xyz) {
    if (xyz != null && xyz.length > 2) {
      for (int i = 0; i < 3; i++) {
        origin[i] = xyz[i];
      }
    }
  }
}
