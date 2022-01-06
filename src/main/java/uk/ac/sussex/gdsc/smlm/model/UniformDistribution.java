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

import java.util.function.Supplier;
import org.apache.commons.math3.random.HaltonSequenceGenerator;
import org.apache.commons.math3.random.RandomVectorGenerator;
import org.apache.commons.rng.UniformRandomProvider;

/**
 * Samples uniformly from the specified volume.
 *
 * <p>Uses a Halton sequence generator by default. This can be overriden by passing in a vector
 * sequence generator or by providing a factory for a random generator (which should produce uniform
 * equi-distributed numbers in the domain [0,1]).
 */
public class UniformDistribution implements SpatialDistribution {
  /**
   * Wrap a standard random generator to create a vector generator for 3 dimensions.
   */
  private static class VectorGeneratorWrapper implements RandomVectorGenerator {
    private final UniformRandomProvider rng1;
    private final UniformRandomProvider rng2;
    private final UniformRandomProvider rng3;

    VectorGeneratorWrapper(Supplier<UniformRandomProvider> randomGeneratorFactory) {
      rng1 = create(randomGeneratorFactory);
      rng2 = create(randomGeneratorFactory);
      rng3 = create(randomGeneratorFactory);
    }

    private static UniformRandomProvider
        create(Supplier<UniformRandomProvider> randomGeneratorFactory) {
      final UniformRandomProvider rng = randomGeneratorFactory.get();
      if (rng == null) {
        throw new NullPointerException(
            "UniformRandomProviderFactory created null UniformRandomProvider");
      }
      return rng;
    }

    @Override
    public double[] nextVector() {
      return new double[] {rng1.nextDouble(), rng2.nextDouble(), rng3.nextDouble()};
    }
  }

  private double[] min;
  private double[] max;
  private double[] range;
  private RandomVectorGenerator vectorGenerator;

  /**
   * Create a new uniform distribution using a Halton sequence. The minimum bounds are set to zero.
   *
   * @param max The maximum bounds for the distribution
   */
  public UniformDistribution(double[] max) {
    init(new double[3], max, null);
  }

  /**
   * Create a new uniform distribution using a Halton sequence.
   *
   * @param min The minimum bounds for the distribution
   * @param max The maximum bounds for the distribution
   */
  public UniformDistribution(double[] min, double[] max) {
    init(min, max, null);
  }

  /**
   * Create a new uniform distribution using a Halton sequence.
   *
   * @param min The minimum bounds for the distribution
   * @param max The maximum bounds for the distribution
   * @param seed Start at the i-th point in the Halton sequence
   */
  public UniformDistribution(double[] min, double[] max, int seed) {
    // The Halton sequence based on the prime of 2 does not provide great variety in the
    // lesser significant digits when simulating a 512x512 pixel image. This is not suitable for
    // PSF fitting since we require variation to at least 3 decimal places. So start at the
    // prime of 3.
    final HaltonSequenceGenerator randomVectorGenerator =
        new HaltonSequenceGenerator(3, new int[] {3, 5, 7}, null);
    randomVectorGenerator.skipTo(Math.abs(seed));
    init(min, max, randomVectorGenerator);
  }

  /**
   * Create a new uniform distribution using the given vector generator.
   *
   * @param min The minimum bounds for the distribution
   * @param max The maximum bounds for the distribution
   * @param randomVectorGenerator Must produce vectors with dimension 3 (or above)
   */
  public UniformDistribution(double[] min, double[] max,
      RandomVectorGenerator randomVectorGenerator) {
    init(min, max, randomVectorGenerator);
  }

  /**
   * Create a new uniform distribution using a new random number generator from the factory for each
   * dimension.
   *
   * @param min The minimum bounds for the distribution
   * @param max The maximum bounds for the distribution
   * @param randomGeneratorFactory Must produce random number generators with uniform random numbers
   *        in the domain [0,1]
   */
  public UniformDistribution(double[] min, double[] max,
      Supplier<UniformRandomProvider> randomGeneratorFactory) {
    final RandomVectorGenerator randomVectorGenerator =
        new VectorGeneratorWrapper(randomGeneratorFactory);
    init(min, max, randomVectorGenerator);
  }

  private void init(double[] min, double[] max, RandomVectorGenerator randomVectorGenerator) {
    if (min == null) {
      min = new double[0];
    }
    if (max == null) {
      max = new double[0];
    }

    this.min = new double[3];
    for (int i = 0; i < min.length; i++) {
      this.min[i] = min[i];
    }
    this.max = new double[3];
    for (int i = 0; i < max.length; i++) {
      if (max[i] < this.min[i]) {
        throw new IllegalArgumentException(
            String.format("Max %f must be greater than min %f", max[i], this.min[i]));
      }
      this.max[i] = max[i];
    }

    this.range = new double[3];
    for (int i = 0; i < this.max.length; i++) {
      this.range[i] = this.max[i] - this.min[i];
    }

    if (randomVectorGenerator == null) {
      randomVectorGenerator = new HaltonSequenceGenerator(3);
    }
    this.vectorGenerator = randomVectorGenerator;
  }

  @Override
  public double[] next() {
    final double[] d = vectorGenerator.nextVector();
    for (int i = 0; i < 3; i++) {
      d[i] = min[i] + d[i] * range[i];
    }
    return d;
  }

  /**
   * Return a vector with values in the unit domain ([0,1]).
   *
   * @return the vector populated with values in the unit domain ([0,1])
   */
  public double[] nextUnit() {
    return vectorGenerator.nextVector();
  }

  @Override
  public boolean isWithin(double[] xyz) {
    for (int i = 0; i < xyz.length; i++) {
      if (xyz[i] < min[i] || xyz[i] > max[i]) {
        return false;
      }
    }
    return true;
  }

  @Override
  public boolean isWithinXy(double[] xyz) {
    for (int i = 0; i < 2; i++) {
      if (xyz[i] < min[i] || xyz[i] > max[i]) {
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
