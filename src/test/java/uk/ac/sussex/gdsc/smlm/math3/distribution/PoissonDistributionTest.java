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

package uk.ac.sussex.gdsc.smlm.math3.distribution;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;

@SuppressWarnings({"javadoc"})
public class PoissonDistributionTest {
  @SeededTest
  public void canComputeProbability(RandomSeed seed) {
    final UniformRandomProvider rng = RngUtils.create(seed.getSeed());

    final PoissonDistribution fpd = new PoissonDistribution(1);
    for (int i = 1; i <= 100; i++) {
      final double mean = rng.nextDouble() * i;
      final org.apache.commons.math3.distribution.PoissonDistribution pd =
          createReferencePoissonDistribution(mean);
      fpd.setMean(mean);
      final int lower = (int) Math.floor(Math.max(0, mean - 5 * Math.sqrt(mean)));
      final int upper = (int) Math.ceil((mean + 5 * Math.sqrt(mean)));
      for (int x = lower; x <= upper; x++) {
        Assertions.assertEquals(pd.probability(x), fpd.probability(x), 1e-16);
      }
    }
  }

  @SeededTest
  public void canComputeCumulativeProbability(RandomSeed seed) {
    final UniformRandomProvider rng = RngUtils.create(seed.getSeed());

    final PoissonDistribution fpd = new PoissonDistribution(1);
    for (int i = 1; i <= 100; i++) {
      final double mean = rng.nextDouble() * i;
      final org.apache.commons.math3.distribution.PoissonDistribution pd =
          createReferencePoissonDistribution(mean);
      fpd.setMean(mean);
      final int lower = (int) Math.floor(Math.max(0, mean - 5 * Math.sqrt(mean)));
      final int upper = (int) Math.ceil((mean + 5 * Math.sqrt(mean)));
      for (int x = lower; x <= upper; x++) {
        Assertions.assertEquals(pd.cumulativeProbability(x), fpd.cumulativeProbability(x), 1e-16);
      }
    }
  }

  @SeededTest
  public void canComputeInverseCumulativeProbability(RandomSeed seed) {
    final UniformRandomProvider rng = RngUtils.create(seed.getSeed());

    final PoissonDistribution fpd = new PoissonDistribution(1);
    for (int i = 1; i <= 100; i++) {
      final double mean = rng.nextDouble() * i;
      final org.apache.commons.math3.distribution.PoissonDistribution pd =
          createReferencePoissonDistribution(mean);
      fpd.setMean(mean);
      for (int j = 0; j <= 10; j++) {
        final double p = 0.1 * j;
        Assertions.assertEquals(pd.inverseCumulativeProbability(p),
            fpd.inverseCumulativeProbability(p), 1e-16);
      }
    }
  }

  private static org.apache.commons.math3.distribution.PoissonDistribution
      createReferencePoissonDistribution(final double mean) {
    final org.apache.commons.math3.distribution.PoissonDistribution pd =
        new org.apache.commons.math3.distribution.PoissonDistribution(null, mean,
            org.apache.commons.math3.distribution.PoissonDistribution.DEFAULT_EPSILON,
            org.apache.commons.math3.distribution.PoissonDistribution.DEFAULT_MAX_ITERATIONS);
    return pd;
  }

  @Test
  public void testMeanProperty() {
    final double mean = 1.23;
    final PoissonDistribution fpd = new PoissonDistribution(mean);
    Assertions.assertEquals(mean, fpd.getMean());
    fpd.setMean(mean - 1.0);
    Assertions.assertEquals(mean - 1.0, fpd.getMean());
    // Test this does not throw
    fpd.setMeanUnsafe(mean - 2.0);
    Assertions.assertEquals(mean - 2.0, fpd.getMean());
  }

  @Test
  public void testProbabilityEdgeCases() {
    final double mean = 1.23;
    final PoissonDistribution fpd = new PoissonDistribution(mean);
    Assertions.assertEquals(0, fpd.probability(-1));
    Assertions.assertEquals(0, fpd.probability(Integer.MAX_VALUE));
    Assertions.assertEquals(Math.exp(-mean), fpd.probability(0));
  }

  @Test
  public void testLogProbabilityEdgeCases() {
    final double mean = 1.23;
    final PoissonDistribution fpd = new PoissonDistribution(mean);
    Assertions.assertEquals(Double.NEGATIVE_INFINITY, fpd.logProbability(-1));
    Assertions.assertEquals(Double.NEGATIVE_INFINITY, fpd.logProbability(Integer.MAX_VALUE));
    Assertions.assertEquals(-mean, fpd.logProbability(0));
  }

  @Test
  public void testCumulativeProbabilityEdgeCases() {
    final double mean = 1.23;
    final PoissonDistribution fpd = new PoissonDistribution(mean);
    Assertions.assertEquals(0, fpd.cumulativeProbability(-1));
    Assertions.assertEquals(1.0, fpd.cumulativeProbability(Integer.MAX_VALUE));
  }

  @Test
  public void testInverseCumulativeProbabilityEdgeCases() {
    final double mean = 1.23;
    final PoissonDistribution fpd = new PoissonDistribution(mean);
    Assertions.assertEquals(0, fpd.inverseCumulativeProbability(0));
    Assertions.assertEquals(Integer.MAX_VALUE, fpd.inverseCumulativeProbability(1.0));
    Assertions.assertThrows(IllegalArgumentException.class,
        () -> fpd.inverseCumulativeProbability(-0.1));
    Assertions.assertThrows(IllegalArgumentException.class,
        () -> fpd.inverseCumulativeProbability(1.1));
  }

  @Test
  public void testLogProbabilityAtExtremeValue() {
    final double mean = Double.MIN_VALUE;
    final PoissonDistribution fpd = new PoissonDistribution(mean);
    Assertions.assertEquals(Double.NEGATIVE_INFINITY, fpd.logProbability(Integer.MAX_VALUE - 1));
  }

  @Test
  public void testInverseCumulativeProbabilityUpperLimitAtExtremeValue() {
    Assertions.assertEquals(Integer.MAX_VALUE, (int) Math.floor(Integer.MAX_VALUE + 0.5));
    Assertions.assertEquals(Integer.MAX_VALUE, (int) Math.floor(Double.MAX_VALUE));
    Assertions.assertEquals(Integer.MAX_VALUE, (int) Math.floor(Double.POSITIVE_INFINITY));
  }
}
