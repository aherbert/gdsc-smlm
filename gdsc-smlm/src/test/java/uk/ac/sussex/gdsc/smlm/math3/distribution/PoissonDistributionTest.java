/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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
import uk.ac.sussex.gdsc.test.api.Predicates;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.function.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngFactory;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;

@SuppressWarnings({"javadoc"})
class PoissonDistributionTest {
  @SeededTest
  void canComputeProbability(RandomSeed seed) {
    final UniformRandomProvider rng = RngFactory.create(seed.get());
    final DoubleDoubleBiPredicate tol = Predicates.doublesAreRelativelyClose(1e-12);

    final PoissonDistribution fpd = new PoissonDistribution(1);
    for (int i = 1; i <= 100; i++) {
      final double mean = rng.nextDouble() * i;
      final org.apache.commons.statistics.distribution.PoissonDistribution pd =
          createReferencePoissonDistribution(mean);
      fpd.setMean(mean);
      final int lower = (int) Math.floor(Math.max(0, mean - 5 * Math.sqrt(mean)));
      final int upper = (int) Math.ceil((mean + 5 * Math.sqrt(mean)));
      for (int x = lower; x <= upper; x++) {
        TestAssertions.assertTest(pd.probability(x), fpd.probability(x), tol);
      }
    }
  }

  private static org.apache.commons.statistics.distribution.PoissonDistribution
      createReferencePoissonDistribution(final double mean) {
    return org.apache.commons.statistics.distribution.PoissonDistribution.of(mean);
  }

  @Test
  void testMeanProperty() {
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
  void testProbabilityEdgeCases() {
    final double mean = 1.23;
    final PoissonDistribution fpd = new PoissonDistribution(mean);
    Assertions.assertEquals(0, fpd.probability(-1));
    Assertions.assertEquals(0, fpd.probability(Integer.MAX_VALUE));
    Assertions.assertEquals(Math.exp(-mean), fpd.probability(0));
  }

  @Test
  void testLogProbabilityEdgeCases() {
    final double mean = 1.23;
    final PoissonDistribution fpd = new PoissonDistribution(mean);
    Assertions.assertEquals(Double.NEGATIVE_INFINITY, fpd.logProbability(-1));
    Assertions.assertEquals(-mean, fpd.logProbability(0));
  }

  @Test
  void testLogProbabilityAtExtremeValue() {
    final double mean = Double.MIN_VALUE;
    final PoissonDistribution fpd = new PoissonDistribution(mean);
    Assertions.assertEquals(Double.NEGATIVE_INFINITY, fpd.logProbability(Integer.MAX_VALUE - 1));
  }

  @Test
  void testInverseCumulativeProbabilityUpperLimitAtExtremeValue() {
    Assertions.assertEquals(Integer.MAX_VALUE, (int) Math.floor(Integer.MAX_VALUE + 0.5));
    Assertions.assertEquals(Integer.MAX_VALUE, (int) Math.floor(Double.MAX_VALUE));
    Assertions.assertEquals(Integer.MAX_VALUE, (int) Math.floor(Double.POSITIVE_INFINITY));
  }
}
