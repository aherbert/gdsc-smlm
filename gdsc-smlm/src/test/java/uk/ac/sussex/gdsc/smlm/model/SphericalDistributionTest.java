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

import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestSettings;

@SuppressWarnings({"javadoc"})
class SphericalDistributionTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(SphericalDistributionTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  @SeededTest
  void canSampleUsingTransformationMethod(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeed());
    final double radius = 10 + rg.nextDouble() * 10;
    final SphericalDistribution dist = new SphericalDistribution(radius, rg);
    dist.setUseRejectionMethod(false);
    for (int i = 100; i-- > 0;) {
      dist.next();
    }
  }

  @SeededTest
  void canSampleUsingRejectionMethod(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeed());
    final double radius = 10 + rg.nextDouble() * 10;
    final SphericalDistribution dist = new SphericalDistribution(radius, rg);
    dist.setUseRejectionMethod(true);
    for (int i = 100; i-- > 0;) {
      dist.next();
    }
  }

  @SeededTest
  void rejectionMethodIsFasterThanTransformationMethod(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final UniformRandomProvider rg = RngUtils.create(seed.getSeed());
    final double radius = 10 + rg.nextDouble() * 10;
    final SphericalDistribution dist = new SphericalDistribution(radius, rg);
    dist.setUseRejectionMethod(false);
    for (int i = 100; i-- > 0;) {
      dist.next();
    }
    dist.setUseRejectionMethod(true);
    for (int i = 100; i-- > 0;) {
      dist.next();
    }

    dist.setUseRejectionMethod(false);
    final long time1 = getRunTime(dist);
    dist.setUseRejectionMethod(true);
    final long time2 = getRunTime(dist);
    Assertions.assertTrue(time1 > time2,
        () -> String.format("Rejection = %d, Transformation = %d", time2, time1));
    logger.log(
        TestLogUtils.getRecord(Level.INFO, "Rejection = %d, Transformation = %d", time2, time1));
  }

  private static long getRunTime(SphericalDistribution dist) {
    final long start = System.nanoTime();
    for (int i = 1000000; i-- > 0;) {
      dist.next();
    }
    return System.nanoTime() - start;
  }
}
