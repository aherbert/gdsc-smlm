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

import java.util.Arrays;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;

@SuppressWarnings({"javadoc"})
class SphericalDistributionTest {
  @SeededTest
  void canSample(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.get());
    final double radius = 10 + rg.nextDouble() * 10;
    final SphericalDistribution dist = new SphericalDistribution(radius, rg);
    for (int i = 100; i-- > 0;) {
      final double[] x = dist.next();
      Assertions.assertTrue(MathUtils.distance(x[0], x[1], x[2], 0, 0, 0) <= radius,
          () -> "Bad coords: " + Arrays.toString(x) + " radius=" + radius);
    }
  }

  @SeededTest
  void canSampleWithNoRadius(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.get());
    final double[] expected = {0, 0, 0};
    for (final double radius : new double[] {0.0, -0.0}) {
      final SphericalDistribution dist = new SphericalDistribution(radius, rg);
      for (int i = 5; i-- > 0;) {
        final double[] x = dist.next();
        Assertions.assertArrayEquals(expected, x);
      }
    }
  }

  @ParameterizedTest
  @ValueSource(doubles = {-1, Double.NaN})
  void testIllegalRadiusThrows(double radius) {
    Assertions.assertThrows(IllegalArgumentException.class,
        () -> new SphericalDistribution(radius));
  }
}
