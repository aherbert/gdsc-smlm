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

package uk.ac.sussex.gdsc.smlm.results.filter;

import java.util.logging.Logger;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import uk.ac.sussex.gdsc.core.utils.XmlUtils;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.RandomSeed;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils.TestLevel;

@SuppressWarnings({"javadoc"})
class FilterTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(FilterTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  @SeededTest
  void canCompareMultiFilter(RandomSeed seed) {
    final UniformRandomProvider UniformRandomProvider = RngUtils.create(seed.get());
    final MultiFilter f = new MultiFilter(0, 0, 0, 0, 0, 0, 0, 0, 0);
    for (int i = 1000; i-- > 0;) {
      final MultiFilter f1 =
          (MultiFilter) f.create(random(f.getNumberOfParameters(), UniformRandomProvider));
      final MultiFilter f2 =
          (MultiFilter) f.create(random(f.getNumberOfParameters(), UniformRandomProvider));
      final int e = f1.weakest((Filter) f2);
      final int o = f1.weakest(f2);
      Assertions.assertEquals(e, o);
    }
  }

  @SeededTest
  void canCompareMultiFilter2(RandomSeed seed) {
    final UniformRandomProvider UniformRandomProvider = RngUtils.create(seed.get());
    final MultiFilter2 f = new MultiFilter2(0, 0, 0, 0, 0, 0, 0, 0, 0);
    for (int i = 1000; i-- > 0;) {
      final MultiFilter2 f1 =
          (MultiFilter2) f.create(random(f.getNumberOfParameters(), UniformRandomProvider));
      final MultiFilter2 f2 =
          (MultiFilter2) f.create(random(f.getNumberOfParameters(), UniformRandomProvider));
      final int e = f1.weakest((Filter) f2);
      final int o = f1.weakest(f2);
      Assertions.assertEquals(e, o);
    }
  }

  private static double[] random(int np, UniformRandomProvider rng) {
    final double[] params = new double[np];
    while (np-- > 0) {
      params[np] = rng.nextInt(3);
    }
    return params;
  }

  @SeededTest
  void canSerialiseMultiFilter(RandomSeed seed) {
    // Check the XStream serialisation supports inheritance
    final UniformRandomProvider rng = RngUtils.create(seed.get());
    testSerialisation(new MultiFilter(0, 0, 0, 0, 0, 0, 0, 0, 0), rng);
    testSerialisation(new MultiFilter2(0, 0, 0, 0, 0, 0, 0, 0, 0), rng);
    testSerialisation(new MultiFilterCrlb(0, 0, 0, 0, 0, 0, 0, 0, 0), rng);
  }

  private static void testSerialisation(MultiFilter filter, UniformRandomProvider rng) {
    for (int i = 10; i-- > 0;) {
      final MultiFilter f1 =
          (MultiFilter) filter.create(random(filter.getNumberOfParameters(), rng));
      final String xml = f1.toXml();
      logger.log(TestLogUtils.getRecord(TestLevel.TEST_DEBUG, XmlUtils.prettyPrintXml(xml)));
      final MultiFilter f2 = (MultiFilter) Filter.fromXml(xml);
      Assertions.assertTrue(f1.getClass().equals(f2.getClass()));
      Assertions.assertEquals(f1, f2);
    }
  }
}
