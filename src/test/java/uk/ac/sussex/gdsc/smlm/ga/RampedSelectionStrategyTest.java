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

package uk.ac.sussex.gdsc.smlm.ga;

import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.utils.BaseTimingTask;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.TimingResult;
import uk.ac.sussex.gdsc.test.utils.TimingService;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.logging.Level;
import java.util.logging.Logger;

@SuppressWarnings({"javadoc"})
public class RampedSelectionStrategyTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(RampedSelectionStrategyTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  @Test
  public void canSearchUsingActualKey() {
    final long[] sum = RampedSelectionStrategy.createRampedSum(10);

    for (int i = 0; i < sum.length - 1; i++) {
      final long key = sum[i];
      final int j = RampedSelectionStrategy.search(sum, key);
      Assertions.assertEquals(i + 1, j);
    }
  }

  @Test
  public void canBinarySearchUsingActualKey() {
    final long[] sum = RampedSelectionStrategy.createRampedSum(10);

    for (int i = 0; i < sum.length - 1; i++) {
      final long key = sum[i];
      final int j = RampedSelectionStrategy.binarySearch(sum, key);
      Assertions.assertEquals(i + 1, j);
    }
  }

  @Test
  public void canSearchUsingNotActualKey() {
    final long[] sum = RampedSelectionStrategy.createRampedSum(10);

    for (int i = 0; i < sum.length; i++) {
      final long key = sum[i] - 1;
      final int j = RampedSelectionStrategy.search(sum, key);
      Assertions.assertEquals(i, j);
    }
  }

  @Test
  public void canBinarySearchUsingNotActualKey() {
    final long[] sum = RampedSelectionStrategy.createRampedSum(10);

    for (int i = 0; i < sum.length; i++) {
      final long key = sum[i] - 1;
      final int j = RampedSelectionStrategy.binarySearch(sum, key);
      Assertions.assertEquals(i, j);
    }
  }

  @Test
  public void binarySearchEqualsSearch() {
    final long[] sum = RampedSelectionStrategy.createRampedSum(100);
    for (int key = (int) sum[sum.length - 1]; key-- > 0;) {
      final int i = RampedSelectionStrategy.search(sum, key);
      final int j = RampedSelectionStrategy.binarySearch(sum, key);
      Assertions.assertEquals(i, j);
    }
  }

  @SpeedTag
  @Test
  public void speedTest50() {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.LOW));
    speedTest(50, false, 10);
  }

  @SpeedTag
  @Test
  public void speedTest200() {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    speedTest(200, true, 5);
  }

  @SpeedTag
  @Test
  public void speedTest1000() {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));
    speedTest(1000, true, 2);
  }

  // Too slow for common use
  @SpeedTag
  @Test
  public void speedTest5000() {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.HIGH));
    speedTest(5000, true, 1);
  }

  private static void speedTest(final int size, boolean faster, int runs) {
    final long[] sum = RampedSelectionStrategy.createRampedSum(size);

    final TimingService ts = new TimingService(runs);

    ts.execute(new BaseTimingTask("search" + size) {
      @Override
      public Object getData(int index) {
        return sum;
      }

      @Override
      public Object run(Object data) {
        for (int key = (int) sum[sum.length - 1]; key-- > 0;) {
          RampedSelectionStrategy.search(sum, key);
        }
        return null;
      }

      @Override
      public int getSize() {
        return 1;
      }
    });

    ts.execute(new BaseTimingTask("binarySearch" + size) {
      @Override
      public Object getData(int index) {
        return sum[index];
      }

      @Override
      public Object run(Object data) {
        for (int key = (int) sum[sum.length - 1]; key-- > 0;) {
          RampedSelectionStrategy.binarySearch(sum, key);
        }
        return null;
      }

      @Override
      public int getSize() {
        return 1;
      }
    });

    final int n = ts.repeat();
    ts.repeat(n);
    if (logger.isLoggable(Level.INFO)) {
      logger.info(ts.getReport());
    }

    final TimingResult slow = ts.get((faster) ? ts.getSize() - 2 : ts.getSize() - 1);
    final TimingResult fast = ts.get((faster) ? ts.getSize() - 1 : ts.getSize() - 2);
    logger.log(TestLogUtils.getTimingRecord(slow, fast));
  }
}
