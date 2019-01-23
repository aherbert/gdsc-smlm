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

package uk.ac.sussex.gdsc.smlm.math3.distribution;

import uk.ac.sussex.gdsc.core.utils.RandomGeneratorAdapter;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.BaseTimingTask;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TimingResult;
import uk.ac.sussex.gdsc.test.utils.TimingService;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;

import java.util.logging.Level;
import java.util.logging.Logger;

@SuppressWarnings({"javadoc"})
public class CustomPoissonDistributionTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(CustomPoissonDistributionTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  private abstract static class MyTimingTask extends BaseTimingTask {
    static final int SAMPLES = 10;
    RandomSeed seed;
    UniformRandomProvider rng;
    double mean;
    double min;
    int numberOfMeans;

    public MyTimingTask(String name, RandomSeed seed, double min, double max) {
      super(String.format("%s %.1f - %.1f", name, min, max));
      this.seed = seed;
      this.min = min;
      mean = min;
      numberOfMeans = 0;
      while (mean < max) {
        numberOfMeans++;
        mean += 1;
      }
    }

    @Override
    public int getSize() {
      return 1;
    }

    @Override
    public Object getData(int index) {
      rng = RngUtils.create(seed.getSeedAsLong());
      mean = min;
      return null;
    }
  }

  private static class StaticTimingTask extends MyTimingTask {
    public StaticTimingTask(RandomSeed seed, double min, double max) {
      super("RandomDataGenerator", seed, min, max);
    }

    @Override
    public Object run(Object data) {
      final RandomDataGenerator rdg = new RandomDataGenerator(new RandomGeneratorAdapter(rng));
      final long[] e = new long[numberOfMeans * SAMPLES];
      for (int i = 0, k = 0; i < numberOfMeans; i++) {
        for (int j = 0; j < SAMPLES; j++, k++) {
          e[k] = rdg.nextPoisson(mean);
        }
        mean += 1;
      }
      return e;
    }
  }

  private static class InstanceTimingTask extends MyTimingTask {
    public InstanceTimingTask(RandomSeed seed, double min, double max) {
      super("CustomPoissonDistribution", seed, min, max);
    }

    @Override
    public Object run(Object data) {
      final CustomPoissonDistribution dist =
          new CustomPoissonDistribution(new RandomGeneratorAdapter(rng), 1);
      final long[] e = new long[numberOfMeans * SAMPLES];
      for (int i = 0, k = 0; i < numberOfMeans; i++) {
        dist.setMean(mean);
        for (int j = 0; j < SAMPLES; j++, k++) {
          e[k] = dist.sample();
        }
        mean += 1;
      }
      return e;
    }
  }

  @SeededTest
  public void canCreateSamples(RandomSeed seed) {
    final StaticTimingTask t1 = new StaticTimingTask(seed, 0.5, 60);
    t1.getData(0);
    final long[] e = (long[]) t1.run(null);

    final InstanceTimingTask t2 = new InstanceTimingTask(seed, 0.5, 60);
    t2.getData(0);
    final long[] o = (long[]) t2.run(null);

    Assertions.assertArrayEquals(e, o);
  }

  @SpeedTag
  @SeededTest
  public void customDistributionIsFasterWithTinyMean(RandomSeed seed) {
    final TimingService ts = new TimingService(5);
    ts.execute(new StaticTimingTask(seed, 0.5, 10));
    ts.execute(new InstanceTimingTask(seed, 0.5, 10));

    final int size = ts.getSize();
    ts.repeat(size);
    if (logger.isLoggable(Level.INFO)) {
      logger.info(ts.getReport(size));
    }

    // Assertions.assertTrue(ts.get(-1).getMean() < ts.get(-2).getMean());
    final TimingResult slow = ts.get(-2);
    final TimingResult fast = ts.get(-1);
    logger.log(TestLogUtils.getTimingRecord(slow, fast));
  }

  @SpeedTag
  @SeededTest
  public void customDistributionIsFasterWithSmallMean(RandomSeed seed) {
    final TimingService ts = new TimingService(5);
    ts.execute(new StaticTimingTask(seed, 10, 38));
    ts.execute(new InstanceTimingTask(seed, 10, 38));

    final int size = ts.getSize();
    ts.repeat(size);
    if (logger.isLoggable(Level.INFO)) {
      logger.info(ts.getReport(size));
    }

    // Assertions.assertTrue(ts.get(-1).getMean() < ts.get(-2).getMean());
    final TimingResult slow = ts.get(-2);
    final TimingResult fast = ts.get(-1);
    logger.log(TestLogUtils.getTimingRecord(slow, fast));
  }

  @SpeedTag
  @SeededTest
  public void customDistributionIsFasterWithBigMean(RandomSeed seed) {
    // When the mean is above 40 the PoissonDistribution switches to a different
    // sampling method and this is so slow that the speed increase from using
    // the instance class is negligible. However test it is still faster. If this fails
    // then Apache commons may have changed their implementation and the custom
    // class should be updated.

    final TimingService ts = new TimingService(5);
    ts.execute(new StaticTimingTask(seed, 40.5, 60));
    ts.execute(new InstanceTimingTask(seed, 40.5, 60));

    final int size = ts.getSize();
    ts.repeat(size);
    if (logger.isLoggable(Level.INFO)) {
      logger.info(ts.getReport(size));
    }

    // Assertions.assertTrue(ts.get(-1).getMean() < ts.get(-2).getMean());
    final TimingResult slow = ts.get(-2);
    final TimingResult fast = ts.get(-1);
    logger.log(TestLogUtils.getTimingRecord(slow, fast));
  }
}
