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
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestSettings;
import uk.ac.sussex.gdsc.test.utils.TimingService;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeAll;

import java.util.logging.Level;
import java.util.logging.Logger;

@SuppressWarnings({"javadoc"})
public class CustomGammaDistributionTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(CustomGammaDistributionTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  private abstract static class MyTimingTask extends BaseTimingTask {
    static final int MAX_SHAPE = 1000;
    static final int SAMPLES = 10;

    RandomSeed seed;
    UniformRandomProvider rng;
    double shape = 0.5;
    double scale = 300;

    public MyTimingTask(String name, RandomSeed seed) {
      super(name);
      this.seed = seed;
    }

    @Override
    public int getSize() {
      return 1;
    }

    @Override
    public Object getData(int index) {
      rng = RngUtils.create(seed.getSeedAsLong());
      shape = 0.5;
      return null;
    }
  }

  private static class StaticTimingTask extends MyTimingTask {
    public StaticTimingTask(RandomSeed seed) {
      super("RandomDataGenerator", seed);
    }

    @Override
    public Object run(Object data) {
      final RandomDataGenerator rdg = new RandomDataGenerator(new RandomGeneratorAdapter(rng));
      final double[] e = new double[MAX_SHAPE * SAMPLES];
      for (int i = 0, k = 0; i < MAX_SHAPE; i++) {
        for (int j = 0; j < SAMPLES; j++, k++) {
          e[k] = rdg.nextGamma(shape, scale);
        }
        shape += 1;
      }
      return e;
    }
  }

  private static class InstanceTimingTask extends MyTimingTask {
    public InstanceTimingTask(RandomSeed seed) {
      super("Instance", seed);
    }

    @Override
    public Object run(Object data) {
      final CustomGammaDistribution dist =
          new CustomGammaDistribution(new RandomGeneratorAdapter(rng), 1, scale);
      final double[] e = new double[MAX_SHAPE * SAMPLES];
      for (int i = 0, k = 0; i < MAX_SHAPE; i++) {
        dist.setShape(shape);
        for (int j = 0; j < SAMPLES; j++, k++) {
          e[k] = dist.sample();
        }
        shape += 1;
      }
      return e;
    }
  }

  @SeededTest
  public void canCreateSamples(RandomSeed seed) {
    final StaticTimingTask t1 = new StaticTimingTask(seed);
    t1.getData(0);
    final double[] e = (double[]) t1.run(null);

    final InstanceTimingTask t2 = new InstanceTimingTask(seed);
    t2.getData(0);
    final double[] o = (double[]) t2.run(null);

    Assertions.assertArrayEquals(e, o);
  }

  @SpeedTag
  @SeededTest
  public void customDistributionIsFaster(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final TimingService ts = new TimingService(5);
    ts.execute(new StaticTimingTask(seed));
    ts.execute(new InstanceTimingTask(seed));

    final int size = ts.getSize();
    ts.repeat(size);
    if (logger.isLoggable(Level.INFO)) {
      logger.info(ts.getReport(size));
    }

    Assertions.assertTrue(ts.get(-1).getMean() < ts.get(-2).getMean());
  }
}
