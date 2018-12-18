package uk.ac.sussex.gdsc.smlm.filters;

import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.FloatFloatBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestSettings;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Assumptions;

import java.util.ArrayList;

@SuppressWarnings({"deprecation", "javadoc"})
public class SumFilterTest extends AbstractFilterTest {
  private static int InternalITER3 = 500;
  private static int InternalITER = 50;
  private static int ITER3 = 200;
  private static int ITER = 20;
  private static FloatFloatBiPredicate equality = TestHelper.floatsAreClose(1e-4, 0);

  /**
   * Check the float arrays are equal, else fail with a formatted message.
   *
   * @param data1 the data 1
   * @param data2 the data 2
   * @param boxSize the box size used for filtering
   * @param format the format
   * @param args the args
   */
  private static void floatArrayEquals(float[] data1, float[] data2, int boxSize, String format,
      Object... args) {
    TestAssertions.assertArrayTest(data1, data2, equality, () -> String.format(format, args));
  }

  /**
   * Check the int arrays are equal, else fail with a formatted message.
   *
   * @param data1 the data 1
   * @param data2 the data 2
   * @param boxSize the box size used for filtering
   * @param format the format
   * @param args the args
   */
  private static void intArrayEquals(int[] data1, int[] data2, int boxSize, String format,
      Object... args) {
    Assertions.assertArrayEquals(data1, data2, () -> String.format(format, args));
  }

  private static float[] floatCreateData(UniformRandomProvider rg, int width, int height) {
    return createData(rg, width, height);
  }

  private static int[] intCreateData(UniformRandomProvider rg, int width, int height) {
    return createIntData(rg, width, height);
  }

  private ArrayList<float[]> floatCreateSpeedData(RandomSeed seed, int iter) {
    return getSpeedData(seed, iter);
  }

  private ArrayList<int[]> intCreateSpeedData(RandomSeed seed, int iter) {
    return getIntSpeedData(seed, iter);
  }

  // COPY CODE FROM HERE...
  @SeededTest
  public void
      floatBlockSumNxNInternalAndRollingBlockSumNxNInternalReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final SumFilter filter = new SumFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        for (final int boxSize : boxSizes) {
          floatCompareBlockSumNxNInternalAndRollingBlockSumNxNInternal(rg, filter, width, height,
              boxSize);
        }
      }
    }
  }

  private static void floatCompareBlockSumNxNInternalAndRollingBlockSumNxNInternal(
      UniformRandomProvider rg, SumFilter filter, int width, int height, int boxSize) {
    final float[] data1 = floatCreateData(rg, width, height);
    final float[] data2 = floatClone(data1);

    filter.blockSumNxNInternal(data1, width, height, boxSize);
    filter.rollingBlockSumNxNInternal(data2, width, height, boxSize);

    floatArrayEquals(data1, data2, boxSize, "Internal arrays do not match: [%dx%d] @ %d", width,
        height, boxSize);
  }

  @SeededTest
  public void
      floatBlockSumNxNInternalAndStripedBlockSumNxNInternalReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final SumFilter filter = new SumFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        for (final int boxSize : boxSizes) {
          floatCompareBlockSumNxNInternalAndStripedBlockSumNxNInternal(rg, filter, width, height,
              boxSize);
        }
      }
    }
  }

  private static void floatCompareBlockSumNxNInternalAndStripedBlockSumNxNInternal(
      UniformRandomProvider rg, SumFilter filter, int width, int height, int boxSize) {
    final float[] data1 = floatCreateData(rg, width, height);
    final float[] data2 = floatClone(data1);

    filter.blockSumNxNInternal(data1, width, height, boxSize);
    filter.stripedBlockSumNxNInternal(data2, width, height, boxSize);

    floatArrayEquals(data1, data2, boxSize, "Internal arrays do not match: [%dx%d] @ %d", width,
        height, boxSize);
  }

  @SeededTest
  public void
      floatBlockSum3x3InternalAndRollingBlockSumNxNInternalReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final SumFilter filter = new SumFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        floatCompareBlockSum3x3InternalAndRollingBlockSumNxNInternal(rg, filter, width, height);
      }
    }
  }

  private static void floatCompareBlockSum3x3InternalAndRollingBlockSumNxNInternal(
      UniformRandomProvider rg, SumFilter filter, int width, int height) {
    final float[] data1 = floatCreateData(rg, width, height);
    final float[] data2 = floatClone(data1);

    filter.blockSum3x3Internal(data1, width, height);
    filter.rollingBlockSumNxNInternal(data2, width, height, 1);

    floatArrayEquals(data1, data2, 1, "Internal arrays do not match: [%dx%d] @ %d", width, height,
        1);
  }

  @SeededTest
  public void
      floatRollingBlockSumNxNInternalAndRollingBlockSumNxNInternalTransposedReturnSameResult(
          RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final SumFilter filter = new SumFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        for (final int boxSize : boxSizes) {
          floatCompareRollingBlockSumNxNInternalAndRollingBlockSumNxNInternalTransposed(rg, filter,
              width, height, boxSize);
        }
      }
    }
  }

  private static void floatCompareRollingBlockSumNxNInternalAndRollingBlockSumNxNInternalTransposed(
      UniformRandomProvider rg, SumFilter filter, int width, int height, int boxSize) {
    final float[] data1 = floatCreateData(rg, width, height);
    final float[] data2 = floatClone(data1);

    filter.rollingBlockSumNxNInternal(data1, width, height, boxSize);
    filter.rollingBlockSumNxNInternalTransposed(data2, width, height, boxSize);

    floatArrayEquals(data1, data2, boxSize, "Internal arrays do not match: [%dx%d] @ %d", width,
        height, boxSize);
  }

  @SpeedTag
  @SeededTest
  public void floatRollingBlockSumNxNInternalIsFasterThanBlockSumNxNInternal(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<float[]> dataSet = floatCreateSpeedData(seed, InternalITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.blockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
    filter.rollingBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0],
        boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final int boxSize : boxSizes) {
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.blockSumNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "float blockSumNxNInternal [%dx%d] @ %d : %d => "
                    + "rollingBlockSumNxNInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS)
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d",
            // width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("float blockSumNxNInternal " + boxSize,
          boxSlowTotal, "rollingBlockSumNxNInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("float blockSumNxNInternal", slowTotal,
        "rollingBlockSumNxNInternal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void floatStripedBlockSumNxNInternalIsFasterThanBlockSumNxNInternal(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<float[]> dataSet = floatCreateSpeedData(seed, InternalITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.blockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
    filter.stripedBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0],
        boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.stripedBlockSumNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final int boxSize : boxSizes) {
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.blockSumNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "float blockSumNxNInternal [%dx%d] @ %d : %d => "
                    + "stripedBlockSumNxNInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS)
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d",
            // width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("float blockSumNxNInternal " + boxSize,
          boxSlowTotal, "stripedBlockSumNxNInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("float blockSumNxNInternal", slowTotal,
        "stripedBlockSumNxNInternal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void
      floatRollingBlockSumNxNInternalIsFasterThanStripedBlockSumNxNInternal(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<float[]> dataSet = floatCreateSpeedData(seed, InternalITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0],
        boxSizes[0]);
    filter.stripedBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0],
        boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final int boxSize : boxSizes) {
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.stripedBlockSumNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "float stripedBlockSumNxNInternal [%dx%d] @ %d : %d => "
                    + "rollingBlockSumNxNInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS)
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d",
            // width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("float stripedBlockSumNxNInternal " + boxSize,
          boxSlowTotal, "rollingBlockSumNxNInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("float stripedBlockSumNxNInternal", slowTotal,
        "rollingBlockSumNxNInternal", fastTotal));
  }

  @SeededTest
  public void floatBlockSum3x3InternalAndBlockSumNxNInternalReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final SumFilter filter = new SumFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        floatCompareBlockSum3x3InternalAndBlockSumNxNInternal(rg, filter, width, height);
      }
    }
  }

  private static void floatCompareBlockSum3x3InternalAndBlockSumNxNInternal(
      UniformRandomProvider rg, SumFilter filter, int width, int height) {
    final float[] data1 = floatCreateData(rg, width, height);
    final float[] data2 = floatClone(data1);

    filter.blockSum3x3Internal(data1, width, height);
    filter.blockSumNxNInternal(data2, width, height, 1);

    floatArrayEquals(data1, data2, 1, "Internal arrays do not match: [%dx%d] @ %d", width, height,
        1);
  }

  @SeededTest
  public void floatBlockSum3x3InternalIsFasterThanBlockSumNxNInternal(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<float[]> dataSet = floatCreateSpeedData(seed, InternalITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.blockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
    filter.blockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.blockSum3x3Internal(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.blockSumNxNInternal(data, width, height, 1);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format(
              "float blockSumNxNInternal [%dx%d] %d => blockSum3x3Internal %d = %.2fx", width,
              height, time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("float blockSumNxNInternal", slowTotal,
        "blockSum3x3Internal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void floatRollingBlockSum3x3InternalIsFasterThanBlockSum3x3Internal(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<float[]> dataSet = floatCreateSpeedData(seed, InternalITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingBlockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);
    filter.blockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.rollingBlockSum3x3Internal(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.blockSum3x3Internal(data, width, height);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format(
              "float blockSum3x3Internal [%dx%d] %d => rollingBlockSum3x3Internal %d = %.2fx",
              width, height, time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("float blockSum3x3Internal", slowTotal,
        "rollingBlockSum3x3Internal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void floatStripedBlockSum3x3InternalIsFasterThanBlockSum3x3Internal(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<float[]> dataSet = floatCreateSpeedData(seed, InternalITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.stripedBlockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);
    filter.blockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.stripedBlockSum3x3Internal(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.blockSum3x3Internal(data, width, height);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format(
              "float blockSum3x3Internal [%dx%d] %d => stripedBlockSum3x3Internal %d = %.2fx",
              width, height, time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("float blockSum3x3Internal", slowTotal,
        "stripedBlockSum3x3Internal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void
      floatRollingBlockSum3x3InternalIsFasterThanStripedBlockSum3x3Internal(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<float[]> dataSet = floatCreateSpeedData(seed, InternalITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingBlockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);
    filter.stripedBlockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.rollingBlockSum3x3Internal(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.stripedBlockSum3x3Internal(data, width, height);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format(
              "float stripedBlockSum3x3Internal [%dx%d] %d => "
                  + "rollingBlockSum3x3Internal %d = %.2fx",
              width, height, time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("float stripedBlockSum3x3Internal", slowTotal,
        "rollingBlockSum3x3Internal", fastTotal));
  }

  @SeededTest
  public void floatRollingBlockSum3x3InternalAndRollingBlockSumNxNInternalReturnSameResult(
      RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final SumFilter filter = new SumFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        floatCompareRollingBlockSum3x3InternalAndRollingBlockSumNxNInternal(rg, filter, width,
            height);
      }
    }
  }

  private static void floatCompareRollingBlockSum3x3InternalAndRollingBlockSumNxNInternal(
      UniformRandomProvider rg, SumFilter filter, int width, int height) {
    final float[] data1 = floatCreateData(rg, width, height);
    final float[] data2 = floatClone(data1);

    filter.rollingBlockSum3x3Internal(data1, width, height);
    filter.rollingBlockSumNxNInternal(data2, width, height, 1);

    floatArrayEquals(data1, data2, 1, "Internal arrays do not match: [%dx%d] @ %d", width, height,
        1);
  }

  @SpeedTag
  @SeededTest
  public void
      floatRollingBlockSum3x3InternalIsFasterThanRollingBlockSumNxNInternal(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<float[]> dataSet = floatCreateSpeedData(seed, InternalITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingBlockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);
    filter.rollingBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], 1);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.rollingBlockSum3x3Internal(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.rollingBlockSumNxNInternal(data, width, height, 1);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format(
              "float rollingBlockSumNxNInternal [%dx%d] %d => "
                  + "rollingBlockSum3x3Internal %d = %.2fx",
              width, height, time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("float rollingBlockSumNxNInternal", slowTotal,
        "rollingBlockSum3x3Internal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void
      floatStripedBlockSum3x3InternalIsFasterThanStripedBlockSumNxNInternal(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<float[]> dataSet = floatCreateSpeedData(seed, InternalITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.stripedBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
    filter.stripedBlockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.stripedBlockSum3x3Internal(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.stripedBlockSumNxNInternal(data, width, height, 1);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format(
              "float stripedBlockSumNxNInternal [%dx%d] %d => "
                  + "stripedBlockSum3x3Internal %d = %.2fx",
              width, height, time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("float stripedBlockSumNxNInternal", slowTotal,
        "stripedBlockSum3x3Internal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void floatRollingBlockSumNxNInternalIsFasterThanRollingBlockSumNxNInternalTransposed(
      RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<float[]> dataSet = floatCreateSpeedData(seed, InternalITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0],
        boxSizes[0]);
    filter.rollingBlockSumNxNInternalTransposed(floatClone(dataSet.get(0)), primes[0], primes[0],
        boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final int boxSize : boxSizes) {
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.rollingBlockSumNxNInternalTransposed(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "float rollingBlockSumNxNInternalTransposed [%dx%d] @ %d : %d => "
                    + "rollingBlockSumNxNInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS)
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d",
            // width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      if (debug) {
        final long boxSlowTotal2 = boxSlowTotal;
        final long boxFastTotal2 = boxFastTotal;
        logger.fine(() -> String.format(
            "float rollingBlockSumNxNInternalTransposed %d : %d => "
                + "rollingBlockSumNxNInternal %d = %.2fx",
            boxSize, boxSlowTotal2, boxFastTotal2, speedUpFactor(boxSlowTotal2, boxFastTotal2)));
      }
      // if (ExtraAssertions.assert_SPEED_TESTS)
      // Assertions.assertTrue(String.format("Not faster: Block %d : %d > %d",
      // boxSize, boxFastTotal, boxSlowTotal),
      // boxFastTotal < boxSlowTotal);
    }
    logger.log(TestLogUtils.getTimingRecord("float rollingBlockSumNxNInternalTransposed", slowTotal,
        "rollingBlockSumNxNInternal", fastTotal));
  }

  @SeededTest
  public void floatBlockSumNxNAndStripedBlockSumNxNReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final SumFilter filter = new SumFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        for (final int boxSize : boxSizes) {
          floatCompareBlockSumNxNAndStripedBlockSumNxN(rg, filter, width, height, boxSize);
        }
      }
    }
  }

  private static void floatCompareBlockSumNxNAndStripedBlockSumNxN(UniformRandomProvider rg,
      SumFilter filter, int width, int height, int boxSize) {
    final float[] data1 = floatCreateData(rg, width, height);
    final float[] data2 = floatClone(data1);

    filter.blockSumNxN(data1, width, height, boxSize);
    filter.rollingBlockSumNxN(data2, width, height, boxSize);

    floatArrayEquals(data1, data2, boxSize, "Arrays do not match: [%dx%d] @ %d", width, height,
        boxSize);
  }

  @SeededTest
  public void floatBlockSumNxNAndRollingBlockSumNxNReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final SumFilter filter = new SumFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        for (final int boxSize : boxSizes) {
          floatCompareBlockSumNxNAndRollingBlockSumNxN(rg, filter, width, height, boxSize);
        }
      }
    }
  }

  private static void floatCompareBlockSumNxNAndRollingBlockSumNxN(UniformRandomProvider rg,
      SumFilter filter, int width, int height, int boxSize) {
    final float[] data1 = floatCreateData(rg, width, height);
    final float[] data2 = floatClone(data1);

    filter.blockSumNxN(data1, width, height, boxSize);
    filter.rollingBlockSumNxN(data2, width, height, boxSize);

    floatArrayEquals(data1, data2, boxSize, "Arrays do not match: [%dx%d] @ %d", width, height,
        boxSize);
  }

  @SpeedTag
  @SeededTest
  public void floatBlockSumInternalNxNIsFasterThanBlockSumNxN(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<float[]> dataSet = floatCreateSpeedData(seed, ITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.blockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
    filter.blockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.blockSumNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final int boxSize : boxSizes) {
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.blockSumNxN(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "float blockSumNxN [%dx%d] @ %d : %d => blockSumNxNInternal %d = %.2fx", width,
                height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS)
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d",
            // width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("float blockSumNxN " + boxSize, boxSlowTotal,
          "blockSumNxNInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("float blockSumNxN", slowTotal, "blockSumNxNInternal",
        fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void floatStripedBlockSumNxNIsFasterThanBlockSumNxN(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<float[]> dataSet = floatCreateSpeedData(seed, ITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.blockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
    filter.stripedBlockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.stripedBlockSumNxN(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final int boxSize : boxSizes) {
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.blockSumNxN(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "float blockSumNxN [%dx%d] @ %d : %d => stripedBlockSumNxN %d = %.2fx", width,
                height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS)
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d",
            // width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("float blockSumNxN " + boxSize, boxSlowTotal,
          "stripedBlockSumNxN", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("float blockSumNxN", slowTotal, "stripedBlockSumNxN",
        fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void floatStripedBlockSumInternalNxNIsFasterThanStripedBlockSumNxN(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<float[]> dataSet = floatCreateSpeedData(seed, ITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.stripedBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0],
        boxSizes[0]);
    filter.stripedBlockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.stripedBlockSumNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final int boxSize : boxSizes) {
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.stripedBlockSumNxN(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "float stripedBlockSumNxN [%dx%d] @ %d : %d => "
                    + "stripedBlockSumNxNInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS)
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d",
            // width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("float stripedBlockSumNxN " + boxSize,
          boxSlowTotal, "stripedBlockSumNxNInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("float stripedBlockSumNxN", slowTotal,
        "stripedBlockSumNxNInternal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void floatRollingBlockSumNxNIsFasterThanBlockSumNxN(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<float[]> dataSet = floatCreateSpeedData(seed, ITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.blockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
    filter.rollingBlockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.rollingBlockSumNxN(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final int boxSize : boxSizes) {
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.blockSumNxN(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "float blockSumNxN [%dx%d] @ %d : %d => rollingBlockSumNxN %d = %.2fx", width,
                height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS)
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d",
            // width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("float blockSumNxN " + boxSize, boxSlowTotal,
          "rollingBlockSumNxN", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("float blockSumNxN", slowTotal, "rollingBlockSumNxN",
        fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void floatRollingBlockSumInternalNxNIsFasterThanRollingBlockSumNxN(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<float[]> dataSet = floatCreateSpeedData(seed, ITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0],
        boxSizes[0]);
    filter.rollingBlockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final int boxSize : boxSizes) {
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.rollingBlockSumNxN(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "float rollingBlockSumNxN [%dx%d] @ %d : %d => "
                    + "rollingBlockSumNxNInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS)
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d",
            // width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("float rollingBlockSumNxN " + boxSize,
          boxSlowTotal, "rollingBlockSumNxNInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("float rollingBlockSumNxN", slowTotal,
        "rollingBlockSumNxNInternal", fastTotal));
  }

  @SeededTest
  public void floatBlockSum3x3AndBlockSumNxNReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final SumFilter filter = new SumFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        floatCompareBlockSum3x3AndBlockSumNxN(rg, filter, width, height);
      }
    }
  }

  private static void floatCompareBlockSum3x3AndBlockSumNxN(UniformRandomProvider rg,
      SumFilter filter, int width, int height) {
    final float[] data1 = floatCreateData(rg, width, height);
    final float[] data2 = floatClone(data1);

    filter.blockSum3x3(data1, width, height);
    filter.blockSumNxN(data2, width, height, 1);

    floatArrayEquals(data1, data2, 1, "Arrays do not match: [%dx%d] @ %d", width, height, 1);
  }

  @SeededTest
  public void floatBlockSum3x3IsFasterThanBlockSumNxN(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<float[]> dataSet = floatCreateSpeedData(seed, ITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.blockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
    filter.blockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.blockSum3x3(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.blockSumNxN(data, width, height, 1);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format("float blockSumNxN [%dx%d] %d => blockSum3x3 %d = %.2fx",
              width, height, time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(
        TestLogUtils.getTimingRecord("float blockSumNxN", slowTotal, "blockSum3x3", fastTotal));
  }

  @SeededTest
  public void floatStripedBlockSum3x3AndStripedBlockSumNxNReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final SumFilter filter = new SumFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        floatCompareStripedBlockSum3x3AndStripedBlockSumNxN(rg, filter, width, height);
      }
    }
  }

  private static void floatCompareStripedBlockSum3x3AndStripedBlockSumNxN(UniformRandomProvider rg,
      SumFilter filter, int width, int height) {
    final float[] data1 = floatCreateData(rg, width, height);
    final float[] data2 = floatClone(data1);

    filter.stripedBlockSum3x3(data1, width, height);
    filter.stripedBlockSumNxN(data2, width, height, 1);

    floatArrayEquals(data1, data2, 1, "Arrays do not match: [%dx%d] @ %d", width, height, 1);
  }

  @SpeedTag
  @SeededTest
  public void floatStripedBlockSum3x3IsFasterThanStripedBlockSumNxN(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<float[]> dataSet = floatCreateSpeedData(seed, ITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.stripedBlockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
    filter.stripedBlockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.stripedBlockSum3x3(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.stripedBlockSumNxN(data, width, height, 1);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format(
              "float stripedBlockSumNxN [%dx%d] %d => stripedBlockSum3x3 %d = %.2fx", width, height,
              time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // stripedBlockTime, time), stripedBlockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("float stripedBlockSumNxN", slowTotal,
        "stripedBlockSum3x3", fastTotal));
  }

  @SeededTest
  public void floatRollingBlockSum3x3AndRollingBlockSumNxNReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final SumFilter filter = new SumFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        floatCompareRollingBlockSum3x3AndRollingBlockSumNxN(rg, filter, width, height);
      }
    }
  }

  private static void floatCompareRollingBlockSum3x3AndRollingBlockSumNxN(UniformRandomProvider rg,
      SumFilter filter, int width, int height) {
    final float[] data1 = floatCreateData(rg, width, height);
    final float[] data2 = floatClone(data1);

    filter.rollingBlockSum3x3(data1, width, height);
    filter.rollingBlockSumNxN(data2, width, height, 1);

    floatArrayEquals(data1, data2, 1, "Arrays do not match: [%dx%d] @ %d", width, height, 1);
  }

  @SpeedTag
  @SeededTest
  public void floatRollingBlockSum3x3IsFasterThanRollingBlockSumNxN(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<float[]> dataSet = floatCreateSpeedData(seed, ITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingBlockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
    filter.rollingBlockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.rollingBlockSum3x3(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.rollingBlockSumNxN(data, width, height, 1);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format(
              "float rollingBlockSumNxN [%dx%d] %d => rollingBlockSum3x3 %d = %.2fx", width, height,
              time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // rollingBlockTime, time), rollingBlockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("float rollingBlockSumNxN", slowTotal,
        "rollingBlockSum3x3", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void floatRollingBlockSum3x3IsFasterThanBlockSum3x3(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<float[]> dataSet = floatCreateSpeedData(seed, ITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingBlockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);
    filter.blockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.rollingBlockSum3x3(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.blockSum3x3(data, width, height);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(
              () -> String.format("float blockSum3x3 [%dx%d] %d => rollingBlockSum3x3 %d = %.2fx",
                  width, height, time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("float blockSum3x3", slowTotal, "rollingBlockSum3x3",
        fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void floatStripedBlockSum3x3IsFasterThanBlockSum3x3(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<float[]> dataSet = floatCreateSpeedData(seed, ITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.stripedBlockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);
    filter.blockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.stripedBlockSum3x3(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.blockSum3x3(data, width, height);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(
              () -> String.format("float blockSum3x3 [%dx%d] %d => stripedBlockSum3x3 %d = %.2fx",
                  width, height, time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("float blockSum3x3", slowTotal, "stripedBlockSum3x3",
        fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void floatRollingBlockSum3x3IsFasterThanStripedBlockSum3x3(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<float[]> dataSet = floatCreateSpeedData(seed, ITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingBlockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);
    filter.stripedBlockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.rollingBlockSum3x3(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.stripedBlockSum3x3(data, width, height);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format(
              "float stripedBlockSum3x3 [%dx%d] %d => rollingBlockSum3x3 %d = %.2fx", width, height,
              time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("float stripedBlockSum3x3", slowTotal,
        "rollingBlockSum3x3", fastTotal));
  }

  @SeededTest
  public void intBlockSumNxNInternalAndRollingBlockSumNxNInternalReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final SumFilter filter = new SumFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        for (final int boxSize : boxSizes) {
          intCompareBlockSumNxNInternalAndRollingBlockSumNxNInternal(rg, filter, width, height,
              boxSize);
        }
      }
    }
  }

  private static void intCompareBlockSumNxNInternalAndRollingBlockSumNxNInternal(
      UniformRandomProvider rg, SumFilter filter, int width, int height, int boxSize) {
    final int[] data1 = intCreateData(rg, width, height);
    final int[] data2 = intClone(data1);

    filter.blockSumNxNInternal(data1, width, height, boxSize);
    filter.rollingBlockSumNxNInternal(data2, width, height, boxSize);

    intArrayEquals(data1, data2, boxSize, "Internal arrays do not match: [%dx%d] @ %d", width,
        height, boxSize);
  }

  @SeededTest
  public void intBlockSumNxNInternalAndStripedBlockSumNxNInternalReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final SumFilter filter = new SumFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        for (final int boxSize : boxSizes) {
          intCompareBlockSumNxNInternalAndStripedBlockSumNxNInternal(rg, filter, width, height,
              boxSize);
        }
      }
    }
  }

  private static void intCompareBlockSumNxNInternalAndStripedBlockSumNxNInternal(
      UniformRandomProvider rg, SumFilter filter, int width, int height, int boxSize) {
    final int[] data1 = intCreateData(rg, width, height);
    final int[] data2 = intClone(data1);

    filter.blockSumNxNInternal(data1, width, height, boxSize);
    filter.stripedBlockSumNxNInternal(data2, width, height, boxSize);

    intArrayEquals(data1, data2, boxSize, "Internal arrays do not match: [%dx%d] @ %d", width,
        height, boxSize);
  }

  @SeededTest
  public void intBlockSum3x3InternalAndRollingBlockSumNxNInternalReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final SumFilter filter = new SumFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        intCompareBlockSum3x3InternalAndRollingBlockSumNxNInternal(rg, filter, width, height);
      }
    }
  }

  private static void intCompareBlockSum3x3InternalAndRollingBlockSumNxNInternal(
      UniformRandomProvider rg, SumFilter filter, int width, int height) {
    final int[] data1 = intCreateData(rg, width, height);
    final int[] data2 = intClone(data1);

    filter.blockSum3x3Internal(data1, width, height);
    filter.rollingBlockSumNxNInternal(data2, width, height, 1);

    intArrayEquals(data1, data2, 1, "Internal arrays do not match: [%dx%d] @ %d", width, height, 1);
  }

  @SeededTest
  public void intRollingBlockSumNxNInternalAndRollingBlockSumNxNInternalTransposedReturnSameResult(
      RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final SumFilter filter = new SumFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        for (final int boxSize : boxSizes) {
          intCompareRollingBlockSumNxNInternalAndRollingBlockSumNxNInternalTransposed(rg, filter,
              width, height, boxSize);
        }
      }
    }
  }

  private static void intCompareRollingBlockSumNxNInternalAndRollingBlockSumNxNInternalTransposed(
      UniformRandomProvider rg, SumFilter filter, int width, int height, int boxSize) {
    final int[] data1 = intCreateData(rg, width, height);
    final int[] data2 = intClone(data1);

    filter.rollingBlockSumNxNInternal(data1, width, height, boxSize);
    filter.rollingBlockSumNxNInternalTransposed(data2, width, height, boxSize);

    intArrayEquals(data1, data2, boxSize, "Internal arrays do not match: [%dx%d] @ %d", width,
        height, boxSize);
  }

  @SpeedTag
  @SeededTest
  public void intRollingBlockSumNxNInternalIsFasterThanBlockSumNxNInternal(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<int[]> dataSet = intCreateSpeedData(seed, InternalITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.blockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
    filter.rollingBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final int[] data : dataSet) {
            dataSet2.add(intClone(data));
          }

          final long start = System.nanoTime();
          for (final int[] data : dataSet2) {
            filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final int boxSize : boxSizes) {
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final int[] data : dataSet) {
            dataSet2.add(intClone(data));
          }

          final long start = System.nanoTime();
          for (final int[] data : dataSet2) {
            filter.blockSumNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "int blockSumNxNInternal [%dx%d] @ %d : %d => "
                    + "rollingBlockSumNxNInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS)
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d",
            // width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("int blockSumNxNInternal " + boxSize,
          boxSlowTotal, "rollingBlockSumNxNInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("int blockSumNxNInternal", slowTotal,
        "rollingBlockSumNxNInternal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void intStripedBlockSumNxNInternalIsFasterThanBlockSumNxNInternal(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<int[]> dataSet = intCreateSpeedData(seed, InternalITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.blockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
    filter.stripedBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final int[] data : dataSet) {
            dataSet2.add(intClone(data));
          }

          final long start = System.nanoTime();
          for (final int[] data : dataSet2) {
            filter.stripedBlockSumNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final int boxSize : boxSizes) {
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final int[] data : dataSet) {
            dataSet2.add(intClone(data));
          }

          final long start = System.nanoTime();
          for (final int[] data : dataSet2) {
            filter.blockSumNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "int blockSumNxNInternal [%dx%d] @ %d : %d => "
                    + "stripedBlockSumNxNInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS)
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d",
            // width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("int blockSumNxNInternal " + boxSize,
          boxSlowTotal, "stripedBlockSumNxNInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("int blockSumNxNInternal", slowTotal,
        "stripedBlockSumNxNInternal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void intRollingBlockSumNxNInternalIsFasterThanStripedBlockSumNxNInternal(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<int[]> dataSet = intCreateSpeedData(seed, InternalITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
    filter.stripedBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final int[] data : dataSet) {
            dataSet2.add(intClone(data));
          }

          final long start = System.nanoTime();
          for (final int[] data : dataSet2) {
            filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final int boxSize : boxSizes) {
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final int[] data : dataSet) {
            dataSet2.add(intClone(data));
          }

          final long start = System.nanoTime();
          for (final int[] data : dataSet2) {
            filter.stripedBlockSumNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "int stripedBlockSumNxNInternal [%dx%d] @ %d : %d => "
                    + "rollingBlockSumNxNInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS)
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d",
            // width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("int stripedBlockSumNxNInternal " + boxSize,
          boxSlowTotal, "rollingBlockSumNxNInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("int stripedBlockSumNxNInternal", slowTotal,
        "rollingBlockSumNxNInternal", fastTotal));
  }

  @SeededTest
  public void intBlockSum3x3InternalAndBlockSumNxNInternalReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final SumFilter filter = new SumFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        intCompareBlockSum3x3InternalAndBlockSumNxNInternal(rg, filter, width, height);
      }
    }
  }

  private static void intCompareBlockSum3x3InternalAndBlockSumNxNInternal(UniformRandomProvider rg,
      SumFilter filter, int width, int height) {
    final int[] data1 = intCreateData(rg, width, height);
    final int[] data2 = intClone(data1);

    filter.blockSum3x3Internal(data1, width, height);
    filter.blockSumNxNInternal(data2, width, height, 1);

    intArrayEquals(data1, data2, 1, "Internal arrays do not match: [%dx%d] @ %d", width, height, 1);
  }

  @SeededTest
  public void intBlockSum3x3InternalIsFasterThanBlockSumNxNInternal(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<int[]> dataSet = intCreateSpeedData(seed, InternalITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.blockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], 1);
    filter.blockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.blockSum3x3Internal(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.blockSumNxNInternal(data, width, height, 1);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format(
              "int blockSumNxNInternal [%dx%d] %d => blockSum3x3Internal %d = %.2fx", width, height,
              time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("int blockSumNxNInternal", slowTotal,
        "blockSum3x3Internal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void intRollingBlockSum3x3InternalIsFasterThanBlockSum3x3Internal(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<int[]> dataSet = intCreateSpeedData(seed, InternalITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingBlockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);
    filter.blockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.rollingBlockSum3x3Internal(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.blockSum3x3Internal(data, width, height);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format(
              "int blockSum3x3Internal [%dx%d] %d => rollingBlockSum3x3Internal %d = %.2fx", width,
              height, time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("int blockSum3x3Internal", slowTotal,
        "rollingBlockSum3x3Internal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void intStripedBlockSum3x3InternalIsFasterThanBlockSum3x3Internal(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<int[]> dataSet = intCreateSpeedData(seed, InternalITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.stripedBlockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);
    filter.blockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.stripedBlockSum3x3Internal(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.blockSum3x3Internal(data, width, height);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format(
              "int blockSum3x3Internal [%dx%d] %d => stripedBlockSum3x3Internal %d = %.2fx", width,
              height, time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("int blockSum3x3Internal", slowTotal,
        "stripedBlockSum3x3Internal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void intRollingBlockSum3x3InternalIsFasterThanStripedBlockSum3x3Internal(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<int[]> dataSet = intCreateSpeedData(seed, InternalITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingBlockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);
    filter.stripedBlockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.rollingBlockSum3x3Internal(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.stripedBlockSum3x3Internal(data, width, height);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format(
              "int stripedBlockSum3x3Internal [%dx%d] %d => rollingBlockSum3x3Internal %d = %.2fx",
              width, height, time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("int stripedBlockSum3x3Internal", slowTotal,
        "rollingBlockSum3x3Internal", fastTotal));
  }

  @SeededTest
  public void
      intRollingBlockSum3x3InternalAndRollingBlockSumNxNInternalReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final SumFilter filter = new SumFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        intCompareRollingBlockSum3x3InternalAndRollingBlockSumNxNInternal(rg, filter, width,
            height);
      }
    }
  }

  private static void intCompareRollingBlockSum3x3InternalAndRollingBlockSumNxNInternal(
      UniformRandomProvider rg, SumFilter filter, int width, int height) {
    final int[] data1 = intCreateData(rg, width, height);
    final int[] data2 = intClone(data1);

    filter.rollingBlockSum3x3Internal(data1, width, height);
    filter.rollingBlockSumNxNInternal(data2, width, height, 1);

    intArrayEquals(data1, data2, 1, "Internal arrays do not match: [%dx%d] @ %d", width, height, 1);
  }

  @SpeedTag
  @SeededTest
  public void intRollingBlockSum3x3InternalIsFasterThanRollingBlockSumNxNInternal(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<int[]> dataSet = intCreateSpeedData(seed, InternalITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingBlockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);
    filter.rollingBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], 1);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.rollingBlockSum3x3Internal(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.rollingBlockSumNxNInternal(data, width, height, 1);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format(
              "int rollingBlockSumNxNInternal [%dx%d] %d => rollingBlockSum3x3Internal %d = %.2fx",
              width, height, time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("int rollingBlockSumNxNInternal", slowTotal,
        "rollingBlockSum3x3Internal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void intStripedBlockSum3x3InternalIsFasterThanStripedBlockSumNxNInternal(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<int[]> dataSet = intCreateSpeedData(seed, InternalITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.stripedBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], 1);
    filter.stripedBlockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.stripedBlockSum3x3Internal(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.stripedBlockSumNxNInternal(data, width, height, 1);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format(
              "int stripedBlockSumNxNInternal [%dx%d] %d => stripedBlockSum3x3Internal %d = %.2fx",
              width, height, time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("int stripedBlockSumNxNInternal", slowTotal,
        "stripedBlockSum3x3Internal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void intRollingBlockSumNxNInternalIsFasterThanRollingBlockSumNxNInternalTransposed(
      RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<int[]> dataSet = intCreateSpeedData(seed, InternalITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
    filter.rollingBlockSumNxNInternalTransposed(intClone(dataSet.get(0)), primes[0], primes[0],
        boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final int[] data : dataSet) {
            dataSet2.add(intClone(data));
          }

          final long start = System.nanoTime();
          for (final int[] data : dataSet2) {
            filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final int boxSize : boxSizes) {
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final int[] data : dataSet) {
            dataSet2.add(intClone(data));
          }

          final long start = System.nanoTime();
          for (final int[] data : dataSet2) {
            filter.rollingBlockSumNxNInternalTransposed(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "int rollingBlockSumNxNInternalTransposed [%dx%d] @ %d : %d => "
                    + "rollingBlockSumNxNInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS)
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d",
            // width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      if (debug) {
        final long boxSlowTotal2 = boxSlowTotal;
        final long boxFastTotal2 = boxFastTotal;
        logger.fine(() -> String.format(
            "int rollingBlockSumNxNInternalTransposed %d : %d => "
                + "rollingBlockSumNxNInternal %d = %.2fx",
            boxSize, boxSlowTotal2, boxFastTotal2, speedUpFactor(boxSlowTotal2, boxFastTotal2)));
      }
      // if (ExtraAssertions.assert_SPEED_TESTS)
      // Assertions.assertTrue(String.format("Not faster: Block %d : %d > %d",
      // boxSize, boxFastTotal, boxSlowTotal),
      // boxFastTotal < boxSlowTotal);
    }
    logger.log(TestLogUtils.getTimingRecord("int rollingBlockSumNxNInternalTransposed", slowTotal,
        "rollingBlockSumNxNInternal", fastTotal));
  }

  @SeededTest
  public void intBlockSumNxNAndStripedBlockSumNxNReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final SumFilter filter = new SumFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        for (final int boxSize : boxSizes) {
          intCompareBlockSumNxNAndStripedBlockSumNxN(rg, filter, width, height, boxSize);
        }
      }
    }
  }

  private static void intCompareBlockSumNxNAndStripedBlockSumNxN(UniformRandomProvider rg,
      SumFilter filter, int width, int height, int boxSize) {
    final int[] data1 = intCreateData(rg, width, height);
    final int[] data2 = intClone(data1);

    filter.blockSumNxN(data1, width, height, boxSize);
    filter.rollingBlockSumNxN(data2, width, height, boxSize);

    intArrayEquals(data1, data2, boxSize, "Arrays do not match: [%dx%d] @ %d", width, height,
        boxSize);
  }

  @SeededTest
  public void intBlockSumNxNAndRollingBlockSumNxNReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final SumFilter filter = new SumFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        for (final int boxSize : boxSizes) {
          intCompareBlockSumNxNAndRollingBlockSumNxN(rg, filter, width, height, boxSize);
        }
      }
    }
  }

  private static void intCompareBlockSumNxNAndRollingBlockSumNxN(UniformRandomProvider rg,
      SumFilter filter, int width, int height, int boxSize) {
    final int[] data1 = intCreateData(rg, width, height);
    final int[] data2 = intClone(data1);

    filter.blockSumNxN(data1, width, height, boxSize);
    filter.rollingBlockSumNxN(data2, width, height, boxSize);

    intArrayEquals(data1, data2, boxSize, "Arrays do not match: [%dx%d] @ %d", width, height,
        boxSize);
  }

  @SpeedTag
  @SeededTest
  public void intBlockSumInternalNxNIsFasterThanBlockSumNxN(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<int[]> dataSet = intCreateSpeedData(seed, ITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.blockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
    filter.blockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final int[] data : dataSet) {
            dataSet2.add(intClone(data));
          }

          final long start = System.nanoTime();
          for (final int[] data : dataSet2) {
            filter.blockSumNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final int boxSize : boxSizes) {
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final int[] data : dataSet) {
            dataSet2.add(intClone(data));
          }

          final long start = System.nanoTime();
          for (final int[] data : dataSet2) {
            filter.blockSumNxN(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "int blockSumNxN [%dx%d] @ %d : %d => blockSumNxNInternal %d = %.2fx", width,
                height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS)
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d",
            // width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("int blockSumNxN " + boxSize, boxSlowTotal,
          "blockSumNxNInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("int blockSumNxN", slowTotal, "blockSumNxNInternal",
        fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void intStripedBlockSumNxNIsFasterThanBlockSumNxN(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<int[]> dataSet = intCreateSpeedData(seed, ITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.blockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
    filter.stripedBlockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final int[] data : dataSet) {
            dataSet2.add(intClone(data));
          }

          final long start = System.nanoTime();
          for (final int[] data : dataSet2) {
            filter.stripedBlockSumNxN(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final int boxSize : boxSizes) {
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final int[] data : dataSet) {
            dataSet2.add(intClone(data));
          }

          final long start = System.nanoTime();
          for (final int[] data : dataSet2) {
            filter.blockSumNxN(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "int blockSumNxN [%dx%d] @ %d : %d => stripedBlockSumNxN %d = %.2fx", width, height,
                boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS)
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d",
            // width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("int blockSumNxN " + boxSize, boxSlowTotal,
          "stripedBlockSumNxN", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("int blockSumNxN", slowTotal, "stripedBlockSumNxN",
        fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void intStripedBlockSumInternalNxNIsFasterThanStripedBlockSumNxN(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<int[]> dataSet = intCreateSpeedData(seed, ITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.stripedBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
    filter.stripedBlockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final int[] data : dataSet) {
            dataSet2.add(intClone(data));
          }

          final long start = System.nanoTime();
          for (final int[] data : dataSet2) {
            filter.stripedBlockSumNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final int boxSize : boxSizes) {
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final int[] data : dataSet) {
            dataSet2.add(intClone(data));
          }

          final long start = System.nanoTime();
          for (final int[] data : dataSet2) {
            filter.stripedBlockSumNxN(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "int stripedBlockSumNxN [%dx%d] @ %d : %d => stripedBlockSumNxNInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS)
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d",
            // width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("int stripedBlockSumNxN " + boxSize,
          boxSlowTotal, "stripedBlockSumNxNInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("int stripedBlockSumNxN", slowTotal,
        "stripedBlockSumNxNInternal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void intRollingBlockSumNxNIsFasterThanBlockSumNxN(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<int[]> dataSet = intCreateSpeedData(seed, ITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.blockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
    filter.rollingBlockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final int[] data : dataSet) {
            dataSet2.add(intClone(data));
          }

          final long start = System.nanoTime();
          for (final int[] data : dataSet2) {
            filter.rollingBlockSumNxN(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final int boxSize : boxSizes) {
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final int[] data : dataSet) {
            dataSet2.add(intClone(data));
          }

          final long start = System.nanoTime();
          for (final int[] data : dataSet2) {
            filter.blockSumNxN(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "int blockSumNxN [%dx%d] @ %d : %d => rollingBlockSumNxN %d = %.2fx", width, height,
                boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS)
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d",
            // width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("int blockSumNxN " + boxSize, boxSlowTotal,
          "rollingBlockSumNxN", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("int blockSumNxN", slowTotal, "rollingBlockSumNxN",
        fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void intRollingBlockSumInternalNxNIsFasterThanRollingBlockSumNxN(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<int[]> dataSet = intCreateSpeedData(seed, ITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
    filter.rollingBlockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final int[] data : dataSet) {
            dataSet2.add(intClone(data));
          }

          final long start = System.nanoTime();
          for (final int[] data : dataSet2) {
            filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final int boxSize : boxSizes) {
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final int[] data : dataSet) {
            dataSet2.add(intClone(data));
          }

          final long start = System.nanoTime();
          for (final int[] data : dataSet2) {
            filter.rollingBlockSumNxN(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "int rollingBlockSumNxN [%dx%d] @ %d : %d => rollingBlockSumNxNInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS)
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d",
            // width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("int rollingBlockSumNxN " + boxSize,
          boxSlowTotal, "rollingBlockSumNxNInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("int rollingBlockSumNxN", slowTotal,
        "rollingBlockSumNxNInternal", fastTotal));
  }

  @SeededTest
  public void intBlockSum3x3AndBlockSumNxNReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final SumFilter filter = new SumFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        intCompareBlockSum3x3AndBlockSumNxN(rg, filter, width, height);
      }
    }
  }

  private static void intCompareBlockSum3x3AndBlockSumNxN(UniformRandomProvider rg,
      SumFilter filter, int width, int height) {
    final int[] data1 = intCreateData(rg, width, height);
    final int[] data2 = intClone(data1);

    filter.blockSum3x3(data1, width, height);
    filter.blockSumNxN(data2, width, height, 1);

    intArrayEquals(data1, data2, 1, "Arrays do not match: [%dx%d] @ %d", width, height, 1);
  }

  @SeededTest
  public void intBlockSum3x3IsFasterThanBlockSumNxN(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<int[]> dataSet = intCreateSpeedData(seed, ITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.blockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], 1);
    filter.blockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.blockSum3x3(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.blockSumNxN(data, width, height, 1);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format("int blockSumNxN [%dx%d] %d => blockSum3x3 %d = %.2fx",
              width, height, time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger
        .log(TestLogUtils.getTimingRecord("int blockSumNxN", slowTotal, "blockSum3x3", fastTotal));
  }

  @SeededTest
  public void intStripedBlockSum3x3AndStripedBlockSumNxNReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final SumFilter filter = new SumFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        intCompareStripedBlockSum3x3AndStripedBlockSumNxN(rg, filter, width, height);
      }
    }
  }

  private static void intCompareStripedBlockSum3x3AndStripedBlockSumNxN(UniformRandomProvider rg,
      SumFilter filter, int width, int height) {
    final int[] data1 = intCreateData(rg, width, height);
    final int[] data2 = intClone(data1);

    filter.stripedBlockSum3x3(data1, width, height);
    filter.stripedBlockSumNxN(data2, width, height, 1);

    intArrayEquals(data1, data2, 1, "Arrays do not match: [%dx%d] @ %d", width, height, 1);
  }

  @SpeedTag
  @SeededTest
  public void intStripedBlockSum3x3IsFasterThanStripedBlockSumNxN(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<int[]> dataSet = intCreateSpeedData(seed, ITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.stripedBlockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], 1);
    filter.stripedBlockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.stripedBlockSum3x3(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.stripedBlockSumNxN(data, width, height, 1);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format(
              "int stripedBlockSumNxN [%dx%d] %d => stripedBlockSum3x3 %d = %.2fx", width, height,
              time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // stripedBlockTime, time), stripedBlockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("int stripedBlockSumNxN", slowTotal,
        "stripedBlockSum3x3", fastTotal));
  }

  @SeededTest
  public void intRollingBlockSum3x3AndRollingBlockSumNxNReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final SumFilter filter = new SumFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        intCompareRollingBlockSum3x3AndRollingBlockSumNxN(rg, filter, width, height);
      }
    }
  }

  private static void intCompareRollingBlockSum3x3AndRollingBlockSumNxN(UniformRandomProvider rg,
      SumFilter filter, int width, int height) {
    final int[] data1 = intCreateData(rg, width, height);
    final int[] data2 = intClone(data1);

    filter.rollingBlockSum3x3(data1, width, height);
    filter.rollingBlockSumNxN(data2, width, height, 1);

    intArrayEquals(data1, data2, 1, "Arrays do not match: [%dx%d] @ %d", width, height, 1);
  }

  @SpeedTag
  @SeededTest
  public void intRollingBlockSum3x3IsFasterThanRollingBlockSumNxN(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<int[]> dataSet = intCreateSpeedData(seed, ITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingBlockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], 1);
    filter.rollingBlockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.rollingBlockSum3x3(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.rollingBlockSumNxN(data, width, height, 1);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format(
              "int rollingBlockSumNxN [%dx%d] %d => rollingBlockSum3x3 %d = %.2fx", width, height,
              time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // rollingBlockTime, time), rollingBlockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("int rollingBlockSumNxN", slowTotal,
        "rollingBlockSum3x3", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void intRollingBlockSum3x3IsFasterThanBlockSum3x3(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<int[]> dataSet = intCreateSpeedData(seed, ITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingBlockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);
    filter.blockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.rollingBlockSum3x3(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.blockSum3x3(data, width, height);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(
              () -> String.format("int blockSum3x3 [%dx%d] %d => rollingBlockSum3x3 %d = %.2fx",
                  width, height, time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("int blockSum3x3", slowTotal, "rollingBlockSum3x3",
        fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void intStripedBlockSum3x3IsFasterThanBlockSum3x3(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<int[]> dataSet = intCreateSpeedData(seed, ITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.stripedBlockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);
    filter.blockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.stripedBlockSum3x3(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.blockSum3x3(data, width, height);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(
              () -> String.format("int blockSum3x3 [%dx%d] %d => stripedBlockSum3x3 %d = %.2fx",
                  width, height, time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("int blockSum3x3", slowTotal, "stripedBlockSum3x3",
        fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void intRollingBlockSum3x3IsFasterThanStripedBlockSum3x3(RandomSeed seed) {
    // These test a deprecated filter
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.VERY_HIGH));

    final SumFilter filter = new SumFilter();

    final ArrayList<int[]> dataSet = intCreateSpeedData(seed, ITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingBlockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);
    filter.stripedBlockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);

    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.rollingBlockSum3x3(data, width, height);
        }
        final long time = System.nanoTime() - start;
        fastTimes.add(time);
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    @SuppressWarnings("unused")
    long boxSlowTotal = 0;
    long boxFastTotal = 0;
    for (final int width : primes) {
      for (final int height : primes) {
        final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final int[] data : dataSet) {
          dataSet2.add(intClone(data));
        }

        final long start = System.nanoTime();
        for (final int[] data : dataSet2) {
          filter.stripedBlockSum3x3(data, width, height);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format(
              "int stripedBlockSum3x3 [%dx%d] %d => rollingBlockSum3x3 %d = %.2fx", width, height,
              time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS)
          // Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width,
          // height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("int stripedBlockSum3x3", slowTotal,
        "rollingBlockSum3x3", fastTotal));
  }
}
