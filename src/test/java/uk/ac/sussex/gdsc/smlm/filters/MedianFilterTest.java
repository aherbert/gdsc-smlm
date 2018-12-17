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
import org.junit.jupiter.api.Assumptions;

import java.util.ArrayList;

@SuppressWarnings({"javadoc"})
public class MedianFilterTest extends AbstractFilterTest {
  private static int InternalITER3 = 200;
  private static int InternalITER = 20;
  private static int ITER3 = 100;
  private static int ITER = 10;
  private static FloatFloatBiPredicate equality = TestHelper.floatsAreClose(1e-3, 0);

  private static void floatArrayEquals(float[] data1, float[] data2,
      @SuppressWarnings("unused") int boxSize, String format, Object... args) {
    try {
      // TestAssertions.assertArrayTest(data1, data2, TestHelper.almostEqualFloats(boxSize, 0) *
      // boxSize * 1e-3);
      TestAssertions.assertArrayTest(data1, data2, equality);
    } catch (final AssertionError ex) {
      throw new AssertionError(String.format(format, args), ex);
    }
  }

  @SeededTest
  public void
      floatBlockMedianNxNInternalAndRollingMedianNxNInternalReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final MedianFilter filter = new MedianFilter();
    for (final int width : primes) {
      for (final int height : primes) {
        for (final int boxSize : boxSizes) {
          floatCompareBlockMedianNxNInternalAndRollingMedianNxNInternal(rg, filter, width, height,
              boxSize);
        }
      }
    }
  }

  private static void floatCompareBlockMedianNxNInternalAndRollingMedianNxNInternal(
      UniformRandomProvider rg, MedianFilter filter, int width, int height, int boxSize) {
    final float[] data1 = createData(rg, width, height);
    final float[] data2 = floatClone(data1);

    filter.blockMedianNxNInternal(data1, width, height, boxSize);
    filter.rollingMedianNxNInternal(data2, width, height, boxSize);

    floatArrayEquals(data1, data2, boxSize, "Internal arrays do not match: [%dx%d] @ %d", width,
        height, boxSize);
  }

  @SeededTest
  public void
      floatBlockMedian3x3InternalAndRollingMedianNxNInternalReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final MedianFilter filter = new MedianFilter();
    for (final int width : primes) {
      for (final int height : primes) {
        floatCompareBlockMedian3x3InternalAndRollingMedianNxNInternal(rg, filter, width, height);
      }
    }
  }

  private static void floatCompareBlockMedian3x3InternalAndRollingMedianNxNInternal(
      UniformRandomProvider rg, MedianFilter filter, int width, int height) {
    final float[] data1 = createData(rg, width, height);
    final float[] data2 = floatClone(data1);

    filter.blockMedian3x3Internal(data1, width, height);
    filter.rollingMedianNxNInternal(data2, width, height, 1);

    floatArrayEquals(data1, data2, 1, "Internal arrays do not match: [%dx%d] @ %d", width, height,
        1);
  }

  @SpeedTag
  @SeededTest
  public void floatBlockMedianNxNInternalIsFasterThanRollingMedianNxNInternal(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final MedianFilter filter = new MedianFilter();

    final ArrayList<float[]> dataSet = getSpeedData(seed, InternalITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.blockMedianNxNInternal(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0],
        boxSizes[0]);
    filter.rollingMedianNxNInternal(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0],
        boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : speedPrimes) {
        for (final int height : speedPrimes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.blockMedianNxNInternal(data, width, height, boxSize);
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
      for (final int width : speedPrimes) {
        for (final int height : speedPrimes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.rollingMedianNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "float rollingMedianNxNInternal [%dx%d] @ %d : %d => "
                + "blockMedianNxNInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not
            // faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("float rollingMedianNxNInternal " + boxSize,
          boxSlowTotal, "blockMedianNxNInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("float rollingMedianNxNInternal", slowTotal,
        "blockMedianNxNInternal", fastTotal));
  }

  @SeededTest
  public void
      floatBlockMedian3x3InternalAndBlockMedianNxNInternalReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final MedianFilter filter = new MedianFilter();
    for (final int width : primes) {
      for (final int height : primes) {
        floatCompareBlockMedian3x3InternalAndBlockMedianNxNInternal(rg, filter, width, height);
      }
    }
  }

  private static void floatCompareBlockMedian3x3InternalAndBlockMedianNxNInternal(
      UniformRandomProvider rg, MedianFilter filter, int width, int height) {
    final float[] data1 = createData(rg, width, height);
    final float[] data2 = floatClone(data1);

    filter.blockMedian3x3Internal(data1, width, height);
    filter.blockMedianNxNInternal(data2, width, height, 1);

    floatArrayEquals(data1, data2, 1, "Internal arrays do not match: [%dx%d] @ %d", width, height,
        1);
  }

  @SpeedTag
  @SeededTest
  public void floatBlockMedian3x3InternalIsFasterThanBlockMedianNxNInternal(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final MedianFilter filter = new MedianFilter();

    final ArrayList<float[]> dataSet = getSpeedData(seed, InternalITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.blockMedianNxNInternal(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0], 1);
    filter.blockMedian3x3Internal(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0]);

    for (final int width : speedPrimes) {
      for (final int height : speedPrimes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.blockMedian3x3Internal(data, width, height);
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
    for (final int width : speedPrimes) {
      for (final int height : speedPrimes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.blockMedianNxNInternal(data, width, height, 1);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format(
              "float blockMedianNxNInternal [%dx%d] %d => "
              + "blockMedian3x3Internal %d = %.2fx", width,
              height, time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not
          // faster: [%dx%d] %d > %d", width, height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("float blockMedianNxNInternal", slowTotal,
        "blockMedian3x3Internal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void floatBlockMedian3x3InternalIsFasterThanRollingMedian3x3Internal(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final MedianFilter filter = new MedianFilter();

    final ArrayList<float[]> dataSet = getSpeedData(seed, InternalITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingMedian3x3Internal(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0]);
    filter.blockMedian3x3Internal(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0]);

    for (final int width : speedPrimes) {
      for (final int height : speedPrimes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.blockMedian3x3Internal(data, width, height);
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
    for (final int width : speedPrimes) {
      for (final int height : speedPrimes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.rollingMedian3x3Internal(data, width, height);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format(
              "float rollingMedian3x3Internal [%dx%d] %d => blockMedian3x3Internal %d = %.2fx",
              width, height, time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not
          // faster: [%dx%d] %d > %d", width, height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("float rollingMedian3x3Internal", slowTotal,
        "blockMedian3x3Internal", fastTotal));
  }

  @SeededTest
  public void
      floatRollingMedian3x3InternalAndRollingMedianNxNInternalReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final MedianFilter filter = new MedianFilter();
    for (final int width : primes) {
      for (final int height : primes) {
        floatCompareRollingMedian3x3InternalAndRollingMedianNxNInternal(rg, filter, width, height);
      }
    }
  }

  private static void floatCompareRollingMedian3x3InternalAndRollingMedianNxNInternal(
      UniformRandomProvider rg, MedianFilter filter, int width, int height) {
    final float[] data1 = createData(rg, width, height);
    final float[] data2 = floatClone(data1);

    filter.rollingMedian3x3Internal(data1, width, height);
    filter.rollingMedianNxNInternal(data2, width, height, 1);

    floatArrayEquals(data1, data2, 1, "Internal arrays do not match: [%dx%d] @ %d", width, height,
        1);
  }

  @SpeedTag
  @SeededTest
  public void floatRollingMedian3x3InternalIsFasterThanRollingMedianNxNInternal(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final MedianFilter filter = new MedianFilter();

    final ArrayList<float[]> dataSet = getSpeedData(seed, InternalITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingMedian3x3Internal(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0]);
    filter.rollingMedianNxNInternal(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0], 1);

    for (final int width : speedPrimes) {
      for (final int height : speedPrimes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.rollingMedian3x3Internal(data, width, height);
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
    for (final int width : speedPrimes) {
      for (final int height : speedPrimes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.rollingMedianNxNInternal(data, width, height, 1);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format(
              "float rollingMedianNxNInternal [%dx%d] %d => rollingMedian3x3Internal %d = %.2fx",
              width, height, time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not
          // faster: [%dx%d] %d > %d", width, height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("float rollingMedianNxNInternal", slowTotal,
        "rollingMedian3x3Internal", fastTotal));
  }

  @SeededTest
  public void floatBlockMedianNxNAndRollingMedianNxNReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final MedianFilter filter = new MedianFilter();
    for (final int width : primes) {
      for (final int height : primes) {
        for (final int boxSize : boxSizes) {
          floatCompareBlockMedianNxNAndRollingMedianNxN(rg, filter, width, height, boxSize);
        }
      }
    }
  }

  private static void floatCompareBlockMedianNxNAndRollingMedianNxN(UniformRandomProvider rg,
      MedianFilter filter, int width, int height, int boxSize) {
    final float[] data1 = createData(rg, width, height);
    final float[] data2 = floatClone(data1);

    filter.blockMedianNxN(data1, width, height, boxSize);
    filter.rollingMedianNxN(data2, width, height, boxSize);

    floatArrayEquals(data1, data2, boxSize, "Arrays do not match: [%dx%d] @ %d", width, height,
        boxSize);
  }

  @SpeedTag
  @SeededTest
  public void floatBlockMedianInternalNxNIsFasterThanBlockMedianNxN(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final MedianFilter filter = new MedianFilter();

    final ArrayList<float[]> dataSet = getSpeedData(seed, ITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.blockMedianNxNInternal(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0],
        boxSizes[0]);
    filter.blockMedianNxN(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0], boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : speedPrimes) {
        for (final int height : speedPrimes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.blockMedianNxNInternal(data, width, height, boxSize);
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
      for (final int width : speedPrimes) {
        for (final int height : speedPrimes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.blockMedianNxN(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "float blockMedianNxN [%dx%d] @ %d : %d => blockMedianNxNInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not
            // faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("float blockMedianNxN " + boxSize, boxSlowTotal,
          "blockMedianNxNInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("float blockMedianNxN", slowTotal,
        "blockMedianNxNInternal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void floatBlockMedianNxNIsFasterThanRollingMedianNxN(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final MedianFilter filter = new MedianFilter();

    final ArrayList<float[]> dataSet = getSpeedData(seed, ITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.blockMedianNxN(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0], boxSizes[0]);
    filter.rollingMedianNxN(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0],
        boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : speedPrimes) {
        for (final int height : speedPrimes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.blockMedianNxN(data, width, height, boxSize);
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
      for (final int width : speedPrimes) {
        for (final int height : speedPrimes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.rollingMedianNxN(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "float rollingMedianNxN [%dx%d] @ %d : %d => blockMedianNxN %d = %.2fx", width,
                height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not
            // faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
            // rollingTime, time), rollingTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("float rollingMedianNxN " + boxSize,
          boxSlowTotal, "blockMedianNxN", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("float rollingMedianNxN", slowTotal, "blockMedianNxN",
        fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void floatRollingMedianInternalNxNIsFasterThanRollingMedianNxN(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final MedianFilter filter = new MedianFilter();

    final ArrayList<float[]> dataSet = getSpeedData(seed, ITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingMedianNxNInternal(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0],
        boxSizes[0]);
    filter.rollingMedianNxN(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0],
        boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : speedPrimes) {
        for (final int height : speedPrimes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.rollingMedianNxNInternal(data, width, height, boxSize);
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
      for (final int width : speedPrimes) {
        for (final int height : speedPrimes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.rollingMedianNxN(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "float rollingMedianNxN [%dx%d] @ %d : %d => rollingMedianNxNInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not
            // faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("float rollingMedianNxN " + boxSize,
          boxSlowTotal, "rollingMedianNxNInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("float rollingMedianNxN", slowTotal,
        "rollingMedianNxNInternal", fastTotal));
  }

  @SeededTest
  public void floatBlockMedian3x3AndBlockMedianNxNReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final MedianFilter filter = new MedianFilter();
    for (final int width : primes) {
      for (final int height : primes) {
        floatCompareBlockMedian3x3AndBlockMedianNxN(rg, filter, width, height);
      }
    }
  }

  private static void floatCompareBlockMedian3x3AndBlockMedianNxN(UniformRandomProvider rg,
      MedianFilter filter, int width, int height) {
    final float[] data1 = createData(rg, width, height);
    final float[] data2 = floatClone(data1);

    filter.blockMedian3x3(data1, width, height);
    filter.blockMedianNxN(data2, width, height, 1);

    floatArrayEquals(data1, data2, 1, "Arrays do not match: [%dx%d] @ %d", width, height, 1);
  }

  @SpeedTag
  @SeededTest
  public void floatBlockMedian3x3IsFasterThanBlockMedianNxN(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final MedianFilter filter = new MedianFilter();

    final ArrayList<float[]> dataSet = getSpeedData(seed, ITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.blockMedianNxN(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0], 1);
    filter.blockMedian3x3(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0]);

    for (final int width : speedPrimes) {
      for (final int height : speedPrimes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.blockMedian3x3(data, width, height);
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
    for (final int width : speedPrimes) {
      for (final int height : speedPrimes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.blockMedianNxN(data, width, height, 1);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(
              () -> String.format("float blockMedianNxN [%dx%d] %d => blockMedian3x3 %d = %.2fx",
                  width, height, time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not
          // faster: [%dx%d] %d > %d", width, height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("float blockMedianNxN", slowTotal, "blockMedian3x3",
        fastTotal));
  }

  @SeededTest
  public void floatRollingMedian3x3AndRollingMedianNxNReturnSameResult(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final MedianFilter filter = new MedianFilter();

    for (final int width : primes) {
      for (final int height : primes) {
        floatCompareRollingMedian3x3AndRollingMedianNxN(rg, filter, width, height);
      }
    }
  }

  private static void floatCompareRollingMedian3x3AndRollingMedianNxN(UniformRandomProvider rg,
      MedianFilter filter, int width, int height) {
    final float[] data1 = createData(rg, width, height);
    final float[] data2 = floatClone(data1);

    filter.rollingMedian3x3(data1, width, height);
    filter.rollingMedianNxN(data2, width, height, 1);

    floatArrayEquals(data1, data2, 1, "Arrays do not match: [%dx%d] @ %d", width, height, 1);
  }

  @SpeedTag
  @SeededTest
  public void floatRollingMedian3x3IsFasterThanRollingMedianNxN(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final MedianFilter filter = new MedianFilter();

    final ArrayList<float[]> dataSet = getSpeedData(seed, ITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingMedianNxN(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0], 1);
    filter.rollingMedian3x3(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0]);

    for (final int width : speedPrimes) {
      for (final int height : speedPrimes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.rollingMedian3x3(data, width, height);
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
    for (final int width : speedPrimes) {
      for (final int height : speedPrimes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.rollingMedianNxN(data, width, height, 1);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(() -> String.format(
              "float rollingMedianNxN [%dx%d] %d => rollingMedian3x3 %d = %.2fx", width, height,
              time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not
          // faster: [%dx%d] %d > %d", width, height,
          // rollingBlockTime, time), rollingBlockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("float rollingMedianNxN", slowTotal, "rollingMedian3x3",
        fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void floatRollingMedian3x3IsFasterThanBlockMedian3x3(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final MedianFilter filter = new MedianFilter();

    final ArrayList<float[]> dataSet = getSpeedData(seed, ITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingMedian3x3(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0]);
    filter.blockMedian3x3(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0]);

    for (final int width : speedPrimes) {
      for (final int height : speedPrimes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.rollingMedian3x3(data, width, height);
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
    for (final int width : speedPrimes) {
      for (final int height : speedPrimes) {
        final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
        for (final float[] data : dataSet) {
          dataSet2.add(floatClone(data));
        }

        final long start = System.nanoTime();
        for (final float[] data : dataSet2) {
          filter.blockMedian3x3(data, width, height);
        }
        final long time = System.nanoTime() - start;

        final long fastTime = fastTimes.get(index++);
        slowTotal += time;
        fastTotal += fastTime;
        boxSlowTotal += time;
        boxFastTotal += fastTime;
        if (debug) {
          logger.fine(
              () -> String.format("float blockMedian3x3 [%dx%d] %d => rollingMedian3x3 %d = %.2fx",
                  width, height, time, fastTime, speedUpFactor(time, fastTime)));
          // if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not
          // faster: [%dx%d] %d > %d", width, height,
          // blockTime, time), blockTime < time);
        }
      }
    }
    logger.log(TestLogUtils.getTimingRecord("float blockMedian3x3", slowTotal, "rollingMedian3x3",
        fastTotal));
  }
}
