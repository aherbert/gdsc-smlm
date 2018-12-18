package uk.ac.sussex.gdsc.smlm.filters;

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
public class AreaAverageFilterTest extends AbstractFilterTest {
  private static final int ITER = 100;
  private static final int InternalITER = 300;

  @SeededTest
  public void areaAverageUsingSumsNxNInternalIsFasterThanAreaAverageNxNInternal(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final AreaAverageFilter filter = new AreaAverageFilter();

    final ArrayList<float[]> dataSet = getSpeedData(seed, InternalITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    for (final float boxSize : fBoxSizes) {
      filter.areaAverageUsingAveragesInternal(dataSet.get(0).clone(), primes[0], primes[0],
          boxSize);
      filter.areaAverageUsingSumsInternal(dataSet.get(0).clone(), primes[0], primes[0], boxSize);
    }

    for (final float boxSize : fBoxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(data.clone());
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.areaAverageUsingSumsInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final float boxSize : fBoxSizes) {
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(data.clone());
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.areaAverageUsingAveragesInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "float areaAverageInternal [%dx%d] @ %.1f : "
                    + "%d => areaAverageUsingSumsInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not
            // faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("float areaAverageInternal " + boxSize,
          boxSlowTotal, "areaAverageUsingSumsInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("float areaAverageInternal", slowTotal,
        "areaAverageUsingSumsInternal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void stripedBlockAverageIsFasterThanAreaAverage(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final AreaAverageFilter filter = new AreaAverageFilter();
    final AverageFilter filter2 = new AverageFilter();

    final ArrayList<float[]> dataSet = getSpeedData(seed, ITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    for (final float boxSize : fBoxSizes) {
      filter.areaAverageUsingAverages(dataSet.get(0).clone(), primes[0], primes[0], boxSize);
      filter2.stripedBlockAverage(dataSet.get(0).clone(), primes[0], primes[0], boxSize);
    }

    for (final float boxSize : fBoxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(data.clone());
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter2.stripedBlockAverage(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final float boxSize : fBoxSizes) {
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(data.clone());
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.areaAverageUsingAverages(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "float areaAverageUsingAverages [%dx%d] @ %.1f : "
                    + "%d => stripedBlockAverage %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not
            // faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("float areaAverageUsingAverages " + boxSize,
          boxSlowTotal, "stripedBlockAverage", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("float areaAverageUsingAverages", slowTotal,
        "stripedBlockAverage", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void stripedBlockAverageInternalIsFasterThanAreaAverageInternal(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final AreaAverageFilter filter = new AreaAverageFilter();
    final AverageFilter filter2 = new AverageFilter();

    final ArrayList<float[]> dataSet = getSpeedData(seed, InternalITER);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    for (final float boxSize : fBoxSizes) {
      filter.areaAverageUsingAveragesInternal(dataSet.get(0).clone(), primes[0], primes[0],
          boxSize);
      filter2.stripedBlockAverageInternal(dataSet.get(0).clone(), primes[0], primes[0], boxSize);
    }

    for (final float boxSize : fBoxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(data.clone());
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter2.stripedBlockAverageInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final float boxSize : fBoxSizes) {
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
          for (final float[] data : dataSet) {
            dataSet2.add(data.clone());
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.areaAverageUsingAveragesInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "float areaAverageUsingAveragesInternal [%dx%d] @ %.1f : "
                    + "%d => stripedBlockAverageInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not
            // faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(
          TestLogUtils.getStageTimingRecord("float areaAverageUsingAveragesInternal " + boxSize,
              boxSlowTotal, "stripedBlockAverageInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("float areaAverageUsingAveragesInternal", slowTotal,
        "stripedBlockAverageInternal", fastTotal));
  }

  @SeededTest
  public void areaAverageCorrectlyInterpolatesBetweenBlocks(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final int max = 50;
    final float[] data = createData(rg, max, max);
    final AreaAverageFilter filter = new AreaAverageFilter();
    final int n = 30;
    final float[][] results = new float[n + 1][];
    final double[] w = new double[n + 1];
    int count = 0;
    for (int i = 0; i <= n; i++) {
      w[count] = i / 10.0;
      results[count] = data.clone();
      filter.areaAverageUsingAverages(results[count], max, max, w[count]);
      count++;
    }

    checkInterpolation(max, n, results, count);
  }

  private static void checkInterpolation(int max, int n, float[][] results, int count) {
    // Pick some points and see if they are monototically interpolated between integer blocks
    final int[] p = new int[] {10, 20, 30, 40};
    for (final int x : p) {
      for (final int y : p) {
        final int index = y * max + x;
        final double[] yy = new double[count];
        int i1 = 0;
        for (final float[] data1 : results) {
          yy[i1++] = data1[index];
        }

        //// Debugging
        // String title = "AreaAverage";
        // ij.gui.Plot plot = new ij.gui.Plot(title, "width", "Mean", w, yy);
        // uk.ac.sussex.gdsc.core.ij.Utils.display(title, plot);

        for (int i = 0; i < n; i += 10) {
          final boolean up = yy[i + 10] > yy[i];
          for (int j = i + 1; j < i + 10; j++) {
            if (up) {
              Assertions.assertTrue(yy[j] >= yy[j - 1]);
              Assertions.assertTrue(yy[j] <= yy[j + 1]);
            } else {
              Assertions.assertTrue(yy[j] <= yy[j - 1]);
              Assertions.assertTrue(yy[j] >= yy[j + 1]);
            }
          }
        }
      }
    }
  }

  @SeededTest
  public void areaAverageInternalCorrectlyInterpolatesBetweenBlocks(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final int max = 50;
    final float[] data = createData(rg, max, max);
    final AreaAverageFilter filter = new AreaAverageFilter();
    final int n = 30;
    final float[][] results = new float[n + 1][];
    final double[] w = new double[n + 1];
    int count = 0;
    for (int i = 0; i <= n; i++) {
      w[count] = i / 10.0;
      results[count] = data.clone();
      filter.areaAverageUsingAveragesInternal(results[count], max, max, w[count]);
      count++;
    }

    checkInterpolation(max, n, results, count);
  }

  @SeededTest
  public void areaAverageUsingSumsCorrectlyInterpolatesBetweenBlocks(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final int max = 50;
    final float[] data = createData(rg, max, max);
    final AreaAverageFilter filter = new AreaAverageFilter();
    filter.setSimpleInterpolation(false);
    final int n = 30;
    final float[][] results = new float[n + 1][];
    final double[] w = new double[n + 1];
    int count = 0;
    for (int i = 0; i <= n; i++) {
      w[count] = i / 10.0;
      results[count] = data.clone();
      filter.areaAverageUsingSums(results[count], max, max, w[count]);
      count++;
    }

    checkInterpolation(max, n, results, count);
  }

  @SeededTest
  public void areaAverageUsingSumsInternalCorrectlyInterpolatesBetweenBlocks(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final int max = 50;
    final float[] data = createData(rg, max, max);
    final AreaAverageFilter filter = new AreaAverageFilter();
    filter.setSimpleInterpolation(false);
    final int n = 30;
    final float[][] results = new float[n + 1][];
    final double[] w = new double[n + 1];
    int count = 0;
    for (int i = 0; i <= n; i++) {
      w[count] = i / 10.0;
      results[count] = data.clone();
      filter.areaAverageUsingSumsInternal(results[count], max, max, w[count]);
      count++;
    }

    checkInterpolation(max, n, results, count);
  }
}
