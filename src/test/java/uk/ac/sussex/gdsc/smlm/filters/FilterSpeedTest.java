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

package uk.ac.sussex.gdsc.smlm.filters;

import java.util.ArrayList;
import org.junit.jupiter.api.Assumptions;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestSettings;

@SuppressWarnings({"javadoc"})
public class FilterSpeedTest extends AbstractFilterTest {
  @SpeedTag
  @SeededTest
  public void
      floatRollingBlockSumNxNInternalIsFasterThanRollingBlockMeanNxNInternal(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final BlockSumFilter filter = new BlockSumFilter();
    final BlockMeanFilter filter2 = new BlockMeanFilter();

    final int iter = 50;
    final ArrayList<float[]> dataSet = getSpeedData(seed, iter);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingBlockFilterNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0],
        boxSizes[0]);
    filter2.rollingBlockFilterNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0],
        boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(iter);
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter.rollingBlockFilterNxNInternal(data, width, height, boxSize);
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
          final ArrayList<float[]> dataSet2 = new ArrayList<>(iter);
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter2.rollingBlockFilterNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "float rollingBlockMeanNxNInternal [%dx%d] @ %d : %d => "
                    + "rollingBlockSumNxNInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width,
            // height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("float rollingBlockMeanNxNInternal " + boxSize,
          boxSlowTotal, "rollingBlockSumNxNInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("float rollingBlockMeanNxNInternal", slowTotal,
        "rollingBlockSumNxNInternal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void floatRollingBlockMeanNxNInternalIsFasterThanBlockMedianNxNInternal(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final BlockMeanFilter filter1 = new BlockMeanFilter();
    final MedianFilter filter2 = new MedianFilter();

    final int iter = 10;
    final ArrayList<float[]> dataSet = getSpeedData(seed, iter);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter1.rollingBlockFilterNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0],
        boxSizes[0]);
    filter2.blockMedianNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(iter);
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter1.rollingBlockFilterNxNInternal(data, width, height, boxSize);
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
          final ArrayList<float[]> dataSet2 = new ArrayList<>(iter);
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter2.blockMedianNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "float blockMedianNxNInternal [%dx%d] @ %d : %d => "
                    + "rollingBlockMeanNxNInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width,
            // height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("float blockMedianNxNInternal " + boxSize,
          boxSlowTotal, "rollingBlockMeanNxNInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("float blockMedianNxNInternal", slowTotal,
        "rollingBlockMeanNxNInternal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void
      floatRollingBlockMeanNxNInternalIsFasterThanRollingMedianNxNInternal(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final BlockMeanFilter filter1 = new BlockMeanFilter();
    final MedianFilter filter2 = new MedianFilter();

    final int iter = 10;
    final ArrayList<float[]> dataSet = getSpeedData(seed, iter);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter1.rollingBlockFilterNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0],
        boxSizes[0]);
    filter2.rollingMedianNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(iter);
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter1.rollingBlockFilterNxNInternal(data, width, height, boxSize);
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
          final ArrayList<float[]> dataSet2 = new ArrayList<>(iter);
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter2.rollingMedianNxNInternal(data, width, height, boxSize);
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
                    + "rollingBlockMeanNxNInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width,
            // height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("float rollingMedianNxNInternal " + boxSize,
          boxSlowTotal, "rollingBlockMeanNxNInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("float rollingMedianNxNInternal", slowTotal,
        "rollingBlockMeanNxNInternal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void floatRollingBlockMeanNxNInternalIsFasterThanGaussianNxNInternal(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final BlockMeanFilter filter1 = new BlockMeanFilter();
    final GaussianFilter filter2 = new GaussianFilter();

    final int iter = 10;
    final ArrayList<float[]> dataSet = getSpeedData(seed, iter);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter1.rollingBlockFilterNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0],
        boxSizes[0]);
    filter2.convolveInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0] / 3.0);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(iter);
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter1.rollingBlockFilterNxNInternal(data, width, height, boxSize);
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
          final ArrayList<float[]> dataSet2 = new ArrayList<>(iter);
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter2.convolveInternal(data, width, height, boxSize / 3.0);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "float convolveInternal [%dx%d] @ %d : %d => "
                    + "rollingBlockMeanNxNInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width,
            // height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("float convolveInternal " + boxSize,
          boxSlowTotal, "rollingBlockMeanNxNInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("float convolveInternal", slowTotal,
        "rollingBlockMeanNxNInternal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void floatRollingBlockMeanNxNInternalIsFasterThanAreaFilterNxNInternal(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final BlockMeanFilter filter1 = new BlockMeanFilter();
    final AreaAverageFilter filter2 = new AreaAverageFilter();

    final int iter = 10;
    final ArrayList<float[]> dataSet = getSpeedData(seed, iter);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    // filter1.rollingBlockFilterNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0],
    // boxSizes[0]);
    // filter2.areaFilterInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0] -
    // 0.05);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(iter);
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          // Initialise
          for (final float[] data : dataSet2) {
            filter1.rollingBlockFilterNxNInternal(data.clone(), width, height, boxSize);
          }
          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter1.rollingBlockFilterNxNInternal(data, width, height, boxSize);
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
          final ArrayList<float[]> dataSet2 = new ArrayList<>(iter);
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          // Initialise
          for (final float[] data : dataSet2) {
            filter2.areaAverageUsingAveragesInternal(data.clone(), width, height, boxSize - 0.05);
          }
          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter2.areaAverageUsingAveragesInternal(data, width, height, boxSize - 0.05);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "float areaFilterInternal [%dx%d] @ %d : %d => "
                    + "rollingBlockMeanNxNInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width,
            // height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("float areaFilterInternal " + boxSize,
          boxSlowTotal, "rollingBlockMeanNxNInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("float areaFilterInternal", slowTotal,
        "rollingBlockMeanNxNInternal", fastTotal));
  }

  @SpeedTag
  @SeededTest
  public void floatStripedBlockMeanNxNInternalIsFasterThanAreaFilterNxNInternal(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final BlockMeanFilter filter1 = new BlockMeanFilter();
    final AreaAverageFilter filter2 = new AreaAverageFilter();

    final int iter = 10;
    final ArrayList<float[]> dataSet = getSpeedData(seed, iter);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    // filter1.rollingBlockFilterNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0],
    // boxSizes[0]);
    // filter2.areaFilterInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0] -
    // 0.05);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(iter);
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          // Initialise
          final float w = (float) (boxSize - 0.05);
          for (final float[] data : dataSet2) {
            filter1.stripedBlockFilterNxNInternal(data.clone(), width, height, w);
          }
          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter1.stripedBlockFilterNxNInternal(data, width, height, w);
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
          final ArrayList<float[]> dataSet2 = new ArrayList<>(iter);
          for (final float[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          // Initialise
          for (final float[] data : dataSet2) {
            filter2.areaAverageUsingAveragesInternal(data.clone(), width, height, boxSize - 0.05);
          }
          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter2.areaAverageUsingAveragesInternal(data, width, height, boxSize - 0.05);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "float areaFilterInternal [%dx%d] @ %d : %d => "
                    + "stripedBlockMeanNxNInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width,
            // height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("float areaFilterInternal " + boxSize,
          boxSlowTotal, "stripedBlockMeanNxNInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("float areaFilterInternal", slowTotal,
        "stripedBlockMeanNxNInternal", fastTotal));
  }

  @SuppressWarnings("deprecation")
  @SpeedTag
  @SeededTest
  public void
      floatRollingBlockSumNxNInternalIsFasterThanIntRollingBlockSumNxNInternal(RandomSeed seed) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    final SumFilter filter = new SumFilter();
    final BlockSumFilter filter2 = new BlockSumFilter();

    final int iter = 50;
    final ArrayList<int[]> dataSet = getIntSpeedData(seed, iter);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    // Initialise
    filter.rollingBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
    filter2.rollingBlockFilterNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0],
        boxSizes[0]);

    for (final int boxSize : boxSizes) {
      for (final int width : primes) {
        for (final int height : primes) {
          final ArrayList<float[]> dataSet2 = new ArrayList<>(iter);
          for (final int[] data : dataSet) {
            dataSet2.add(floatClone(data));
          }

          final long start = System.nanoTime();
          for (final float[] data : dataSet2) {
            filter2.rollingBlockFilterNxNInternal(data, width, height, boxSize);
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
          final ArrayList<int[]> dataSet2 = new ArrayList<>(iter);
          for (final int[] data : dataSet) {
            dataSet2.add(intClone(data));
          }

          final long start = System.nanoTime();
          for (final int[] data : dataSet2) {
            filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format(
                "int rollingBlockSumNxNInternal [%dx%d] @ %d : %d => "
                    + "float rollingBlockSumNxNInternal %d = %.2fx",
                width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime)));
            // Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width,
            // height, boxSize,
            // blockTime, time), blockTime < time);
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("float rollingBlockSumNxNInternal " + boxSize,
          boxSlowTotal, "rollingBlockSumNxNInternal", boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord("float rollingBlockSumNxNInternal", slowTotal,
        "rollingBlockSumNxNInternal", fastTotal));
  }

  // TODO
  // int sum faster than float sum
  // Internal version vs complete version -> Is the speed hit significant?
}
