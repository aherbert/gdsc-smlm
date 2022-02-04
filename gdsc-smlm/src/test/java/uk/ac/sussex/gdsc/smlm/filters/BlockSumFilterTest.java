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

package uk.ac.sussex.gdsc.smlm.filters;

import java.util.ArrayList;
import java.util.Arrays;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.ContinuousSampler;
import org.junit.jupiter.api.Assumptions;
import uk.ac.sussex.gdsc.core.utils.FloatEquality;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.junit5.SpeedTag;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;
import uk.ac.sussex.gdsc.test.utils.TestSettings;

@SuppressWarnings({"javadoc"})
class BlockSumFilterTest extends AbstractFilterTest {
  private static final int INTERNAL_ITER3 = 500;
  private static final int INTERNAL_ITER = 50;
  private static final int ITER3 = 200;
  private static final int ITER = 20;

  /**
   * Do a simple and stupid sum filter.
   *
   * @param data the data
   * @param maxx the maxx
   * @param maxy the maxy
   * @param boxSize the box size
   */
  public static void sum(float[] data, int maxx, int maxy, float boxSize) {
    if (boxSize <= 0) {
      return;
    }

    final int n = (int) Math.ceil(boxSize);
    final int size = 2 * n + 1;
    final float[] weight = new float[size];
    Arrays.fill(weight, 1);
    if (boxSize != n) {
      weight[0] = weight[weight.length - 1] = boxSize - (n - 1);
    }

    final float[] out = new float[data.length];

    final int[] oy = new int[size];
    final int[] ox = new int[size];

    for (int y = 0; y < maxy; y++) {
      // Cache offset
      for (int yy = 0; yy < size; yy++) {
        int yyy = y + yy - n;
        if (yyy < 0) {
          yyy = 0;
        } else if (yyy >= maxy) {
          yyy = maxy - 1;
        }
        oy[yy] = yyy * maxx;
      }

      for (int x = 0; x < maxx; x++) {
        // Cache offset
        for (int xx = 0; xx < size; xx++) {
          int xxx = x + xx - n;
          if (xxx < 0) {
            xxx = 0;
          } else if (xxx >= maxx) {
            xxx = maxx - 1;
          }
          ox[xx] = xxx;
        }

        double sum = 0;
        for (int yy = 0; yy < size; yy++) {
          // int yyy = y + yy - n;
          // if (yyy < 0)
          // yyy = 0;
          // else if (yyy >= maxy)
          // yyy = maxy - 1;

          final int index = oy[yy];
          final float wy = weight[yy];
          for (int xx = 0; xx < size; xx++) {
            sum += data[index + ox[xx]] * wy * weight[xx];
          }
        }
        out[y * maxx + x] = (float) (sum);
      }
    }
    System.arraycopy(out, 0, data, 0, out.length);
  }

  /**
   * Do a simple and stupid sum filter with weights.
   *
   * @param data the data
   * @param weights the weights
   * @param maxx the maxx
   * @param maxy the maxy
   * @param boxSize the box size
   */
  public static void weightedSum(float[] data, float[] weights, int maxx, int maxy, float boxSize) {
    if (boxSize <= 0) {
      return;
    }

    final int n = (int) Math.ceil(boxSize);
    final int size = 2 * n + 1;
    final float[] weight = new float[size];
    Arrays.fill(weight, 1);
    if (boxSize != n) {
      weight[0] = weight[weight.length - 1] = boxSize - (n - 1);
    }
    final double area = MathUtils.pow2(2 * boxSize + 1);

    final float[] out = new float[data.length];

    for (int y = 0; y < maxy; y++) {
      for (int x = 0; x < maxx; x++) {
        double sum = 0;
        double sumw = 0;
        for (int yy = 0; yy < size; yy++) {
          int yyy = y + yy - n;
          if (yyy < 0) {
            yyy = 0;
          }
          if (yyy >= maxy) {
            yyy = maxy - 1;
          }
          for (int xx = 0; xx < size; xx++) {
            int xxx = x + xx - n;
            if (xxx < 0) {
              xxx = 0;
            }
            if (xxx >= maxx) {
              xxx = maxx - 1;
            }
            final int index = yyy * maxx + xxx;
            final double w2 = weights[index] * weight[yy] * weight[xx];
            sum += data[index] * w2;
            sumw += w2;
          }
        }
        // The sum should not be effected by the weights.
        out[y * maxx + x] = (float) (sum / (sumw / area));
      }
    }
    System.arraycopy(out, 0, data, 0, out.length);
  }

  /**
   * Used to test the filter methods calculate the correct result.
   */
  private abstract class BlockSumDataFilter extends DataFilter {
    public BlockSumDataFilter(String name, boolean isInterpolated) {
      super(name, isInterpolated);
    }

    BlockSumFilter filter = new BlockSumFilter();

    @Override
    public void setWeights(float[] weights, int width, int height) {
      filter.setWeights(weights, width, height);
    }
  }

  private static void sumIsCorrect(float[] data, int width, int height, float boxSize,
      boolean internal, BlockSumDataFilter filter) {
    final float[] data1 = data.clone();
    final float[] data2 = data.clone();
    final FloatEquality eq = new FloatEquality(1e-3f, 1e-10f);

    sum(data1, width, height, boxSize);
    if (internal) {
      filter.filterInternal(data2, width, height, boxSize);
      floatArrayEquals(eq, data1, data2, width, height, boxSize,
          "Internal arrays do not match: [%dx%d] @ %.1f", width, height, boxSize);
    } else {
      filter.filter(data2, width, height, boxSize);
      floatArrayEquals(eq, data1, data2, width, height, 0, "Arrays do not match: [%dx%d] @ %.1f",
          width, height, boxSize);
    }
  }

  private static void weightedSumIsCorrect(float[] data, float[] weights, int width, int height,
      float boxSize, boolean internal, BlockSumDataFilter filter) {
    final float[] data1 = data.clone();
    final float[] data2 = data.clone();
    final FloatEquality eq = new FloatEquality(1e-3f, 1e-10f);

    weightedSum(data1, weights, width, height, boxSize);

    //// Check the weights do not alter the image sum
    // double u1 = uk.ac.sussex.gdsc.core.utils.Maths.sum(sum(data.clone(), width, height,
    //// boxSize));
    // double u2 = uk.ac.sussex.gdsc.core.utils.Maths.sum(data1);
    // logger.fine(() -> String.format("[%dx%d] @ %.1f : %g => %g (%g)", width, height, boxSize, u1,
    //// u2,
    // uk.ac.sussex.gdsc.core.utils.DoubleEquality.relativeError(u1, u2));

    if (internal) {
      filter.filterInternal(data2, width, height, boxSize);
      floatArrayEquals(eq, data1, data2, width, height, boxSize,
          "Internal arrays do not match: [%dx%d] @ %.1f", width, height, boxSize);
    } else {
      filter.filter(data2, width, height, boxSize);
      floatArrayEquals(eq, data1, data2, width, height, 0, "Arrays do not match: [%dx%d] @ %.1f",
          width, height, boxSize);
    }
  }

  private static void checkIsCorrect(RandomSeed seed, BlockSumDataFilter filter) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeed());
    final ContinuousSampler ed = SamplerUtils.createExponentialSampler(rg, 57);

    for (final int width : primes) {
      for (final int height : primes) {
        final float[] data = createData(rg, width, height);

        filter.filter.setWeights(null, 0, 0);
        for (final float boxSize : boxSizes) {
          for (final boolean internal : checkInternal) {
            sumIsCorrect(data, width, height, boxSize, internal, filter);
            if (filter.isInterpolated) {
              sumIsCorrect(data, width, height, boxSize - 0.3f, internal, filter);
              sumIsCorrect(data, width, height, boxSize - 0.6f, internal, filter);
            }
          }
        }

        // Uniform weights
        final float[] w = new float[width * height];
        Arrays.fill(w, 1);
        filter.filter.setWeights(w, width, height);
        for (final float boxSize : boxSizes) {
          for (final boolean internal : checkInternal) {
            weightedSumIsCorrect(data, w, width, height, boxSize, internal, filter);
            if (filter.isInterpolated) {
              weightedSumIsCorrect(data, w, width, height, boxSize - 0.3f, internal, filter);
              weightedSumIsCorrect(data, w, width, height, boxSize - 0.6f, internal, filter);
            }
          }
        }

        // Weights simulating the variance of sCMOS pixels
        for (int i = 0; i < w.length; i++) {
          w[i] = (float) (1.0 / Math.max(0.01, ed.sample()));
        }

        filter.filter.setWeights(w, width, height);
        for (final float boxSize : boxSizes) {
          for (final boolean internal : checkInternal) {
            weightedSumIsCorrect(data, w, width, height, boxSize, internal, filter);
            if (filter.isInterpolated) {
              weightedSumIsCorrect(data, w, width, height, boxSize - 0.3f, internal, filter);
              weightedSumIsCorrect(data, w, width, height, boxSize - 0.6f, internal, filter);
            }
          }
        }
      }
    }
  }

  @SeededTest
  void blockFilterIsCorrect(RandomSeed seed) {
    final BlockSumDataFilter filter = new BlockSumDataFilter("block", true) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.blockFilter(data, width, height, boxSize);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.blockFilterInternal(data, width, height, boxSize);
      }
    };
    checkIsCorrect(seed, filter);
  }

  @SeededTest
  void stripedBlockFilterIsCorrect(RandomSeed seed) {
    final BlockSumDataFilter filter = new BlockSumDataFilter("stripedBlock", true) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilter(data, width, height, boxSize);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilterInternal(data, width, height, boxSize);
      }
    };
    checkIsCorrect(seed, filter);
  }

  @SeededTest
  void rollingBlockFilterIsCorrect(RandomSeed seed) {
    final BlockSumDataFilter filter = new BlockSumDataFilter("rollingBlock", false) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.rollingBlockFilter(data, width, height, (int) boxSize);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.rollingBlockFilterInternal(data, width, height, (int) boxSize);
      }
    };
    checkIsCorrect(seed, filter);
  }

  private void speedTest(RandomSeed seed, BlockSumDataFilter fast, BlockSumDataFilter slow) {
    speedTest(seed, fast, slow, boxSizes);
  }

  private void speedTest(RandomSeed seed, BlockSumDataFilter fast, BlockSumDataFilter slow,
      int[] testBoxSizes) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    ArrayList<float[]> dataSet = getSpeedData(seed, ITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    final float[] boxSizes = new float[testBoxSizes.length];
    final float offset = (fast.isInterpolated && slow.isInterpolated) ? 0.3f : 0;
    for (int i = 0; i < boxSizes.length; i++) {
      boxSizes[i] = testBoxSizes[i] - offset;
    }

    // Initialise
    for (final float boxSize : boxSizes) {
      fast.filter(dataSet.get(0).clone(), speedPrimes[0], speedPrimes[0], boxSize);
      slow.filter(dataSet.get(0).clone(), speedPrimes[0], speedPrimes[0], boxSize);
    }

    for (final float boxSize : boxSizes) {
      final int iter = (boxSize == 1) ? ITER3 : ITER;
      for (final int width : speedPrimes) {
        for (final int height : speedPrimes) {
          dataSet = getSpeedData(seed, iter);

          final long start = System.nanoTime();
          for (final float[] data : dataSet) {
            fast.filter(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final float boxSize : boxSizes) {
      final int iter = (boxSize == 1) ? ITER3 : ITER;
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : speedPrimes) {
        for (final int height : speedPrimes) {
          dataSet = getSpeedData(seed, iter);

          final long start = System.nanoTime();
          for (final float[] data : dataSet) {
            slow.filter(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format("%s [%dx%d] @ %.1f : %d => %s %d = %.2fx", slow.name,
                width, height, boxSize, time, fast.name, fastTime, speedUpFactor(time, fastTime)));
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord(slow.name + " " + boxSize, boxSlowTotal,
          fast.name, boxFastTotal));
    }
    logger.log(TestLogUtils.getTimingRecord(slow.name, slowTotal, fast.name, fastTotal));
  }

  private void speedTestInternal(RandomSeed seed, BlockSumDataFilter fast,
      BlockSumDataFilter slow) {
    speedTestInternal(seed, fast, slow, boxSizes);
  }

  private void speedTestInternal(RandomSeed seed, BlockSumDataFilter fast, BlockSumDataFilter slow,
      int[] testBoxSizes) {
    Assumptions.assumeTrue(TestSettings.allow(TestComplexity.MEDIUM));

    ArrayList<float[]> dataSet = getSpeedData(seed, INTERNAL_ITER3);

    final ArrayList<Long> fastTimes = new ArrayList<>();

    final float[] boxSizes = new float[testBoxSizes.length];
    final float offset = (fast.isInterpolated && slow.isInterpolated) ? 0.3f : 0;
    for (int i = 0; i < boxSizes.length; i++) {
      boxSizes[i] = testBoxSizes[i] - offset;
    }

    // Initialise
    for (final float boxSize : boxSizes) {
      fast.filterInternal(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0], boxSize);
      slow.filterInternal(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0], boxSize);
    }

    for (final float boxSize : boxSizes) {
      final int iter = (boxSize == 1) ? INTERNAL_ITER3 : INTERNAL_ITER;
      for (final int width : speedPrimes) {
        for (final int height : speedPrimes) {
          dataSet = getSpeedData(seed, iter);

          final long start = System.nanoTime();
          for (final float[] data : dataSet) {
            fast.filterInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;
          fastTimes.add(time);
        }
      }
    }

    long slowTotal = 0;
    long fastTotal = 0;
    int index = 0;
    for (final float boxSize : boxSizes) {
      final int iter = (boxSize == 1) ? INTERNAL_ITER3 : INTERNAL_ITER;
      long boxSlowTotal = 0;
      long boxFastTotal = 0;
      for (final int width : speedPrimes) {
        for (final int height : speedPrimes) {
          dataSet = getSpeedData(seed, iter);

          final long start = System.nanoTime();
          for (final float[] data : dataSet) {
            slow.filterInternal(data, width, height, boxSize);
          }
          final long time = System.nanoTime() - start;

          final long fastTime = fastTimes.get(index++);
          slowTotal += time;
          fastTotal += fastTime;
          boxSlowTotal += time;
          boxFastTotal += fastTime;
          if (debug) {
            logger.fine(() -> String.format("Internal %s [%dx%d] @ %.1f : %d => %s %d = %.2fx",
                slow.name, width, height, boxSize, time, fast.name, fastTime,
                speedUpFactor(time, fastTime)));
          }
        }
      }
      // if (debug)
      logger.log(TestLogUtils.getStageTimingRecord("Internal " + slow.name + " " + boxSize,
          boxSlowTotal, fast.name, boxFastTotal));
    }
    logger.log(
        TestLogUtils.getTimingRecord("Internal " + slow.name, slowTotal, fast.name, fastTotal));

  }

  @SpeedTag
  @SeededTest
  void stripedBlockIsFasterThanBlock(RandomSeed seed) {
    final BlockSumDataFilter slow = new BlockSumDataFilter("block", false) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.blockFilter(data, width, height, (int) boxSize);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.blockFilterInternal(data, width, height, (int) boxSize);
      }
    };
    final BlockSumDataFilter fast = new BlockSumDataFilter("stripedBlock", false) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilter(data, width, height, (int) boxSize);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilterInternal(data, width, height, (int) boxSize);
      }
    };

    speedTest(seed, fast, slow);
    speedTestInternal(seed, fast, slow);
  }

  @SpeedTag
  @SeededTest
  void interpolatedStripedBlockIsFasterThanBlock(RandomSeed seed) {
    final BlockSumDataFilter slow = new BlockSumDataFilter("block", true) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.blockFilter(data, width, height, boxSize);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.blockFilterInternal(data, width, height, boxSize);
      }
    };
    final BlockSumDataFilter fast = new BlockSumDataFilter("stripedBlock", true) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilter(data, width, height, boxSize);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilterInternal(data, width, height, boxSize);
      }
    };

    speedTest(seed, fast, slow);
    speedTestInternal(seed, fast, slow);
  }

  @SpeedTag
  @SeededTest
  void rollingBlockIsFasterThanBlock(RandomSeed seed) {
    final BlockSumDataFilter slow = new BlockSumDataFilter("block", false) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.blockFilter(data, width, height, (int) boxSize);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.blockFilterInternal(data, width, height, (int) boxSize);
      }
    };
    final BlockSumDataFilter fast = new BlockSumDataFilter("rollingBlock", false) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.rollingBlockFilter(data, width, height, (int) boxSize);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.rollingBlockFilterInternal(data, width, height, (int) boxSize);
      }
    };

    speedTest(seed, fast, slow);
    speedTestInternal(seed, fast, slow);
  }

  @SpeedTag
  @SeededTest
  void rollingBlockIsFasterThanStripedBlock(RandomSeed seed) {
    final BlockSumDataFilter slow = new BlockSumDataFilter("stripedBlock", false) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilter(data, width, height, (int) boxSize);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilterInternal(data, width, height, (int) boxSize);
      }
    };
    final BlockSumDataFilter fast = new BlockSumDataFilter("rollingBlock", false) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.rollingBlockFilter(data, width, height, (int) boxSize);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.rollingBlockFilterInternal(data, width, height, (int) boxSize);
      }
    };

    speedTest(seed, fast, slow);
    speedTestInternal(seed, fast, slow);
  }

  @SpeedTag
  @SeededTest
  void stripedBlock3x3IsFasterThanStripedBlockNxN(RandomSeed seed) {
    final BlockSumDataFilter slow = new BlockSumDataFilter("stripedBlockNxN", false) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilterNxN(data, width, height, (int) boxSize);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilterNxNInternal(data, width, height, (int) boxSize);
      }
    };
    final BlockSumDataFilter fast = new BlockSumDataFilter("stripedBlock3x3", false) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilter3x3(data, width, height);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilter3x3Internal(data, width, height);
      }
    };

    final int[] testBoxSizes = new int[] {1};
    speedTest(seed, fast, slow, testBoxSizes);
    speedTestInternal(seed, fast, slow, testBoxSizes);
  }

  @SpeedTag
  @SeededTest
  void interpolatedStripedBlock3x3IsFasterThanStripedBlockNxN(RandomSeed seed) {
    final BlockSumDataFilter slow = new BlockSumDataFilter("stripedBlockNxN", true) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilterNxN(data, width, height, boxSize);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilterNxNInternal(data, width, height, boxSize);
      }
    };
    final BlockSumDataFilter fast = new BlockSumDataFilter("stripedBlock3x3", true) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilter3x3(data, width, height, boxSize);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilter3x3Internal(data, width, height, boxSize);
      }
    };

    final int[] testBoxSizes = new int[] {1};
    speedTest(seed, fast, slow, testBoxSizes);
    speedTestInternal(seed, fast, slow, testBoxSizes);
  }

  @SpeedTag
  @SeededTest
  void stripedBlock5x5IsFasterThanStripedBlockNxN(RandomSeed seed) {
    final BlockSumDataFilter slow = new BlockSumDataFilter("stripedBlockNxN", false) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilterNxN(data, width, height, (int) boxSize);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilterNxNInternal(data, width, height, (int) boxSize);
      }
    };
    final BlockSumDataFilter fast = new BlockSumDataFilter("stripedBlock5x5", false) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilter5x5(data, width, height);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilter5x5Internal(data, width, height);
      }
    };

    final int[] testBoxSizes = new int[] {2};
    speedTest(seed, fast, slow, testBoxSizes);
    speedTestInternal(seed, fast, slow, testBoxSizes);
  }

  @SpeedTag
  @SeededTest
  void interpolatedStripedBlock5x5IsFasterThanStripedBlockNxN(RandomSeed seed) {
    final BlockSumDataFilter slow = new BlockSumDataFilter("stripedBlockNxN", true) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilterNxN(data, width, height, boxSize);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilterNxNInternal(data, width, height, boxSize);
      }
    };
    final BlockSumDataFilter fast = new BlockSumDataFilter("stripedBlock5x5", true) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilter5x5(data, width, height, boxSize);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilter5x5Internal(data, width, height, boxSize);
      }
    };

    final int[] testBoxSizes = new int[] {2};
    speedTest(seed, fast, slow, testBoxSizes);
    speedTestInternal(seed, fast, slow, testBoxSizes);
  }

  @SpeedTag
  @SeededTest
  void stripedBlock7x7IsFasterThanStripedBlockNxN(RandomSeed seed) {
    final BlockSumDataFilter slow = new BlockSumDataFilter("stripedBlockNxN", false) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilterNxN(data, width, height, (int) boxSize);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilterNxNInternal(data, width, height, (int) boxSize);
      }
    };
    final BlockSumDataFilter fast = new BlockSumDataFilter("stripedBlock7x7", false) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilter7x7(data, width, height);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilter7x7Internal(data, width, height);
      }
    };

    final int[] testBoxSizes = new int[] {3};
    speedTest(seed, fast, slow, testBoxSizes);
    speedTestInternal(seed, fast, slow, testBoxSizes);
  }

  @SpeedTag
  @SeededTest
  void interpolatedStripedBlock7x7IsFasterThanStripedBlockNxN(RandomSeed seed) {
    final BlockSumDataFilter slow = new BlockSumDataFilter("stripedBlockNxN", true) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilterNxN(data, width, height, boxSize);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilterNxNInternal(data, width, height, boxSize);
      }
    };
    final BlockSumDataFilter fast = new BlockSumDataFilter("stripedBlock7x7", true) {
      @Override
      public void filter(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilter7x7(data, width, height, boxSize);
      }

      @Override
      public void filterInternal(float[] data, int width, int height, float boxSize) {
        filter.stripedBlockFilter7x7Internal(data, width, height, boxSize);
      }
    };

    final int[] testBoxSizes = new int[] {3};
    speedTest(seed, fast, slow, testBoxSizes);
    speedTestInternal(seed, fast, slow, testBoxSizes);
  }
}
