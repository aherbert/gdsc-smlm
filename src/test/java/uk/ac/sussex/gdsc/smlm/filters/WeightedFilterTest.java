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

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.RandomUtils;
import uk.ac.sussex.gdsc.core.utils.rng.GaussianSamplerUtils;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.FloatFloatBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.GaussianSampler;
import org.junit.jupiter.api.Assertions;

import java.util.Arrays;

@SuppressWarnings({"javadoc"})
public abstract class WeightedFilterTest {
  /** The primes used for the width./height of images during filter testing. */
  static int[] primes = new int[] {113, /* 97, 53, */ 29};
  /**
   * The box sizes used during filter testing. 15 is required to make the box larger than the
   * smallest image.
   */
  static int[] boxSizes = new int[] {15, 5, 3, 2, 1};
  static float[] offsets = new float[] {0, 0.3f, 0.6f};
  /** The check internal flags [true,false]. */
  static boolean[] checkInternal = new boolean[] {true, false};

  float[] createData(int width, int height, UniformRandomProvider rg) {
    final float[] data = new float[width * height];
    for (int i = data.length; i-- > 0;) {
      data[i] = i;
    }
    RandomUtils.shuffle(data, rg);
    return data;
  }

  abstract DataFilter createDataFilter();

  @SeededTest
  public void evenWeightsDoesNotAlterFiltering(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());

    final DataFilter filter1 = createDataFilter();
    final DataFilter filter2 = createDataFilter();
    final float[] offsets = getOffsets(filter1);
    final int[] boxSizes = getBoxSizes(filter1);

    final int[] primes =
        Arrays.copyOf(WeightedFilterTest.primes, WeightedFilterTest.primes.length - 1);

    final FloatFloatBiPredicate equality = TestHelper.floatsAreClose(1e-4, 0);

    for (final int width : primes) {
      for (final int height : primes) {
        final float[] data = createData(width, height, rg);

        // Uniform weights
        final float[] w = new float[width * height];
        Arrays.fill(w, 0.5f);
        filter2.setWeights(w, width, height);

        for (final int boxSize : boxSizes) {
          for (final float offset : offsets) {
            for (final boolean internal : checkInternal) {
              final float[] e = filter(data, width, height, boxSize - offset, internal, filter1);
              final float[] o = filter(data, width, height, boxSize - offset, internal, filter2);
              try {
                TestAssertions.assertArrayTest(e, o, equality);
              } catch (final AssertionError ex) {
                final String msg = String.format("%s : [%dx%d] @ %.1f [internal=%b]", filter2.name,
                    width, height, boxSize - offset, internal);
                throw new AssertionError(msg, ex);
              }
            }
          }
        }
      }
    }
  }

  protected static float[] getOffsets(DataFilter filter1) {
    final float[] offsets = (filter1.isInterpolated) ? WeightedFilterTest.offsets : new float[1];
    return offsets;
  }

  protected static int[] getBoxSizes(DataFilter filter1) {
    if (filter1.minBoxSize == 0) {
      return boxSizes;
    }
    final TIntArrayList list = new TIntArrayList();
    for (final int b : boxSizes) {
      if (b >= filter1.minBoxSize) {
        list.add(b);
      }
    }
    return list.toArray();
  }

  @SeededTest
  public void filterDoesNotAlterFilteredImageMean(RandomSeed seed) {
    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    // ExponentialDistribution ed = new ExponentialDistribution(rand, 57,
    // ExponentialDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);

    final DataFilter filter = createDataFilter();
    final float[] offsets = getOffsets(filter);
    final int[] boxSizes = getBoxSizes(filter);

    final TDoubleArrayList l1 = new TDoubleArrayList();

    final GaussianSampler gs = GaussianSamplerUtils.createGaussianSampler(rg, 2, 0.2);

    for (final int width : primes) {
      for (final int height : primes) {
        final float[] data = createData(width, height, rg);
        l1.reset();

        filter.setWeights(null, width, height);
        for (final int boxSize : boxSizes) {
          for (final float offset : offsets) {
            for (final boolean internal : checkInternal) {
              l1.add(getMean(data, width, height, boxSize - offset, internal, filter));
            }
          }
        }

        final double[] e = l1.toArray();
        int ei = 0;

        // Uniform weights
        final float[] w = new float[width * height];

        Arrays.fill(w, 0.5f);
        filter.setWeights(w, width, height);
        for (final int boxSize : boxSizes) {
          for (final float offset : offsets) {
            for (final boolean internal : checkInternal) {
              testMean(data, width, height, boxSize - offset, internal, filter, "w=0.5", e[ei++],
                  1e-5);
            }
          }
        }

        // Weights simulating the variance of sCMOS pixels
        for (int i = 0; i < w.length; i++) {
          // w[i] = (float) (1.0 / Math.max(0.01, ed.sample()));
          w[i] = (float) (1.0 / Math.max(0.01, gs.sample()));
          // w[i] = 0.5f;
        }

        ei = 0;
        filter.setWeights(w, width, height);
        for (final int boxSize : boxSizes) {
          for (final float offset : offsets) {
            for (final boolean internal : checkInternal) {
              testMean(data, width, height, boxSize - offset, internal, filter, "w=?", e[ei++],
                  5e-2);
            }
          }
        }
      }
    }
  }

  protected static float[] filter(float[] data, int width, int height, float boxSize,
      boolean internal, DataFilter filter) {
    final float[] data1 = data.clone();
    if (internal) {
      filter.filterInternal(data1, width, height, boxSize);
    } else {
      filter.filter(data1, width, height, boxSize);
    }
    return data1;
  }

  protected static double getMean(float[] data, int width, int height, float boxSize,
      boolean internal, DataFilter filter) {
    return MathUtils.sum(filter(data, width, height, boxSize, internal, filter)) / data.length;
  }

  protected static double testMean(float[] data, int width, int height, float boxSize,
      boolean internal, DataFilter filter, String title, double u1, double tol) {
    final double u2 = getMean(data, width, height, boxSize, internal, filter);
    final double error = DoubleEquality.relativeError(u1, u2);
    // logger.fine(() -> String.format("%s : %s [%dx%d] @ %.1f [internal=%b] : %g => %g (%g)",
    // filter.name, title,
    // width, height, boxSize, internal, u1, u2, error);
    Assertions.assertEquals(0, error, tol,
        FunctionUtils.getSupplier("%s : %s [%dx%d] @ %.1f [internal=%b] : %g => %g  (%g)",
            filter.name, title, width, height, boxSize, internal, u1, u2, error));
    return u2;
  }
}
