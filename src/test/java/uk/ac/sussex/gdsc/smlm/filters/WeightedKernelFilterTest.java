/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
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

import uk.ac.sussex.gdsc.core.utils.rng.GaussianSamplerUtils;
import uk.ac.sussex.gdsc.test.api.TestAssertions;
import uk.ac.sussex.gdsc.test.api.TestHelper;
import uk.ac.sussex.gdsc.test.api.function.FloatFloatBiPredicate;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.functions.FunctionUtils;

import gnu.trove.list.array.TDoubleArrayList;

import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.GaussianSampler;

import java.util.Arrays;

@SuppressWarnings({"javadoc"})
public abstract class WeightedKernelFilterTest extends WeightedFilterTest {
  @SeededTest
  public void filterPerformsWeightedKernelFiltering(RandomSeed seed) {
    final DataFilter filter = createDataFilter();

    final UniformRandomProvider rg = RngUtils.create(seed.getSeedAsLong());
    final GaussianSampler gs = GaussianSamplerUtils.createGaussianSampler(rg, 2, 0.2);

    final float[] offsets = getOffsets(filter);
    final int[] boxSizes = getBoxSizes(filter);

    final TDoubleArrayList l1 = new TDoubleArrayList();
    final FloatFloatBiPredicate equality = TestHelper.floatsAreClose(1e-6, 0);

    for (final int width : primes) {
      for (final int height : new int[] {29}) {
        final float[] data = createData(width, height, rg);
        l1.reset();

        // Uniform weights
        final float[] w1 = new float[width * height];
        Arrays.fill(w1, 0.5f);

        // Weights simulating the variance of sCMOS pixels
        final float[] w2 = new float[width * height];
        for (int i = 0; i < w2.length; i++) {
          w2[i] = (float) (1.0 / Math.max(0.01, gs.sample()));
        }

        for (final int boxSize : boxSizes) {
          for (final float offset : offsets) {
            for (final boolean internal : checkInternal) {
              // For each pixel over the range around the pixel (vi).
              // kernel filter: sum(vi * ki) / sum(ki)
              // Weighted kernel filter: sum(vi * wi * ki) / sum(ki * wi)
              // Note: The kernel filter is like a weighted filter
              // (New kernel = wi * ki)

              filter.setWeights(null, width, height);

              // Uniform weights
              testfilterPerformsWeightedFiltering(filter, width, height, data, w1, boxSize, offset,
                  internal, equality);

              // Random weights.
              testfilterPerformsWeightedFiltering(filter, width, height, data, w2, boxSize, offset,
                  internal, equality);
            }
          }
        }
      }
    }
  }

  private static void testfilterPerformsWeightedFiltering(DataFilter filter, int width, int height,
      float[] data, float[] w, int boxSize, float offset, boolean internal,
      FloatFloatBiPredicate equality) throws AssertionError {
    // The kernel filter f(x) should compute:
    // sum(vi * wi * ki) / sum(ki * wi)
    // Note: The kernel filter is like a weighted filter
    // (New kernel = wi * ki)
    // If the kernel is normalised to 1 then this is equal to:
    // f(vi * wi) / f(wi)

    filter.setWeights(null, width, height);
    final float[] fWi = filter(w, width, height, boxSize - offset, internal, filter);
    final float[] e = data.clone();
    for (int i = 0; i < e.length; i++) {
      e[i] = data[i] * w[i];
    }
    final float[] fViWi = filter(e, width, height, boxSize - offset, internal, filter);
    for (int i = 0; i < e.length; i++) {
      e[i] = fViWi[i] / fWi[i];
    }

    filter.setWeights(w, width, height);
    final float[] o = filter(data, width, height, boxSize - offset, internal, filter);

    TestAssertions.assertArrayTest(e, o, equality,
        FunctionUtils.getSupplier("%s : [%dx%d] @ %.1f [internal=%b]", filter.name, width, height,
            boxSize - offset, internal));
  }
}
