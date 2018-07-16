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
package gdsc.smlm.filters;

import java.util.Arrays;

import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Test;

import gdsc.test.TestAssert;
import gdsc.test.TestSettings;
import gnu.trove.list.array.TDoubleArrayList;

@SuppressWarnings({ "javadoc" })
public abstract class WeightedKernelFilterTest extends WeightedFilterTest
{
	@Test
	public void filterPerformsWeightedKernelFiltering()
	{
		DataFilter filter = createDataFilter();

		RandomGenerator rg = TestSettings.getRandomGenerator();

		float[] offsets = getOffsets(filter);
		int[] boxSizes = getBoxSizes(filter);

		boolean[] checkInternal = new boolean[] { false };

		TDoubleArrayList l1 = new TDoubleArrayList();

		for (int width : primes)
			for (int height : new int[] { 29 })
			{
				float[] data = createData(width, height, rg);
				l1.reset();

				// Uniform weights
				float[] w1 = new float[width * height];
				Arrays.fill(w1, 0.5f);

				// Weights simulating the variance of sCMOS pixels
				float[] w2 = new float[width * height];
				for (int i = 0; i < w2.length; i++)
				{
					w2[i] = (float) (1.0 / Math.max(0.01, rg.nextGaussian() * 0.2 + 2));
				}

				for (int boxSize : boxSizes)
					for (float offset : offsets)
						for (boolean internal : checkInternal)
						{
							// For each pixel over the range around the pixel (vi).
							// kernel filter: sum(vi * ki) / sum(ki)
							// Weighted kernel filter: sum(vi * wi * ki) / sum(ki * wi)
							// Note: The kernel filter is like a weighted filter 
							// (New kernel = wi * ki)

							filter.setWeights(null, width, height);

							// Uniform weights
							testfilterPerformsWeightedFiltering(filter, width, height, data, w1, boxSize, offset,
									internal);

							// Random weights.
							testfilterPerformsWeightedFiltering(filter, width, height, data, w2, boxSize, offset,
									internal);
						}
			}
	}

	private static void testfilterPerformsWeightedFiltering(DataFilter filter, int width, int height, float[] data,
			float[] w, int boxSize, float offset, boolean internal) throws AssertionError
	{
		// The kernel filter f(x) should compute:
		//    sum(vi * wi * ki) / sum(ki * wi)
		// Note: The kernel filter is like a weighted filter 
		// (New kernel = wi * ki)
		// If the kernel is normalised to 1 then this is equal to:
		//    f(vi * wi) / f(wi)

		filter.setWeights(null, width, height);
		float[] fWi = filter(w, width, height, boxSize - offset, internal, filter);
		float[] e = data.clone();
		for (int i = 0; i < e.length; i++)
			e[i] = data[i] * w[i];
		float[] fViWi = filter(e, width, height, boxSize - offset, internal, filter);
		for (int i = 0; i < e.length; i++)
			e[i] = fViWi[i] / fWi[i];

		filter.setWeights(w, width, height);
		float[] o = filter(data, width, height, boxSize - offset, internal, filter);

		TestAssert.assertArrayEqualsRelative(e, o, 1e-4f, "%s : [%dx%d] @ %.1f [internal=%b]", filter.name, width,
				height, boxSize - offset, internal);
	}
}
