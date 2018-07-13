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

import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.Maths;
import gdsc.core.utils.Random;
import gdsc.test.TestAssert;
import gdsc.test.TestSettings;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;

@SuppressWarnings({ "javadoc" })
public abstract class WeightedFilterTest
{
	/** The primes used for the width/height of images during filter testing. */
	static int[] primes = new int[] { 113, /* 97, 53, */ 29 };
	/**
	 * The box sizes used during filter testing.
	 * 15 is required to make the box larger than the smallest image.
	 */
	static int[] boxSizes = new int[] { 15, 5, 3, 2, 1 };
	static float[] offsets = new float[] { 0, 0.3f, 0.6f };
	/** The check internal flags [true,false]. */
	static boolean[] checkInternal = new boolean[] { true, false };

	float[] createData(int width, int height, RandomGenerator rg)
	{
		float[] data = new float[width * height];
		for (int i = data.length; i-- > 0;)
			data[i] = i;
		Random.shuffle(data, rg);
		return data;
	}

	abstract DataFilter createDataFilter();

	@Test
	public void evenWeightsDoesNotAlterFiltering()
	{
		RandomGenerator rg = TestSettings.getRandomGenerator();

		DataFilter filter1 = createDataFilter();
		DataFilter filter2 = createDataFilter();
		float[] offsets = getOffsets(filter1);
		int[] boxSizes = getBoxSizes(filter1);

		int[] primes = Arrays.copyOf(WeightedFilterTest.primes, WeightedFilterTest.primes.length - 1);

		for (int width : primes)
			for (int height : primes)
			{
				float[] data = createData(width, height, rg);

				// Uniform weights
				float[] w = new float[width * height];
				Arrays.fill(w, 0.5f);
				filter2.setWeights(w, width, height);

				for (int boxSize : boxSizes)
					for (float offset : offsets)
						for (boolean internal : checkInternal)
						{
							float[] e = filter(data, width, height, boxSize - offset, internal, filter1);
							float[] o = filter(data, width, height, boxSize - offset, internal, filter2);
							try
							{
								TestAssert.assertArrayEqualsRelative(e, o, 1e-4f);
							}
							catch (AssertionError ex)
							{
								String msg = String.format("%s : [%dx%d] @ %.1f [internal=%b]", filter2.name, width,
										height, boxSize - offset, internal);
								throw new AssertionError(msg, ex);
							}
						}
			}
	}

	private static float[] getOffsets(DataFilter filter1)
	{
		float[] offsets = (filter1.isInterpolated) ? WeightedFilterTest.offsets : new float[1];
		return offsets;
	}

	private static int[] getBoxSizes(DataFilter filter1)
	{
		if (filter1.minBoxSize == 0)
			return boxSizes;
		TIntArrayList list = new TIntArrayList();
		for (int b : boxSizes)
			if (b >= filter1.minBoxSize)
				list.add(b);
		return list.toArray();
	}

	@Test
	public void filterDoesNotAlterFilteredImageMean()
	{
		RandomGenerator rg = TestSettings.getRandomGenerator();
		//ExponentialDistribution ed = new ExponentialDistribution(rand, 57,
		//		ExponentialDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);

		DataFilter filter = createDataFilter();
		float[] offsets = getOffsets(filter);
		int[] boxSizes = getBoxSizes(filter);

		TDoubleArrayList l1 = new TDoubleArrayList();

		for (int width : primes)
			for (int height : primes)
			{
				float[] data = createData(width, height, rg);
				l1.reset();

				filter.setWeights(null, width, height);
				for (int boxSize : boxSizes)
					for (float offset : offsets)
						for (boolean internal : checkInternal)
						{
							l1.add(getMean(data, width, height, boxSize - offset, internal, filter));
						}

				double[] e = l1.toArray();
				int ei = 0;

				// Uniform weights
				float[] w = new float[width * height];

				Arrays.fill(w, 0.5f);
				filter.setWeights(w, width, height);
				for (int boxSize : boxSizes)
					for (float offset : offsets)
						for (boolean internal : checkInternal)
						{
							testMean(data, width, height, boxSize - offset, internal, filter, "w=0.5", e[ei++], 1e-5);
						}

				// Weights simulating the variance of sCMOS pixels
				for (int i = 0; i < w.length; i++)
				{
					//w[i] = (float) (1.0 / Math.max(0.01, ed.sample()));
					w[i] = (float) (1.0 / Math.max(0.01, rg.nextGaussian() * 0.2 + 2));
					//w[i] = 0.5f;
				}

				ei = 0;
				filter.setWeights(w, width, height);
				for (int boxSize : boxSizes)
					for (float offset : offsets)
						for (boolean internal : checkInternal)
						{
							testMean(data, width, height, boxSize - offset, internal, filter, "w=?", e[ei++], 5e-2);
						}
			}
	}

	@Test
	public void filterPerformsWeightedFiltering()
	{
		DataFilter filter = createDataFilter();
		//org.junit.Assume.assumeFalse(filter.isSumFilter());

		RandomGenerator rg = TestSettings.getRandomGenerator();

		float[] offsets = getOffsets(filter);
		int[] boxSizes = getBoxSizes(filter);

		TDoubleArrayList l1 = new TDoubleArrayList();

		for (int width : primes)
			for (int height : primes)
			{
				float[] data = createData(width, height, rg);
				l1.reset();

				// Uniform weights
				float[] w1 = new float[width * height];
				Arrays.fill(w1, 1f);
				float[] w = new float[width * height];
				Arrays.fill(w, 0.5f);

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
							// The filter f(x) should compute:
							//   f(vi * wi) / (f(wi) * f(1)) 
							// For each pixel over the range around the pixel (vi).
							// Use of f(1) handles sum filters.

							filter.setWeights(null, width, height);

							float[] f1 = filter(w1, width, height, boxSize - offset, internal, filter);
							
							// Uniform weights
							testfilterPerformsWeightedFiltering(filter, width, height, data, w, boxSize, offset, internal, f1);
							
							// Random weights.
							// This fails for the Gaussian and Kernel filters
							testfilterPerformsWeightedFiltering(filter, width, height, data, w2, boxSize, offset, internal, f1);
						}
			}
	}

	private static void testfilterPerformsWeightedFiltering(DataFilter filter, int width, int height, float[] data, float[] w, int boxSize, float offset,
			boolean internal, float[] f1) throws AssertionError
	{
		// The filter f(x) should compute:
		//   f(vi * wi) / (f(wi) * f(1)) 
		// For each pixel over the range around the pixel (vi).
		// Use of f(1) handles sum filters.
		
		filter.setWeights(null, width, height);
		float[] fWi = filter(w, width, height, boxSize - offset, internal, filter);
		float[] e = data.clone();
		for (int i = 0; i < e.length; i++)
			e[i] = data[i] * w[i];
		float[] fViWi = filter(e, width, height, boxSize - offset, internal, filter);
		for (int i = 0; i < e.length; i++)
			e[i] = fViWi[i] / (fWi[i] / f1[i]);

		filter.setWeights(w, width, height);
		float[] o = filter(data, width, height, boxSize - offset, internal, filter);

		try
		{
			TestAssert.assertArrayEqualsRelative(e, o, 1e-4f);
		}
		catch (AssertionError ex)
		{
			String msg = String.format("%s : [%dx%d] @ %.1f [internal=%b]", filter.name, width,
					height, boxSize - offset, internal);
			throw new AssertionError(msg, ex);
		}
	}

	private static float[] filter(float[] data, int width, int height, float boxSize, boolean internal,
			DataFilter filter)
	{
		float[] data1 = data.clone();
		if (internal)
		{
			filter.filterInternal(data1, width, height, boxSize);
		}
		else
		{
			filter.filter(data1, width, height, boxSize);
		}
		return data1;
	}

	private static double getMean(float[] data, int width, int height, float boxSize, boolean internal,
			DataFilter filter)
	{
		return Maths.sum(filter(data, width, height, boxSize, internal, filter)) / data.length;
	}

	private static double testMean(float[] data, int width, int height, float boxSize, boolean internal,
			DataFilter filter, String title, double u1, double tol)
	{
		double u2 = getMean(data, width, height, boxSize, internal, filter);
		double error = DoubleEquality.relativeError(u1, u2);
		//System.out.printf("%s : %s [%dx%d] @ %.1f [internal=%b] : %g => %g  (%g)\n", filter.name, title,
		//		width, height, boxSize, internal, u1, u2, error);
		TestAssert.assertTrue(error <= tol, "%s : %s [%dx%d] @ %.1f [internal=%b] : %g => %g  (%g)", filter.name, title,
				width, height, boxSize, internal, u1, u2, error);
		return u2;
	}
}
