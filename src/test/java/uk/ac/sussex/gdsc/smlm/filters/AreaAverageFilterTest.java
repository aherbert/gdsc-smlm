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

import java.util.ArrayList;

import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.TestSettings;

@SuppressWarnings({ "deprecation", "javadoc" })
public class AreaAverageFilterTest extends AbstractFilterTest
{
	private final int ITER = 100;
	private final int InternalITER = 300;

	@Test
	public void areaAverageUsingSumsNxNInternalIsFasterThanAreaAverageNxNInternal()
	{
		TestSettings.assumeMediumComplexity();

		final AreaAverageFilter filter = new AreaAverageFilter();

		final ArrayList<float[]> dataSet = getSpeedData(InternalITER);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		for (final float boxSize : fBoxSizes)
		{
			filter.areaAverageUsingAveragesInternal(dataSet.get(0).clone(), primes[0], primes[0], boxSize);
			filter.areaAverageUsingSumsInternal(dataSet.get(0).clone(), primes[0], primes[0], boxSize);
		}

		for (final float boxSize : fBoxSizes)
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(data.clone());

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter.areaAverageUsingSumsInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (final float boxSize : fBoxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(data.clone());

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter.areaAverageUsingAveragesInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					final long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"float areaAverageInternal [%dx%d] @ %.1f : %d => areaAverageUsingSumsInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestAssert.assert_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			TestLog.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"float areaAverageInternal %.1f : %d => areaAverageUsingSumsInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float areaAverageInternal %d => areaAverageUsingSumsInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void stripedBlockAverageIsFasterThanAreaAverage()
	{
		TestSettings.assumeSpeedTest();

		final AreaAverageFilter filter = new AreaAverageFilter();
		final AverageFilter filter2 = new AverageFilter();

		final ArrayList<float[]> dataSet = getSpeedData(ITER);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		for (final float boxSize : fBoxSizes)
		{
			filter.areaAverageUsingAverages(dataSet.get(0).clone(), primes[0], primes[0], boxSize);
			filter2.stripedBlockAverage(dataSet.get(0).clone(), primes[0], primes[0], boxSize);
		}

		for (final float boxSize : fBoxSizes)
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(data.clone());

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter2.stripedBlockAverage(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (final float boxSize : fBoxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(data.clone());

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter.areaAverageUsingAverages(data, width, height, boxSize);
					time = System.nanoTime() - time;

					final long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"float areaAverageUsingAverages [%dx%d] @ %.1f : %d => stripedBlockAverage %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestAssert.assert_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			TestLog.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"float areaAverageUsingAverages %.1f : %d => stripedBlockAverage %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float areaAverageUsingAverages %d => stripedBlockAverage %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void stripedBlockAverageInternalIsFasterThanAreaAverageInternal()
	{
		TestSettings.assumeSpeedTest();

		final AreaAverageFilter filter = new AreaAverageFilter();
		final AverageFilter filter2 = new AverageFilter();

		final ArrayList<float[]> dataSet = getSpeedData(InternalITER);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		for (final float boxSize : fBoxSizes)
		{
			filter.areaAverageUsingAveragesInternal(dataSet.get(0).clone(), primes[0], primes[0], boxSize);
			filter2.stripedBlockAverageInternal(dataSet.get(0).clone(), primes[0], primes[0], boxSize);
		}

		for (final float boxSize : fBoxSizes)
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(data.clone());

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter2.stripedBlockAverageInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (final float boxSize : fBoxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(data.clone());

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter.areaAverageUsingAveragesInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					final long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"float areaAverageUsingAveragesInternal [%dx%d] @ %.1f : %d => stripedBlockAverageInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestAssert.assert_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			TestLog.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"float areaAverageUsingAveragesInternal %.1f : %d => stripedBlockAverageInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float areaAverageUsingAveragesInternal %d => stripedBlockAverageInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void areaAverageCorrectlyInterpolatesBetweenBlocks()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final int max = 50;
		final float[] data = createData(rg, max, max);
		final AreaAverageFilter filter = new AreaAverageFilter();
		final int n = 30;
		final float[][] results = new float[n + 1][];
		final double[] w = new double[n + 1];
		int count = 0;
		for (int i = 0; i <= n; i++)
		{
			w[count] = i / 10.0;
			results[count] = data.clone();
			filter.areaAverageUsingAverages(results[count], max, max, w[count]);
			count++;
		}

		checkInterpolation(max, n, results, count);
	}

	public void checkInterpolation(int max, int n, float[][] results, int count)
	{
		// Pick some points and see if they are monototically interpolated between integer blocks
		final int[] p = new int[] { 10, 20, 30, 40 };
		for (final int x : p)
			for (final int y : p)
			{
				final int index = y * max + x;
				final double[] yy = new double[count];
				int c = 0;
				for (final float[] data1 : results)
					yy[c++] = data1[index];

				//// Debugging
				//String title = "AreaAverage";
				//ij.gui.Plot plot = new ij.gui.Plot(title, "width", "Mean", w, yy);
				//uk.ac.sussex.gdsc.core.ij.Utils.display(title, plot);

				for (int i = 0; i < n; i += 10)
				{
					final boolean up = yy[i + 10] > yy[i];
					for (int j = i + 1; j < i + 10; j++)
						if (up)
						{
							Assert.assertTrue(yy[j] >= yy[j - 1]);
							Assert.assertTrue(yy[j] <= yy[j + 1]);
						}
						else
						{
							Assert.assertTrue(yy[j] <= yy[j - 1]);
							Assert.assertTrue(yy[j] >= yy[j + 1]);
						}
				}

			}
	}

	@Test
	public void areaAverageInternalCorrectlyInterpolatesBetweenBlocks()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final int max = 50;
		final float[] data = createData(rg, max, max);
		final AreaAverageFilter filter = new AreaAverageFilter();
		final int n = 30;
		final float[][] results = new float[n + 1][];
		final double[] w = new double[n + 1];
		int count = 0;
		for (int i = 0; i <= n; i++)
		{
			w[count] = i / 10.0;
			results[count] = data.clone();
			filter.areaAverageUsingAveragesInternal(results[count], max, max, w[count]);
			count++;
		}

		checkInterpolation(max, n, results, count);
	}

	@Test
	public void areaAverageUsingSumsCorrectlyInterpolatesBetweenBlocks()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final int max = 50;
		final float[] data = createData(rg, max, max);
		final AreaAverageFilter filter = new AreaAverageFilter();
		filter.setSimpleInterpolation(false);
		final int n = 30;
		final float[][] results = new float[n + 1][];
		final double[] w = new double[n + 1];
		int count = 0;
		for (int i = 0; i <= n; i++)
		{
			w[count] = i / 10.0;
			results[count] = data.clone();
			filter.areaAverageUsingSums(results[count], max, max, w[count]);
			count++;
		}

		checkInterpolation(max, n, results, count);
	}

	@Test
	public void areaAverageUsingSumsInternalCorrectlyInterpolatesBetweenBlocks()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final int max = 50;
		final float[] data = createData(rg, max, max);
		final AreaAverageFilter filter = new AreaAverageFilter();
		filter.setSimpleInterpolation(false);
		final int n = 30;
		final float[][] results = new float[n + 1][];
		final double[] w = new double[n + 1];
		int count = 0;
		for (int i = 0; i <= n; i++)
		{
			w[count] = i / 10.0;
			results[count] = data.clone();
			filter.areaAverageUsingSumsInternal(results[count], max, max, w[count]);
			count++;
		}

		checkInterpolation(max, n, results, count);
	}
}
