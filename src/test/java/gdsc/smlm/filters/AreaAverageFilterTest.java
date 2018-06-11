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

import gdsc.smlm.TestSettings;

import java.util.ArrayList;

import org.junit.Assert;
import org.junit.Test;

@SuppressWarnings("deprecation")
public class AreaAverageFilterTest
{
	private gdsc.core.utils.Random rand;

	private boolean debug = false;
	private int ITER = 100;
	private int InternalITER = 300;

	// TODO - The test data should be representative of the final use case
	int[] primes = new int[] { 113, 97, 53, 29 };
	//int[] primes = new int[] { 1024 };
	float[] boxSizes = new float[] { 15.5f, 9.5f, 5.5f, 3.5f, 2.5f, 1.5f, 0.5f };

	private float[] createData(int width, int height)
	{
		float[] data = new float[width * height];
		for (int i = data.length; i-- > 0;)
			data[i] = rand.next();

		return data;
	}

	private double speedUpFactor(long slowTotal, long fastTotal)
	{
		return (1.0 * slowTotal) / fastTotal;
	}

	@Test
	public void areaAverageUsingSumsNxNInternalIsFasterThanAreaAverageNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		AreaAverageFilter filter = new AreaAverageFilter();

		ArrayList<float[]> dataSet = createSpeedData(InternalITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		for (float boxSize : boxSizes)
		{
			filter.areaAverageUsingAveragesInternal(dataSet.get(0).clone(), primes[0], primes[0], boxSize);
			filter.areaAverageUsingSumsInternal(dataSet.get(0).clone(), primes[0], primes[0], boxSize);
		}

		for (float boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(data.clone());

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.areaAverageUsingSumsInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (float boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(data.clone());

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.areaAverageUsingAveragesInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("float areaAverageInternal [%dx%d] @ %.1f : %d => areaAverageUsingSumsInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float areaAverageInternal %.1f : %d => areaAverageUsingSumsInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float areaAverageInternal %d => areaAverageUsingSumsInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void stripedBlockAverageIsFasterThanAreaAverage()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		AreaAverageFilter filter = new AreaAverageFilter();
		AverageFilter filter2 = new AverageFilter();

		ArrayList<float[]> dataSet = createSpeedData(ITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		for (float boxSize : boxSizes)
		{
			filter.areaAverageUsingAverages(dataSet.get(0).clone(), primes[0], primes[0], boxSize);
			filter2.stripedBlockAverage(dataSet.get(0).clone(), primes[0], primes[0], boxSize);
		}

		for (float boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(data.clone());

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter2.stripedBlockAverage(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (float boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(data.clone());

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.areaAverageUsingAverages(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("float areaAverageUsingAverages [%dx%d] @ %.1f : %d => stripedBlockAverage %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float areaAverageUsingAverages %.1f : %d => stripedBlockAverage %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float areaAverageUsingAverages %d => stripedBlockAverage %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void stripedBlockAverageInternalIsFasterThanAreaAverageInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		AreaAverageFilter filter = new AreaAverageFilter();
		AverageFilter filter2 = new AverageFilter();

		ArrayList<float[]> dataSet = createSpeedData(InternalITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		for (float boxSize : boxSizes)
		{
			filter.areaAverageUsingAveragesInternal(dataSet.get(0).clone(), primes[0], primes[0], boxSize);
			filter2.stripedBlockAverageInternal(dataSet.get(0).clone(), primes[0], primes[0], boxSize);
		}

		for (float boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(data.clone());

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter2.stripedBlockAverageInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (float boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(data.clone());

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.areaAverageUsingAveragesInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("float areaAverageUsingAveragesInternal [%dx%d] @ %.1f : %d => stripedBlockAverageInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float areaAverageUsingAveragesInternal %.1f : %d => stripedBlockAverageInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float areaAverageUsingAveragesInternal %d => stripedBlockAverageInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}
	
	private ArrayList<float[]> createSpeedData(int iter)
	{
		ArrayList<float[]> dataSet = new ArrayList<float[]>(iter);
		for (int i = iter; i-- > 0;)
		{
			dataSet.add(createData(primes[0], primes[0]));
		}
		return dataSet;
	}

	@Test
	public void areaAverageCorrectlyInterpolatesBetweenBlocks()
	{
		rand = new gdsc.core.utils.Random(-30051976);
		int max = 50;
		float[] data = createData(max, max);
		AreaAverageFilter filter = new AreaAverageFilter();
		int n = 30;
		float[][] results = new float[n + 1][];
		double[] w = new double[n + 1];
		int count = 0;
		for (int i = 0; i <= n; i++)
		{
			w[count] = i / 10.0;
			results[count] = data.clone();
			filter.areaAverageUsingAverages(results[count], max, max, w[count]);
			count++;
		}

		checkInterpolation(max, n, w, results, count);
	}

	public void checkInterpolation(int max, int n, double[] w, float[][] results, int count)
	{
		// Pick some points and see if they are monototically interpolated between integer blocks
		int[] p = new int[] { 10, 20, 30, 40 };
		for (int x : p)
		{
			for (int y : p)
			{
				int index = y * max + x;
				double[] yy = new double[count];
				int c = 0;
				for (float[] data1 : results)
				{
					yy[c++] = data1[index];
				}

				//// Debugging
				//String title = "AreaAverage";
				//ij.gui.Plot plot = new ij.gui.Plot(title, "width", "Mean", w, yy);
				//gdsc.core.ij.Utils.display(title, plot);

				for (int i = 0; i < n; i += 10)
				{
					boolean up = yy[i + 10] > yy[i];
					for (int j = i + 1; j < i + 10; j++)
					{
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
		}
	}

	@Test
	public void areaAverageInternalCorrectlyInterpolatesBetweenBlocks()
	{
		rand = new gdsc.core.utils.Random(-30051976);
		int max = 50;
		float[] data = createData(max, max);
		AreaAverageFilter filter = new AreaAverageFilter();
		int n = 30;
		float[][] results = new float[n + 1][];
		double[] w = new double[n + 1];
		int count = 0;
		for (int i = 0; i <= n; i++)
		{
			w[count] = i / 10.0;
			results[count] = data.clone();
			filter.areaAverageUsingAveragesInternal(results[count], max, max, w[count]);
			count++;
		}

		checkInterpolation(max, n, w, results, count);
	}

	@Test
	public void areaAverageUsingSumsCorrectlyInterpolatesBetweenBlocks()
	{
		rand = new gdsc.core.utils.Random(-30051976);
		int max = 50;
		float[] data = createData(max, max);
		AreaAverageFilter filter = new AreaAverageFilter();
		filter.setSimpleInterpolation(false);
		int n = 30;
		float[][] results = new float[n + 1][];
		double[] w = new double[n + 1];
		int count = 0;
		for (int i = 0; i <= n; i++)
		{
			w[count] = i / 10.0;
			results[count] = data.clone();
			filter.areaAverageUsingSums(results[count], max, max, w[count]);
			count++;
		}

		checkInterpolation(max, n, w, results, count);
	}

	@Test
	public void areaAverageUsingSumsInternalCorrectlyInterpolatesBetweenBlocks()
	{
		rand = new gdsc.core.utils.Random(-30051976);
		int max = 50;
		float[] data = createData(max, max);
		AreaAverageFilter filter = new AreaAverageFilter();
		filter.setSimpleInterpolation(false);
		int n = 30;
		float[][] results = new float[n + 1][];
		double[] w = new double[n + 1];
		int count = 0;
		for (int i = 0; i <= n; i++)
		{
			w[count] = i / 10.0;
			results[count] = data.clone();
			filter.areaAverageUsingSumsInternal(results[count], max, max, w[count]);
			count++;
		}

		checkInterpolation(max, n, w, results, count);
	}
}
