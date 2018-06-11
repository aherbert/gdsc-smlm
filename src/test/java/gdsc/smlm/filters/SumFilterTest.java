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
import gdsc.smlm.filters.SumFilter;

import java.util.ArrayList;
import java.util.Arrays;

import org.junit.Assert;
import org.junit.Test;
import org.junit.internal.ArrayComparisonFailure;

@SuppressWarnings("deprecation")
public class SumFilterTest
{
	private gdsc.core.utils.Random rand;

	private boolean debug = false;
	private int InternalITER3 = 500;
	private int InternalITER = 50;
	private int ITER3 = 200;
	private int ITER = 20;

	// TODO - The test data should be representative of the final use case
	int[] primes = new int[] { 113, 97, 53, 29 };
	//int[] primes = new int[] { 1024 };
	int[] boxSizes = new int[] { 15, 9, 5, 3, 2 };

	private void floatArrayEquals(String message, float[] data1, float[] data2, int boxSize)
	{
		Assert.assertArrayEquals(message, data1, data2, boxSize * boxSize * 1e-3f);
	}

	private void intArrayEquals(String message, int[] data1, int[] data2, int boxSize)
	{
		Assert.assertArrayEquals(message, data1, data2);
	}

	private double speedUpFactor(long slowTotal, long fastTotal)
	{
		return (1.0 * slowTotal) / fastTotal;
	}

	// COPY CODE FROM HERE...
	@Test
	public void floatBlockSumNxNInternalAndRollingBlockSumNxNInternalReturnSameResult()
	{
		SumFilter filter = new SumFilter();

		for (int width : primes)
			for (int height : primes)
				for (int boxSize : boxSizes)
					floatCompareBlockSumNxNInternalAndRollingBlockSumNxNInternal(filter, width, height, boxSize);
	}

	private void floatCompareBlockSumNxNInternalAndRollingBlockSumNxNInternal(SumFilter filter, int width,
			int height, int boxSize) throws ArrayComparisonFailure
	{
		rand = new gdsc.core.utils.Random(-30051976);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.blockSumNxNInternal(data1, width, height, boxSize);
		filter.rollingBlockSumNxNInternal(data2, width, height, boxSize);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d] @ %d", width, height, boxSize), data1, data2,
				boxSize);
	}

	@Test
	public void floatBlockSumNxNInternalAndStripedBlockSumNxNInternalReturnSameResult()
	{
		SumFilter filter = new SumFilter();

		for (int width : primes)
			for (int height : primes)
				for (int boxSize : boxSizes)
					floatCompareBlockSumNxNInternalAndStripedBlockSumNxNInternal(filter, width, height, boxSize);
	}

	private void floatCompareBlockSumNxNInternalAndStripedBlockSumNxNInternal(SumFilter filter, int width,
			int height, int boxSize) throws ArrayComparisonFailure
	{
		rand = new gdsc.core.utils.Random(-30051976);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.blockSumNxNInternal(data1, width, height, boxSize);
		filter.stripedBlockSumNxNInternal(data2, width, height, boxSize);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d] @ %d", width, height, boxSize), data1, data2, 0);
	}

	@Test
	public void floatBlockSum3x3InternalAndRollingBlockSumNxNInternalReturnSameResult()
	{
		SumFilter filter = new SumFilter();

		for (int width : primes)
			for (int height : primes)
				floatCompareBlockSum3x3InternalAndRollingBlockSumNxNInternal(filter, width, height);
	}

	private void floatCompareBlockSum3x3InternalAndRollingBlockSumNxNInternal(SumFilter filter, int width,
			int height) throws ArrayComparisonFailure
	{
		rand = new gdsc.core.utils.Random(-30051977);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.blockSum3x3Internal(data1, width, height);
		filter.rollingBlockSumNxNInternal(data2, width, height, 1);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);
	}

	@Test
	public void floatRollingBlockSumNxNInternalAndRollingBlockSumNxNInternalTransposedReturnSameResult()
	{
		SumFilter filter = new SumFilter();

		for (int width : primes)
			for (int height : primes)
				for (int boxSize : boxSizes)
					floatCompareRollingBlockSumNxNInternalAndRollingBlockSumNxNInternalTransposed(filter,
							width, height, boxSize);
	}

	private void floatCompareRollingBlockSumNxNInternalAndRollingBlockSumNxNInternalTransposed(
			SumFilter filter, int width, int height, int boxSize) throws ArrayComparisonFailure
	{
		rand = new gdsc.core.utils.Random(-30051977);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.rollingBlockSumNxNInternal(data1, width, height, boxSize);
		filter.rollingBlockSumNxNInternalTransposed(data2, width, height, boxSize);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d] @ %d", width, height, boxSize), data1, data2,
				boxSize);
	}

	private float[] floatClone(float[] data1)
	{
		float[] data2 = Arrays.copyOf(data1, data1.length);
		return data2;
	}

	@Test
	public void floatRollingBlockSumNxNInternalIsFasterThanBlockSumNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		SumFilter filter = new SumFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.blockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("float blockSumNxNInternal [%dx%d] @ %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float blockSumNxNInternal %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float blockSumNxNInternal %d => rollingBlockSumNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	private ArrayList<float[]> floatCreateSpeedData(int iter)
	{
		ArrayList<float[]> dataSet = new ArrayList<float[]>(iter);
		for (int i = iter; i-- > 0;)
		{
			dataSet.add(floatCreateData(primes[0], primes[0]));
		}
		return dataSet;
	}

	@Test
	public void floatStripedBlockSumNxNInternalIsFasterThanBlockSumNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		SumFilter filter = new SumFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.stripedBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.stripedBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.blockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("float blockSumNxNInternal [%dx%d] @ %d : %d => stripedBlockSumNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float blockSumNxNInternal %d : %d => stripedBlockSumNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float blockSumNxNInternal %d => stripedBlockSumNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockSumNxNInternalIsFasterThanStripedBlockSumNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		SumFilter filter = new SumFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.stripedBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.stripedBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("float stripedBlockSumNxNInternal [%dx%d] @ %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf(
					"float stripedBlockSumNxNInternal %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float stripedBlockSumNxNInternal %d => rollingBlockSumNxNInternal %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatBlockSum3x3InternalAndBlockSumNxNInternalReturnSameResult()
	{
		SumFilter filter = new SumFilter();

		for (int width : primes)
			for (int height : primes)
				floatCompareBlockSum3x3InternalAndBlockSumNxNInternal(filter, width, height);
	}

	private void floatCompareBlockSum3x3InternalAndBlockSumNxNInternal(SumFilter filter, int width,
			int height) throws ArrayComparisonFailure
	{
		rand = new gdsc.core.utils.Random(-30051977);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.blockSum3x3Internal(data1, width, height);
		filter.blockSumNxNInternal(data2, width, height, 1);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);
	}

	@Test
	public void floatBlockSum3x3InternalIsFasterThanBlockSumNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.blockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.blockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.blockSumNxNInternal(data, width, height, 1);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf(
							"float blockSumNxNInternal [%dx%d] %d => blockSum3x3Internal %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float blockSumNxNInternal %d => blockSum3x3Internal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockSum3x3InternalIsFasterThanBlockSum3x3Internal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.blockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.rollingBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.blockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf(
							"float blockSum3x3Internal [%dx%d] %d => rollingBlockSum3x3Internal %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float blockSum3x3Internal %d => rollingBlockSum3x3Internal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatStripedBlockSum3x3InternalIsFasterThanBlockSum3x3Internal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.stripedBlockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.blockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.stripedBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.blockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf(
							"float blockSum3x3Internal [%dx%d] %d => stripedBlockSum3x3Internal %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float blockSum3x3Internal %d => stripedBlockSum3x3Internal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockSum3x3InternalIsFasterThanStripedBlockSum3x3Internal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.stripedBlockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.rollingBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.stripedBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out
							.printf("float stripedBlockSum3x3Internal [%dx%d] %d => rollingBlockSum3x3Internal %d = %.2fx\n",
									width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float stripedBlockSum3x3Internal %d => rollingBlockSum3x3Internal %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockSum3x3InternalAndRollingBlockSumNxNInternalReturnSameResult()
	{
		SumFilter filter = new SumFilter();

		for (int width : primes)
			for (int height : primes)
				floatCompareRollingBlockSum3x3InternalAndRollingBlockSumNxNInternal(filter, width, height);
	}

	private void floatCompareRollingBlockSum3x3InternalAndRollingBlockSumNxNInternal(SumFilter filter,
			int width, int height) throws ArrayComparisonFailure
	{
		rand = new gdsc.core.utils.Random(-30051977);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.rollingBlockSum3x3Internal(data1, width, height);
		filter.rollingBlockSumNxNInternal(data2, width, height, 1);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);
	}

	@Test
	public void floatRollingBlockSum3x3InternalIsFasterThanRollingBlockSumNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.rollingBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], 1);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.rollingBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.rollingBlockSumNxNInternal(data, width, height, 1);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out
							.printf("float rollingBlockSumNxNInternal [%dx%d] %d => rollingBlockSum3x3Internal %d = %.2fx\n",
									width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float rollingBlockSumNxNInternal %d => rollingBlockSum3x3Internal %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatStripedBlockSum3x3InternalIsFasterThanStripedBlockSumNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.stripedBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.stripedBlockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.stripedBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.stripedBlockSumNxNInternal(data, width, height, 1);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out
							.printf("float stripedBlockSumNxNInternal [%dx%d] %d => stripedBlockSum3x3Internal %d = %.2fx\n",
									width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float stripedBlockSumNxNInternal %d => stripedBlockSum3x3Internal %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockSumNxNInternalIsFasterThanRollingBlockSumNxNInternalTransposed()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		SumFilter filter = new SumFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingBlockSumNxNInternalTransposed(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.rollingBlockSumNxNInternalTransposed(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("float rollingBlockSumNxNInternalTransposed [%dx%d] @ %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			if (debug)
				System.out
						.printf("float rollingBlockSumNxNInternalTransposed %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
								boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		System.out.printf(
				"float rollingBlockSumNxNInternalTransposed %d => rollingBlockSumNxNInternal %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	private float[] floatCreateData(int width, int height)
	{
		float[] data = new float[width * height];
		for (int i = data.length; i-- > 0;)
			data[i] = i;

		rand.shuffle(data);

		return data;
	}

	@Test
	public void floatBlockSumNxNAndStripedBlockSumNxNReturnSameResult()
	{
		SumFilter filter = new SumFilter();

		for (int width : primes)
			for (int height : primes)
				for (int boxSize : boxSizes)
					floatCompareBlockSumNxNAndStripedBlockSumNxN(filter, width, height, boxSize);
	}

	private void floatCompareBlockSumNxNAndStripedBlockSumNxN(SumFilter filter, int width, int height,
			int boxSize) throws ArrayComparisonFailure
	{
		rand = new gdsc.core.utils.Random(-30051976);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.blockSumNxN(data1, width, height, boxSize);
		filter.rollingBlockSumNxN(data2, width, height, boxSize);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d] @ %d", width, height, boxSize), data1, data2, 0);
	}

	@Test
	public void floatBlockSumNxNAndRollingBlockSumNxNReturnSameResult()
	{
		SumFilter filter = new SumFilter();

		for (int width : primes)
			for (int height : primes)
				for (int boxSize : boxSizes)
					floatCompareBlockSumNxNAndRollingBlockSumNxN(filter, width, height, boxSize);
	}

	private void floatCompareBlockSumNxNAndRollingBlockSumNxN(SumFilter filter, int width, int height,
			int boxSize) throws ArrayComparisonFailure
	{
		rand = new gdsc.core.utils.Random(-30051976);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.blockSumNxN(data1, width, height, boxSize);
		filter.rollingBlockSumNxN(data2, width, height, boxSize);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d] @ %d", width, height, boxSize), data1, data2, 0);
	}

	@Test
	public void floatBlockSumInternalNxNIsFasterThanBlockSumNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		SumFilter filter = new SumFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.blockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.blockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.blockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"float blockSumNxN [%dx%d] @ %d : %d => blockSumNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float blockSumNxN %d : %d => blockSumNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float blockSumNxN %d => blockSumNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatStripedBlockSumNxNIsFasterThanBlockSumNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		SumFilter filter = new SumFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.stripedBlockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.stripedBlockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.blockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"float blockSumNxN [%dx%d] @ %d : %d => stripedBlockSumNxN %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float blockSumNxN %d : %d => stripedBlockSumNxN %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float blockSumNxN %d => stripedBlockSumNxN %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatStripedBlockSumInternalNxNIsFasterThanStripedBlockSumNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		SumFilter filter = new SumFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.stripedBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.stripedBlockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.stripedBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.stripedBlockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("float stripedBlockSumNxN [%dx%d] @ %d : %d => stripedBlockSumNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float stripedBlockSumNxN %d : %d => stripedBlockSumNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float stripedBlockSumNxN %d => stripedBlockSumNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockSumNxNIsFasterThanBlockSumNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		SumFilter filter = new SumFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingBlockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.rollingBlockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.blockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"float blockSumNxN [%dx%d] @ %d : %d => rollingBlockSumNxN %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float blockSumNxN %d : %d => rollingBlockSumNxN %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float blockSumNxN %d => rollingBlockSumNxN %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockSumInternalNxNIsFasterThanRollingBlockSumNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		SumFilter filter = new SumFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingBlockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.rollingBlockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("float rollingBlockSumNxN [%dx%d] @ %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float rollingBlockSumNxN %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float rollingBlockSumNxN %d => rollingBlockSumNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatBlockSum3x3AndBlockSumNxNReturnSameResult()
	{
		SumFilter filter = new SumFilter();

		for (int width : primes)
			for (int height : primes)
				floatCompareBlockSum3x3AndBlockSumNxN(filter, width, height);
	}

	private void floatCompareBlockSum3x3AndBlockSumNxN(SumFilter filter, int width, int height)
			throws ArrayComparisonFailure
	{
		rand = new gdsc.core.utils.Random(-30051977);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.blockSum3x3(data1, width, height);
		filter.blockSumNxN(data2, width, height, 1);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);
	}

	@Test
	public void floatBlockSum3x3IsFasterThanBlockSumNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.blockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.blockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.blockSumNxN(data, width, height, 1);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float blockSumNxN [%dx%d] %d => blockSum3x3 %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float blockSumNxN %d => blockSum3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatStripedBlockSum3x3AndStripedBlockSumNxNReturnSameResult()
	{
		SumFilter filter = new SumFilter();

		for (int width : primes)
			for (int height : primes)
				floatCompareStripedBlockSum3x3AndStripedBlockSumNxN(filter, width, height);
	}

	private void floatCompareStripedBlockSum3x3AndStripedBlockSumNxN(SumFilter filter, int width, int height)
			throws ArrayComparisonFailure
	{
		rand = new gdsc.core.utils.Random(-30051977);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.stripedBlockSum3x3(data1, width, height);
		filter.stripedBlockSumNxN(data2, width, height, 1);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);
	}

	@Test
	public void floatStripedBlockSum3x3IsFasterThanStripedBlockSumNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.stripedBlockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.stripedBlockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.stripedBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.stripedBlockSumNxN(data, width, height, 1);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float stripedBlockSumNxN [%dx%d] %d => stripedBlockSum3x3 %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		stripedBlockTime, time), stripedBlockTime < time);
			}
		System.out.printf("float stripedBlockSumNxN %d => stripedBlockSum3x3 %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockSum3x3AndRollingBlockSumNxNReturnSameResult()
	{
		SumFilter filter = new SumFilter();

		for (int width : primes)
			for (int height : primes)
				floatCompareRollingBlockSum3x3AndRollingBlockSumNxN(filter, width, height);
	}

	private void floatCompareRollingBlockSum3x3AndRollingBlockSumNxN(SumFilter filter, int width, int height)
			throws ArrayComparisonFailure
	{
		rand = new gdsc.core.utils.Random(-30051977);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.rollingBlockSum3x3(data1, width, height);
		filter.rollingBlockSumNxN(data2, width, height, 1);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);
	}

	@Test
	public void floatRollingBlockSum3x3IsFasterThanRollingBlockSumNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.rollingBlockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.rollingBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.rollingBlockSumNxN(data, width, height, 1);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float rollingBlockSumNxN [%dx%d] %d => rollingBlockSum3x3 %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		rollingBlockTime, time), rollingBlockTime < time);
			}
		System.out.printf("float rollingBlockSumNxN %d => rollingBlockSum3x3 %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockSum3x3IsFasterThanBlockSum3x3()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.blockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.rollingBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.blockSum3x3(data, width, height);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float blockSum3x3 [%dx%d] %d => rollingBlockSum3x3 %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float blockSum3x3 %d => rollingBlockSum3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatStripedBlockSum3x3IsFasterThanBlockSum3x3()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.stripedBlockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.blockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.stripedBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.blockSum3x3(data, width, height);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float blockSum3x3 [%dx%d] %d => stripedBlockSum3x3 %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float blockSum3x3 %d => stripedBlockSum3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockSum3x3IsFasterThanStripedBlockSum3x3()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.stripedBlockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.rollingBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.stripedBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float stripedBlockSum3x3 [%dx%d] %d => rollingBlockSum3x3 %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float stripedBlockSum3x3 %d => rollingBlockSum3x3 %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void intBlockSumNxNInternalAndRollingBlockSumNxNReturnSameResult()
	{
		SumFilter filter = new SumFilter();

		for (int width : primes)
			for (int height : primes)
				for (int boxSize : boxSizes)
					intCompareBlockSumNxNInternalAndRollingBlockSumNxNInternal(filter, width, height, boxSize);
	}

	private void intCompareBlockSumNxNInternalAndRollingBlockSumNxNInternal(SumFilter filter, int width,
			int height, int boxSize) throws ArrayComparisonFailure
	{
		rand = new gdsc.core.utils.Random(-30051976);
		int[] data1 = intCreateData(width, height);
		int[] data2 = intClone(data1);

		filter.blockSumNxNInternal(data1, width, height, boxSize);
		filter.rollingBlockSumNxNInternal(data2, width, height, boxSize);

		intArrayEquals(String.format("Arrays do not match: [%dx%d] @ %d", width, height, boxSize), data1, data2,
				boxSize);
	}

	@Test
	public void intBlockSumNxNInternalAndStripedBlockSumNxNReturnSameResult()
	{
		SumFilter filter = new SumFilter();

		for (int width : primes)
			for (int height : primes)
				for (int boxSize : boxSizes)
					intCompareBlockSumNxNInternalAndStripedBlockSumNxNInternal(filter, width, height, boxSize);
	}

	private void intCompareBlockSumNxNInternalAndStripedBlockSumNxNInternal(SumFilter filter, int width,
			int height, int boxSize) throws ArrayComparisonFailure
	{
		rand = new gdsc.core.utils.Random(-30051976);
		int[] data1 = intCreateData(width, height);
		int[] data2 = intClone(data1);

		filter.blockSumNxNInternal(data1, width, height, boxSize);
		filter.stripedBlockSumNxNInternal(data2, width, height, boxSize);

		intArrayEquals(String.format("Arrays do not match: [%dx%d] @ %d", width, height, boxSize), data1, data2, 0);
	}

	@Test
	public void intBlockSum3x3InternalAndRollingBlockSumNxNReturnSameResult()
	{
		SumFilter filter = new SumFilter();

		for (int width : primes)
			for (int height : primes)
				intCompareBlockSum3x3InternalAndRollingBlockSumNxNInternal(filter, width, height);
	}

	private void intCompareBlockSum3x3InternalAndRollingBlockSumNxNInternal(SumFilter filter, int width,
			int height) throws ArrayComparisonFailure
	{
		rand = new gdsc.core.utils.Random(-30051977);
		int[] data1 = intCreateData(width, height);
		int[] data2 = intClone(data1);

		filter.blockSum3x3Internal(data1, width, height);
		filter.rollingBlockSumNxNInternal(data2, width, height, 1);

		intArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);
	}

	@Test
	public void intRollingBlockSumNxNInternalAndRollingBlockSumNxNInternalTransposedReturnSameResult()
	{
		SumFilter filter = new SumFilter();

		for (int width : primes)
			for (int height : primes)
				for (int boxSize : boxSizes)
					intCompareRollingBlockSumNxNInternalAndRollingBlockSumNxNInternalTransposed(filter, width,
							height, boxSize);
	}

	private void intCompareRollingBlockSumNxNInternalAndRollingBlockSumNxNInternalTransposed(
			SumFilter filter, int width, int height, int boxSize) throws ArrayComparisonFailure
	{
		rand = new gdsc.core.utils.Random(-30051977);
		int[] data1 = intCreateData(width, height);
		int[] data2 = intClone(data1);

		filter.rollingBlockSumNxNInternal(data1, width, height, boxSize);
		filter.rollingBlockSumNxNInternalTransposed(data2, width, height, boxSize);

		intArrayEquals(String.format("Arrays do not match: [%dx%d] @ %d", width, height, boxSize), data1, data2,
				boxSize);
	}

	private int[] intClone(int[] data1)
	{
		int[] data2 = Arrays.copyOf(data1, data1.length);
		return data2;
	}

	@Test
	public void intRollingBlockSumNxNInternalIsFasterThanBlockSumNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		SumFilter filter = new SumFilter();

		ArrayList<int[]> dataSet = intCreateSpeedData(InternalITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
					for (int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (int[] data : dataSet2)
						filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
					for (int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (int[] data : dataSet2)
						filter.blockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("int blockSumNxNInternal [%dx%d] @ %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("int blockSumNxNInternal %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("int blockSumNxNInternal %d => rollingBlockSumNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	private ArrayList<int[]> intCreateSpeedData(int iter)
	{
		ArrayList<int[]> dataSet = new ArrayList<int[]>(iter);
		for (int i = iter; i-- > 0;)
		{
			dataSet.add(intCreateData(primes[0], primes[0]));
		}
		return dataSet;
	}

	@Test
	public void intStripedBlockSumNxNInternalIsFasterThanBlockSumNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		SumFilter filter = new SumFilter();

		ArrayList<int[]> dataSet = intCreateSpeedData(InternalITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.stripedBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
					for (int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (int[] data : dataSet2)
						filter.stripedBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
					for (int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (int[] data : dataSet2)
						filter.blockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("int blockSumNxNInternal [%dx%d] @ %d : %d => stripedBlockSumNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("int blockSumNxNInternal %d : %d => stripedBlockSumNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("int blockSumNxNInternal %d => stripedBlockSumNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void intRollingBlockSumNxNInternalIsFasterThanStripedBlockSumNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		SumFilter filter = new SumFilter();

		ArrayList<int[]> dataSet = intCreateSpeedData(InternalITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.stripedBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
					for (int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (int[] data : dataSet2)
						filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
					for (int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (int[] data : dataSet2)
						filter.stripedBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("int stripedBlockSumNxNInternal [%dx%d] @ %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf(
					"int stripedBlockSumNxNInternal %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("int stripedBlockSumNxNInternal %d => rollingBlockSumNxNInternal %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void intBlockSum3x3InternalAndBlockSumNxNInternalReturnSameResult()
	{
		SumFilter filter = new SumFilter();

		for (int width : primes)
			for (int height : primes)
				intCompareBlockSum3x3InternalAndBlockSumNxNInternal(filter, width, height);
	}

	private void intCompareBlockSum3x3InternalAndBlockSumNxNInternal(SumFilter filter, int width, int height)
			throws ArrayComparisonFailure
	{
		rand = new gdsc.core.utils.Random(-30051977);
		int[] data1 = intCreateData(width, height);
		int[] data2 = intClone(data1);

		filter.blockSum3x3Internal(data1, width, height);
		filter.blockSumNxNInternal(data2, width, height, 1);

		intArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);
	}

	@Test
	public void intBlockSum3x3InternalIsFasterThanBlockSumNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<int[]> dataSet = intCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.blockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.blockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.blockSumNxNInternal(data, width, height, 1);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("int blockSumNxNInternal [%dx%d] %d => blockSum3x3Internal %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("int blockSumNxNInternal %d => blockSum3x3Internal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void intRollingBlockSum3x3InternalIsFasterThanBlockSum3x3Internal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<int[]> dataSet = intCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);
		filter.blockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.rollingBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.blockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf(
							"int blockSum3x3Internal [%dx%d] %d => rollingBlockSum3x3Internal %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("int blockSum3x3Internal %d => rollingBlockSum3x3Internal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void intStripedBlockSum3x3InternalIsFasterThanBlockSum3x3Internal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<int[]> dataSet = intCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.stripedBlockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);
		filter.blockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.stripedBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.blockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf(
							"int blockSum3x3Internal [%dx%d] %d => stripedBlockSum3x3Internal %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("int blockSum3x3Internal %d => stripedBlockSum3x3Internal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void intRollingBlockSum3x3InternalIsFasterThanStripedBlockSum3x3Internal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<int[]> dataSet = intCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);
		filter.stripedBlockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.rollingBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.stripedBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out
							.printf("int stripedBlockSum3x3Internal [%dx%d] %d => rollingBlockSum3x3Internal %d = %.2fx\n",
									width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("int stripedBlockSum3x3Internal %d => rollingBlockSum3x3Internal %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void intRollingBlockSum3x3InternalAndRollingBlockSumNxNInternalReturnSameResult()
	{
		SumFilter filter = new SumFilter();

		for (int width : primes)
			for (int height : primes)
				intCompareRollingBlockSum3x3InternalAndRollingBlockSumNxNInternal(filter, width, height);
	}

	private void intCompareRollingBlockSum3x3InternalAndRollingBlockSumNxNInternal(SumFilter filter,
			int width, int height) throws ArrayComparisonFailure
	{
		rand = new gdsc.core.utils.Random(-30051977);
		int[] data1 = intCreateData(width, height);
		int[] data2 = intClone(data1);

		filter.rollingBlockSum3x3Internal(data1, width, height);
		filter.rollingBlockSumNxNInternal(data2, width, height, 1);

		intArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);
	}

	@Test
	public void intRollingBlockSum3x3InternalIsFasterThanRollingBlockSumNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<int[]> dataSet = intCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);
		filter.rollingBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], 1);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.rollingBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.rollingBlockSumNxNInternal(data, width, height, 1);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out
							.printf("int rollingBlockSumNxNInternal [%dx%d] %d => rollingBlockSum3x3Internal %d = %.2fx\n",
									width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("int rollingBlockSumNxNInternal %d => rollingBlockSum3x3Internal %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void intStripedBlockSum3x3InternalIsFasterThanStripedBlockSumNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<int[]> dataSet = intCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.stripedBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.stripedBlockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.stripedBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.stripedBlockSumNxNInternal(data, width, height, 1);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out
							.printf("int stripedBlockSumNxNInternal [%dx%d] %d => stripedBlockSum3x3Internal %d = %.2fx\n",
									width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("int stripedBlockSumNxNInternal %d => stripedBlockSum3x3Internal %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void intRollingBlockSumNxNInternalIsFasterThanRollingBlockSumNxNInternalTransposed()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		SumFilter filter = new SumFilter();

		ArrayList<int[]> dataSet = intCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingBlockSumNxNInternalTransposed(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
					for (int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (int[] data : dataSet2)
						filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
					for (int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (int[] data : dataSet2)
						filter.rollingBlockSumNxNInternalTransposed(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("int rollingBlockSumNxNInternalTransposed [%dx%d] @ %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			if (debug)
				System.out
						.printf("int rollingBlockSumNxNInternalTransposed %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
								boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		System.out.printf(
				"int rollingBlockSumNxNInternalTransposed %d => rollingBlockSumNxNInternal %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	private int[] intCreateData(int width, int height)
	{
		int[] data = new int[width * height];
		for (int i = data.length; i-- > 0;)
			data[i] = i;

		rand.shuffle(data);

		return data;
	}

	@Test
	public void intBlockSumNxNAndStripedBlockSumNxNReturnSameResult()
	{
		SumFilter filter = new SumFilter();

		for (int width : primes)
			for (int height : primes)
				for (int boxSize : boxSizes)
					intCompareBlockSumNxNAndStripedBlockSumNxN(filter, width, height, boxSize);
	}

	private void intCompareBlockSumNxNAndStripedBlockSumNxN(SumFilter filter, int width, int height,
			int boxSize) throws ArrayComparisonFailure
	{
		rand = new gdsc.core.utils.Random(-30051976);
		int[] data1 = intCreateData(width, height);
		int[] data2 = intClone(data1);

		filter.blockSumNxN(data1, width, height, boxSize);
		filter.rollingBlockSumNxN(data2, width, height, boxSize);

		intArrayEquals(String.format("Arrays do not match: [%dx%d] @ %d", width, height, boxSize), data1, data2, 0);
	}

	@Test
	public void intBlockSumNxNAndRollingBlockSumNxNReturnSameResult()
	{
		SumFilter filter = new SumFilter();

		for (int width : primes)
			for (int height : primes)
				for (int boxSize : boxSizes)
					intCompareBlockSumNxNAndRollingBlockSumNxN(filter, width, height, boxSize);
	}

	private void intCompareBlockSumNxNAndRollingBlockSumNxN(SumFilter filter, int width, int height,
			int boxSize) throws ArrayComparisonFailure
	{
		rand = new gdsc.core.utils.Random(-30051976);
		int[] data1 = intCreateData(width, height);
		int[] data2 = intClone(data1);

		filter.blockSumNxN(data1, width, height, boxSize);
		filter.rollingBlockSumNxN(data2, width, height, boxSize);

		intArrayEquals(String.format("Arrays do not match: [%dx%d] @ %d", width, height, boxSize), data1, data2, 0);
	}

	@Test
	public void intBlockSumInternalNxNIsFasterThanBlockSumNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		SumFilter filter = new SumFilter();

		ArrayList<int[]> dataSet = intCreateSpeedData(ITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.blockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
					for (int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (int[] data : dataSet2)
						filter.blockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
					for (int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (int[] data : dataSet2)
						filter.blockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"int blockSumNxN [%dx%d] @ %d : %d => blockSumNxNInternal %d = %.2fx\n", width,
								height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("int blockSumNxN %d : %d => blockSumNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("int blockSumNxN %d => blockSumNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void intStripedBlockSumNxNIsFasterThanBlockSumNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		SumFilter filter = new SumFilter();

		ArrayList<int[]> dataSet = intCreateSpeedData(ITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.stripedBlockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
					for (int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (int[] data : dataSet2)
						filter.stripedBlockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
					for (int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (int[] data : dataSet2)
						filter.blockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"int blockSumNxN [%dx%d] @ %d : %d => stripedBlockSumNxN %d = %.2fx\n", width,
								height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("int blockSumNxN %d : %d => stripedBlockSumNxN %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("int blockSumNxN %d => stripedBlockSumNxN %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void intStripedBlockSumInternalNxNIsFasterThanStripedBlockSumNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		SumFilter filter = new SumFilter();

		ArrayList<int[]> dataSet = intCreateSpeedData(ITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.stripedBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.stripedBlockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
					for (int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (int[] data : dataSet2)
						filter.stripedBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
					for (int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (int[] data : dataSet2)
						filter.stripedBlockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("int stripedBlockSumNxN [%dx%d] @ %d : %d => stripedBlockSumNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("int stripedBlockSumNxN %d : %d => stripedBlockSumNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("int stripedBlockSumNxN %d => stripedBlockSumNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void intRollingBlockSumNxNIsFasterThanBlockSumNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		SumFilter filter = new SumFilter();

		ArrayList<int[]> dataSet = intCreateSpeedData(ITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingBlockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
					for (int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (int[] data : dataSet2)
						filter.rollingBlockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
					for (int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (int[] data : dataSet2)
						filter.blockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"int blockSumNxN [%dx%d] @ %d : %d => rollingBlockSumNxN %d = %.2fx\n", width,
								height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("int blockSumNxN %d : %d => rollingBlockSumNxN %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("int blockSumNxN %d => rollingBlockSumNxN %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void intRollingBlockSumInternalNxNIsFasterThanRollingBlockSumNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		SumFilter filter = new SumFilter();

		ArrayList<int[]> dataSet = intCreateSpeedData(ITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingBlockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
					for (int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (int[] data : dataSet2)
						filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
					for (int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (int[] data : dataSet2)
						filter.rollingBlockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("int rollingBlockSumNxN [%dx%d] @ %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("int rollingBlockSumNxN %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("int rollingBlockSumNxN %d => rollingBlockSumNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void intBlockSum3x3AndBlockSumNxNReturnSameResult()
	{
		SumFilter filter = new SumFilter();

		for (int width : primes)
			for (int height : primes)
				intCompareBlockSum3x3AndBlockSumNxN(filter, width, height);
	}

	private void intCompareBlockSum3x3AndBlockSumNxN(SumFilter filter, int width, int height)
			throws ArrayComparisonFailure
	{
		rand = new gdsc.core.utils.Random(-30051977);
		int[] data1 = intCreateData(width, height);
		int[] data2 = intClone(data1);

		filter.blockSum3x3(data1, width, height);
		filter.blockSumNxN(data2, width, height, 1);

		intArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);
	}

	@Test
	public void intBlockSum3x3IsFasterThanBlockSumNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<int[]> dataSet = intCreateSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.blockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.blockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.blockSumNxN(data, width, height, 1);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("int blockSumNxN [%dx%d] %d => blockSum3x3 %d = %.2fx\n", width, height,
							time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("int blockSumNxN %d => blockSum3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void intStripedBlockSum3x3AndStripedBlockSumNxNReturnSameResult()
	{
		SumFilter filter = new SumFilter();

		for (int width : primes)
			for (int height : primes)
				intCompareStripedBlockSum3x3AndStripedBlockSumNxN(filter, width, height);
	}

	private void intCompareStripedBlockSum3x3AndStripedBlockSumNxN(SumFilter filter, int width, int height)
			throws ArrayComparisonFailure
	{
		rand = new gdsc.core.utils.Random(-30051977);
		int[] data1 = intCreateData(width, height);
		int[] data2 = intClone(data1);

		filter.stripedBlockSum3x3(data1, width, height);
		filter.stripedBlockSumNxN(data2, width, height, 1);

		intArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);
	}

	@Test
	public void intStripedBlockSum3x3IsFasterThanStripedBlockSumNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<int[]> dataSet = intCreateSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.stripedBlockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.stripedBlockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.stripedBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.stripedBlockSumNxN(data, width, height, 1);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("int stripedBlockSumNxN [%dx%d] %d => stripedBlockSum3x3 %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		stripedBlockTime, time), stripedBlockTime < time);
			}
		System.out.printf("int stripedBlockSumNxN %d => stripedBlockSum3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void intRollingBlockSum3x3AndRollingBlockSumNxNReturnSameResult()
	{
		SumFilter filter = new SumFilter();

		for (int width : primes)
			for (int height : primes)
				intCompareRollingBlockSum3x3AndRollingBlockSumNxN(filter, width, height);
	}

	private void intCompareRollingBlockSum3x3AndRollingBlockSumNxN(SumFilter filter, int width, int height)
			throws ArrayComparisonFailure
	{
		rand = new gdsc.core.utils.Random(-30051977);
		int[] data1 = intCreateData(width, height);
		int[] data2 = intClone(data1);

		filter.rollingBlockSum3x3(data1, width, height);
		filter.rollingBlockSumNxN(data2, width, height, 1);

		intArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);
	}

	@Test
	public void intRollingBlockSum3x3IsFasterThanRollingBlockSumNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<int[]> dataSet = intCreateSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.rollingBlockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.rollingBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.rollingBlockSumNxN(data, width, height, 1);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("int rollingBlockSumNxN [%dx%d] %d => rollingBlockSum3x3 %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		rollingBlockTime, time), rollingBlockTime < time);
			}
		System.out.printf("int rollingBlockSumNxN %d => rollingBlockSum3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void intRollingBlockSum3x3IsFasterThanBlockSum3x3()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<int[]> dataSet = intCreateSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);
		filter.blockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.rollingBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.blockSum3x3(data, width, height);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("int blockSum3x3 [%dx%d] %d => rollingBlockSum3x3 %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("int blockSum3x3 %d => rollingBlockSum3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void intStripedBlockSum3x3IsFasterThanBlockSum3x3()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<int[]> dataSet = intCreateSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.stripedBlockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);
		filter.blockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.stripedBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.blockSum3x3(data, width, height);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("int blockSum3x3 [%dx%d] %d => stripedBlockSum3x3 %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("int blockSum3x3 %d => stripedBlockSum3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void intRollingBlockSum3x3IsFasterThanStripedBlockSum3x3()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-30051977);

		SumFilter filter = new SumFilter();

		ArrayList<int[]> dataSet = intCreateSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);
		filter.stripedBlockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.rollingBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<int[]> dataSet2 = new ArrayList<int[]>(dataSet.size());
				for (int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (int[] data : dataSet2)
					filter.stripedBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("int stripedBlockSum3x3 [%dx%d] %d => rollingBlockSum3x3 %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("int stripedBlockSum3x3 %d => rollingBlockSum3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}
}
