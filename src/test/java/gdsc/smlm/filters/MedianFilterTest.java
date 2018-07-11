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

import java.util.ArrayList;

import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Test;
import org.junit.internal.ArrayComparisonFailure;

import gdsc.test.TestAssert;
import gdsc.test.TestSettings;

@SuppressWarnings({ "javadoc" })
public class MedianFilterTest extends AbstractFilterTest
{
	private int InternalITER3 = 200;
	private int InternalITER = 20;
	private int ITER3 = 100;
	private int ITER = 10;

	private void floatArrayEquals(float[] data1, float[] data2, int boxSize, String format, Object... args)
	{
		try
		{
			TestAssert.assertArrayEqualsRelative(null, data1, data2, boxSize * boxSize * 1e-3);
		}
		catch (AssertionError e)
		{
			throw new AssertionError(String.format(format, args), e);
		}
	}

	@Test
	public void floatBlockMedianNxNInternalAndRollingMedianNxNInternalReturnSameResult()
	{
		RandomGenerator rg = TestSettings.getRandomGenerator();
		MedianFilter filter = new MedianFilter();
		for (int width : primes)
			for (int height : primes)
				for (int boxSize : boxSizes)
					floatCompareBlockMedianNxNInternalAndRollingMedianNxNInternal(rg, filter, width, height, boxSize);
	}

	private void floatCompareBlockMedianNxNInternalAndRollingMedianNxNInternal(RandomGenerator rg, MedianFilter filter,
			int width, int height, int boxSize) throws ArrayComparisonFailure
	{
		float[] data1 = createData(rg, width, height);
		float[] data2 = floatClone(data1);

		filter.blockMedianNxNInternal(data1, width, height, boxSize);
		filter.rollingMedianNxNInternal(data2, width, height, boxSize);

		floatArrayEquals(data1, data2, boxSize, "Internal arrays do not match: [%dx%d] @ %d", width, height, boxSize);
	}

	@Test
	public void floatBlockMedian3x3InternalAndRollingMedianNxNInternalReturnSameResult()
	{
		RandomGenerator rg = TestSettings.getRandomGenerator();
		MedianFilter filter = new MedianFilter();
		for (int width : primes)
			for (int height : primes)
				floatCompareBlockMedian3x3InternalAndRollingMedianNxNInternal(rg, filter, width, height);
	}

	private void floatCompareBlockMedian3x3InternalAndRollingMedianNxNInternal(RandomGenerator rg, MedianFilter filter,
			int width, int height) throws ArrayComparisonFailure
	{
		float[] data1 = createData(rg, width, height);
		float[] data2 = floatClone(data1);

		filter.blockMedian3x3Internal(data1, width, height);
		filter.rollingMedianNxNInternal(data2, width, height, 1);

		floatArrayEquals(data1, data2, 1, "Internal arrays do not match: [%dx%d] @ %d", width, height, 1);
	}

	@Test
	public void floatBlockMedianNxNInternalIsFasterThanRollingMedianNxNInternal()
	{
		TestSettings.assumeSpeedTest();

		MedianFilter filter = new MedianFilter();

		ArrayList<float[]> dataSet = getSpeedData(InternalITER);

		ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.blockMedianNxNInternal(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0], boxSizes[0]);
		filter.rollingMedianNxNInternal(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : speedPrimes)
				for (int height : speedPrimes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.blockMedianNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : speedPrimes)
				for (int height : speedPrimes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.rollingMedianNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"float rollingMedianNxNInternal [%dx%d] @ %d : %d => blockMedianNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestAssert.assert_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			TestSettings.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"float rollingMedianNxNInternal %d : %d => blockMedianNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestSettings.logSpeedTestResult(fastTotal < slowTotal,
				"float rollingMedianNxNInternal %d => blockMedianNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatBlockMedian3x3InternalAndBlockMedianNxNInternalReturnSameResult()
	{
		RandomGenerator rg = TestSettings.getRandomGenerator();
		MedianFilter filter = new MedianFilter();
		for (int width : primes)
			for (int height : primes)
				floatCompareBlockMedian3x3InternalAndBlockMedianNxNInternal(rg, filter, width, height);
	}

	private void floatCompareBlockMedian3x3InternalAndBlockMedianNxNInternal(RandomGenerator rg, MedianFilter filter,
			int width, int height) throws ArrayComparisonFailure
	{
		float[] data1 = createData(rg, width, height);
		float[] data2 = floatClone(data1);

		filter.blockMedian3x3Internal(data1, width, height);
		filter.blockMedianNxNInternal(data2, width, height, 1);

		floatArrayEquals(data1, data2, 1, "Internal arrays do not match: [%dx%d] @ %d", width, height, 1);
	}

	@Test
	public void floatBlockMedian3x3InternalIsFasterThanBlockMedianNxNInternal()
	{
		TestSettings.assumeSpeedTest();

		MedianFilter filter = new MedianFilter();

		ArrayList<float[]> dataSet = getSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.blockMedianNxNInternal(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0], 1);
		filter.blockMedian3x3Internal(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0]);

		for (int width : speedPrimes)
			for (int height : speedPrimes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.blockMedian3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : speedPrimes)
			for (int height : speedPrimes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.blockMedianNxNInternal(data, width, height, 1);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float blockMedianNxNInternal [%dx%d] %d => blockMedian3x3Internal %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestAssert.assert_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestSettings.logSpeedTestResult(fastTotal < slowTotal,
				"float blockMedianNxNInternal %d => blockMedian3x3Internal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatBlockMedian3x3InternalIsFasterThanRollingMedian3x3Internal()
	{
		TestSettings.assumeSpeedTest();

		MedianFilter filter = new MedianFilter();

		ArrayList<float[]> dataSet = getSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingMedian3x3Internal(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0]);
		filter.blockMedian3x3Internal(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0]);

		for (int width : speedPrimes)
			for (int height : speedPrimes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.blockMedian3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : speedPrimes)
			for (int height : speedPrimes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.rollingMedian3x3Internal(data, width, height);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf(
							"float rollingMedian3x3Internal [%dx%d] %d => blockMedian3x3Internal %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestAssert.assert_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestSettings.logSpeedTestResult(fastTotal < slowTotal,
				"float rollingMedian3x3Internal %d => blockMedian3x3Internal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatRollingMedian3x3InternalAndRollingMedianNxNInternalReturnSameResult()
	{
		RandomGenerator rg = TestSettings.getRandomGenerator();
		MedianFilter filter = new MedianFilter();
		for (int width : primes)
			for (int height : primes)
				floatCompareRollingMedian3x3InternalAndRollingMedianNxNInternal(rg, filter, width, height);
	}

	private void floatCompareRollingMedian3x3InternalAndRollingMedianNxNInternal(RandomGenerator rg,
			MedianFilter filter, int width, int height) throws ArrayComparisonFailure
	{
		float[] data1 = createData(rg, width, height);
		float[] data2 = floatClone(data1);

		filter.rollingMedian3x3Internal(data1, width, height);
		filter.rollingMedianNxNInternal(data2, width, height, 1);

		floatArrayEquals(data1, data2, 1, "Internal arrays do not match: [%dx%d] @ %d", width, height, 1);
	}

	@Test
	public void floatRollingMedian3x3InternalIsFasterThanRollingMedianNxNInternal()
	{
		TestSettings.assumeSpeedTest();

		MedianFilter filter = new MedianFilter();

		ArrayList<float[]> dataSet = getSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingMedian3x3Internal(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0]);
		filter.rollingMedianNxNInternal(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0], 1);

		for (int width : speedPrimes)
			for (int height : speedPrimes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.rollingMedian3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : speedPrimes)
			for (int height : speedPrimes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.rollingMedianNxNInternal(data, width, height, 1);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf(
							"float rollingMedianNxNInternal [%dx%d] %d => rollingMedian3x3Internal %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestAssert.assert_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestSettings.logSpeedTestResult(fastTotal < slowTotal,
				"float rollingMedianNxNInternal %d => rollingMedian3x3Internal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatBlockMedianNxNAndRollingMedianNxNReturnSameResult()
	{
		RandomGenerator rg = TestSettings.getRandomGenerator();
		MedianFilter filter = new MedianFilter();
		for (int width : primes)
		{
			for (int height : primes)
				for (int boxSize : boxSizes)
					floatCompareBlockMedianNxNAndRollingMedianNxN(rg, filter, width, height, boxSize);
		}
	}

	private void floatCompareBlockMedianNxNAndRollingMedianNxN(RandomGenerator rg, MedianFilter filter, int width,
			int height, int boxSize) throws ArrayComparisonFailure
	{
		float[] data1 = createData(rg, width, height);
		float[] data2 = floatClone(data1);

		filter.blockMedianNxN(data1, width, height, boxSize);
		filter.rollingMedianNxN(data2, width, height, boxSize);

		floatArrayEquals(data1, data2, boxSize, "Arrays do not match: [%dx%d] @ %d", width, height, boxSize);
	}

	@Test
	public void floatBlockMedianInternalNxNIsFasterThanBlockMedianNxN()
	{
		TestSettings.assumeSpeedTest();

		MedianFilter filter = new MedianFilter();

		ArrayList<float[]> dataSet = getSpeedData(ITER);

		ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.blockMedianNxNInternal(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0], boxSizes[0]);
		filter.blockMedianNxN(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : speedPrimes)
				for (int height : speedPrimes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.blockMedianNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : speedPrimes)
				for (int height : speedPrimes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.blockMedianNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"float blockMedianNxN [%dx%d] @ %d : %d => blockMedianNxNInternal %d = %.2fx\n", width,
								height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestAssert.assert_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			TestSettings.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"float blockMedianNxN %d : %d => blockMedianNxNInternal %d = %.2fx\n", boxSize, boxSlowTotal,
					boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestSettings.logSpeedTestResult(fastTotal < slowTotal,
				"float blockMedianNxN %d => blockMedianNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatBlockMedianNxNIsFasterThanRollingMedianNxN()
	{
		TestSettings.assumeSpeedTest();

		MedianFilter filter = new MedianFilter();

		ArrayList<float[]> dataSet = getSpeedData(ITER);

		ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.blockMedianNxN(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0], boxSizes[0]);
		filter.rollingMedianNxN(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : speedPrimes)
				for (int height : speedPrimes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.blockMedianNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : speedPrimes)
				for (int height : speedPrimes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.rollingMedianNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf("float rollingMedianNxN [%dx%d] @ %d : %d => blockMedianNxN %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestAssert.assert_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		rollingTime, time), rollingTime < time);
				}
			//if (debug)
			TestSettings.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"float rollingMedianNxN %d : %d => blockMedianNxN %d = %.2fx\n", boxSize, boxSlowTotal,
					boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestSettings.logSpeedTestResult(fastTotal < slowTotal,
				"float rollingMedianNxN %d => blockMedianNxN %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatRollingMedianInternalNxNIsFasterThanRollingMedianNxN()
	{
		TestSettings.assumeSpeedTest();

		MedianFilter filter = new MedianFilter();

		ArrayList<float[]> dataSet = getSpeedData(ITER);

		ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingMedianNxNInternal(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0], boxSizes[0]);
		filter.rollingMedianNxN(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : speedPrimes)
				for (int height : speedPrimes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.rollingMedianNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : speedPrimes)
				for (int height : speedPrimes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.rollingMedianNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"float rollingMedianNxN [%dx%d] @ %d : %d => rollingMedianNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestAssert.assert_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			TestSettings.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"float rollingMedianNxN %d : %d => rollingMedianNxNInternal %d = %.2fx\n", boxSize, boxSlowTotal,
					boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestSettings.logSpeedTestResult(fastTotal < slowTotal,
				"float rollingMedianNxN %d => rollingMedianNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatBlockMedian3x3AndBlockMedianNxNReturnSameResult()
	{
		RandomGenerator rg = TestSettings.getRandomGenerator();
		MedianFilter filter = new MedianFilter();
		for (int width : primes)
			for (int height : primes)
				floatCompareBlockMedian3x3AndBlockMedianNxN(rg, filter, width, height);
	}

	private void floatCompareBlockMedian3x3AndBlockMedianNxN(RandomGenerator rg, MedianFilter filter, int width,
			int height) throws ArrayComparisonFailure
	{
		float[] data1 = createData(rg, width, height);
		float[] data2 = floatClone(data1);

		filter.blockMedian3x3(data1, width, height);
		filter.blockMedianNxN(data2, width, height, 1);

		floatArrayEquals(data1, data2, 1, "Arrays do not match: [%dx%d] @ %d", width, height, 1);
	}

	@Test
	public void floatBlockMedian3x3IsFasterThanBlockMedianNxN()
	{
		TestSettings.assumeSpeedTest();

		MedianFilter filter = new MedianFilter();

		ArrayList<float[]> dataSet = getSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.blockMedianNxN(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0], 1);
		filter.blockMedian3x3(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0]);

		for (int width : speedPrimes)
			for (int height : speedPrimes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.blockMedian3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : speedPrimes)
			for (int height : speedPrimes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.blockMedianNxN(data, width, height, 1);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float blockMedianNxN [%dx%d] %d => blockMedian3x3 %d = %.2fx\n", width, height,
							time, fastTime, speedUpFactor(time, fastTime));
				//if (TestAssert.assert_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestSettings.logSpeedTestResult(fastTotal < slowTotal, "float blockMedianNxN %d => blockMedian3x3 %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatRollingMedian3x3AndRollingMedianNxNReturnSameResult()
	{
		RandomGenerator rg = TestSettings.getRandomGenerator();
		MedianFilter filter = new MedianFilter();

		for (int width : primes)
			for (int height : primes)
				floatCompareRollingMedian3x3AndRollingMedianNxN(rg, filter, width, height);
	}

	private void floatCompareRollingMedian3x3AndRollingMedianNxN(RandomGenerator rg, MedianFilter filter, int width,
			int height) throws ArrayComparisonFailure
	{
		float[] data1 = createData(rg, width, height);
		float[] data2 = floatClone(data1);

		filter.rollingMedian3x3(data1, width, height);
		filter.rollingMedianNxN(data2, width, height, 1);

		floatArrayEquals(data1, data2, 1, "Arrays do not match: [%dx%d] @ %d", width, height, 1);
	}

	@Test
	public void floatRollingMedian3x3IsFasterThanRollingMedianNxN()
	{
		TestSettings.assumeSpeedTest();

		MedianFilter filter = new MedianFilter();

		ArrayList<float[]> dataSet = getSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingMedianNxN(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0], 1);
		filter.rollingMedian3x3(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0]);

		for (int width : speedPrimes)
			for (int height : speedPrimes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.rollingMedian3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : speedPrimes)
			for (int height : speedPrimes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.rollingMedianNxN(data, width, height, 1);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float rollingMedianNxN [%dx%d] %d => rollingMedian3x3 %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestAssert.assert_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		rollingBlockTime, time), rollingBlockTime < time);
			}
		TestSettings.logSpeedTestResult(fastTotal < slowTotal,
				"float rollingMedianNxN %d => rollingMedian3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatRollingMedian3x3IsFasterThanBlockMedian3x3()
	{
		TestSettings.assumeSpeedTest();

		MedianFilter filter = new MedianFilter();

		ArrayList<float[]> dataSet = getSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingMedian3x3(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0]);
		filter.blockMedian3x3(floatClone(dataSet.get(0)), speedPrimes[0], speedPrimes[0]);

		for (int width : speedPrimes)
			for (int height : speedPrimes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.rollingMedian3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (int width : speedPrimes)
			for (int height : speedPrimes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.blockMedian3x3(data, width, height);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float blockMedian3x3 [%dx%d] %d => rollingMedian3x3 %d = %.2fx\n", width, height,
							time, fastTime, speedUpFactor(time, fastTime));
				//if (TestAssert.assert_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestSettings.logSpeedTestResult(fastTotal < slowTotal,
				"float blockMedian3x3 %d => rollingMedian3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}
}
