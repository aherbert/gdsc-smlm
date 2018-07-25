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
import org.junit.jupiter.api.Test;
import org.junit.internal.ArrayComparisonFailure;

import uk.ac.sussex.gdsc.test.TestComplexity;
import uk.ac.sussex.gdsc.test.TestLog;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssertions;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssumptions;

@SuppressWarnings({ "deprecation", "javadoc" })
public class SumFilterTest extends AbstractFilterTest
{
	private static int InternalITER3 = 500;
	private static int InternalITER = 50;
	private static int ITER3 = 200;
	private static int ITER = 20;

	/**
	 * Check the float arrays are equal, else fail with a formatted message
	 *
	 * @param data1
	 *            the data 1
	 * @param data2
	 *            the data 2
	 * @param boxSize
	 *            the box size used for filtering
	 * @param format
	 *            the format
	 * @param args
	 *            the args
	 */
	private static void floatArrayEquals(float[] data1, float[] data2, int boxSize, String format, Object... args)
	{
		ExtraAssertions.assertArrayEqualsRelative(data1, data2, 1e-4, format, args);
	}

	/**
	 * Check the int arrays are equal, else fail with a formatted message
	 *
	 * @param data1
	 *            the data 1
	 * @param data2
	 *            the data 2
	 * @param boxSize
	 *            the box size used for filtering
	 * @param format
	 *            the format
	 * @param args
	 *            the args
	 */
	private static void intArrayEquals(int[] data1, int[] data2, int boxSize, String format, Object... args)
	{
		ExtraAssertions.assertArrayEquals(data1, data2, format, args);
	}

	private static float[] floatCreateData(RandomGenerator rg, int width, int height)
	{
		return createData(rg, width, height);
	}

	private static int[] intCreateData(RandomGenerator rg, int width, int height)
	{
		return createIntData(rg, width, height);
	}

	private static ArrayList<float[]> floatCreateSpeedData(int iter)
	{
		return getSpeedData(iter);
	}

	private static ArrayList<int[]> intCreateSpeedData(int iter)
	{
		return getIntSpeedData(iter);
	}

	// COPY CODE FROM HERE...
	@Test
	public void floatBlockSumNxNInternalAndRollingBlockSumNxNInternalReturnSameResult()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SumFilter filter = new SumFilter();

		for (final int width : primes)
			for (final int height : primes)
				for (final int boxSize : boxSizes)
					floatCompareBlockSumNxNInternalAndRollingBlockSumNxNInternal(rg, filter, width, height, boxSize);
	}

	private static void floatCompareBlockSumNxNInternalAndRollingBlockSumNxNInternal(RandomGenerator rg,
			SumFilter filter, int width, int height, int boxSize) throws ArrayComparisonFailure
	{
		final float[] data1 = floatCreateData(rg, width, height);
		final float[] data2 = floatClone(data1);

		filter.blockSumNxNInternal(data1, width, height, boxSize);
		filter.rollingBlockSumNxNInternal(data2, width, height, boxSize);

		floatArrayEquals(data1, data2, boxSize, "Internal arrays do not match: [%dx%d] @ %d", width, height, boxSize);
	}

	@Test
	public void floatBlockSumNxNInternalAndStripedBlockSumNxNInternalReturnSameResult()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SumFilter filter = new SumFilter();

		for (final int width : primes)
			for (final int height : primes)
				for (final int boxSize : boxSizes)
					floatCompareBlockSumNxNInternalAndStripedBlockSumNxNInternal(rg, filter, width, height, boxSize);
	}

	private static void floatCompareBlockSumNxNInternalAndStripedBlockSumNxNInternal(RandomGenerator rg,
			SumFilter filter, int width, int height, int boxSize) throws ArrayComparisonFailure
	{
		final float[] data1 = floatCreateData(rg, width, height);
		final float[] data2 = floatClone(data1);

		filter.blockSumNxNInternal(data1, width, height, boxSize);
		filter.stripedBlockSumNxNInternal(data2, width, height, boxSize);

		floatArrayEquals(data1, data2, boxSize, "Internal arrays do not match: [%dx%d] @ %d", width, height, boxSize);
	}

	@Test
	public void floatBlockSum3x3InternalAndRollingBlockSumNxNInternalReturnSameResult()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SumFilter filter = new SumFilter();

		for (final int width : primes)
			for (final int height : primes)
				floatCompareBlockSum3x3InternalAndRollingBlockSumNxNInternal(rg, filter, width, height);
	}

	private static void floatCompareBlockSum3x3InternalAndRollingBlockSumNxNInternal(RandomGenerator rg,
			SumFilter filter, int width, int height) throws ArrayComparisonFailure
	{
		final float[] data1 = floatCreateData(rg, width, height);
		final float[] data2 = floatClone(data1);

		filter.blockSum3x3Internal(data1, width, height);
		filter.rollingBlockSumNxNInternal(data2, width, height, 1);

		floatArrayEquals(data1, data2, 1, "Internal arrays do not match: [%dx%d] @ %d", width, height, 1);
	}

	@Test
	public void floatRollingBlockSumNxNInternalAndRollingBlockSumNxNInternalTransposedReturnSameResult()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SumFilter filter = new SumFilter();

		for (final int width : primes)
			for (final int height : primes)
				for (final int boxSize : boxSizes)
					floatCompareRollingBlockSumNxNInternalAndRollingBlockSumNxNInternalTransposed(rg, filter, width,
							height, boxSize);
	}

	private static void floatCompareRollingBlockSumNxNInternalAndRollingBlockSumNxNInternalTransposed(
			RandomGenerator rg, SumFilter filter, int width, int height, int boxSize) throws ArrayComparisonFailure
	{
		final float[] data1 = floatCreateData(rg, width, height);
		final float[] data2 = floatClone(data1);

		filter.rollingBlockSumNxNInternal(data1, width, height, boxSize);
		filter.rollingBlockSumNxNInternalTransposed(data2, width, height, boxSize);

		floatArrayEquals(data1, data2, boxSize, "Internal arrays do not match: [%dx%d] @ %d", width, height, boxSize);
	}

	@Test
	public void floatRollingBlockSumNxNInternalIsFasterThanBlockSumNxNInternal()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.blockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (final int boxSize : boxSizes)
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (final int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter.blockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					final long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"float blockSumNxNInternal [%dx%d] @ %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			TestLog.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"float blockSumNxNInternal %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float blockSumNxNInternal %d => rollingBlockSumNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatStripedBlockSumNxNInternalIsFasterThanBlockSumNxNInternal()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.blockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.stripedBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (final int boxSize : boxSizes)
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter.stripedBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (final int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter.blockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					final long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"float blockSumNxNInternal [%dx%d] @ %d : %d => stripedBlockSumNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			TestLog.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"float blockSumNxNInternal %d : %d => stripedBlockSumNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float blockSumNxNInternal %d => stripedBlockSumNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatRollingBlockSumNxNInternalIsFasterThanStripedBlockSumNxNInternal()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.stripedBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (final int boxSize : boxSizes)
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (final int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter.stripedBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					final long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"float stripedBlockSumNxNInternal [%dx%d] @ %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			TestLog.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"float stripedBlockSumNxNInternal %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float stripedBlockSumNxNInternal %d => rollingBlockSumNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatBlockSum3x3InternalAndBlockSumNxNInternalReturnSameResult()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SumFilter filter = new SumFilter();

		for (final int width : primes)
			for (final int height : primes)
				floatCompareBlockSum3x3InternalAndBlockSumNxNInternal(rg, filter, width, height);
	}

	private static void floatCompareBlockSum3x3InternalAndBlockSumNxNInternal(RandomGenerator rg, SumFilter filter,
			int width, int height) throws ArrayComparisonFailure
	{
		final float[] data1 = floatCreateData(rg, width, height);
		final float[] data2 = floatClone(data1);

		filter.blockSum3x3Internal(data1, width, height);
		filter.blockSumNxNInternal(data2, width, height, 1);

		floatArrayEquals(data1, data2, 1, "Internal arrays do not match: [%dx%d] @ %d", width, height, 1);
	}

	@Test
	public void floatBlockSum3x3InternalIsFasterThanBlockSumNxNInternal()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.blockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.blockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.blockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.blockSumNxNInternal(data, width, height, 1);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float blockSumNxNInternal [%dx%d] %d => blockSum3x3Internal %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float blockSumNxNInternal %d => blockSum3x3Internal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatRollingBlockSum3x3InternalIsFasterThanBlockSum3x3Internal()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingBlockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.blockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.rollingBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.blockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float blockSum3x3Internal [%dx%d] %d => rollingBlockSum3x3Internal %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float blockSum3x3Internal %d => rollingBlockSum3x3Internal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatStripedBlockSum3x3InternalIsFasterThanBlockSum3x3Internal()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.stripedBlockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.blockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.stripedBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.blockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float blockSum3x3Internal [%dx%d] %d => stripedBlockSum3x3Internal %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float blockSum3x3Internal %d => stripedBlockSum3x3Internal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatRollingBlockSum3x3InternalIsFasterThanStripedBlockSum3x3Internal()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingBlockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.stripedBlockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.rollingBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.stripedBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf(
							"float stripedBlockSum3x3Internal [%dx%d] %d => rollingBlockSum3x3Internal %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float stripedBlockSum3x3Internal %d => rollingBlockSum3x3Internal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatRollingBlockSum3x3InternalAndRollingBlockSumNxNInternalReturnSameResult()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SumFilter filter = new SumFilter();

		for (final int width : primes)
			for (final int height : primes)
				floatCompareRollingBlockSum3x3InternalAndRollingBlockSumNxNInternal(rg, filter, width, height);
	}

	private static void floatCompareRollingBlockSum3x3InternalAndRollingBlockSumNxNInternal(RandomGenerator rg,
			SumFilter filter, int width, int height) throws ArrayComparisonFailure
	{
		final float[] data1 = floatCreateData(rg, width, height);
		final float[] data2 = floatClone(data1);

		filter.rollingBlockSum3x3Internal(data1, width, height);
		filter.rollingBlockSumNxNInternal(data2, width, height, 1);

		floatArrayEquals(data1, data2, 1, "Internal arrays do not match: [%dx%d] @ %d", width, height, 1);
	}

	@Test
	public void floatRollingBlockSum3x3InternalIsFasterThanRollingBlockSumNxNInternal()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingBlockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.rollingBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], 1);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.rollingBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.rollingBlockSumNxNInternal(data, width, height, 1);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf(
							"float rollingBlockSumNxNInternal [%dx%d] %d => rollingBlockSum3x3Internal %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float rollingBlockSumNxNInternal %d => rollingBlockSum3x3Internal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatStripedBlockSum3x3InternalIsFasterThanStripedBlockSumNxNInternal()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.stripedBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.stripedBlockSum3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.stripedBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.stripedBlockSumNxNInternal(data, width, height, 1);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf(
							"float stripedBlockSumNxNInternal [%dx%d] %d => stripedBlockSum3x3Internal %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float stripedBlockSumNxNInternal %d => stripedBlockSum3x3Internal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatRollingBlockSumNxNInternalIsFasterThanRollingBlockSumNxNInternalTransposed()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingBlockSumNxNInternalTransposed(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (final int boxSize : boxSizes)
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (final int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter.rollingBlockSumNxNInternalTransposed(data, width, height, boxSize);
					time = System.nanoTime() - time;

					final long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"float rollingBlockSumNxNInternalTransposed [%dx%d] @ %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			if (debug)
				System.out.printf(
						"float rollingBlockSumNxNInternalTransposed %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
						boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float rollingBlockSumNxNInternalTransposed %d => rollingBlockSumNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatBlockSumNxNAndStripedBlockSumNxNReturnSameResult()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SumFilter filter = new SumFilter();

		for (final int width : primes)
			for (final int height : primes)
				for (final int boxSize : boxSizes)
					floatCompareBlockSumNxNAndStripedBlockSumNxN(rg, filter, width, height, boxSize);
	}

	private static void floatCompareBlockSumNxNAndStripedBlockSumNxN(RandomGenerator rg, SumFilter filter, int width,
			int height, int boxSize) throws ArrayComparisonFailure
	{
		final float[] data1 = floatCreateData(rg, width, height);
		final float[] data2 = floatClone(data1);

		filter.blockSumNxN(data1, width, height, boxSize);
		filter.rollingBlockSumNxN(data2, width, height, boxSize);

		floatArrayEquals(data1, data2, boxSize, "Arrays do not match: [%dx%d] @ %d", width, height, boxSize);
	}

	@Test
	public void floatBlockSumNxNAndRollingBlockSumNxNReturnSameResult()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SumFilter filter = new SumFilter();

		for (final int width : primes)
			for (final int height : primes)
				for (final int boxSize : boxSizes)
					floatCompareBlockSumNxNAndRollingBlockSumNxN(rg, filter, width, height, boxSize);
	}

	private static void floatCompareBlockSumNxNAndRollingBlockSumNxN(RandomGenerator rg, SumFilter filter, int width,
			int height, int boxSize) throws ArrayComparisonFailure
	{
		final float[] data1 = floatCreateData(rg, width, height);
		final float[] data2 = floatClone(data1);

		filter.blockSumNxN(data1, width, height, boxSize);
		filter.rollingBlockSumNxN(data2, width, height, boxSize);

		floatArrayEquals(data1, data2, boxSize, "Arrays do not match: [%dx%d] @ %d", width, height, boxSize);
	}

	@Test
	public void floatBlockSumInternalNxNIsFasterThanBlockSumNxN()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<float[]> dataSet = floatCreateSpeedData(ITER);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.blockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.blockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (final int boxSize : boxSizes)
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter.blockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (final int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter.blockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					final long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf("float blockSumNxN [%dx%d] @ %d : %d => blockSumNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			TestLog.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"float blockSumNxN %d : %d => blockSumNxNInternal %d = %.2fx\n", boxSize, boxSlowTotal,
					boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float blockSumNxN %d => blockSumNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatStripedBlockSumNxNIsFasterThanBlockSumNxN()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<float[]> dataSet = floatCreateSpeedData(ITER);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.blockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.stripedBlockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (final int boxSize : boxSizes)
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter.stripedBlockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (final int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter.blockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					final long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf("float blockSumNxN [%dx%d] @ %d : %d => stripedBlockSumNxN %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			TestLog.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"float blockSumNxN %d : %d => stripedBlockSumNxN %d = %.2fx\n", boxSize, boxSlowTotal, boxFastTotal,
					speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float blockSumNxN %d => stripedBlockSumNxN %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatStripedBlockSumInternalNxNIsFasterThanStripedBlockSumNxN()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<float[]> dataSet = floatCreateSpeedData(ITER);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.stripedBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.stripedBlockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (final int boxSize : boxSizes)
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter.stripedBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (final int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter.stripedBlockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					final long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"float stripedBlockSumNxN [%dx%d] @ %d : %d => stripedBlockSumNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			TestLog.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"float stripedBlockSumNxN %d : %d => stripedBlockSumNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float stripedBlockSumNxN %d => stripedBlockSumNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatRollingBlockSumNxNIsFasterThanBlockSumNxN()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<float[]> dataSet = floatCreateSpeedData(ITER);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.blockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingBlockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (final int boxSize : boxSizes)
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter.rollingBlockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (final int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter.blockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					final long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf("float blockSumNxN [%dx%d] @ %d : %d => rollingBlockSumNxN %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			TestLog.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"float blockSumNxN %d : %d => rollingBlockSumNxN %d = %.2fx\n", boxSize, boxSlowTotal, boxFastTotal,
					speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float blockSumNxN %d => rollingBlockSumNxN %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatRollingBlockSumInternalNxNIsFasterThanRollingBlockSumNxN()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<float[]> dataSet = floatCreateSpeedData(ITER);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingBlockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (final int boxSize : boxSizes)
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (final int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (final float[] data : dataSet2)
						filter.rollingBlockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					final long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"float rollingBlockSumNxN [%dx%d] @ %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			TestLog.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"float rollingBlockSumNxN %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float rollingBlockSumNxN %d => rollingBlockSumNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatBlockSum3x3AndBlockSumNxNReturnSameResult()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SumFilter filter = new SumFilter();

		for (final int width : primes)
			for (final int height : primes)
				floatCompareBlockSum3x3AndBlockSumNxN(rg, filter, width, height);
	}

	private static void floatCompareBlockSum3x3AndBlockSumNxN(RandomGenerator rg, SumFilter filter, int width,
			int height) throws ArrayComparisonFailure
	{
		final float[] data1 = floatCreateData(rg, width, height);
		final float[] data2 = floatClone(data1);

		filter.blockSum3x3(data1, width, height);
		filter.blockSumNxN(data2, width, height, 1);

		floatArrayEquals(data1, data2, 1, "Arrays do not match: [%dx%d] @ %d", width, height, 1);
	}

	@Test
	public void floatBlockSum3x3IsFasterThanBlockSumNxN()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<float[]> dataSet = floatCreateSpeedData(ITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.blockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.blockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.blockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.blockSumNxN(data, width, height, 1);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float blockSumNxN [%dx%d] %d => blockSum3x3 %d = %.2fx\n", width, height, time,
							fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal, "float blockSumNxN %d => blockSum3x3 %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatStripedBlockSum3x3AndStripedBlockSumNxNReturnSameResult()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SumFilter filter = new SumFilter();

		for (final int width : primes)
			for (final int height : primes)
				floatCompareStripedBlockSum3x3AndStripedBlockSumNxN(rg, filter, width, height);
	}

	private static void floatCompareStripedBlockSum3x3AndStripedBlockSumNxN(RandomGenerator rg, SumFilter filter,
			int width, int height) throws ArrayComparisonFailure
	{
		final float[] data1 = floatCreateData(rg, width, height);
		final float[] data2 = floatClone(data1);

		filter.stripedBlockSum3x3(data1, width, height);
		filter.stripedBlockSumNxN(data2, width, height, 1);

		floatArrayEquals(data1, data2, 1, "Arrays do not match: [%dx%d] @ %d", width, height, 1);
	}

	@Test
	public void floatStripedBlockSum3x3IsFasterThanStripedBlockSumNxN()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<float[]> dataSet = floatCreateSpeedData(ITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.stripedBlockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.stripedBlockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.stripedBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.stripedBlockSumNxN(data, width, height, 1);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float stripedBlockSumNxN [%dx%d] %d => stripedBlockSum3x3 %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		stripedBlockTime, time), stripedBlockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float stripedBlockSumNxN %d => stripedBlockSum3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatRollingBlockSum3x3AndRollingBlockSumNxNReturnSameResult()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SumFilter filter = new SumFilter();

		for (final int width : primes)
			for (final int height : primes)
				floatCompareRollingBlockSum3x3AndRollingBlockSumNxN(rg, filter, width, height);
	}

	private static void floatCompareRollingBlockSum3x3AndRollingBlockSumNxN(RandomGenerator rg, SumFilter filter,
			int width, int height) throws ArrayComparisonFailure
	{
		final float[] data1 = floatCreateData(rg, width, height);
		final float[] data2 = floatClone(data1);

		filter.rollingBlockSum3x3(data1, width, height);
		filter.rollingBlockSumNxN(data2, width, height, 1);

		floatArrayEquals(data1, data2, 1, "Arrays do not match: [%dx%d] @ %d", width, height, 1);
	}

	@Test
	public void floatRollingBlockSum3x3IsFasterThanRollingBlockSumNxN()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<float[]> dataSet = floatCreateSpeedData(ITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingBlockSumNxN(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.rollingBlockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.rollingBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.rollingBlockSumNxN(data, width, height, 1);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float rollingBlockSumNxN [%dx%d] %d => rollingBlockSum3x3 %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		rollingBlockTime, time), rollingBlockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float rollingBlockSumNxN %d => rollingBlockSum3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatRollingBlockSum3x3IsFasterThanBlockSum3x3()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<float[]> dataSet = floatCreateSpeedData(ITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingBlockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.blockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.rollingBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.blockSum3x3(data, width, height);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float blockSum3x3 [%dx%d] %d => rollingBlockSum3x3 %d = %.2fx\n", width, height,
							time, fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float blockSum3x3 %d => rollingBlockSum3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatStripedBlockSum3x3IsFasterThanBlockSum3x3()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<float[]> dataSet = floatCreateSpeedData(ITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.stripedBlockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.blockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.stripedBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.blockSum3x3(data, width, height);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float blockSum3x3 [%dx%d] %d => stripedBlockSum3x3 %d = %.2fx\n", width, height,
							time, fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float blockSum3x3 %d => stripedBlockSum3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatRollingBlockSum3x3IsFasterThanStripedBlockSum3x3()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<float[]> dataSet = floatCreateSpeedData(ITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingBlockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.stripedBlockSum3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.rollingBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<float[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (final float[] data : dataSet2)
					filter.stripedBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float stripedBlockSum3x3 [%dx%d] %d => rollingBlockSum3x3 %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"float stripedBlockSum3x3 %d => rollingBlockSum3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void intBlockSumNxNInternalAndRollingBlockSumNxNInternalReturnSameResult()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SumFilter filter = new SumFilter();

		for (final int width : primes)
			for (final int height : primes)
				for (final int boxSize : boxSizes)
					intCompareBlockSumNxNInternalAndRollingBlockSumNxNInternal(rg, filter, width, height, boxSize);
	}

	private static void intCompareBlockSumNxNInternalAndRollingBlockSumNxNInternal(RandomGenerator rg, SumFilter filter,
			int width, int height, int boxSize) throws ArrayComparisonFailure
	{
		final int[] data1 = intCreateData(rg, width, height);
		final int[] data2 = intClone(data1);

		filter.blockSumNxNInternal(data1, width, height, boxSize);
		filter.rollingBlockSumNxNInternal(data2, width, height, boxSize);

		intArrayEquals(data1, data2, boxSize, "Internal arrays do not match: [%dx%d] @ %d", width, height, boxSize);
	}

	@Test
	public void intBlockSumNxNInternalAndStripedBlockSumNxNInternalReturnSameResult()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SumFilter filter = new SumFilter();

		for (final int width : primes)
			for (final int height : primes)
				for (final int boxSize : boxSizes)
					intCompareBlockSumNxNInternalAndStripedBlockSumNxNInternal(rg, filter, width, height, boxSize);
	}

	private static void intCompareBlockSumNxNInternalAndStripedBlockSumNxNInternal(RandomGenerator rg, SumFilter filter,
			int width, int height, int boxSize) throws ArrayComparisonFailure
	{
		final int[] data1 = intCreateData(rg, width, height);
		final int[] data2 = intClone(data1);

		filter.blockSumNxNInternal(data1, width, height, boxSize);
		filter.stripedBlockSumNxNInternal(data2, width, height, boxSize);

		intArrayEquals(data1, data2, boxSize, "Internal arrays do not match: [%dx%d] @ %d", width, height, boxSize);
	}

	@Test
	public void intBlockSum3x3InternalAndRollingBlockSumNxNInternalReturnSameResult()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SumFilter filter = new SumFilter();

		for (final int width : primes)
			for (final int height : primes)
				intCompareBlockSum3x3InternalAndRollingBlockSumNxNInternal(rg, filter, width, height);
	}

	private static void intCompareBlockSum3x3InternalAndRollingBlockSumNxNInternal(RandomGenerator rg, SumFilter filter,
			int width, int height) throws ArrayComparisonFailure
	{
		final int[] data1 = intCreateData(rg, width, height);
		final int[] data2 = intClone(data1);

		filter.blockSum3x3Internal(data1, width, height);
		filter.rollingBlockSumNxNInternal(data2, width, height, 1);

		intArrayEquals(data1, data2, 1, "Internal arrays do not match: [%dx%d] @ %d", width, height, 1);
	}

	@Test
	public void intRollingBlockSumNxNInternalAndRollingBlockSumNxNInternalTransposedReturnSameResult()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SumFilter filter = new SumFilter();

		for (final int width : primes)
			for (final int height : primes)
				for (final int boxSize : boxSizes)
					intCompareRollingBlockSumNxNInternalAndRollingBlockSumNxNInternalTransposed(rg, filter, width,
							height, boxSize);
	}

	private static void intCompareRollingBlockSumNxNInternalAndRollingBlockSumNxNInternalTransposed(RandomGenerator rg,
			SumFilter filter, int width, int height, int boxSize) throws ArrayComparisonFailure
	{
		final int[] data1 = intCreateData(rg, width, height);
		final int[] data2 = intClone(data1);

		filter.rollingBlockSumNxNInternal(data1, width, height, boxSize);
		filter.rollingBlockSumNxNInternalTransposed(data2, width, height, boxSize);

		intArrayEquals(data1, data2, boxSize, "Internal arrays do not match: [%dx%d] @ %d", width, height, boxSize);
	}

	@Test
	public void intRollingBlockSumNxNInternalIsFasterThanBlockSumNxNInternal()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<int[]> dataSet = intCreateSpeedData(InternalITER);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.blockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (final int boxSize : boxSizes)
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (final int[] data : dataSet2)
						filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (final int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (final int[] data : dataSet2)
						filter.blockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					final long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"int blockSumNxNInternal [%dx%d] @ %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			TestLog.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"int blockSumNxNInternal %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n", boxSize, boxSlowTotal,
					boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"int blockSumNxNInternal %d => rollingBlockSumNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void intStripedBlockSumNxNInternalIsFasterThanBlockSumNxNInternal()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<int[]> dataSet = intCreateSpeedData(InternalITER);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.blockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.stripedBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (final int boxSize : boxSizes)
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (final int[] data : dataSet2)
						filter.stripedBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (final int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (final int[] data : dataSet2)
						filter.blockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					final long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"int blockSumNxNInternal [%dx%d] @ %d : %d => stripedBlockSumNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			TestLog.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"int blockSumNxNInternal %d : %d => stripedBlockSumNxNInternal %d = %.2fx\n", boxSize, boxSlowTotal,
					boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"int blockSumNxNInternal %d => stripedBlockSumNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void intRollingBlockSumNxNInternalIsFasterThanStripedBlockSumNxNInternal()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<int[]> dataSet = intCreateSpeedData(InternalITER);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.stripedBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (final int boxSize : boxSizes)
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (final int[] data : dataSet2)
						filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (final int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (final int[] data : dataSet2)
						filter.stripedBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					final long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"int stripedBlockSumNxNInternal [%dx%d] @ %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			TestLog.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"int stripedBlockSumNxNInternal %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"int stripedBlockSumNxNInternal %d => rollingBlockSumNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void intBlockSum3x3InternalAndBlockSumNxNInternalReturnSameResult()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SumFilter filter = new SumFilter();

		for (final int width : primes)
			for (final int height : primes)
				intCompareBlockSum3x3InternalAndBlockSumNxNInternal(rg, filter, width, height);
	}

	private static void intCompareBlockSum3x3InternalAndBlockSumNxNInternal(RandomGenerator rg, SumFilter filter,
			int width, int height) throws ArrayComparisonFailure
	{
		final int[] data1 = intCreateData(rg, width, height);
		final int[] data2 = intClone(data1);

		filter.blockSum3x3Internal(data1, width, height);
		filter.blockSumNxNInternal(data2, width, height, 1);

		intArrayEquals(data1, data2, 1, "Internal arrays do not match: [%dx%d] @ %d", width, height, 1);
	}

	@Test
	public void intBlockSum3x3InternalIsFasterThanBlockSumNxNInternal()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<int[]> dataSet = intCreateSpeedData(InternalITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.blockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.blockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.blockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.blockSumNxNInternal(data, width, height, 1);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("int blockSumNxNInternal [%dx%d] %d => blockSum3x3Internal %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"int blockSumNxNInternal %d => blockSum3x3Internal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void intRollingBlockSum3x3InternalIsFasterThanBlockSum3x3Internal()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<int[]> dataSet = intCreateSpeedData(InternalITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingBlockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);
		filter.blockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.rollingBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.blockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("int blockSum3x3Internal [%dx%d] %d => rollingBlockSum3x3Internal %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"int blockSum3x3Internal %d => rollingBlockSum3x3Internal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void intStripedBlockSum3x3InternalIsFasterThanBlockSum3x3Internal()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<int[]> dataSet = intCreateSpeedData(InternalITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.stripedBlockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);
		filter.blockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.stripedBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.blockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("int blockSum3x3Internal [%dx%d] %d => stripedBlockSum3x3Internal %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"int blockSum3x3Internal %d => stripedBlockSum3x3Internal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void intRollingBlockSum3x3InternalIsFasterThanStripedBlockSum3x3Internal()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<int[]> dataSet = intCreateSpeedData(InternalITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingBlockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);
		filter.stripedBlockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.rollingBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.stripedBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf(
							"int stripedBlockSum3x3Internal [%dx%d] %d => rollingBlockSum3x3Internal %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"int stripedBlockSum3x3Internal %d => rollingBlockSum3x3Internal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void intRollingBlockSum3x3InternalAndRollingBlockSumNxNInternalReturnSameResult()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SumFilter filter = new SumFilter();

		for (final int width : primes)
			for (final int height : primes)
				intCompareRollingBlockSum3x3InternalAndRollingBlockSumNxNInternal(rg, filter, width, height);
	}

	private static void intCompareRollingBlockSum3x3InternalAndRollingBlockSumNxNInternal(RandomGenerator rg,
			SumFilter filter, int width, int height) throws ArrayComparisonFailure
	{
		final int[] data1 = intCreateData(rg, width, height);
		final int[] data2 = intClone(data1);

		filter.rollingBlockSum3x3Internal(data1, width, height);
		filter.rollingBlockSumNxNInternal(data2, width, height, 1);

		intArrayEquals(data1, data2, 1, "Internal arrays do not match: [%dx%d] @ %d", width, height, 1);
	}

	@Test
	public void intRollingBlockSum3x3InternalIsFasterThanRollingBlockSumNxNInternal()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<int[]> dataSet = intCreateSpeedData(InternalITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingBlockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);
		filter.rollingBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], 1);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.rollingBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.rollingBlockSumNxNInternal(data, width, height, 1);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf(
							"int rollingBlockSumNxNInternal [%dx%d] %d => rollingBlockSum3x3Internal %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"int rollingBlockSumNxNInternal %d => rollingBlockSum3x3Internal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void intStripedBlockSum3x3InternalIsFasterThanStripedBlockSumNxNInternal()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<int[]> dataSet = intCreateSpeedData(InternalITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.stripedBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.stripedBlockSum3x3Internal(intClone(dataSet.get(0)), primes[0], primes[0]);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.stripedBlockSum3x3Internal(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.stripedBlockSumNxNInternal(data, width, height, 1);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf(
							"int stripedBlockSumNxNInternal [%dx%d] %d => stripedBlockSum3x3Internal %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"int stripedBlockSumNxNInternal %d => stripedBlockSum3x3Internal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void intRollingBlockSumNxNInternalIsFasterThanRollingBlockSumNxNInternalTransposed()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<int[]> dataSet = intCreateSpeedData(InternalITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingBlockSumNxNInternalTransposed(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (final int boxSize : boxSizes)
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (final int[] data : dataSet2)
						filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (final int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (final int[] data : dataSet2)
						filter.rollingBlockSumNxNInternalTransposed(data, width, height, boxSize);
					time = System.nanoTime() - time;

					final long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"int rollingBlockSumNxNInternalTransposed [%dx%d] @ %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			if (debug)
				System.out.printf(
						"int rollingBlockSumNxNInternalTransposed %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
						boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"int rollingBlockSumNxNInternalTransposed %d => rollingBlockSumNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void intBlockSumNxNAndStripedBlockSumNxNReturnSameResult()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SumFilter filter = new SumFilter();

		for (final int width : primes)
			for (final int height : primes)
				for (final int boxSize : boxSizes)
					intCompareBlockSumNxNAndStripedBlockSumNxN(rg, filter, width, height, boxSize);
	}

	private static void intCompareBlockSumNxNAndStripedBlockSumNxN(RandomGenerator rg, SumFilter filter, int width,
			int height, int boxSize) throws ArrayComparisonFailure
	{
		final int[] data1 = intCreateData(rg, width, height);
		final int[] data2 = intClone(data1);

		filter.blockSumNxN(data1, width, height, boxSize);
		filter.rollingBlockSumNxN(data2, width, height, boxSize);

		intArrayEquals(data1, data2, boxSize, "Arrays do not match: [%dx%d] @ %d", width, height, boxSize);
	}

	@Test
	public void intBlockSumNxNAndRollingBlockSumNxNReturnSameResult()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SumFilter filter = new SumFilter();

		for (final int width : primes)
			for (final int height : primes)
				for (final int boxSize : boxSizes)
					intCompareBlockSumNxNAndRollingBlockSumNxN(rg, filter, width, height, boxSize);
	}

	private static void intCompareBlockSumNxNAndRollingBlockSumNxN(RandomGenerator rg, SumFilter filter, int width,
			int height, int boxSize) throws ArrayComparisonFailure
	{
		final int[] data1 = intCreateData(rg, width, height);
		final int[] data2 = intClone(data1);

		filter.blockSumNxN(data1, width, height, boxSize);
		filter.rollingBlockSumNxN(data2, width, height, boxSize);

		intArrayEquals(data1, data2, boxSize, "Arrays do not match: [%dx%d] @ %d", width, height, boxSize);
	}

	@Test
	public void intBlockSumInternalNxNIsFasterThanBlockSumNxN()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<int[]> dataSet = intCreateSpeedData(ITER);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.blockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.blockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (final int boxSize : boxSizes)
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (final int[] data : dataSet2)
						filter.blockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (final int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (final int[] data : dataSet2)
						filter.blockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					final long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf("int blockSumNxN [%dx%d] @ %d : %d => blockSumNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			TestLog.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"int blockSumNxN %d : %d => blockSumNxNInternal %d = %.2fx\n", boxSize, boxSlowTotal, boxFastTotal,
					speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestLog.logSpeedTestResult(fastTotal < slowTotal, "int blockSumNxN %d => blockSumNxNInternal %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void intStripedBlockSumNxNIsFasterThanBlockSumNxN()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<int[]> dataSet = intCreateSpeedData(ITER);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.blockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.stripedBlockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (final int boxSize : boxSizes)
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (final int[] data : dataSet2)
						filter.stripedBlockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (final int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (final int[] data : dataSet2)
						filter.blockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					final long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf("int blockSumNxN [%dx%d] @ %d : %d => stripedBlockSumNxN %d = %.2fx\n", width,
								height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			TestLog.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"int blockSumNxN %d : %d => stripedBlockSumNxN %d = %.2fx\n", boxSize, boxSlowTotal, boxFastTotal,
					speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestLog.logSpeedTestResult(fastTotal < slowTotal, "int blockSumNxN %d => stripedBlockSumNxN %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void intStripedBlockSumInternalNxNIsFasterThanStripedBlockSumNxN()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<int[]> dataSet = intCreateSpeedData(ITER);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.stripedBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.stripedBlockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (final int boxSize : boxSizes)
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (final int[] data : dataSet2)
						filter.stripedBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (final int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (final int[] data : dataSet2)
						filter.stripedBlockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					final long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"int stripedBlockSumNxN [%dx%d] @ %d : %d => stripedBlockSumNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			TestLog.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"int stripedBlockSumNxN %d : %d => stripedBlockSumNxNInternal %d = %.2fx\n", boxSize, boxSlowTotal,
					boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"int stripedBlockSumNxN %d => stripedBlockSumNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void intRollingBlockSumNxNIsFasterThanBlockSumNxN()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<int[]> dataSet = intCreateSpeedData(ITER);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.blockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingBlockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (final int boxSize : boxSizes)
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (final int[] data : dataSet2)
						filter.rollingBlockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (final int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (final int[] data : dataSet2)
						filter.blockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					final long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf("int blockSumNxN [%dx%d] @ %d : %d => rollingBlockSumNxN %d = %.2fx\n", width,
								height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			TestLog.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"int blockSumNxN %d : %d => rollingBlockSumNxN %d = %.2fx\n", boxSize, boxSlowTotal, boxFastTotal,
					speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestLog.logSpeedTestResult(fastTotal < slowTotal, "int blockSumNxN %d => rollingBlockSumNxN %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void intRollingBlockSumInternalNxNIsFasterThanRollingBlockSumNxN()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<int[]> dataSet = intCreateSpeedData(ITER);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingBlockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (final int boxSize : boxSizes)
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (final int[] data : dataSet2)
						filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (final int boxSize : boxSizes)
		{
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (final int width : primes)
				for (final int height : primes)
				{
					final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
					for (final int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (final int[] data : dataSet2)
						filter.rollingBlockSumNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					final long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"int rollingBlockSumNxN [%dx%d] @ %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			TestLog.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"int rollingBlockSumNxN %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n", boxSize, boxSlowTotal,
					boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"int rollingBlockSumNxN %d => rollingBlockSumNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void intBlockSum3x3AndBlockSumNxNReturnSameResult()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SumFilter filter = new SumFilter();

		for (final int width : primes)
			for (final int height : primes)
				intCompareBlockSum3x3AndBlockSumNxN(rg, filter, width, height);
	}

	private static void intCompareBlockSum3x3AndBlockSumNxN(RandomGenerator rg, SumFilter filter, int width, int height)
			throws ArrayComparisonFailure
	{
		final int[] data1 = intCreateData(rg, width, height);
		final int[] data2 = intClone(data1);

		filter.blockSum3x3(data1, width, height);
		filter.blockSumNxN(data2, width, height, 1);

		intArrayEquals(data1, data2, 1, "Arrays do not match: [%dx%d] @ %d", width, height, 1);
	}

	@Test
	public void intBlockSum3x3IsFasterThanBlockSumNxN()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<int[]> dataSet = intCreateSpeedData(ITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.blockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.blockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.blockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.blockSumNxN(data, width, height, 1);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("int blockSumNxN [%dx%d] %d => blockSum3x3 %d = %.2fx\n", width, height, time,
							fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal, "int blockSumNxN %d => blockSum3x3 %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void intStripedBlockSum3x3AndStripedBlockSumNxNReturnSameResult()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SumFilter filter = new SumFilter();

		for (final int width : primes)
			for (final int height : primes)
				intCompareStripedBlockSum3x3AndStripedBlockSumNxN(rg, filter, width, height);
	}

	private static void intCompareStripedBlockSum3x3AndStripedBlockSumNxN(RandomGenerator rg, SumFilter filter,
			int width, int height) throws ArrayComparisonFailure
	{
		final int[] data1 = intCreateData(rg, width, height);
		final int[] data2 = intClone(data1);

		filter.stripedBlockSum3x3(data1, width, height);
		filter.stripedBlockSumNxN(data2, width, height, 1);

		intArrayEquals(data1, data2, 1, "Arrays do not match: [%dx%d] @ %d", width, height, 1);
	}

	@Test
	public void intStripedBlockSum3x3IsFasterThanStripedBlockSumNxN()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<int[]> dataSet = intCreateSpeedData(ITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.stripedBlockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.stripedBlockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.stripedBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.stripedBlockSumNxN(data, width, height, 1);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("int stripedBlockSumNxN [%dx%d] %d => stripedBlockSum3x3 %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		stripedBlockTime, time), stripedBlockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"int stripedBlockSumNxN %d => stripedBlockSum3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void intRollingBlockSum3x3AndRollingBlockSumNxNReturnSameResult()
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		final SumFilter filter = new SumFilter();

		for (final int width : primes)
			for (final int height : primes)
				intCompareRollingBlockSum3x3AndRollingBlockSumNxN(rg, filter, width, height);
	}

	private static void intCompareRollingBlockSum3x3AndRollingBlockSumNxN(RandomGenerator rg, SumFilter filter,
			int width, int height) throws ArrayComparisonFailure
	{
		final int[] data1 = intCreateData(rg, width, height);
		final int[] data2 = intClone(data1);

		filter.rollingBlockSum3x3(data1, width, height);
		filter.rollingBlockSumNxN(data2, width, height, 1);

		intArrayEquals(data1, data2, 1, "Arrays do not match: [%dx%d] @ %d", width, height, 1);
	}

	@Test
	public void intRollingBlockSum3x3IsFasterThanRollingBlockSumNxN()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<int[]> dataSet = intCreateSpeedData(ITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingBlockSumNxN(intClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.rollingBlockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.rollingBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.rollingBlockSumNxN(data, width, height, 1);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("int rollingBlockSumNxN [%dx%d] %d => rollingBlockSum3x3 %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		rollingBlockTime, time), rollingBlockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"int rollingBlockSumNxN %d => rollingBlockSum3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void intRollingBlockSum3x3IsFasterThanBlockSum3x3()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<int[]> dataSet = intCreateSpeedData(ITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingBlockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);
		filter.blockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.rollingBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.blockSum3x3(data, width, height);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("int blockSum3x3 [%dx%d] %d => rollingBlockSum3x3 %d = %.2fx\n", width, height,
							time, fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal, "int blockSum3x3 %d => rollingBlockSum3x3 %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void intStripedBlockSum3x3IsFasterThanBlockSum3x3()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<int[]> dataSet = intCreateSpeedData(ITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.stripedBlockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);
		filter.blockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.stripedBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.blockSum3x3(data, width, height);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("int blockSum3x3 [%dx%d] %d => stripedBlockSum3x3 %d = %.2fx\n", width, height,
							time, fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal, "int blockSum3x3 %d => stripedBlockSum3x3 %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void intRollingBlockSum3x3IsFasterThanStripedBlockSum3x3()
	{
		// These test a deprecated filter
		ExtraAssumptions.assumeSpeedTest(TestComplexity.VERY_HIGH);

		final SumFilter filter = new SumFilter();

		final ArrayList<int[]> dataSet = intCreateSpeedData(ITER3);

		final ArrayList<Long> fastTimes = new ArrayList<>();

		// Initialise
		filter.rollingBlockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);
		filter.stripedBlockSum3x3(intClone(dataSet.get(0)), primes[0], primes[0]);

		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.rollingBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;
				fastTimes.add(time);
			}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		@SuppressWarnings("unused")
		long boxSlowTotal = 0, boxFastTotal = 0;
		for (final int width : primes)
			for (final int height : primes)
			{
				final ArrayList<int[]> dataSet2 = new ArrayList<>(dataSet.size());
				for (final int[] data : dataSet)
					dataSet2.add(intClone(data));

				long time = System.nanoTime();
				for (final int[] data : dataSet2)
					filter.stripedBlockSum3x3(data, width, height);
				time = System.nanoTime() - time;

				final long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("int stripedBlockSum3x3 [%dx%d] %d => rollingBlockSum3x3 %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (ExtraAssertions.assert_SPEED_TESTS) Assertions.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		TestLog.logSpeedTestResult(fastTotal < slowTotal,
				"int stripedBlockSum3x3 %d => rollingBlockSum3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}
}
