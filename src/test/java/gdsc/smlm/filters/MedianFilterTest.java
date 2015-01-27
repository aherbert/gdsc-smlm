package gdsc.smlm.filters;

import gdsc.smlm.TestSettings;
import gdsc.smlm.filters.MedianFilter;

import java.util.ArrayList;
import java.util.Arrays;

import org.junit.Assert;
import org.junit.Test;
import org.junit.internal.ArrayComparisonFailure;

public class MedianFilterTest
{
	private gdsc.smlm.utils.Random rand;

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

	private double speedUpFactor(long slowTotal, long fastTotal)
	{
		return (1.0 * slowTotal) / fastTotal;
	}

	@Test
	public void floatBlockMedianNxNInternalAndRollingMedianNxNInternalReturnSameResult()
	{
		MedianFilter filter = new MedianFilter();

		for (int width : primes)
			for (int height : primes)
				for (int boxSize : boxSizes)
					floatCompareBlockMedianNxNInternalAndRollingMedianNxNInternal(filter, width, height, boxSize);
	}

	private void floatCompareBlockMedianNxNInternalAndRollingMedianNxNInternal(MedianFilter filter, int width,
			int height, int boxSize) throws ArrayComparisonFailure
	{
		rand = new gdsc.smlm.utils.Random(-30051976);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.blockMedianNxNInternal(data1, width, height, boxSize);
		filter.rollingMedianNxNInternal(data2, width, height, boxSize);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d] @ %d", width, height, boxSize), data1, data2,
				boxSize);
	}

	@Test
	public void floatBlockMedian3x3InternalAndRollingMedianNxNInternalReturnSameResult()
	{
		MedianFilter filter = new MedianFilter();

		for (int width : primes)
			for (int height : primes)
				floatCompareBlockMedian3x3InternalAndRollingMedianNxNInternal(filter, width, height);
	}

	private void floatCompareBlockMedian3x3InternalAndRollingMedianNxNInternal(MedianFilter filter, int width,
			int height) throws ArrayComparisonFailure
	{
		rand = new gdsc.smlm.utils.Random(-30051977);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.blockMedian3x3Internal(data1, width, height);
		filter.rollingMedianNxNInternal(data2, width, height, 1);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);
	}

	private float[] floatClone(float[] data1)
	{
		float[] data2 = Arrays.copyOf(data1, data1.length);
		return data2;
	}

	@Test
	public void floatRollingMedianNxNInternalIsFasterThanBlockMedianNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		MedianFilter filter = new MedianFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockMedianNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingMedianNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
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
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.blockMedianNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("float blockMedianNxNInternal [%dx%d] @ %d : %d => rollingMedianNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float blockMedianNxNInternal %d : %d => rollingMedianNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float blockMedianNxNInternal %d => rollingMedianNxNInternal %d = %.2fx\n", slowTotal,
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
	public void floatBlockMedian3x3InternalAndBlockMedianNxNInternalReturnSameResult()
	{
		MedianFilter filter = new MedianFilter();

		for (int width : primes)
			for (int height : primes)
				floatCompareBlockMedian3x3InternalAndBlockMedianNxNInternal(filter, width, height);
	}

	private void floatCompareBlockMedian3x3InternalAndBlockMedianNxNInternal(MedianFilter filter, int width,
			int height) throws ArrayComparisonFailure
	{
		rand = new gdsc.smlm.utils.Random(-30051977);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.blockMedian3x3Internal(data1, width, height);
		filter.blockMedianNxNInternal(data2, width, height, 1);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);
	}

	@Test
	public void floatBlockMedian3x3InternalIsFasterThanBlockMedianNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-30051977);

		MedianFilter filter = new MedianFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockMedianNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.blockMedian3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
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
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
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
					System.out.printf(
							"float blockMedianNxNInternal [%dx%d] %d => blockMedian3x3Internal %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float blockMedianNxNInternal %d => blockMedian3x3Internal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingMedian3x3InternalIsFasterThanBlockMedian3x3Internal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-30051977);

		MedianFilter filter = new MedianFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingMedian3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.blockMedian3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
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
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.blockMedian3x3Internal(data, width, height);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf(
							"float blockMedian3x3Internal [%dx%d] %d => rollingMedian3x3Internal %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float blockMedian3x3Internal %d => rollingMedian3x3Internal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}


	@Test
	public void floatRollingMedian3x3InternalAndRollingMedianNxNInternalReturnSameResult()
	{
		MedianFilter filter = new MedianFilter();

		for (int width : primes)
			for (int height : primes)
				floatCompareRollingMedian3x3InternalAndRollingMedianNxNInternal(filter, width, height);
	}

	private void floatCompareRollingMedian3x3InternalAndRollingMedianNxNInternal(MedianFilter filter,
			int width, int height) throws ArrayComparisonFailure
	{
		rand = new gdsc.smlm.utils.Random(-30051977);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.rollingMedian3x3Internal(data1, width, height);
		filter.rollingMedianNxNInternal(data2, width, height, 1);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);
	}

	@Test
	public void floatRollingMedian3x3InternalIsFasterThanRollingMedianNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-30051977);

		MedianFilter filter = new MedianFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingMedian3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.rollingMedianNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], 1);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
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
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
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
					System.out
							.printf("float rollingMedianNxNInternal [%dx%d] %d => rollingMedian3x3Internal %d = %.2fx\n",
									width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float rollingMedianNxNInternal %d => rollingMedian3x3Internal %d = %.2fx\n",
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
	public void floatBlockMedianNxNAndRollingMedianNxNReturnSameResult()
	{
		MedianFilter filter = new MedianFilter();

		for (int width : primes)
		{
			for (int height : primes)
				for (int boxSize : boxSizes)
					floatCompareBlockMedianNxNAndRollingMedianNxN(filter, width, height, boxSize);
		}
	}

	private void floatCompareBlockMedianNxNAndRollingMedianNxN(MedianFilter filter, int width, int height,
			int boxSize) throws ArrayComparisonFailure
	{
		rand = new gdsc.smlm.utils.Random(-30051976);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.blockMedianNxN(data1, width, height, boxSize);
		filter.rollingMedianNxN(data2, width, height, boxSize);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d] @ %d", width, height, boxSize), data1, data2, 0);
	}

	@Test
	public void floatBlockMedianInternalNxNIsFasterThanBlockMedianNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		MedianFilter filter = new MedianFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockMedianNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.blockMedianNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
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
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
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
								"float blockMedianNxN [%dx%d] @ %d : %d => blockMedianNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float blockMedianNxN %d : %d => blockMedianNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float blockMedianNxN %d => blockMedianNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingMedianNxNIsFasterThanBlockMedianNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		MedianFilter filter = new MedianFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockMedianNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingMedianNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.rollingMedianNxN(data, width, height, boxSize);
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
						filter.blockMedianNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"float blockMedianNxN [%dx%d] @ %d : %d => rollingMedianNxN %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float blockMedianNxN %d : %d => rollingMedianNxN %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float blockMedianNxN %d => rollingMedianNxN %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingMedianInternalNxNIsFasterThanRollingMedianNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		MedianFilter filter = new MedianFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingMedianNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingMedianNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
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
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
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
						System.out
								.printf("float rollingMedianNxN [%dx%d] @ %d : %d => rollingMedianNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float rollingMedianNxN %d : %d => rollingMedianNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float rollingMedianNxN %d => rollingMedianNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatBlockMedian3x3AndBlockMedianNxNReturnSameResult()
	{
		MedianFilter filter = new MedianFilter();

		for (int width : primes)
			for (int height : primes)
				floatCompareBlockMedian3x3AndBlockMedianNxN(filter, width, height);
	}

	private void floatCompareBlockMedian3x3AndBlockMedianNxN(MedianFilter filter, int width, int height)
			throws ArrayComparisonFailure
	{
		rand = new gdsc.smlm.utils.Random(-30051977);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.blockMedian3x3(data1, width, height);
		filter.blockMedianNxN(data2, width, height, 1);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);
	}

	@Test
	public void floatBlockMedian3x3IsFasterThanBlockMedianNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-30051977);

		MedianFilter filter = new MedianFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockMedianNxN(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.blockMedian3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
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
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
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
					System.out.printf("float blockMedianNxN [%dx%d] %d => blockMedian3x3 %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float blockMedianNxN %d => blockMedian3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingMedian3x3AndRollingMedianNxNReturnSameResult()
	{
		MedianFilter filter = new MedianFilter();

		for (int width : primes)
			for (int height : primes)
				floatCompareRollingMedian3x3AndRollingMedianNxN(filter, width, height);
	}

	private void floatCompareRollingMedian3x3AndRollingMedianNxN(MedianFilter filter, int width, int height)
			throws ArrayComparisonFailure
	{
		rand = new gdsc.smlm.utils.Random(-30051977);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.rollingMedian3x3(data1, width, height);
		filter.rollingMedianNxN(data2, width, height, 1);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);
	}

	@Test
	public void floatRollingMedian3x3IsFasterThanRollingMedianNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-30051977);

		MedianFilter filter = new MedianFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingMedianNxN(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.rollingMedian3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
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
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
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
					System.out.printf("float rollingMedianNxN [%dx%d] %d => rollingMedian3x3 %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		rollingBlockTime, time), rollingBlockTime < time);
			}
		System.out.printf("float rollingMedianNxN %d => rollingMedian3x3 %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingMedian3x3IsFasterThanBlockMedian3x3()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-30051977);

		MedianFilter filter = new MedianFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingMedian3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.blockMedian3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
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
		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
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
					System.out.printf("float blockMedian3x3 [%dx%d] %d => rollingMedian3x3 %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float blockMedian3x3 %d => rollingMedian3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}
}
