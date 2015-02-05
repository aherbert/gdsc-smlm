package gdsc.smlm.filters;

import gdsc.smlm.TestSettings;
import gdsc.smlm.filters.AverageFilter;
import gdsc.smlm.utils.FloatEquality;

import java.util.ArrayList;
import java.util.Arrays;

import org.junit.Assert;
import org.junit.Test;
import org.junit.internal.ArrayComparisonFailure;

public class AverageFilterTest
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
		floatArrayEquals(message, data1, data2, 0, 0, boxSize);
	}

	private void floatArrayEquals(String message, float[] data1, float[] data2, int maxx, int maxy, float boxSize)
	{
		FloatEquality eq = new FloatEquality(1e-5f, 1e-10f);
		//gdsc.smlm.ij.utils.Utils.display("data1", new ij.process.FloatProcessor(maxx, maxy, data1));
		//gdsc.smlm.ij.utils.Utils.display("data2", new ij.process.FloatProcessor(maxx, maxy, data2));
		
		// Ignore the border
		int border = (int) Math.ceil(boxSize);
		for (int y = border; y < maxy - border - 1; y++)
		{
			int index = y * maxx + border;
			for (int x = border; x < maxx - border - 1; x++, index++)
			{
				if (!eq.almostEqualRelativeOrAbsolute(data1[index], data2[index]))
				{
					Assert.assertTrue(message + " [" + x + "," + y + "]", false);
				}
			}
		}
	}

	private double speedUpFactor(long slowTotal, long fastTotal)
	{
		return (1.0 * slowTotal) / fastTotal;
	}

	// COPY CODE FROM HERE...
	@Test
	public void floatBlockAverageNxNInternalAndRollingBlockAverageNxNInternalReturnSameResult()
	{
		AverageFilter filter = new AverageFilter();

		for (int width : primes)
			for (int height : primes)
				for (int boxSize : boxSizes)
					floatCompareBlockAverageNxNInternalAndRollingBlockAverageNxNInternal(filter, width, height, boxSize);
	}

	private void floatCompareBlockAverageNxNInternalAndRollingBlockAverageNxNInternal(AverageFilter filter, int width,
			int height, int boxSize) throws ArrayComparisonFailure
	{
		rand = new gdsc.smlm.utils.Random(-30051976);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.blockAverageNxNInternal(data1, width, height, boxSize);
		filter.rollingBlockAverageNxNInternal(data2, width, height, boxSize);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d] @ %d", width, height, boxSize), data1, data2,
				boxSize);
	}

	@Test
	public void floatBlockAverageNxNInternalAndStripedBlockAverageNxNInternalReturnSameResult()
	{
		AverageFilter filter = new AverageFilter();

		for (int width : primes)
			for (int height : primes)
				for (int boxSize : boxSizes)
					floatCompareBlockAverageNxNInternalAndStripedBlockAverageNxNInternal(filter, width, height, boxSize);
	}

	private void floatCompareBlockAverageNxNInternalAndStripedBlockAverageNxNInternal(AverageFilter filter, int width,
			int height, int boxSize) throws ArrayComparisonFailure
	{
		rand = new gdsc.smlm.utils.Random(-30051976);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.blockAverageNxNInternal(data1, width, height, boxSize);
		filter.stripedBlockAverageNxNInternal(data2, width, height, boxSize);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d] @ %d", width, height, boxSize), data1, data2, 0);
	}

	@Test
	public void floatBlockAverage3x3InternalAndRollingBlockAverageNxNInternalReturnSameResult()
	{
		AverageFilter filter = new AverageFilter();

		for (int width : primes)
			for (int height : primes)
				floatCompareBlockAverage3x3InternalAndRollingBlockAverageNxNInternal(filter, width, height);
	}

	private void floatCompareBlockAverage3x3InternalAndRollingBlockAverageNxNInternal(AverageFilter filter, int width,
			int height) throws ArrayComparisonFailure
	{
		rand = new gdsc.smlm.utils.Random(-30051977);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.blockAverage3x3Internal(data1, width, height);
		filter.rollingBlockAverageNxNInternal(data2, width, height, 1);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);
	}

	@Test
	public void floatRollingBlockAverageNxNInternalAndRollingBlockAverageNxNInternalTransposedReturnSameResult()
	{
		AverageFilter filter = new AverageFilter();

		for (int width : primes)
			for (int height : primes)
				for (int boxSize : boxSizes)
					floatCompareRollingBlockAverageNxNInternalAndRollingBlockAverageNxNInternalTransposed(filter,
							width, height, boxSize);
	}

	@SuppressWarnings("deprecation")
	private void floatCompareRollingBlockAverageNxNInternalAndRollingBlockAverageNxNInternalTransposed(
			AverageFilter filter, int width, int height, int boxSize) throws ArrayComparisonFailure
	{
		rand = new gdsc.smlm.utils.Random(-30051977);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.rollingBlockAverageNxNInternal(data1, width, height, boxSize);
		filter.rollingBlockAverageNxNInternalTransposed(data2, width, height, boxSize);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d] @ %d", width, height, boxSize), data1, data2,
				boxSize);
	}

	private float[] floatClone(float[] data1)
	{
		float[] data2 = Arrays.copyOf(data1, data1.length);
		return data2;
	}

	@Test
	public void floatRollingBlockAverageNxNInternalIsFasterThanBlockAverageNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		AverageFilter filter = new AverageFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockAverageNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingBlockAverageNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.rollingBlockAverageNxNInternal(data, width, height, boxSize);
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
						filter.blockAverageNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("float blockAverageNxNInternal [%dx%d] @ %d : %d => rollingBlockAverageNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float blockAverageNxNInternal %d : %d => rollingBlockAverageNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float blockAverageNxNInternal %d => rollingBlockAverageNxNInternal %d = %.2fx\n", slowTotal,
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
	public void floatStripedBlockAverageNxNInternalIsFasterThanBlockAverageNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		AverageFilter filter = new AverageFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockAverageNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.stripedBlockAverageNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.stripedBlockAverageNxNInternal(data, width, height, boxSize);
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
						filter.blockAverageNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("float blockAverageNxNInternal [%dx%d] @ %d : %d => stripedBlockAverageNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float blockAverageNxNInternal %d : %d => stripedBlockAverageNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float blockAverageNxNInternal %d => stripedBlockAverageNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockAverageNxNInternalIsFasterThanStripedBlockAverageNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		AverageFilter filter = new AverageFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockAverageNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.stripedBlockAverageNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.rollingBlockAverageNxNInternal(data, width, height, boxSize);
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
						filter.stripedBlockAverageNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("float stripedBlockAverageNxNInternal [%dx%d] @ %d : %d => rollingBlockAverageNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf(
					"float stripedBlockAverageNxNInternal %d : %d => rollingBlockAverageNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float stripedBlockAverageNxNInternal %d => rollingBlockAverageNxNInternal %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatBlockAverage3x3InternalAndBlockAverageNxNInternalReturnSameResult()
	{
		AverageFilter filter = new AverageFilter();

		for (int width : primes)
			for (int height : primes)
				floatCompareBlockAverage3x3InternalAndBlockAverageNxNInternal(filter, width, height);
	}

	private void floatCompareBlockAverage3x3InternalAndBlockAverageNxNInternal(AverageFilter filter, int width,
			int height) throws ArrayComparisonFailure
	{
		rand = new gdsc.smlm.utils.Random(-30051977);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.blockAverage3x3Internal(data1, width, height);
		filter.blockAverageNxNInternal(data2, width, height, 1);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);
	}

	@Test
	public void floatBlockAverage3x3InternalIsFasterThanBlockAverageNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-30051977);

		AverageFilter filter = new AverageFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockAverageNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.blockAverage3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.blockAverage3x3Internal(data, width, height);
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
					filter.blockAverageNxNInternal(data, width, height, 1);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf(
							"float blockAverageNxNInternal [%dx%d] %d => blockAverage3x3Internal %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float blockAverageNxNInternal %d => blockAverage3x3Internal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockAverage3x3InternalIsFasterThanBlockAverage3x3Internal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-30051977);

		AverageFilter filter = new AverageFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockAverage3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.blockAverage3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.rollingBlockAverage3x3Internal(data, width, height);
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
					filter.blockAverage3x3Internal(data, width, height);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf(
							"float blockAverage3x3Internal [%dx%d] %d => rollingBlockAverage3x3Internal %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float blockAverage3x3Internal %d => rollingBlockAverage3x3Internal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatStripedBlockAverage3x3InternalIsFasterThanBlockAverage3x3Internal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-30051977);

		AverageFilter filter = new AverageFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.stripedBlockAverage3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.blockAverage3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.stripedBlockAverage3x3Internal(data, width, height);
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
					filter.blockAverage3x3Internal(data, width, height);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf(
							"float blockAverage3x3Internal [%dx%d] %d => stripedBlockAverage3x3Internal %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float blockAverage3x3Internal %d => stripedBlockAverage3x3Internal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockAverage3x3InternalIsFasterThanStripedBlockAverage3x3Internal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-30051977);

		AverageFilter filter = new AverageFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockAverage3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.stripedBlockAverage3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.rollingBlockAverage3x3Internal(data, width, height);
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
					filter.stripedBlockAverage3x3Internal(data, width, height);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out
							.printf("float stripedBlockAverage3x3Internal [%dx%d] %d => rollingBlockAverage3x3Internal %d = %.2fx\n",
									width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float stripedBlockAverage3x3Internal %d => rollingBlockAverage3x3Internal %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockAverage3x3InternalAndRollingBlockAverageNxNInternalReturnSameResult()
	{
		AverageFilter filter = new AverageFilter();

		for (int width : primes)
			for (int height : primes)
				floatCompareRollingBlockAverage3x3InternalAndRollingBlockAverageNxNInternal(filter, width, height);
	}

	private void floatCompareRollingBlockAverage3x3InternalAndRollingBlockAverageNxNInternal(AverageFilter filter,
			int width, int height) throws ArrayComparisonFailure
	{
		rand = new gdsc.smlm.utils.Random(-30051977);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.rollingBlockAverage3x3Internal(data1, width, height);
		filter.rollingBlockAverageNxNInternal(data2, width, height, 1);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);
	}

	@Test
	public void floatRollingBlockAverage3x3InternalIsFasterThanRollingBlockAverageNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-30051977);

		AverageFilter filter = new AverageFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockAverage3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.rollingBlockAverageNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], 1);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.rollingBlockAverage3x3Internal(data, width, height);
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
					filter.rollingBlockAverageNxNInternal(data, width, height, 1);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out
							.printf("float rollingBlockAverageNxNInternal [%dx%d] %d => rollingBlockAverage3x3Internal %d = %.2fx\n",
									width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float rollingBlockAverageNxNInternal %d => rollingBlockAverage3x3Internal %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatStripedBlockAverage3x3InternalIsFasterThanStripedBlockAverageNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-30051977);

		AverageFilter filter = new AverageFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.stripedBlockAverageNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.stripedBlockAverage3x3Internal(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.stripedBlockAverage3x3Internal(data, width, height);
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
					filter.stripedBlockAverageNxNInternal(data, width, height, 1);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out
							.printf("float stripedBlockAverageNxNInternal [%dx%d] %d => stripedBlockAverage3x3Internal %d = %.2fx\n",
									width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float stripedBlockAverageNxNInternal %d => stripedBlockAverage3x3Internal %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@SuppressWarnings("deprecation")
	@Test
	public void floatRollingBlockAverageNxNInternalIsFasterThanRollingBlockAverageNxNInternalTransposed()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		AverageFilter filter = new AverageFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockAverageNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingBlockAverageNxNInternalTransposed(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.rollingBlockAverageNxNInternal(data, width, height, boxSize);
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
						filter.rollingBlockAverageNxNInternalTransposed(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("float rollingBlockAverageNxNInternalTransposed [%dx%d] @ %d : %d => rollingBlockAverageNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			if (debug)
				System.out
						.printf("float rollingBlockAverageNxNInternalTransposed %d : %d => rollingBlockAverageNxNInternal %d = %.2fx\n",
								boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		System.out.printf(
				"float rollingBlockAverageNxNInternalTransposed %d => rollingBlockAverageNxNInternal %d = %.2fx\n",
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
	public void floatBlockAverageNxNAndStripedBlockAverageNxNReturnSameResult()
	{
		AverageFilter filter = new AverageFilter();

		for (int width : primes)
			for (int height : primes)
				for (int boxSize : boxSizes)
					floatCompareBlockAverageNxNAndStripedBlockAverageNxN(filter, width, height, boxSize);
	}

	private void floatCompareBlockAverageNxNAndStripedBlockAverageNxN(AverageFilter filter, int width, int height,
			int boxSize) throws ArrayComparisonFailure
	{
		rand = new gdsc.smlm.utils.Random(-30051976);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.blockAverageNxN(data1, width, height, boxSize);
		filter.stripedBlockAverageNxN(data2, width, height, boxSize);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d] @ %d", width, height, boxSize), data1, data2, 0);
	}

	@Test
	public void floatBlockAverageNxNAndRollingBlockAverageNxNReturnSameResult()
	{
		AverageFilter filter = new AverageFilter();

		for (int width : primes)
			for (int height : primes)
				for (int boxSize : boxSizes)
					floatCompareBlockAverageNxNAndRollingBlockAverageNxN(filter, width, height, boxSize);
	}

	private void floatCompareBlockAverageNxNAndRollingBlockAverageNxN(AverageFilter filter, int width, int height,
			int boxSize) throws ArrayComparisonFailure
	{
		rand = new gdsc.smlm.utils.Random(-30051976);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.blockAverageNxN(data1, width, height, boxSize);
		filter.rollingBlockAverageNxN(data2, width, height, boxSize);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d] @ %d", width, height, boxSize), data1, data2, 0);
	}

	@Test
	public void floatBlockAverageInternalNxNIsFasterThanBlockAverageNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		AverageFilter filter = new AverageFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockAverageNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.blockAverageNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.blockAverageNxNInternal(data, width, height, boxSize);
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
						filter.blockAverageNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"float blockAverageNxN [%dx%d] @ %d : %d => blockAverageNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float blockAverageNxN %d : %d => blockAverageNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float blockAverageNxN %d => blockAverageNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatStripedBlockAverageNxNIsFasterThanBlockAverageNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		AverageFilter filter = new AverageFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockAverageNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.stripedBlockAverageNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.stripedBlockAverageNxN(data, width, height, boxSize);
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
						filter.blockAverageNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"float blockAverageNxN [%dx%d] @ %d : %d => stripedBlockAverageNxN %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float blockAverageNxN %d : %d => stripedBlockAverageNxN %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float blockAverageNxN %d => stripedBlockAverageNxN %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatStripedBlockAverageInternalNxNIsFasterThanStripedBlockAverageNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		AverageFilter filter = new AverageFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.stripedBlockAverageNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.stripedBlockAverageNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.stripedBlockAverageNxNInternal(data, width, height, boxSize);
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
						filter.stripedBlockAverageNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("float stripedBlockAverageNxN [%dx%d] @ %d : %d => stripedBlockAverageNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float stripedBlockAverageNxN %d : %d => stripedBlockAverageNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float stripedBlockAverageNxN %d => stripedBlockAverageNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockAverageNxNIsFasterThanBlockAverageNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		AverageFilter filter = new AverageFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockAverageNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingBlockAverageNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.rollingBlockAverageNxN(data, width, height, boxSize);
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
						filter.blockAverageNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"float blockAverageNxN [%dx%d] @ %d : %d => rollingBlockAverageNxN %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float blockAverageNxN %d : %d => rollingBlockAverageNxN %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float blockAverageNxN %d => rollingBlockAverageNxN %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockAverageInternalNxNIsFasterThanRollingBlockAverageNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		AverageFilter filter = new AverageFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockAverageNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingBlockAverageNxN(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.rollingBlockAverageNxNInternal(data, width, height, boxSize);
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
						filter.rollingBlockAverageNxN(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("float rollingBlockAverageNxN [%dx%d] @ %d : %d => rollingBlockAverageNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float rollingBlockAverageNxN %d : %d => rollingBlockAverageNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float rollingBlockAverageNxN %d => rollingBlockAverageNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatBlockAverage3x3AndBlockAverageNxNReturnSameResult()
	{
		AverageFilter filter = new AverageFilter();

		for (int width : primes)
			for (int height : primes)
				floatCompareBlockAverage3x3AndBlockAverageNxN(filter, width, height);
	}

	private void floatCompareBlockAverage3x3AndBlockAverageNxN(AverageFilter filter, int width, int height)
			throws ArrayComparisonFailure
	{
		rand = new gdsc.smlm.utils.Random(-30051977);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.blockAverage3x3(data1, width, height);
		filter.blockAverageNxN(data2, width, height, 1);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);

		for (float f : new float[] { 0.3f, 0.6f })
		{
			filter.blockAverage3x3(data1, width, height, f);
			filter.blockAverageNxN(data2, width, height, f);
			floatArrayEquals(String.format("Arrays do not match: [%dx%d] w="+f, width, height), data1, data2, 1);
		}
	}

	@Test
	public void floatBlockAverage3x3IsFasterThanBlockAverageNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-30051977);

		AverageFilter filter = new AverageFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.blockAverageNxN(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.blockAverage3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.blockAverage3x3(data, width, height);
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
					filter.blockAverageNxN(data, width, height, 1);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float blockAverageNxN [%dx%d] %d => blockAverage3x3 %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float blockAverageNxN %d => blockAverage3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatStripedBlockAverage3x3AndStripedBlockAverageNxNReturnSameResult()
	{
		AverageFilter filter = new AverageFilter();

		for (int width : primes)
			for (int height : primes)
				floatCompareStripedBlockAverage3x3AndStripedBlockAverageNxN(filter, width, height);
	}

	private void floatCompareStripedBlockAverage3x3AndStripedBlockAverageNxN(AverageFilter filter, int width, int height)
			throws ArrayComparisonFailure
	{
		rand = new gdsc.smlm.utils.Random(-30051977);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.stripedBlockAverage3x3(data1, width, height);
		filter.stripedBlockAverageNxN(data2, width, height, 1);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);
	}

	@Test
	public void floatStripedBlockAverage3x3IsFasterThanStripedBlockAverageNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-30051977);

		AverageFilter filter = new AverageFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.stripedBlockAverageNxN(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.stripedBlockAverage3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.stripedBlockAverage3x3(data, width, height);
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
					filter.stripedBlockAverageNxN(data, width, height, 1);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float stripedBlockAverageNxN [%dx%d] %d => stripedBlockAverage3x3 %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		stripedBlockTime, time), stripedBlockTime < time);
			}
		System.out.printf("float stripedBlockAverageNxN %d => stripedBlockAverage3x3 %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockAverage3x3AndRollingBlockAverageNxNReturnSameResult()
	{
		AverageFilter filter = new AverageFilter();

		for (int width : primes)
			for (int height : primes)
				floatCompareRollingBlockAverage3x3AndRollingBlockAverageNxN(filter, width, height);
	}

	private void floatCompareRollingBlockAverage3x3AndRollingBlockAverageNxN(AverageFilter filter, int width, int height)
			throws ArrayComparisonFailure
	{
		rand = new gdsc.smlm.utils.Random(-30051977);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		filter.rollingBlockAverage3x3(data1, width, height);
		filter.rollingBlockAverageNxN(data2, width, height, 1);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d]", width, height), data1, data2, 1);
	}

	@Test
	public void floatRollingBlockAverage3x3IsFasterThanRollingBlockAverageNxN()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-30051977);

		AverageFilter filter = new AverageFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockAverageNxN(floatClone(dataSet.get(0)), primes[0], primes[0], 1);
		filter.rollingBlockAverage3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.rollingBlockAverage3x3(data, width, height);
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
					filter.rollingBlockAverageNxN(data, width, height, 1);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float rollingBlockAverageNxN [%dx%d] %d => rollingBlockAverage3x3 %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		rollingBlockTime, time), rollingBlockTime < time);
			}
		System.out.printf("float rollingBlockAverageNxN %d => rollingBlockAverage3x3 %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockAverage3x3IsFasterThanBlockAverage3x3()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-30051977);

		AverageFilter filter = new AverageFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockAverage3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.blockAverage3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.rollingBlockAverage3x3(data, width, height);
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
					filter.blockAverage3x3(data, width, height);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float blockAverage3x3 [%dx%d] %d => rollingBlockAverage3x3 %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float blockAverage3x3 %d => rollingBlockAverage3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatStripedBlockAverage3x3IsFasterThanBlockAverage3x3()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-30051977);

		AverageFilter filter = new AverageFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.stripedBlockAverage3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.blockAverage3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.stripedBlockAverage3x3(data, width, height);
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
					filter.blockAverage3x3(data, width, height);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float blockAverage3x3 [%dx%d] %d => stripedBlockAverage3x3 %d = %.2fx\n", width,
							height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float blockAverage3x3 %d => stripedBlockAverage3x3 %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockAverage3x3IsFasterThanStripedBlockAverage3x3()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-30051977);

		AverageFilter filter = new AverageFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockAverage3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);
		filter.stripedBlockAverage3x3(floatClone(dataSet.get(0)), primes[0], primes[0]);

		for (int width : primes)
			for (int height : primes)
			{
				ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
				for (float[] data : dataSet)
					dataSet2.add(floatClone(data));

				long time = System.nanoTime();
				for (float[] data : dataSet2)
					filter.rollingBlockAverage3x3(data, width, height);
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
					filter.stripedBlockAverage3x3(data, width, height);
				time = System.nanoTime() - time;

				long fastTime = fastTimes.get(index++);
				slowTotal += time;
				fastTotal += fastTime;
				boxSlowTotal += time;
				boxFastTotal += fastTime;
				if (debug)
					System.out.printf("float stripedBlockAverage3x3 [%dx%d] %d => rollingBlockAverage3x3 %d = %.2fx\n",
							width, height, time, fastTime, speedUpFactor(time, fastTime));
				//if (TestSettings.ASSERT_SPEED_TESTS) Assert.assertTrue(String.format("Not faster: [%dx%d] %d > %d", width, height,
				//		blockTime, time), blockTime < time);
			}
		System.out.printf("float stripedBlockAverage3x3 %d => rollingBlockAverage3x3 %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	/**
	 * Do a simple and stupid mean filter
	 * 
	 * @param data
	 * @param maxx
	 * @param maxy
	 * @param boxSize
	 */
	public static void average(float[] data, int maxx, int maxy, float boxSize)
	{
		if (boxSize <= 0)
			return;

		int n = (int) Math.ceil(boxSize);
		int size = 2 * n + 1;
		float[] weight = new float[size];
		Arrays.fill(weight, 1);
		if (boxSize != n)
			weight[0] = weight[weight.length - 1] = boxSize - (n - 1);

		float norm = 0;
		for (int yy = 0; yy < size; yy++)
			for (int xx = 0; xx < size; xx++)
				norm += weight[yy] * weight[xx];
		norm = (float) (1.0 / norm);

		float[] out = new float[data.length];

		for (int y = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++)
			{
				float sum = 0;
				for (int yy = 0; yy < size; yy++)
				{
					int yyy = y + yy - n;
					if (yyy < 0)
						yyy = 0;
					if (yyy >= maxy)
						yyy = maxy - 1;
					for (int xx = 0; xx < size; xx++)
					{
						int xxx = x + xx - n;
						if (xxx < 0)
							xxx = 0;
						if (xxx >= maxx)
							xxx = maxx - 1;
						int index = yyy * maxx + xxx;
						sum += data[index] * weight[yy] * weight[xx];
					}
				}
				out[y * maxx + x] = sum * norm;
			}
		}
		System.arraycopy(out, 0, data, 0, out.length);
	}

	@Test
	public void floatBlockAverageInternalInterpolatedIsCorrect()
	{
		for (int width : primes)
			for (int height : primes)
				for (float boxSize : boxSizes)
				{
					floatBlockAverageInternalInterpolatedIsCorrect(width, height, boxSize + 0.3f);
					floatBlockAverageInternalInterpolatedIsCorrect(width, height, boxSize + 0.6f);
				}
	}

	private void floatBlockAverageInternalInterpolatedIsCorrect(int width, int height, float boxSize)
			throws ArrayComparisonFailure
	{
		rand = new gdsc.smlm.utils.Random(-30051976);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = data1.clone();

		AverageFilter filter = new AverageFilter();
		AverageFilterTest.average(data1, width, height, boxSize);
		filter.blockAverageInternal(data2, width, height, boxSize);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d] @ %.1f", width, height, boxSize), data1, data2,
				width, height, boxSize);
	}

	@Test
	public void floatBlockAverageInterpolatedIsCorrect()
	{
		for (int width : primes)
			for (int height : primes)
				for (float boxSize : boxSizes)
				{
					floatBlockAverageInterpolatedIsCorrect(width, height, boxSize + 0.3f);
					floatBlockAverageInterpolatedIsCorrect(width, height, boxSize + 0.6f);
				}
	}

	private void floatBlockAverageInterpolatedIsCorrect(int width, int height, float boxSize)
			throws ArrayComparisonFailure
	{
		rand = new gdsc.smlm.utils.Random(-30051976);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = data1.clone();

		AverageFilter filter = new AverageFilter();
		AverageFilterTest.average(data1, width, height, boxSize);
		filter.blockAverage(data2, width, height, boxSize);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d] @ %.1f", width, height, boxSize), data1, data2, 0,
				0, boxSize);
	}
}
