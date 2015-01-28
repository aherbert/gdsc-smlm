package gdsc.smlm.filters;

import gdsc.smlm.TestSettings;
import gdsc.smlm.filters.AverageFilter;
import gdsc.smlm.filters.SumFilter;

import java.util.ArrayList;
import java.util.Arrays;

import org.junit.Assert;
import org.junit.Test;

public class FilterTest
{
	private gdsc.smlm.utils.Random rand;

	private boolean debug = false;

	int[] primes = new int[] { 113, 97, 53, 29 };
	int[] boxSizes = new int[] { 15, 9, 5, 3, 2, 1 };

	private float[] floatClone(int[] data1)
	{
		float[] data2 = new float[data1.length];
		for (int i = data2.length; i-- > 0;)
			data2[i] = data1[i];
		return data2;
	}

	private float[] floatClone(float[] data1)
	{
		float[] data2 = Arrays.copyOf(data1, data1.length);
		return data2;
	}

	private float[] floatCreateData(int width, int height)
	{
		float[] data = new float[width * height];
		for (int i = data.length; i-- > 0;)
			//data[i] = i;
			data[i] = rand.next();
		//rand.shuffle(data);

		return data;
	}

	private int[] intClone(int[] data1)
	{
		int[] data2 = Arrays.copyOf(data1, data1.length);
		return data2;
	}

	private int[] intCreateData(int width, int height)
	{
		int[] data = new int[width * height];
		for (int i = data.length; i-- > 0;)
			data[i] = i;

		rand.shuffle(data);

		return data;
	}

	private double speedUpFactor(long slowTotal, long fastTotal)
	{
		return (1.0 * slowTotal) / fastTotal;
	}

	@Test
	public void intRollingBlockSumNxNInternalIsFasterThanRollingBlockAverageNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		SumFilter filter = new SumFilter();
		AverageFilter filter2 = new AverageFilter();

		int iter = 50;

		ArrayList<int[]> dataSet = new ArrayList<int[]>(iter);
		for (int i = iter; i-- > 0;)
		{
			dataSet.add(intCreateData(primes[0], primes[0]));
		}

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter2.rollingBlockAverageNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<int[]> dataSet2 = new ArrayList<int[]>(iter);
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
					ArrayList<int[]> dataSet2 = new ArrayList<int[]>(iter);
					for (int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (int[] data : dataSet2)
						filter2.rollingBlockAverageNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("int rollingBlockAverageNxNInternal [%dx%d] @ %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			if (debug)
				System.out.printf(
						"int rollingBlockAverageNxNInternal %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
						boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		System.out.printf("int rollingBlockAverageNxNInternal %d => rollingBlockSumNxNInternal %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockSumNxNInternalIsFasterThanRollingBlockAverageNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		SumFilter filter = new SumFilter();
		AverageFilter filter2 = new AverageFilter();

		int iter = 50;

		ArrayList<float[]> dataSet = new ArrayList<float[]>(iter);
		for (int i = iter; i-- > 0;)
		{
			dataSet.add(floatCreateData(primes[0], primes[0]));
		}

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter2.rollingBlockAverageNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(iter);
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
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(iter);
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter2.rollingBlockAverageNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("float rollingBlockAverageNxNInternal [%dx%d] @ %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf(
					"float rollingBlockAverageNxNInternal %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float rollingBlockAverageNxNInternal %d => rollingBlockSumNxNInternal %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockAverageNxNInternalIsFasterThanBlockMedianNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		AverageFilter filter1 = new AverageFilter();
		MedianFilter filter2 = new MedianFilter();

		int iter = 10;

		ArrayList<float[]> dataSet = new ArrayList<float[]>(iter);
		for (int i = iter; i-- > 0;)
		{
			dataSet.add(floatCreateData(primes[0], primes[0]));
		}

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter1.rollingBlockAverageNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter2.blockMedianNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(iter);
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter1.rollingBlockAverageNxNInternal(data, width, height, boxSize);
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
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(iter);
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter2.blockMedianNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("float blockMedianNxNInternal [%dx%d] @ %d : %d => rollingBlockAverageNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float blockMedianNxNInternal %d : %d => rollingBlockAverageNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float blockMedianNxNInternal %d => rollingBlockAverageNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockAverageNxNInternalIsFasterThanRollingMedianNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		AverageFilter filter1 = new AverageFilter();
		MedianFilter filter2 = new MedianFilter();

		int iter = 10;

		ArrayList<float[]> dataSet = new ArrayList<float[]>(iter);
		for (int i = iter; i-- > 0;)
		{
			dataSet.add(floatCreateData(primes[0], primes[0]));
		}

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter1.rollingBlockAverageNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter2.rollingMedianNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(iter);
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter1.rollingBlockAverageNxNInternal(data, width, height, boxSize);
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
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(iter);
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter2.rollingMedianNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("float rollingMedianNxNInternal [%dx%d] @ %d : %d => rollingBlockAverageNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float rollingMedianNxNInternal %d : %d => rollingBlockAverageNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float rollingMedianNxNInternal %d => rollingBlockAverageNxNInternal %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockAverageNxNInternalIsFasterThanGaussianNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		AverageFilter filter1 = new AverageFilter();
		GaussianFilter filter2 = new GaussianFilter();

		int iter = 10;

		ArrayList<float[]> dataSet = new ArrayList<float[]>(iter);
		for (int i = iter; i-- > 0;)
		{
			dataSet.add(floatCreateData(primes[0], primes[0]));
		}

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter1.rollingBlockAverageNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter2.convolveInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0] / 3.0);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(iter);
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter1.rollingBlockAverageNxNInternal(data, width, height, boxSize);
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
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(iter);
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter2.convolveInternal(data, width, height, boxSize / 3.0);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("float convolveInternal [%dx%d] @ %d : %d => rollingBlockAverageNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float convolveInternal %d : %d => rollingBlockAverageNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float convolveInternal %d => rollingBlockAverageNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockAverageNxNInternalIsFasterThanAreaAverageNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		AverageFilter filter1 = new AverageFilter();
		AreaAverageFilter filter2 = new AreaAverageFilter();

		int iter = 10;

		ArrayList<float[]> dataSet = new ArrayList<float[]>(iter);
		for (int i = iter; i-- > 0;)
		{
			dataSet.add(floatCreateData(primes[0], primes[0]));
		}

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		//filter1.rollingBlockAverageNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		//filter2.areaAverageInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0] - 0.05);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(iter);
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					// Initialise
					for (float[] data : dataSet2)
						filter1.rollingBlockAverageNxNInternal(data.clone(), width, height, boxSize);
					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter1.rollingBlockAverageNxNInternal(data, width, height, boxSize);
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
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(iter);
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					// Initialise
					for (float[] data : dataSet2)
						filter2.areaAverageInternal(data.clone(), width, height, boxSize - 0.05);
					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter2.areaAverageInternal(data, width, height, boxSize - 0.05);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("float areaAverageInternal [%dx%d] @ %d : %d => rollingBlockAverageNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float areaAverageInternal %d : %d => rollingBlockAverageNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float areaAverageInternal %d => rollingBlockAverageNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}
	
	@Test
	public void floatRollingBlockSumNxNInternalIsFasterThanIntRollingBlockSumNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		SumFilter filter = new SumFilter();

		int iter = 50;

		ArrayList<int[]> dataSet = new ArrayList<int[]>(iter);
		for (int i = iter; i-- > 0;)
		{
			dataSet.add(intCreateData(primes[0], primes[0]));
		}

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.rollingBlockSumNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(iter);
					for (int[] data : dataSet)
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
					ArrayList<int[]> dataSet2 = new ArrayList<int[]>(iter);
					for (int[] data : dataSet)
						dataSet2.add(intClone(data));

					long time = System.nanoTime();
					for (int[] data : dataSet2)
						filter.rollingBlockSumNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out
								.printf("int rollingBlockSumNxNInternal [%dx%d] @ %d : %d => float rollingBlockSumNxNInternal %d = %.2fx\n",
										width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf(
					"int rollingBlockSumNxNInternal %d : %d => float rollingBlockSumNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
					boxFastTotal < boxSlowTotal);
		}
		System.out.printf("int rollingBlockSumNxNInternal %d => float rollingBlockSumNxNInternal %d = %.2fx\n",
				slowTotal, fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	// TODO
	// int sum faster than float sum
	// Internal version vs complete version -> Is the speed hit significant?
}
