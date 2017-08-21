package gdsc.smlm.filters;

import gdsc.smlm.TestSettings;
import gdsc.smlm.filters.BlockMeanFilter;
import gdsc.smlm.filters.BlockSumFilter;

import java.util.ArrayList;
import java.util.Arrays;

import org.junit.Assert;
import org.junit.Test;

public class FilterTest
{
	private gdsc.core.utils.Random rand;

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
	public void floatRollingBlockSumNxNInternalIsFasterThanRollingBlockMeanNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		BlockSumFilter filter = new BlockSumFilter();
		BlockMeanFilter filter2 = new BlockMeanFilter();

		int iter = 50;

		ArrayList<float[]> dataSet = new ArrayList<float[]>(iter);
		for (int i = iter; i-- > 0;)
		{
			dataSet.add(floatCreateData(primes[0], primes[0]));
		}

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockFilterNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter2.rollingBlockFilterNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(iter);
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.rollingBlockFilterNxNInternal(data, width, height, boxSize);
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
						filter2.rollingBlockFilterNxNInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"float rollingBlockMeanNxNInternal [%dx%d] @ %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float rollingBlockMeanNxNInternal %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float rollingBlockMeanNxNInternal %d => rollingBlockSumNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockMeanNxNInternalIsFasterThanBlockMedianNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		BlockMeanFilter filter1 = new BlockMeanFilter();
		MedianFilter filter2 = new MedianFilter();

		int iter = 10;

		ArrayList<float[]> dataSet = new ArrayList<float[]>(iter);
		for (int i = iter; i-- > 0;)
		{
			dataSet.add(floatCreateData(primes[0], primes[0]));
		}

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter1.rollingBlockFilterNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
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
						filter1.rollingBlockFilterNxNInternal(data, width, height, boxSize);
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
						System.out.printf(
								"float blockMedianNxNInternal [%dx%d] @ %d : %d => rollingBlockMeanNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float blockMedianNxNInternal %d : %d => rollingBlockMeanNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float blockMedianNxNInternal %d => rollingBlockMeanNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockMeanNxNInternalIsFasterThanRollingMedianNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		BlockMeanFilter filter1 = new BlockMeanFilter();
		MedianFilter filter2 = new MedianFilter();

		int iter = 10;

		ArrayList<float[]> dataSet = new ArrayList<float[]>(iter);
		for (int i = iter; i-- > 0;)
		{
			dataSet.add(floatCreateData(primes[0], primes[0]));
		}

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter1.rollingBlockFilterNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
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
						filter1.rollingBlockFilterNxNInternal(data, width, height, boxSize);
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
						System.out.printf(
								"float rollingMedianNxNInternal [%dx%d] @ %d : %d => rollingBlockMeanNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float rollingMedianNxNInternal %d : %d => rollingBlockMeanNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float rollingMedianNxNInternal %d => rollingBlockMeanNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockMeanNxNInternalIsFasterThanGaussianNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		BlockMeanFilter filter1 = new BlockMeanFilter();
		GaussianFilter filter2 = new GaussianFilter();

		int iter = 10;

		ArrayList<float[]> dataSet = new ArrayList<float[]>(iter);
		for (int i = iter; i-- > 0;)
		{
			dataSet.add(floatCreateData(primes[0], primes[0]));
		}

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter1.rollingBlockFilterNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
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
						filter1.rollingBlockFilterNxNInternal(data, width, height, boxSize);
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
						System.out.printf(
								"float convolveInternal [%dx%d] @ %d : %d => rollingBlockMeanNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float convolveInternal %d : %d => rollingBlockMeanNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float convolveInternal %d => rollingBlockMeanNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatRollingBlockMeanNxNInternalIsFasterThanAreaFilterNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		BlockMeanFilter filter1 = new BlockMeanFilter();
		AreaAverageFilter filter2 = new AreaAverageFilter();

		int iter = 10;

		ArrayList<float[]> dataSet = new ArrayList<float[]>(iter);
		for (int i = iter; i-- > 0;)
		{
			dataSet.add(floatCreateData(primes[0], primes[0]));
		}

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		//filter1.rollingBlockFilterNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		//filter2.areaFilterInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0] - 0.05);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(iter);
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					// Initialise
					for (float[] data : dataSet2)
						filter1.rollingBlockFilterNxNInternal(data.clone(), width, height, boxSize);
					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter1.rollingBlockFilterNxNInternal(data, width, height, boxSize);
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
						filter2.areaAverageUsingAveragesInternal(data.clone(), width, height, boxSize - 0.05);
					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter2.areaAverageUsingAveragesInternal(data, width, height, boxSize - 0.05);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"float areaFilterInternal [%dx%d] @ %d : %d => rollingBlockMeanNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float areaFilterInternal %d : %d => rollingBlockMeanNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float areaFilterInternal %d => rollingBlockMeanNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void floatStripedBlockMeanNxNInternalIsFasterThanAreaFilterNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		BlockMeanFilter filter1 = new BlockMeanFilter();
		AreaAverageFilter filter2 = new AreaAverageFilter();

		int iter = 10;

		ArrayList<float[]> dataSet = new ArrayList<float[]>(iter);
		for (int i = iter; i-- > 0;)
		{
			dataSet.add(floatCreateData(primes[0], primes[0]));
		}

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		//filter1.rollingBlockFilterNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		//filter2.areaFilterInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0] - 0.05);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(iter);
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

					// Initialise
					float w = (float) (boxSize - 0.05);
					for (float[] data : dataSet2)
						filter1.stripedBlockFilterNxNInternal(data.clone(), width, height, w);
					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter1.stripedBlockFilterNxNInternal(data, width, height, w);
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
						filter2.areaAverageUsingAveragesInternal(data.clone(), width, height, boxSize - 0.05);
					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter2.areaAverageUsingAveragesInternal(data, width, height, boxSize - 0.05);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf(
								"float areaFilterInternal [%dx%d] @ %d : %d => stripedBlockMeanNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("float areaFilterInternal %d : %d => stripedBlockMeanNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		System.out.printf("float areaFilterInternal %d => stripedBlockMeanNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@SuppressWarnings("deprecation")
	@Test
	public void floatRollingBlockSumNxNInternalIsFasterThanIntRollingBlockSumNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

		SumFilter filter = new SumFilter();
		BlockSumFilter filter2 = new BlockSumFilter();

		int iter = 50;

		ArrayList<int[]> dataSet = new ArrayList<int[]>(iter);
		for (int i = iter; i-- > 0;)
		{
			dataSet.add(intCreateData(primes[0], primes[0]));
		}

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.rollingBlockSumNxNInternal(intClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter2.rollingBlockFilterNxNInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (int boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(iter);
					for (int[] data : dataSet)
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter2.rollingBlockFilterNxNInternal(data, width, height, boxSize);
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
						System.out.printf(
								"int rollingBlockSumNxNInternal [%dx%d] @ %d : %d => float rollingBlockSumNxNInternal %d = %.2fx\n",
								width, height, boxSize, time, fastTime, speedUpFactor(time, fastTime));
					//Assert.assertTrue(String.format("Not faster: [%dx%d] @ %d : %d > %d", width, height, boxSize,
					//		blockTime, time), blockTime < time);
				}
			//if (debug)
			System.out.printf("int rollingBlockSumNxNInternal %d : %d => float rollingBlockSumNxNInternal %d = %.2fx\n",
					boxSize, boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
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
