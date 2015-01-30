package gdsc.smlm.filters;

import gdsc.smlm.TestSettings;

import java.util.ArrayList;
import java.util.Arrays;

import org.junit.Assert;
import org.junit.Test;
import org.junit.internal.ArrayComparisonFailure;

public class AreaAverageFilterTest
{
	private gdsc.smlm.utils.Random rand;

	private boolean debug = false;
	private int InternalITER = 50;

	// TODO - The test data should be representative of the final use case
	int[] primes = new int[] { 113, 97, 53, 29 };
	//int[] primes = new int[] { 1024 };
	float[] boxSizes = new float[] { 15.5f, 9.5f, 5.5f, 3.5f, 2.5f };


	private float[] floatCreateData(int width, int height)
	{
		float[] data = new float[width * height];
		for (int i = data.length; i-- > 0;)
			data[i] = rand.next();

		return data;
	}
	
	private float[] floatClone(float[] data1)
	{
		float[] data2 = Arrays.copyOf(data1, data1.length);
		return data2;
	}
	
	private void floatArrayEquals(String message, float[] data1, float[] data2, float boxSize)
	{
		Assert.assertArrayEquals(message, data1, data2, 1e-2f);
	}

	private double speedUpFactor(long slowTotal, long fastTotal)
	{
		return (1.0 * slowTotal) / fastTotal;
	}

	// TODO - This currently fails 
	//@Test
	public void floatAreaAverageNxNInternalAndAreaAverageUsingSumsNxNInternalReturnSameResult()
	{
		for (int width : primes)
			for (int height : primes)
				for (float boxSize : boxSizes)
					floatCompareAreaAverageNxNInternalAndAreaAverageUsingSumsNxNInternalReturnSameResult(width, height, boxSize);
	}

	private void floatCompareAreaAverageNxNInternalAndAreaAverageUsingSumsNxNInternalReturnSameResult(int width,
			int height, float boxSize) throws ArrayComparisonFailure
	{
		rand = new gdsc.smlm.utils.Random(-30051976);
		float[] data1 = floatCreateData(width, height);
		float[] data2 = floatClone(data1);

		AreaAverageFilter filter = new AreaAverageFilter();
		filter.areaAverageInternal(data1, width, height, boxSize);
		filter.areaAverageUsingSumsInternal(data2, width, height, boxSize);

		floatArrayEquals(String.format("Arrays do not match: [%dx%d] @ %.1f", width, height, boxSize), data1, data2,
				boxSize);
	}

	@Test
	public void floatAreaAverageUsingSumsNxNInternalIsFasterThanAreaAverageNxNInternal()
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		AreaAverageFilter filter = new AreaAverageFilter();

		ArrayList<float[]> dataSet = floatCreateSpeedData(InternalITER);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		// Initialise
		filter.areaAverageInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);
		filter.areaAverageUsingSumsInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSizes[0]);

		for (float boxSize : boxSizes)
			for (int width : primes)
				for (int height : primes)
				{
					ArrayList<float[]> dataSet2 = new ArrayList<float[]>(dataSet.size());
					for (float[] data : dataSet)
						dataSet2.add(floatClone(data));

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
						dataSet2.add(floatClone(data));

					long time = System.nanoTime();
					for (float[] data : dataSet2)
						filter.areaAverageInternal(data, width, height, boxSize);
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

	private ArrayList<float[]> floatCreateSpeedData(int iter)
	{
		ArrayList<float[]> dataSet = new ArrayList<float[]>(iter);
		for (int i = iter; i-- > 0;)
		{
			dataSet.add(floatCreateData(primes[0], primes[0]));
		}
		return dataSet;
	}
}