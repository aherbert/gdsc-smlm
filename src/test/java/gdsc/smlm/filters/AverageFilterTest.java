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
	int[] boxSizes = new int[] { 15, 9, 5, 3, 2, 1 };
	boolean[] checkInternal = new boolean[] { true, false };

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

	private void floatArrayEquals(String message, float[] data1, float[] data2, int maxx, int maxy, float boxSize)
	{
		FloatEquality eq = new FloatEquality(1e-5f, 1e-10f);
		// Debug: show the images
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
					Assert.assertTrue(String.format("%s [%d,%d] %f != %f", message, x, y, data1[index], data2[index]),
							false);
				}
			}
		}
	}

	/**
	 * Used to test the filter methods calculate the correct result
	 */
	private abstract class DataFilter
	{
		final String name;
		final boolean isInterpolated;

		public DataFilter(String name, boolean isInterpolated)
		{
			this.name = name;
			this.isInterpolated = isInterpolated;
		}

		AverageFilter f = new AverageFilter();

		public abstract void filter(float[] data, int width, int height, float boxSize);

		public abstract void filterInternal(float[] data, int width, int height, float boxSize);
	}

	private void averageIsCorrect(int width, int height, float boxSize, boolean internal, DataFilter filter)
			throws ArrayComparisonFailure
	{
		rand = new gdsc.smlm.utils.Random(-30051976);
		float[] data1 = createData(width, height);
		float[] data2 = data1.clone();

		AverageFilterTest.average(data1, width, height, boxSize);
		if (internal)
		{
			filter.filterInternal(data2, width, height, boxSize);
			floatArrayEquals(String.format("Internal arrays do not match: [%dx%d] @ %.1f", width, height, boxSize),
					data1, data2, width, height, boxSize);
		}
		else
		{
			filter.filter(data2, width, height, boxSize);
			floatArrayEquals(String.format("Arrays do not match: [%dx%d] @ %.1f", width, height, boxSize), data1,
					data2, width, height, 0);
		}
	}

	private void checkIsCorrect(DataFilter filter)
	{
		for (int width : primes)
			for (int height : primes)
				for (float boxSize : boxSizes)
					for (boolean internal : checkInternal)
					{
						averageIsCorrect(width, height, boxSize, internal, filter);
						if (filter.isInterpolated)
						{
							averageIsCorrect(width, height, boxSize - 0.3f, internal, filter);
							averageIsCorrect(width, height, boxSize - 0.6f, internal, filter);
						}
					}
	}

	@Test
	public void blockAverageIsCorrect()
	{
		DataFilter filter = new DataFilter("block", true)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.blockAverage(data, width, height, boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.blockAverageInternal(data, width, height, boxSize);
			}
		};
		checkIsCorrect(filter);
	}

	@Test
	public void stripedBlockAverageIsCorrect()
	{
		DataFilter filter = new DataFilter("stripedBlock", true)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverage(data, width, height, boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverageInternal(data, width, height, boxSize);
			}
		};
		checkIsCorrect(filter);
	}

	@Test
	public void rollingBlockAverageIsCorrect()
	{
		DataFilter filter = new DataFilter("rollingBlock", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.rollingBlockAverage(data, width, height, (int) boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.rollingBlockAverageInternal(data, width, height, (int) boxSize);
			}
		};
		checkIsCorrect(filter);
	}

	private double speedUpFactor(long slowTotal, long fastTotal)
	{
		return (1.0 * slowTotal) / fastTotal;
	}

	private float[] floatClone(float[] data1)
	{
		float[] data2 = Arrays.copyOf(data1, data1.length);
		return data2;
	}

	private float[] createData(int width, int height)
	{
		float[] data = new float[width * height];
		for (int i = data.length; i-- > 0;)
			data[i] = i;

		rand.shuffle(data);

		return data;
	}

	private ArrayList<float[]> floatCreateSpeedData(int iter)
	{
		ArrayList<float[]> dataSet = new ArrayList<float[]>(iter);
		for (int i = iter; i-- > 0;)
		{
			dataSet.add(createData(primes[0], primes[0]));
		}
		return dataSet;
	}

	static ArrayList<float[]> dataSet = null;

	private ArrayList<float[]> getSpeedData(int iter)
	{
		if (dataSet == null || dataSet.size() < iter)
		{
			dataSet = floatCreateSpeedData(iter);
		}
		ArrayList<float[]> dataSet2 = new ArrayList<float[]>(iter);
		for (int i = 0; i < iter; i++)
			dataSet2.add(dataSet.get(i).clone());
		return dataSet2;
	}

	private void speedTest(DataFilter fast, DataFilter slow)
	{
		speedTest(fast, slow, boxSizes);
	}

	private void speedTest(DataFilter fast, DataFilter slow, int[] testBoxSizes)
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		ArrayList<float[]> dataSet = getSpeedData(ITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		float[] boxSizes = new float[testBoxSizes.length];
		float offset = (fast.isInterpolated && slow.isInterpolated) ? 0.3f : 0;
		for (int i = 0; i < boxSizes.length; i++)
			boxSizes[i] = testBoxSizes[i] - offset;

		// Initialise
		for (float boxSize : boxSizes)
		{
			fast.filter(dataSet.get(0).clone(), primes[0], primes[0], boxSize);
			slow.filter(dataSet.get(0).clone(), primes[0], primes[0], boxSize);
		}

		for (float boxSize : boxSizes)
		{
			int iter = (boxSize == 1) ? ITER3 : ITER;
			for (int width : primes)
				for (int height : primes)
				{
					dataSet = getSpeedData(iter);

					long time = System.nanoTime();
					for (float[] data : dataSet)
						fast.filter(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}
		}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (float boxSize : boxSizes)
		{
			int iter = (boxSize == 1) ? ITER3 : ITER;
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					dataSet = getSpeedData(iter);

					long time = System.nanoTime();
					for (float[] data : dataSet)
						slow.filter(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf("%s [%dx%d] @ %.1f : %d => %s %d = %.2fx\n", fast.name, width, height,
								boxSize, time, slow.name, fastTime, speedUpFactor(time, fastTime));
				}
			//if (debug)
			System.out.printf("%s %.1f : %d => %s %d = %.2fx\n", fast.name, boxSize, boxSlowTotal, slow.name,
					boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(
						String.format("Not faster: Block %.1f : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("%s %d => %s %d = %.2fx\n", fast.name, slowTotal, slow.name, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	private void speedTestInternal(DataFilter fast, DataFilter slow)
	{
		speedTestInternal(fast, slow, boxSizes);
	}

	private void speedTestInternal(DataFilter fast, DataFilter slow, int[] testBoxSizes)
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.smlm.utils.Random(-300519);

		ArrayList<float[]> dataSet = getSpeedData(InternalITER3);

		ArrayList<Long> fastTimes = new ArrayList<Long>();

		float[] boxSizes = new float[testBoxSizes.length];
		float offset = (fast.isInterpolated && slow.isInterpolated) ? 0.3f : 0;
		for (int i = 0; i < boxSizes.length; i++)
			boxSizes[i] = testBoxSizes[i] - offset;

		// Initialise
		for (float boxSize : boxSizes)
		{
			fast.filterInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSize);
			slow.filterInternal(floatClone(dataSet.get(0)), primes[0], primes[0], boxSize);
		}

		for (float boxSize : boxSizes)
		{
			int iter = (boxSize == 1) ? InternalITER3 : InternalITER;
			for (int width : primes)
				for (int height : primes)
				{
					dataSet = getSpeedData(iter);

					long time = System.nanoTime();
					for (float[] data : dataSet)
						fast.filterInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;
					fastTimes.add(time);
				}
		}

		long slowTotal = 0, fastTotal = 0;
		int index = 0;
		for (float boxSize : boxSizes)
		{
			int iter = (boxSize == 1) ? InternalITER3 : InternalITER;
			long boxSlowTotal = 0, boxFastTotal = 0;
			for (int width : primes)
				for (int height : primes)
				{
					dataSet = getSpeedData(iter);

					long time = System.nanoTime();
					for (float[] data : dataSet)
						slow.filterInternal(data, width, height, boxSize);
					time = System.nanoTime() - time;

					long fastTime = fastTimes.get(index++);
					slowTotal += time;
					fastTotal += fastTime;
					boxSlowTotal += time;
					boxFastTotal += fastTime;
					if (debug)
						System.out.printf("Internal %s [%dx%d] @ %.1f : %d => %s %d = %.2fx\n", fast.name, width,
								height, boxSize, time, slow.name, fastTime, speedUpFactor(time, fastTime));
				}
			//if (debug)
			System.out.printf("Internal %s %.1f : %d => %s %d = %.2fx\n", fast.name, boxSize, boxSlowTotal, slow.name,
					boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			if (TestSettings.ASSERT_SPEED_TESTS)
				Assert.assertTrue(
						String.format("Not faster: Block %.1f : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
						boxFastTotal < boxSlowTotal);
		}
		System.out.printf("Internal %s %d => %s %d = %.2fx\n", fast.name, slowTotal, slow.name, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
		if (TestSettings.ASSERT_SPEED_TESTS)
			Assert.assertTrue(String.format("Not faster: %d > %d", fastTotal, slowTotal), fastTotal < slowTotal);
	}

	@Test
	public void stripedBlockIsFasterThanBlock()
	{
		DataFilter slow = new DataFilter("block", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.blockAverage(data, width, height, (int) boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.blockAverageInternal(data, width, height, (int) boxSize);
			}
		};
		DataFilter fast = new DataFilter("stripedBlock", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverage(data, width, height, (int) boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverageInternal(data, width, height, (int) boxSize);
			}
		};

		speedTest(fast, slow);
		speedTestInternal(fast, slow);
	}

	@Test
	public void interpolatedStripedBlockIsFasterThanBlock()
	{
		DataFilter slow = new DataFilter("block", true)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.blockAverage(data, width, height, boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.blockAverageInternal(data, width, height, boxSize);
			}
		};
		DataFilter fast = new DataFilter("stripedBlock", true)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverage(data, width, height, boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverageInternal(data, width, height, boxSize);
			}
		};

		speedTest(fast, slow);
		speedTestInternal(fast, slow);
	}

	@Test
	public void rollingBlockIsFasterThanBlock()
	{
		DataFilter slow = new DataFilter("block", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.blockAverage(data, width, height, (int) boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.blockAverageInternal(data, width, height, (int) boxSize);
			}
		};
		DataFilter fast = new DataFilter("rollingBlock", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.rollingBlockAverage(data, width, height, (int) boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.rollingBlockAverageInternal(data, width, height, (int) boxSize);
			}
		};

		speedTest(fast, slow);
		speedTestInternal(fast, slow);
	}

	@Test
	public void rollingBlockIsFasterThanStripedBlock()
	{
		DataFilter slow = new DataFilter("stripedBlock", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverage(data, width, height, (int) boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverageInternal(data, width, height, (int) boxSize);
			}
		};
		DataFilter fast = new DataFilter("rollingBlock", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.rollingBlockAverage(data, width, height, (int) boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.rollingBlockAverageInternal(data, width, height, (int) boxSize);
			}
		};

		speedTest(fast, slow);
		speedTestInternal(fast, slow);
	}

	@Test
	public void stripedBlock3x3IsFasterThanStripedBlockNxN()
	{
		DataFilter slow = new DataFilter("stripedBlockNxN", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverageNxN(data, width, height, (int) boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverageNxNInternal(data, width, height, (int) boxSize);
			}
		};
		DataFilter fast = new DataFilter("stripedBlock3x3", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverage3x3(data, width, height);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverage3x3Internal(data, width, height);
			}
		};

		int[] testBoxSizes = new int[] { 1 };
		speedTest(fast, slow, testBoxSizes);
		speedTestInternal(fast, slow, testBoxSizes);
	}

	@Test
	public void interpolatedStripedBlock3x3IsFasterThanStripedBlockNxN()
	{
		DataFilter slow = new DataFilter("stripedBlockNxN", true)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverageNxN(data, width, height, boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverageNxNInternal(data, width, height, boxSize);
			}
		};
		DataFilter fast = new DataFilter("stripedBlock3x3", true)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverage3x3(data, width, height, boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverage3x3Internal(data, width, height, boxSize);
			}
		};

		int[] testBoxSizes = new int[] { 1 };
		speedTest(fast, slow, testBoxSizes);
		speedTestInternal(fast, slow, testBoxSizes);
	}

	@Test
	public void stripedBlock5x5IsFasterThanStripedBlockNxN()
	{
		DataFilter slow = new DataFilter("stripedBlockNxN", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverageNxN(data, width, height, (int) boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverageNxNInternal(data, width, height, (int) boxSize);
			}
		};
		DataFilter fast = new DataFilter("stripedBlock5x5", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverage5x5(data, width, height);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverage5x5Internal(data, width, height);
			}
		};

		int[] testBoxSizes = new int[] { 2 };
		speedTest(fast, slow, testBoxSizes);
		speedTestInternal(fast, slow, testBoxSizes);
	}

	@Test
	public void interpolatedStripedBlock5x5IsFasterThanStripedBlockNxN()
	{
		DataFilter slow = new DataFilter("stripedBlockNxN", true)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverageNxN(data, width, height, boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverageNxNInternal(data, width, height, boxSize);
			}
		};
		DataFilter fast = new DataFilter("stripedBlock5x5", true)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverage5x5(data, width, height, boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockAverage5x5Internal(data, width, height, boxSize);
			}
		};

		int[] testBoxSizes = new int[] { 2 };
		speedTest(fast, slow, testBoxSizes);
		speedTestInternal(fast, slow, testBoxSizes);
	}
}
