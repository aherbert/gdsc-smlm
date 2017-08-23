package gdsc.smlm.filters;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.junit.Assert;
import org.junit.Test;
import org.junit.internal.ArrayComparisonFailure;

import gdsc.core.utils.DoubleEquality;
import gdsc.core.utils.FloatEquality;
import gdsc.core.utils.Maths;
import gdsc.smlm.TestSettings;
import gnu.trove.list.array.TDoubleArrayList;

public class BlockMeanFilterTest
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
	int[] boxSizes = new int[] { 15, 9, 5, 3, 2, 1 };
	boolean[] checkInternal = new boolean[] { true, false };

	/**
	 * Do a simple and stupid mean filter.
	 *
	 * @param data
	 *            the data
	 * @param maxx
	 *            the maxx
	 * @param maxy
	 *            the maxy
	 * @param boxSize
	 *            the box size
	 */
	public static void mean(float[] data, int maxx, int maxy, float boxSize)
	{
		if (boxSize <= 0)
			return;

		int n = (int) Math.ceil(boxSize);
		int size = 2 * n + 1;
		float[] weight = new float[size];
		Arrays.fill(weight, 1);
		if (boxSize != n)
			weight[0] = weight[weight.length - 1] = boxSize - (n - 1);

		double norm = 0;
		for (int yy = 0; yy < size; yy++)
			for (int xx = 0; xx < size; xx++)
				norm += weight[yy] * weight[xx];
		norm = 1.0 / norm;

		float[] out = new float[data.length];

		for (int y = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++)
			{
				double sum = 0;
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
				out[y * maxx + x] = (float) (sum * norm);
			}
		}
		System.arraycopy(out, 0, data, 0, out.length);
	}

	/**
	 * Do a simple and stupid mean filter with weights.
	 *
	 * @param data
	 *            the data
	 * @param w
	 *            the weights
	 * @param maxx
	 *            the maxx
	 * @param maxy
	 *            the maxy
	 * @param boxSize
	 *            the box size
	 */
	public static void weightedMean(float[] data, float[] w, int maxx, int maxy, float boxSize)
	{
		if (boxSize <= 0)
			return;

		int n = (int) Math.ceil(boxSize);
		int size = 2 * n + 1;
		float[] weight = new float[size];
		Arrays.fill(weight, 1);
		if (boxSize != n)
			weight[0] = weight[weight.length - 1] = boxSize - (n - 1);

		float[] out = new float[data.length];

		for (int y = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++)
			{
				double sum = 0, norm = 0;
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
						double w2 = w[index] * weight[yy] * weight[xx];
						sum += data[index] * w2;
						norm += w2;
					}
				}
				out[y * maxx + x] = (float) (sum / norm);
			}
		}
		System.arraycopy(out, 0, data, 0, out.length);
	}

	private void floatArrayEquals(String message, float[] data1, float[] data2, int maxx, int maxy, float boxSize)
	{
		FloatEquality eq = new FloatEquality(2e-4f, 1e-10f);
		// Debug: show the images
		//gdsc.core.ij.Utils.display("data1", new ij.process.FloatProcessor(maxx, maxy, data1));
		//gdsc.core.ij.Utils.display("data2", new ij.process.FloatProcessor(maxx, maxy, data2));

		// Ignore the border
		int border = (int) Math.ceil(boxSize);
		for (int y = border; y < maxy - border - 1; y++)
		{
			int index = y * maxx + border;
			for (int x = border; x < maxx - border - 1; x++, index++)
			{
				if (!eq.almostEqualRelativeOrAbsolute(data1[index], data2[index]))
				{
					Assert.fail(String.format("%s [%d,%d] %f != %f  (%g)", message, x, y, data1[index], data2[index],
							FloatEquality.relativeError(data1[index], data2[index])));
				}
			}
		}
	}

	/**
	 * Used to test the filter methods calculate the correct result
	 */
	private abstract class BlockMeanDataFilter extends DataFilter
	{
		public BlockMeanDataFilter(String name, boolean isInterpolated)
		{
			super(name, isInterpolated);
		}

		BlockMeanFilter f = new BlockMeanFilter();

		@Override
		public void setWeights(float[] w, int width, int height)
		{
			f.setWeights(w, width, height);
		}
	}

	private void meanIsCorrect(float[] data, int width, int height, float boxSize, boolean internal, BlockMeanDataFilter filter)
			throws ArrayComparisonFailure
	{
		float[] data1 = data.clone();
		float[] data2 = data.clone();

		mean(data1, width, height, boxSize);
		if (internal)
		{
			filter.filterInternal(data2, width, height, boxSize);
			floatArrayEquals(String.format("Internal arrays do not match: [%dx%d] @ %.1f", width, height, boxSize),
					data1, data2, width, height, boxSize);
		}
		else
		{
			filter.filter(data2, width, height, boxSize);
			floatArrayEquals(String.format("Arrays do not match: [%dx%d] @ %.1f", width, height, boxSize), data1, data2,
					width, height, 0);
		}
	}

	private void weightedMeanIsCorrect(float[] data, float[] w, int width, int height, float boxSize, boolean internal,
			BlockMeanDataFilter filter) throws ArrayComparisonFailure
	{
		float[] data1 = data.clone();
		float[] data2 = data.clone();

		weightedMean(data1, w, width, height, boxSize);

		//// Check the weights do not alter the image mean
		//double u1 = Maths.sum(data) / data.length;
		//double u2 = Maths.sum(data1) / data.length;
		//System.out.printf("[%dx%d] @ %.1f : %g => %g  (%g)\n", width, height, boxSize, u1, u2,
		//		DoubleEquality.relativeError(u1, u2));

		if (internal)
		{
			filter.filterInternal(data2, width, height, boxSize);
			floatArrayEquals(String.format("Internal arrays do not match: [%dx%d] @ %.1f", width, height, boxSize),
					data1, data2, width, height, boxSize);
		}
		else
		{
			filter.filter(data2, width, height, boxSize);
			floatArrayEquals(String.format("Arrays do not match: [%dx%d] @ %.1f", width, height, boxSize), data1, data2,
					width, height, 0);
		}
	}

	private void checkIsCorrect(BlockMeanDataFilter filter)
	{
		rand = new gdsc.core.utils.Random(-30051976);
		ExponentialDistribution ed = new ExponentialDistribution(rand, 57,
				ExponentialDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);

		for (int width : primes)
			for (int height : primes)
			{
				float[] data = createData(width, height);

				filter.f.setWeights(null, 0, 0);
				for (float boxSize : boxSizes)
					for (boolean internal : checkInternal)
					{
						meanIsCorrect(data, width, height, boxSize, internal, filter);
						if (filter.isInterpolated)
						{
							meanIsCorrect(data, width, height, boxSize - 0.3f, internal, filter);
							meanIsCorrect(data, width, height, boxSize - 0.6f, internal, filter);
						}
					}

				// Uniform weights
				float[] w = new float[width * height];
				Arrays.fill(w, 1);
				filter.f.setWeights(w, width, height);
				for (float boxSize : boxSizes)
					for (boolean internal : checkInternal)
					{
						weightedMeanIsCorrect(data, w, width, height, boxSize, internal, filter);
						if (filter.isInterpolated)
						{
							weightedMeanIsCorrect(data, w, width, height, boxSize - 0.3f, internal, filter);
							weightedMeanIsCorrect(data, w, width, height, boxSize - 0.6f, internal, filter);
						}
					}

				// Weights simulating the variance of sCMOS pixels
				for (int i = 0; i < w.length; i++)
				{
					w[i] = (float) (1.0 / Math.max(0.01, ed.sample()));
				}

				filter.f.setWeights(w, width, height);
				for (float boxSize : boxSizes)
					for (boolean internal : checkInternal)
					{
						weightedMeanIsCorrect(data, w, width, height, boxSize, internal, filter);
						if (filter.isInterpolated)
						{
							weightedMeanIsCorrect(data, w, width, height, boxSize - 0.3f, internal, filter);
							weightedMeanIsCorrect(data, w, width, height, boxSize - 0.6f, internal, filter);
						}
					}
			}
	}

	@Test
	public void blockFilterIsCorrect()
	{
		BlockMeanDataFilter filter = new BlockMeanDataFilter("block", true)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.blockFilter(data, width, height, boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.blockFilterInternal(data, width, height, boxSize);
			}
		};
		checkIsCorrect(filter);
	}

	@Test
	public void stripedBlockFilterIsCorrect()
	{
		BlockMeanDataFilter filter = new BlockMeanDataFilter("stripedBlock", true)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter(data, width, height, boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterInternal(data, width, height, boxSize);
			}
		};
		checkIsCorrect(filter);
	}

	@Test
	public void rollingBlockFilterIsCorrect()
	{
		BlockMeanDataFilter filter = new BlockMeanDataFilter("rollingBlock", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.rollingBlockFilter(data, width, height, (int) boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.rollingBlockFilterInternal(data, width, height, (int) boxSize);
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

	private void speedTest(BlockMeanDataFilter fast, BlockMeanDataFilter slow)
	{
		speedTest(fast, slow, boxSizes);
	}

	private void speedTest(BlockMeanDataFilter fast, BlockMeanDataFilter slow, int[] testBoxSizes)
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

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

	private void speedTestInternal(BlockMeanDataFilter fast, BlockMeanDataFilter slow)
	{
		speedTestInternal(fast, slow, boxSizes);
	}

	private void speedTestInternal(BlockMeanDataFilter fast, BlockMeanDataFilter slow, int[] testBoxSizes)
	{
		org.junit.Assume.assumeTrue(TestSettings.RUN_SPEED_TESTS);

		rand = new gdsc.core.utils.Random(-300519);

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
		BlockMeanDataFilter slow = new BlockMeanDataFilter("block", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.blockFilter(data, width, height, (int) boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.blockFilterInternal(data, width, height, (int) boxSize);
			}
		};
		BlockMeanDataFilter fast = new BlockMeanDataFilter("stripedBlock", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter(data, width, height, (int) boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterInternal(data, width, height, (int) boxSize);
			}
		};

		speedTest(fast, slow);
		speedTestInternal(fast, slow);
	}

	@Test
	public void interpolatedStripedBlockIsFasterThanBlock()
	{
		BlockMeanDataFilter slow = new BlockMeanDataFilter("block", true)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.blockFilter(data, width, height, boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.blockFilterInternal(data, width, height, boxSize);
			}
		};
		BlockMeanDataFilter fast = new BlockMeanDataFilter("stripedBlock", true)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter(data, width, height, boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterInternal(data, width, height, boxSize);
			}
		};

		speedTest(fast, slow);
		speedTestInternal(fast, slow);
	}

	@Test
	public void rollingBlockIsFasterThanBlock()
	{
		BlockMeanDataFilter slow = new BlockMeanDataFilter("block", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.blockFilter(data, width, height, (int) boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.blockFilterInternal(data, width, height, (int) boxSize);
			}
		};
		BlockMeanDataFilter fast = new BlockMeanDataFilter("rollingBlock", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.rollingBlockFilter(data, width, height, (int) boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.rollingBlockFilterInternal(data, width, height, (int) boxSize);
			}
		};

		speedTest(fast, slow);
		speedTestInternal(fast, slow);
	}

	@Test
	public void rollingBlockIsFasterThanStripedBlock()
	{
		BlockMeanDataFilter slow = new BlockMeanDataFilter("stripedBlock", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter(data, width, height, (int) boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterInternal(data, width, height, (int) boxSize);
			}
		};
		BlockMeanDataFilter fast = new BlockMeanDataFilter("rollingBlock", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.rollingBlockFilter(data, width, height, (int) boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.rollingBlockFilterInternal(data, width, height, (int) boxSize);
			}
		};

		speedTest(fast, slow);
		speedTestInternal(fast, slow);
	}

	@Test
	public void stripedBlock3x3IsFasterThanStripedBlockNxN()
	{
		BlockMeanDataFilter slow = new BlockMeanDataFilter("stripedBlockNxN", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxN(data, width, height, (int) boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxNInternal(data, width, height, (int) boxSize);
			}
		};
		BlockMeanDataFilter fast = new BlockMeanDataFilter("stripedBlock3x3", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter3x3(data, width, height);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter3x3Internal(data, width, height);
			}
		};

		int[] testBoxSizes = new int[] { 1 };
		speedTest(fast, slow, testBoxSizes);
		speedTestInternal(fast, slow, testBoxSizes);
	}

	@Test
	public void interpolatedStripedBlock3x3IsFasterThanStripedBlockNxN()
	{
		BlockMeanDataFilter slow = new BlockMeanDataFilter("stripedBlockNxN", true)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxN(data, width, height, boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxNInternal(data, width, height, boxSize);
			}
		};
		BlockMeanDataFilter fast = new BlockMeanDataFilter("stripedBlock3x3", true)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter3x3(data, width, height, boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter3x3Internal(data, width, height, boxSize);
			}
		};

		int[] testBoxSizes = new int[] { 1 };
		speedTest(fast, slow, testBoxSizes);
		speedTestInternal(fast, slow, testBoxSizes);
	}

	@Test
	public void stripedBlock5x5IsFasterThanStripedBlockNxN()
	{
		BlockMeanDataFilter slow = new BlockMeanDataFilter("stripedBlockNxN", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxN(data, width, height, (int) boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxNInternal(data, width, height, (int) boxSize);
			}
		};
		BlockMeanDataFilter fast = new BlockMeanDataFilter("stripedBlock5x5", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter5x5(data, width, height);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter5x5Internal(data, width, height);
			}
		};

		int[] testBoxSizes = new int[] { 2 };
		speedTest(fast, slow, testBoxSizes);
		speedTestInternal(fast, slow, testBoxSizes);
	}

	@Test
	public void interpolatedStripedBlock5x5IsFasterThanStripedBlockNxN()
	{
		BlockMeanDataFilter slow = new BlockMeanDataFilter("stripedBlockNxN", true)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxN(data, width, height, boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxNInternal(data, width, height, boxSize);
			}
		};
		BlockMeanDataFilter fast = new BlockMeanDataFilter("stripedBlock5x5", true)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter5x5(data, width, height, boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter5x5Internal(data, width, height, boxSize);
			}
		};

		int[] testBoxSizes = new int[] { 2 };
		speedTest(fast, slow, testBoxSizes);
		speedTestInternal(fast, slow, testBoxSizes);
	}

	@Test
	public void stripedBlock7x7IsFasterThanStripedBlockNxN()
	{
		BlockMeanDataFilter slow = new BlockMeanDataFilter("stripedBlockNxN", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxN(data, width, height, (int) boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxNInternal(data, width, height, (int) boxSize);
			}
		};
		BlockMeanDataFilter fast = new BlockMeanDataFilter("stripedBlock7x7", false)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter7x7(data, width, height);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter7x7Internal(data, width, height);
			}
		};

		int[] testBoxSizes = new int[] { 3 };
		speedTest(fast, slow, testBoxSizes);
		speedTestInternal(fast, slow, testBoxSizes);
	}

	@Test
	public void interpolatedStripedBlock7x7IsFasterThanStripedBlockNxN()
	{
		BlockMeanDataFilter slow = new BlockMeanDataFilter("stripedBlockNxN", true)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxN(data, width, height, boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxNInternal(data, width, height, boxSize);
			}
		};
		BlockMeanDataFilter fast = new BlockMeanDataFilter("stripedBlock7x7", true)
		{
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter7x7(data, width, height, boxSize);
			}

			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter7x7Internal(data, width, height, boxSize);
			}
		};

		int[] testBoxSizes = new int[] { 3 };
		speedTest(fast, slow, testBoxSizes);
		speedTestInternal(fast, slow, testBoxSizes);
	}
}
