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
import java.util.Arrays;

import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Test;
import org.junit.internal.ArrayComparisonFailure;

import gdsc.core.utils.FloatEquality;
import gdsc.core.utils.Maths;
import gdsc.test.TestSettings;
import gdsc.test.TestSettings.LogLevel;
import gdsc.test.TestSettings.TestComplexity;

public class BlockSumFilterTest extends AbstractFilterTest
{

	private int InternalITER3 = 500;
	private int InternalITER = 50;
	private int ITER3 = 200;
	private int ITER = 20;

	/**
	 * Do a simple and stupid sum filter.
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
	public static void sum(float[] data, int maxx, int maxy, float boxSize)
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

		int[] oy = new int[size];
		int[] ox = new int[size];

		for (int y = 0; y < maxy; y++)
		{
			// Cache offset
			for (int yy = 0; yy < size; yy++)
			{
				int yyy = y + yy - n;
				if (yyy < 0)
					yyy = 0;
				else if (yyy >= maxy)
					yyy = maxy - 1;
				oy[yy] = yyy * maxx;
			}

			for (int x = 0; x < maxx; x++)
			{
				// Cache offset
				for (int xx = 0; xx < size; xx++)
				{
					int xxx = x + xx - n;
					if (xxx < 0)
						xxx = 0;
					else if (xxx >= maxx)
						xxx = maxx - 1;
					ox[xx] = xxx;
				}

				double sum = 0;
				for (int yy = 0; yy < size; yy++)
				{
					//int yyy = y + yy - n;
					//if (yyy < 0)
					//	yyy = 0;
					//else if (yyy >= maxy)
					//	yyy = maxy - 1;

					final int index = oy[yy];
					final float wy = weight[yy];
					for (int xx = 0; xx < size; xx++)
					{
						//int xxx = x + xx - n;
						//if (xxx < 0)
						//	xxx = 0;
						//else if (xxx >= maxx)
						//	xxx = maxx - 1;
						//int index = yyy * maxx + xxx;
						//sum += data[index] * weight[yy] * weight[xx];

						sum += data[index + ox[xx]] * wy * weight[xx];
					}
				}
				out[y * maxx + x] = (float) (sum);
			}
		}
		System.arraycopy(out, 0, data, 0, out.length);
	}

	/**
	 * Do a simple and stupid sum filter with weights.
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
	public static void weightedSum(float[] data, float[] w, int maxx, int maxy, float boxSize)
	{
		if (boxSize <= 0)
			return;

		int n = (int) Math.ceil(boxSize);
		int size = 2 * n + 1;
		float[] weight = new float[size];
		Arrays.fill(weight, 1);
		if (boxSize != n)
			weight[0] = weight[weight.length - 1] = boxSize - (n - 1);
		double area = Maths.pow2(2 * boxSize + 1);

		float[] out = new float[data.length];

		for (int y = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++)
			{
				double sum = 0, sumw = 0;
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
						sumw += w2;
					}
				}
				// The sum should not be effected by the weights.
				out[y * maxx + x] = (float) (sum / (sumw / area));
			}
		}
		System.arraycopy(out, 0, data, 0, out.length);
	}

	/**
	 * Used to test the filter methods calculate the correct result
	 */
	private abstract class BlockSumDataFilter extends DataFilter
	{
		public BlockSumDataFilter(String name, boolean isInterpolated)
		{
			super(name, isInterpolated);
		}

		BlockSumFilter f = new BlockSumFilter();

		@Override
		public void setWeights(float[] w, int width, int height)
		{
			f.setWeights(w, width, height);
		}
	}

	private void sumIsCorrect(float[] data, int width, int height, float boxSize, boolean internal,
			BlockSumDataFilter filter) throws ArrayComparisonFailure
	{
		float[] data1 = data.clone();
		float[] data2 = data.clone();
		FloatEquality eq = new FloatEquality(2e-4f, 1e-10f);

		sum(data1, width, height, boxSize);
		if (internal)
		{
			filter.filterInternal(data2, width, height, boxSize);
			floatArrayEquals(eq, data1, data2, width, height, boxSize, "Internal arrays do not match: [%dx%d] @ %.1f",
					width, height, boxSize);
		}
		else
		{
			filter.filter(data2, width, height, boxSize);
			floatArrayEquals(eq, data1, data2, width, height, 0, "Arrays do not match: [%dx%d] @ %.1f", width, height,
					boxSize);
		}
	}

	private void weightedSumIsCorrect(float[] data, float[] w, int width, int height, float boxSize, boolean internal,
			BlockSumDataFilter filter) throws ArrayComparisonFailure
	{
		float[] data1 = data.clone();
		float[] data2 = data.clone();
		FloatEquality eq = new FloatEquality(2e-4f, 1e-10f);

		weightedSum(data1, w, width, height, boxSize);

		//// Check the weights do not alter the image sum
		//double u1 = gdsc.core.utils.Maths.sum(sum(data.clone(), width, height, boxSize));
		//double u2 = gdsc.core.utils.Maths.sum(data1);
		//System.out.printf("[%dx%d] @ %.1f : %g => %g  (%g)\n", width, height, boxSize, u1, u2,
		//		gdsc.core.utils.DoubleEquality.relativeError(u1, u2));

		if (internal)
		{
			filter.filterInternal(data2, width, height, boxSize);
			floatArrayEquals(eq, data1, data2, width, height, boxSize, "Internal arrays do not match: [%dx%d] @ %.1f",
					width, height, boxSize);
		}
		else
		{
			filter.filter(data2, width, height, boxSize);
			floatArrayEquals(eq, data1, data2, width, height, 0, "Arrays do not match: [%dx%d] @ %.1f", width, height,
					boxSize);
		}
	}

	private void checkIsCorrect(BlockSumDataFilter filter)
	{
		RandomGenerator rg = TestSettings.getRandomGenerator();
		ExponentialDistribution ed = new ExponentialDistribution(rg, 57,
				ExponentialDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);

		for (int width : primes)
			for (int height : primes)
			{
				float[] data = createData(rg, width, height);

				filter.f.setWeights(null, 0, 0);
				for (float boxSize : boxSizes)
					for (boolean internal : checkInternal)
					{
						sumIsCorrect(data, width, height, boxSize, internal, filter);
						if (filter.isInterpolated)
						{
							sumIsCorrect(data, width, height, boxSize - 0.3f, internal, filter);
							sumIsCorrect(data, width, height, boxSize - 0.6f, internal, filter);
						}
					}

				// Uniform weights
				float[] w = new float[width * height];
				Arrays.fill(w, 1);
				filter.f.setWeights(w, width, height);
				for (float boxSize : boxSizes)
					for (boolean internal : checkInternal)
					{
						weightedSumIsCorrect(data, w, width, height, boxSize, internal, filter);
						if (filter.isInterpolated)
						{
							weightedSumIsCorrect(data, w, width, height, boxSize - 0.3f, internal, filter);
							weightedSumIsCorrect(data, w, width, height, boxSize - 0.6f, internal, filter);
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
						weightedSumIsCorrect(data, w, width, height, boxSize, internal, filter);
						if (filter.isInterpolated)
						{
							weightedSumIsCorrect(data, w, width, height, boxSize - 0.3f, internal, filter);
							weightedSumIsCorrect(data, w, width, height, boxSize - 0.6f, internal, filter);
						}
					}
			}
	}

	@Test
	public void blockFilterIsCorrect()
	{
		BlockSumDataFilter filter = new BlockSumDataFilter("block", true)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.blockFilter(data, width, height, boxSize);
			}

			@Override
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
		BlockSumDataFilter filter = new BlockSumDataFilter("stripedBlock", true)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter(data, width, height, boxSize);
			}

			@Override
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
		BlockSumDataFilter filter = new BlockSumDataFilter("rollingBlock", false)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.rollingBlockFilter(data, width, height, (int) boxSize);
			}

			@Override
			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.rollingBlockFilterInternal(data, width, height, (int) boxSize);
			}
		};
		checkIsCorrect(filter);
	}

	private void speedTest(BlockSumDataFilter fast, BlockSumDataFilter slow)
	{
		speedTest(fast, slow, boxSizes);
	}

	private void speedTest(BlockSumDataFilter fast, BlockSumDataFilter slow, int[] testBoxSizes)
	{
		TestSettings.assumeSpeedTest();

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
			TestSettings.logSpeedTestResult(boxFastTotal < boxSlowTotal, "%s %.1f : %d => %s %d = %.2fx\n", fast.name,
					boxSize, boxSlowTotal, slow.name, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestSettings.logSpeedTestResult(fastTotal < slowTotal, "%s %d => %s %d = %.2fx\n", fast.name, slowTotal,
				slow.name, fastTotal, speedUpFactor(slowTotal, fastTotal));
	}

	private void speedTestInternal(BlockSumDataFilter fast, BlockSumDataFilter slow)
	{
		speedTestInternal(fast, slow, boxSizes);
	}

	private void speedTestInternal(BlockSumDataFilter fast, BlockSumDataFilter slow, int[] testBoxSizes)
	{
		TestSettings.assumeSpeedTest();

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
			TestSettings.logSpeedTestResult(boxFastTotal < boxSlowTotal, "Internal %s %.1f : %d => %s %d = %.2fx\n",
					fast.name, boxSize, boxSlowTotal, slow.name, boxFastTotal,
					speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestSettings.logSpeedTestResult(fastTotal < slowTotal, "Internal %s %d => %s %d = %.2fx\n", fast.name,
				slowTotal, slow.name, fastTotal, speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void stripedBlockIsFasterThanBlock()
	{
		BlockSumDataFilter slow = new BlockSumDataFilter("block", false)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.blockFilter(data, width, height, (int) boxSize);
			}

			@Override
			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.blockFilterInternal(data, width, height, (int) boxSize);
			}
		};
		BlockSumDataFilter fast = new BlockSumDataFilter("stripedBlock", false)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter(data, width, height, (int) boxSize);
			}

			@Override
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
		BlockSumDataFilter slow = new BlockSumDataFilter("block", true)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.blockFilter(data, width, height, boxSize);
			}

			@Override
			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.blockFilterInternal(data, width, height, boxSize);
			}
		};
		BlockSumDataFilter fast = new BlockSumDataFilter("stripedBlock", true)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter(data, width, height, boxSize);
			}

			@Override
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
		BlockSumDataFilter slow = new BlockSumDataFilter("block", false)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.blockFilter(data, width, height, (int) boxSize);
			}

			@Override
			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.blockFilterInternal(data, width, height, (int) boxSize);
			}
		};
		BlockSumDataFilter fast = new BlockSumDataFilter("rollingBlock", false)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.rollingBlockFilter(data, width, height, (int) boxSize);
			}

			@Override
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
		BlockSumDataFilter slow = new BlockSumDataFilter("stripedBlock", false)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter(data, width, height, (int) boxSize);
			}

			@Override
			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterInternal(data, width, height, (int) boxSize);
			}
		};
		BlockSumDataFilter fast = new BlockSumDataFilter("rollingBlock", false)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.rollingBlockFilter(data, width, height, (int) boxSize);
			}

			@Override
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
		BlockSumDataFilter slow = new BlockSumDataFilter("stripedBlockNxN", false)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxN(data, width, height, (int) boxSize);
			}

			@Override
			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxNInternal(data, width, height, (int) boxSize);
			}
		};
		BlockSumDataFilter fast = new BlockSumDataFilter("stripedBlock3x3", false)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter3x3(data, width, height);
			}

			@Override
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
		BlockSumDataFilter slow = new BlockSumDataFilter("stripedBlockNxN", true)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxN(data, width, height, boxSize);
			}

			@Override
			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxNInternal(data, width, height, boxSize);
			}
		};
		BlockSumDataFilter fast = new BlockSumDataFilter("stripedBlock3x3", true)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter3x3(data, width, height, boxSize);
			}

			@Override
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
		BlockSumDataFilter slow = new BlockSumDataFilter("stripedBlockNxN", false)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxN(data, width, height, (int) boxSize);
			}

			@Override
			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxNInternal(data, width, height, (int) boxSize);
			}
		};
		BlockSumDataFilter fast = new BlockSumDataFilter("stripedBlock5x5", false)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter5x5(data, width, height);
			}

			@Override
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
		BlockSumDataFilter slow = new BlockSumDataFilter("stripedBlockNxN", true)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxN(data, width, height, boxSize);
			}

			@Override
			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxNInternal(data, width, height, boxSize);
			}
		};
		BlockSumDataFilter fast = new BlockSumDataFilter("stripedBlock5x5", true)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter5x5(data, width, height, boxSize);
			}

			@Override
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
		BlockSumDataFilter slow = new BlockSumDataFilter("stripedBlockNxN", false)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxN(data, width, height, (int) boxSize);
			}

			@Override
			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxNInternal(data, width, height, (int) boxSize);
			}
		};
		BlockSumDataFilter fast = new BlockSumDataFilter("stripedBlock7x7", false)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter7x7(data, width, height);
			}

			@Override
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
		BlockSumDataFilter slow = new BlockSumDataFilter("stripedBlockNxN", true)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxN(data, width, height, boxSize);
			}

			@Override
			public void filterInternal(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilterNxNInternal(data, width, height, boxSize);
			}
		};
		BlockSumDataFilter fast = new BlockSumDataFilter("stripedBlock7x7", true)
		{
			@Override
			public void filter(float[] data, int width, int height, float boxSize)
			{
				f.stripedBlockFilter7x7(data, width, height, boxSize);
			}

			@Override
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
