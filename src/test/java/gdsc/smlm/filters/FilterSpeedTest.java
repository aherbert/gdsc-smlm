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

import org.junit.Test;

import gdsc.test.TestSettings;

@SuppressWarnings({ "javadoc" })
public class FilterSpeedTest extends AbstractFilterTest
{

	@Test
	public void floatRollingBlockSumNxNInternalIsFasterThanRollingBlockMeanNxNInternal()
	{
		TestSettings.assumeSpeedTest();

		BlockSumFilter filter = new BlockSumFilter();
		BlockMeanFilter filter2 = new BlockMeanFilter();

		int iter = 50;
		ArrayList<float[]> dataSet = getSpeedData(iter);

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
			TestSettings.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"float rollingBlockMeanNxNInternal %d : %d => rollingBlockSumNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		TestSettings.logSpeedTestResult(fastTotal < slowTotal,
				"float rollingBlockMeanNxNInternal %d => rollingBlockSumNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatRollingBlockMeanNxNInternalIsFasterThanBlockMedianNxNInternal()
	{
		TestSettings.assumeSpeedTest();

		BlockMeanFilter filter1 = new BlockMeanFilter();
		MedianFilter filter2 = new MedianFilter();

		int iter = 10;
		ArrayList<float[]> dataSet = getSpeedData(iter);

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
			TestSettings.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"float blockMedianNxNInternal %d : %d => rollingBlockMeanNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		TestSettings.logSpeedTestResult(fastTotal < slowTotal,
				"float blockMedianNxNInternal %d => rollingBlockMeanNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatRollingBlockMeanNxNInternalIsFasterThanRollingMedianNxNInternal()
	{
		TestSettings.assumeSpeedTest();

		BlockMeanFilter filter1 = new BlockMeanFilter();
		MedianFilter filter2 = new MedianFilter();

		int iter = 10;
		ArrayList<float[]> dataSet = getSpeedData(iter);

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
			TestSettings.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"float rollingMedianNxNInternal %d : %d => rollingBlockMeanNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		TestSettings.logSpeedTestResult(fastTotal < slowTotal,
				"float rollingMedianNxNInternal %d => rollingBlockMeanNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatRollingBlockMeanNxNInternalIsFasterThanGaussianNxNInternal()
	{
		TestSettings.assumeSpeedTest();

		BlockMeanFilter filter1 = new BlockMeanFilter();
		GaussianFilter filter2 = new GaussianFilter();

		int iter = 10;
		ArrayList<float[]> dataSet = getSpeedData(iter);

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
			TestSettings.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"float convolveInternal %d : %d => rollingBlockMeanNxNInternal %d = %.2fx\n", boxSize, boxSlowTotal,
					boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		TestSettings.logSpeedTestResult(fastTotal < slowTotal,
				"float convolveInternal %d => rollingBlockMeanNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatRollingBlockMeanNxNInternalIsFasterThanAreaFilterNxNInternal()
	{
		TestSettings.assumeSpeedTest();

		BlockMeanFilter filter1 = new BlockMeanFilter();
		AreaAverageFilter filter2 = new AreaAverageFilter();

		int iter = 10;
		ArrayList<float[]> dataSet = getSpeedData(iter);

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
			TestSettings.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"float areaFilterInternal %d : %d => rollingBlockMeanNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		TestSettings.logSpeedTestResult(fastTotal < slowTotal,
				"float areaFilterInternal %d => rollingBlockMeanNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@Test
	public void floatStripedBlockMeanNxNInternalIsFasterThanAreaFilterNxNInternal()
	{
		TestSettings.assumeSpeedTest();

		BlockMeanFilter filter1 = new BlockMeanFilter();
		AreaAverageFilter filter2 = new AreaAverageFilter();

		int iter = 10;
		ArrayList<float[]> dataSet = getSpeedData(iter);

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
			TestSettings.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"float areaFilterInternal %d : %d => stripedBlockMeanNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
			//			Assert.assertTrue(String.format("Not faster: Block %d : %d > %d", boxSize, boxFastTotal, boxSlowTotal),
			//					boxFastTotal < boxSlowTotal);
		}
		TestSettings.logSpeedTestResult(fastTotal < slowTotal,
				"float areaFilterInternal %d => stripedBlockMeanNxNInternal %d = %.2fx\n", slowTotal, fastTotal,
				speedUpFactor(slowTotal, fastTotal));
	}

	@SuppressWarnings("deprecation")
	@Test
	public void floatRollingBlockSumNxNInternalIsFasterThanIntRollingBlockSumNxNInternal()
	{
		TestSettings.assumeSpeedTest();

		SumFilter filter = new SumFilter();
		BlockSumFilter filter2 = new BlockSumFilter();

		int iter = 50;
		ArrayList<int[]> dataSet = getIntSpeedData(iter);

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
			TestSettings.logSpeedTestStageResult(boxFastTotal < boxSlowTotal,
					"int rollingBlockSumNxNInternal %d : %d => float rollingBlockSumNxNInternal %d = %.2fx\n", boxSize,
					boxSlowTotal, boxFastTotal, speedUpFactor(boxSlowTotal, boxFastTotal));
		}
		TestSettings.logSpeedTestResult(fastTotal < slowTotal,
				"int rollingBlockSumNxNInternal %d => float rollingBlockSumNxNInternal %d = %.2fx\n", slowTotal,
				fastTotal, speedUpFactor(slowTotal, fastTotal));
	}

	// TODO
	// int sum faster than float sum
	// Internal version vs complete version -> Is the speed hit significant?
}
