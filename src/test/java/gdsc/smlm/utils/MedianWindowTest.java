package gdsc.smlm.utils;

import java.util.Arrays;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public class MedianWindowTest
{
	Random random = new Random(30051977);
	int dataSize = 2000;
	int[] radii = new int[] { 0, 1, 2, 4, 8, 16 };
	double[] values = new double[] { 0, -1.1, 2.2 };
	int[] speedRadii = new int[] { 16, 32, 64 };
	int[] speedIncrement = new int[] { 1, 2, 4, 8, 16 };

	@Test
	public void testClassCanComputeActualMedian()
	{
		// Verify the internal median method using the Apache commons maths library

		double[] data = createRandomData(dataSize);
		for (int radius : radii)
		{
			for (int i = 0; i < data.length; i++)
			{
				double median = calculateMedian(data, i, radius);
				double median2 = calculateMedian2(data, i, radius);
				Assert.assertEquals(String.format("Position %d, Radius %d", i, radius), median2, median, 1e-6);
			}
		}
		data = createRandomData(dataSize + 1);
		for (int radius : radii)
		{
			for (int i = 0; i < data.length; i++)
			{
				double median = calculateMedian(data, i, radius);
				double median2 = calculateMedian2(data, i, radius);
				Assert.assertEquals(String.format("Position %d, Radius %d", i, radius), median2, median, 1e-6);
			}
		}
	}

	@Test
	public void canComputeMedianForRandomDataUsingDynamicLinkedListIfNewDataIsAboveMedian()
	{
		double[] data = new double[] {1,2,3,4,5};

		MedianWindowDLL mw = new MedianWindowDLL(data);
		double median = mw.getMedian();
		double median2 = calculateMedian(data, 2, 2);
		Assert.assertEquals("Before insert", median2, median, 1e-6);

		mw.add(6);
		median = mw.getMedian();
		data[0] = 6;
		median2 = calculateMedian(data, 2, 2);
		Assert.assertEquals("After insert", median2, median, 1e-6);
	}

	@Test
	public void canComputeMedianForRandomDataUsingDynamicLinkedList()
	{
		double[] data = createRandomData(dataSize);
		for (int radius : radii)
		{
			double[] startData = Arrays.copyOf(data, 2 * radius + 1);
			MedianWindowDLL mw = new MedianWindowDLL(startData);
			for (int i = radius, j = startData.length; j < data.length; i++, j++)
			{
				double median = mw.getMedian();
				mw.add(data[j]);
				double median2 = calculateMedian(data, i, radius);
				System.out.printf("Position %d, Radius %d : %g vs %g\n", i, radius, median2, median);
				Assert.assertEquals(String.format("Position %d, Radius %d", i, radius), median2, median, 1e-6);
			}
		}
	}

	@Test
	public void canComputeMedianForRandomDataUsingSingleIncrement()
	{
		canComputeMedianForDataUsingSingleIncrement(createRandomData(dataSize));
	}

	@Test
	public void canComputeMedianForRandomDataUsingSetPosition()
	{
		canComputeMedianForDataUsingSetPosition(createRandomData(dataSize));
	}

	@Test
	public void canComputeMedianForRandomDataUsingBigIncrement()
	{
		canComputeMedianForDataUsingBigIncrement(createRandomData(dataSize));
	}

	@Test
	public void canComputeMedianForDuplicateDataUsingSingleIncrement()
	{
		for (double value : values)
			canComputeMedianForDataUsingSingleIncrement(createDuplicateData(dataSize, value));
	}

	@Test
	public void canComputeMedianForDuplicateDataUsingSetPosition()
	{
		for (double value : values)
			canComputeMedianForDataUsingSetPosition(createDuplicateData(dataSize, value));
	}

	@Test
	public void canComputeMedianForDuplicateDataUsingBigIncrement()
	{
		for (double value : values)
			canComputeMedianForDataUsingBigIncrement(createDuplicateData(dataSize, value));
	}

	@Test
	public void canComputeMedianForSparseDataUsingSingleIncrement()
	{
		for (double value : values)
			canComputeMedianForDataUsingSingleIncrement(createSparseData(dataSize, value));
	}

	@Test
	public void canComputeMedianForSparseDataUsingSetPosition()
	{
		for (double value : values)
			canComputeMedianForDataUsingSetPosition(createSparseData(dataSize, value));
	}

	@Test
	public void canComputeMedianForSparseDataUsingBigIncrement()
	{
		for (double value : values)
			canComputeMedianForDataUsingBigIncrement(createSparseData(dataSize, value));
	}

	private void canComputeMedianForDataUsingSingleIncrement(double[] data)
	{
		for (int radius : radii)
		{
			MedianWindow mw = new MedianWindow(data, radius);
			for (int i = 0; i < data.length; i++)
			{
				double median = mw.getMedian();
				mw.increment();
				double median2 = calculateMedian(data, i, radius);
				Assert.assertEquals(String.format("Position %d, Radius %d", i, radius), median2, median, 1e-6);
			}
		}
	}

	private void canComputeMedianForDataUsingSetPosition(double[] data)
	{
		for (int radius : radii)
		{
			MedianWindow mw = new MedianWindow(data, radius);
			for (int i = 0; i < data.length; i += 10)
			{
				mw.setPosition(i);
				double median = mw.getMedian();
				double median2 = calculateMedian(data, i, radius);
				Assert.assertEquals(String.format("Position %d, Radius %d", i, radius), median2, median, 1e-6);
			}
		}
	}

	private void canComputeMedianForDataUsingBigIncrement(double[] data)
	{
		int increment = 10;
		for (int radius : radii)
		{
			MedianWindow mw = new MedianWindow(data, radius);
			for (int i = 0; i < data.length; i += increment)
			{
				double median = mw.getMedian();
				mw.increment(increment);
				double median2 = calculateMedian(data, i, radius);
				Assert.assertEquals(String.format("Position %d, Radius %d", i, radius), median2, median, 1e-6);
			}
		}
	}

	@Test
	public void canIncrementThroughTheDataArray()
	{
		double[] data = createRandomData(300);
		for (int radius : radii)
		{
			MedianWindow mw = new MedianWindow(data, radius);
			int i = 0;
			while (mw.isValidPosition())
			{
				double median = mw.getMedian();
				double median2 = calculateMedian(data, i, radius);
				Assert.assertEquals(String.format("Position %d, Radius %d", i, radius), median2, median, 1e-6);

				mw.increment();
				i++;
			}
			Assert.assertEquals("Not all data interated", i, data.length);

			mw = new MedianWindow(data, radius);
			i = 0;
			do
			{
				double median = mw.getMedian();
				double median2 = calculateMedian(data, i, radius);
				Assert.assertEquals(String.format("Position %d, Radius %d", i, radius), median2, median, 1e-6);

				i++;
			} while (mw.increment());
			Assert.assertEquals("Not all data interated", i, data.length);
		}
	}

	@Test
	public void canIncrementThroughTheDataArrayUsingBigIncrement()
	{
		double[] data = createRandomData(300);
		int increment = 10;
		for (int radius : radii)
		{
			MedianWindow mw = new MedianWindow(data, radius);
			int i = 0;
			while (mw.isValidPosition())
			{
				double median = mw.getMedian();
				double median2 = calculateMedian(data, i, radius);
				Assert.assertEquals(String.format("Position %d, Radius %d", i, radius), median2, median, 1e-6);

				mw.increment(increment);
				i += increment;
			}
		}
	}

	@Test
	public void returnNaNForInvalidPositions()
	{
		double[] data = createRandomData(300);
		for (int radius : radii)
		{
			MedianWindow mw = new MedianWindow(data, radius);
			for (int i = 0; i < data.length; i++)
			{
				mw.increment();
			}
			Assert.assertEquals(Double.NaN, mw.getMedian(), 1e-6);

			mw = new MedianWindow(data, radius);
			while (mw.isValidPosition())
			{
				mw.increment();
			}
			Assert.assertEquals(Double.NaN, mw.getMedian(), 1e-6);

			mw = new MedianWindow(data, radius);
			mw.setPosition(data.length + 10);
			Assert.assertEquals(Double.NaN, mw.getMedian(), 1e-6);
		}
	}

	@Test
	public void isFaster()
	{
		for (int radius : speedRadii)
		{
			for (int increment : speedIncrement)
			{
				isFaster(radius, increment);
			}
		}
	}

	private void isFaster(int radius, int increment)
	{
		int iterations = 20;
		double[][] data = new double[iterations][];
		for (int i = 0; i < iterations; i++)
			data[i] = createRandomData(dataSize);

		double[] m1 = new double[dataSize];
		// Initialise class
		MedianWindow mw = new MedianWindow(data[0], radius);
		long t1;
		if (increment == 1)
		{
			do
			{
				mw.getMedian();
			} while (mw.increment());

			long s1 = System.nanoTime();
			for (int iter = 0; iter < iterations; iter++)
			{
				mw = new MedianWindow(data[iter], radius);
				int j = 0;
				do
				{
					m1[j++] = mw.getMedian();
				} while (mw.increment());
			}
			t1 = System.nanoTime() - s1;
		}
		else
		{
			while (mw.isValidPosition())
			{
				mw.getMedian();
				mw.increment(increment);
			}

			long s1 = System.nanoTime();
			for (int iter = 0; iter < iterations; iter++)
			{
				mw = new MedianWindow(data[iter], radius);
				int j = 0;
				while (mw.isValidPosition())
				{
					m1[j++] = mw.getMedian();
					mw.increment(increment);
				}
			}
			t1 = System.nanoTime() - s1;
		}

		double[] m2 = new double[dataSize];
		// Initialise
		for (int i = 0; i < dataSize; i += increment)
		{
			calculateMedian(data[0], i, radius);
		}
		long s2 = System.nanoTime();
		for (int iter = 0; iter < iterations; iter++)
		{
			for (int i = 0, j = 0; i < dataSize; i += increment)
			{
				m2[j++] = calculateMedian(data[iter], i, radius);
			}
		}
		long t2 = System.nanoTime() - s2;

		Assert.assertArrayEquals(m1, m2, 1e-6);
		System.out.printf("Radius %d, Increment %d : window %d : standard %d = %fx faster\n", radius, increment, t1,
				t2, (double) t2 / t1);

		// Only test the largest radii 
		if (radius == speedRadii[speedRadii.length - 1])
			Assert.assertTrue(String.format("Radius %d, Increment %d", radius, increment), t1 < t2);
	}

	@Test
	public void floatVersionIsFasterThanDoubleVersion()
	{
		for (int radius : speedRadii)
		{
			for (int increment : speedIncrement)
			{
				floatVersionIsFasterThanDoubleVersion(radius, increment);
			}
		}
	}

	private void floatVersionIsFasterThanDoubleVersion(int radius, int increment)
	{
		int iterations = 20;
		double[][] data = new double[iterations][];
		float[][] data2 = new float[iterations][];
		for (int i = 0; i < iterations; i++)
		{
			data[i] = createRandomData(dataSize);
			data2[i] = copyData(data[i]);
		}

		double[] m1 = new double[dataSize];
		// Initialise class
		MedianWindow mw = new MedianWindow(data[0], radius);
		long t1;
		if (increment == 1)
		{
			do
			{
				mw.getMedian();
			} while (mw.increment());

			long s1 = System.nanoTime();
			for (int iter = 0; iter < iterations; iter++)
			{
				mw = new MedianWindow(data[iter], radius);
				int j = 0;
				do
				{
					m1[j++] = mw.getMedian();
				} while (mw.increment());
			}
			t1 = System.nanoTime() - s1;
		}
		else
		{
			while (mw.isValidPosition())
			{
				mw.getMedian();
				mw.increment(increment);
			}

			long s1 = System.nanoTime();
			for (int iter = 0; iter < iterations; iter++)
			{
				mw = new MedianWindow(data[iter], radius);
				int j = 0;
				while (mw.isValidPosition())
				{
					m1[j++] = mw.getMedian();
					mw.increment(increment);
				}
			}
			t1 = System.nanoTime() - s1;
		}

		double[] m2 = new double[dataSize];
		// Initialise
		MedianWindowFloat mw2 = new MedianWindowFloat(data2[0], radius);
		long t2;
		if (increment == 1)
		{
			do
			{
				mw2.getMedian();
			} while (mw2.increment());

			long s2 = System.nanoTime();
			for (int iter = 0; iter < iterations; iter++)
			{
				mw2 = new MedianWindowFloat(data2[iter], radius);
				int j = 0;
				do
				{
					m2[j++] = mw2.getMedian();
				} while (mw2.increment());
			}
			t2 = System.nanoTime() - s2;
		}
		else
		{
			while (mw2.isValidPosition())
			{
				mw2.getMedian();
				mw2.increment(increment);
			}

			long s2 = System.nanoTime();
			for (int iter = 0; iter < iterations; iter++)
			{
				mw2 = new MedianWindowFloat(data2[iter], radius);
				int j = 0;
				while (mw2.isValidPosition())
				{
					m2[j++] = mw2.getMedian();
					mw2.increment(increment);
				}
			}
			t2 = System.nanoTime() - s2;
		}

		Assert.assertArrayEquals(m1, m2, 1e-3);
		System.out.printf("Radius %d, Increment %d : double %d : float %d = %fx faster\n", radius, increment, t1, t2,
				(double) t1 / t2);

		// Only test the largest radii 
		if (radius == speedRadii[speedRadii.length - 1])
			Assert.assertTrue(String.format("Radius %d, Increment %d", radius, increment), t2 < t1);
	}

	@Test
	public void intVersionIsFasterThanDoubleVersion()
	{
		for (int radius : speedRadii)
		{
			for (int increment : speedIncrement)
			{
				intVersionIsFasterThanDoubleVersion(radius, increment);
			}
		}
	}

	private void intVersionIsFasterThanDoubleVersion(int radius, int increment)
	{
		int iterations = 20;
		double[][] data = new double[iterations][];
		int[][] data2 = new int[iterations][];
		for (int i = 0; i < iterations; i++)
		{
			data[i] = createRandomData(dataSize);
			data2[i] = copyDataInt(data[i]);
		}

		// Initialise class
		MedianWindow mw = new MedianWindow(data[0], radius);
		long t1;
		if (increment == 1)
		{
			do
			{
				mw.getMedian();
			} while (mw.increment());

			long s1 = System.nanoTime();
			for (int iter = 0; iter < iterations; iter++)
			{
				mw = new MedianWindow(data[iter], radius);
				do
				{
					mw.getMedian();
				} while (mw.increment());
			}
			t1 = System.nanoTime() - s1;
		}
		else
		{
			while (mw.isValidPosition())
			{
				mw.getMedian();
				mw.increment(increment);
			}

			long s1 = System.nanoTime();
			for (int iter = 0; iter < iterations; iter++)
			{
				mw = new MedianWindow(data[iter], radius);
				while (mw.isValidPosition())
				{
					mw.getMedian();
					mw.increment(increment);
				}
			}
			t1 = System.nanoTime() - s1;
		}

		// Initialise
		MedianWindowInt mw2 = new MedianWindowInt(data2[0], radius);
		long t2;
		if (increment == 1)
		{
			do
			{
				mw2.getMedian();
			} while (mw2.increment());

			long s2 = System.nanoTime();
			for (int iter = 0; iter < iterations; iter++)
			{
				mw2 = new MedianWindowInt(data2[iter], radius);
				do
				{
					mw2.getMedian();
				} while (mw2.increment());
			}
			t2 = System.nanoTime() - s2;
		}
		else
		{
			while (mw2.isValidPosition())
			{
				mw2.getMedian();
				mw2.increment(increment);
			}

			long s2 = System.nanoTime();
			for (int iter = 0; iter < iterations; iter++)
			{
				mw2 = new MedianWindowInt(data2[iter], radius);
				while (mw2.isValidPosition())
				{
					mw2.getMedian();
					mw2.increment(increment);
				}
			}
			t2 = System.nanoTime() - s2;
		}

		System.out.printf("Radius %d, Increment %d : double %d : int %d = %fx faster\n", radius, increment, t1, t2,
				(double) t1 / t2);

		// Only test the largest radii 
		if (radius == speedRadii[speedRadii.length - 1])
			Assert.assertTrue(String.format("Radius %d, Increment %d", radius, increment), t2 < t1);
	}

	static double calculateMedian(double[] data, int position, int radius)
	{
		final int start = FastMath.max(0, position - radius);
		final int end = FastMath.min(position + radius + 1, data.length);
		double[] cache = new double[end - start];
		for (int i = start, j = 0; i < end; i++, j++)
			cache[j] = data[i];
		//System.out.println(Arrays.toString(cache));
		Arrays.sort(cache);
		return (cache[(cache.length - 1) / 2] + cache[cache.length / 2]) * 0.5;
	}

	static double calculateMedian2(double[] data, int position, int radius)
	{
		final int start = FastMath.max(0, position - radius);
		final int end = FastMath.min(position + radius + 1, data.length);
		double[] cache = new double[end - start];
		for (int i = start, j = 0; i < end; i++, j++)
			cache[j] = data[i];
		DescriptiveStatistics stats = new DescriptiveStatistics(cache);
		return stats.getPercentile(50);
	}

	double[] createRandomData(int size)
	{
		double[] data = new double[size];
		for (int i = 0; i < data.length; i++)
			data[i] = random.next() * size;
		return data;
	}

	double[] createDuplicateData(int size, double value)
	{
		double[] data = new double[size];
		Arrays.fill(data, value);
		return data;
	}

	double[] createSparseData(int size, double value)
	{
		double[] data = new double[size];
		for (int i = 0; i < data.length; i++)
		{
			data[i] = value;
			if (i % 32 == 0)
				value++;
		}
		random.shuffle(data);
		return data;
	}

	static float[] copyData(double[] data)
	{
		float[] data2 = new float[data.length];
		for (int i = 0; i < data.length; i++)
		{
			data2[i] = (float) data[i];
		}
		return data2;
	}

	static int[] copyDataInt(double[] data)
	{
		int[] data2 = new int[data.length];
		for (int i = 0; i < data.length; i++)
		{
			data2[i] = (int) data[i];
		}
		return data2;
	}
}
