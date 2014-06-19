package gdsc.smlm.utils;

import java.util.Arrays;

import org.junit.Assert;
import org.junit.Test;

public class MedianWindowDLLTest
{
	MedianWindowTest mwt = new MedianWindowTest();
	Random random = new Random();
	int dataSize = 2000;
	int[] radii = new int[] { 0, 1, 2, 4, 8, 16 };
	double[] values = new double[] { 0, -1.1, 2.2 };
	int[] speedRadii = new int[] { 16, 32, 64 };
	int[] speedIncrement = new int[] { 1, 2, 4, 6, 8, 12, 16, 24, 32, 48 };

	@Test
	public void canComputeMedianForRandomDataUsingDynamicLinkedListIfNewDataIsAboveMedian()
	{
		double[] data = new double[] { 1, 2, 3, 4, 5 };

		MedianWindowDLL mw = new MedianWindowDLL(data);
		double median = mw.getMedian();
		double median2 = MedianWindowTest.calculateMedian(data, 2, 2);
		Assert.assertEquals("Before insert", median2, median, 1e-6);

		double[] insert = new double[] { 6, 7, 6, 7 };
		for (int i = 0; i < insert.length; i++)
		{
			mw.add(insert[i]);
			median = mw.getMedian();
			data[i] = insert[i];
			median2 = MedianWindowTest.calculateMedian(data, 2, 2);
			Assert.assertEquals("After insert", median2, median, 1e-6);
		}
	}

	@Test
	public void canComputeMedianForRandomDataUsingDynamicLinkedListIfNewDataIsBelowMedian()
	{
		double[] data = new double[] { 4, 5, 6, 7, 8 };

		MedianWindowDLL mw = new MedianWindowDLL(data);
		double median = mw.getMedian();
		double median2 = MedianWindowTest.calculateMedian(data, 2, 2);
		Assert.assertEquals("Before insert", median2, median, 1e-6);

		double[] insert = new double[] { 3, 2, 3, 2 };
		for (int i = 0; i < insert.length; i++)
		{
			mw.add(insert[i]);
			median = mw.getMedian();
			data[i] = insert[i];
			median2 = MedianWindowTest.calculateMedian(data, 2, 2);
			Assert.assertEquals("After insert", median2, median, 1e-6);
		}
	}

	@Test
	public void canComputeMedianForRandomDataUsingDynamicLinkedListIfNewDataIsMedianOrAbove()
	{
		double[] data = new double[] { 1, 2, 3, 4, 5 };

		MedianWindowDLL mw = new MedianWindowDLL(data);
		double median = mw.getMedian();
		double median2 = MedianWindowTest.calculateMedian(data, 2, 2);
		Assert.assertEquals("Before insert", median2, median, 1e-6);

		double[] insert = new double[] { 3, 6, 3, 6 };
		for (int i = 0; i < insert.length; i++)
		{
			mw.add(insert[i]);
			median = mw.getMedian();
			data[i] = insert[i];
			median2 = MedianWindowTest.calculateMedian(data, 2, 2);
			Assert.assertEquals("After insert", median2, median, 1e-6);
		}
	}

	@Test
	public void canComputeMedianForRandomDataUsingDynamicLinkedListIfNewDataIsMedianOrBelow()
	{
		double[] data = new double[] { 1, 2, 3, 4, 5 };

		MedianWindowDLL mw = new MedianWindowDLL(data);
		double median = mw.getMedian();
		double median2 = MedianWindowTest.calculateMedian(data, 2, 2);
		Assert.assertEquals("Before insert", median2, median, 1e-6);

		double[] insert = new double[] { 3, 0, 3, 0 };
		for (int i = 0; i < insert.length; i++)
		{
			mw.add(insert[i]);
			median = mw.getMedian();
			data[i] = insert[i];
			median2 = MedianWindowTest.calculateMedian(data, 2, 2);
			Assert.assertEquals("After insert", median2, median, 1e-6);
		}
	}

	@Test
	public void canComputeMedianForRandomDataUsingDynamicLinkedList()
	{
		double[] data = mwt.createRandomData(dataSize);
		for (int radius : radii)
		{
			double[] startData = Arrays.copyOf(data, 2 * radius + 1);
			MedianWindowDLL mw = new MedianWindowDLL(startData);
			int p = 0;
			for (int i = 0; i < radius; i++, p++)
			{
				double median = mw.getMedianOldest(i + 1 + radius);
				double median2 = MedianWindowTest.calculateMedian(data, p, radius);
				//System.out.printf("Position %d, Radius %d : %g vs %g\n", p, radius, median2, median);
				Assert.assertEquals(String.format("Position %d, Radius %d", p, radius), median2, median, 1e-6);
			}
			for (int j = startData.length; j < data.length; j++, p++)
			{
				double median = mw.getMedian();
				mw.add(data[j]);
				double median2 = MedianWindowTest.calculateMedian(data, p, radius);
				//System.out.printf("Position %d, Radius %d : %g vs %g\n", p, radius, median2, median);
				Assert.assertEquals(String.format("Position %d, Radius %d", p, radius), median2, median, 1e-6);
			}
			for (int i = 2*radius+1; i-- > 0; p++)
			{
				double median = mw.getMedianYoungest(i + 1);
				double median2 = MedianWindowTest.calculateMedian(data, p, radius);
				//System.out.printf("Position %d, Radius %d : %g vs %g\n", p, radius, median2, median);
				Assert.assertEquals(String.format("Position %d, Radius %d", p, radius), median2, median, 1e-6);
			}
		}
	}

	@Test
	public void canComputeMedianForSparseDataUsingDynamicLinkedList()
	{
		for (double value : values)
		{
			double[] data = mwt.createSparseData(dataSize, value);
			for (int radius : radii)
			{
				double[] startData = Arrays.copyOf(data, 2 * radius + 1);
				MedianWindowDLL mw = new MedianWindowDLL(startData);
				int p = 0;
				for (int i = 0; i < radius; i++, p++)
				{
					double median = mw.getMedianOldest(i + 1 + radius);
					double median2 = MedianWindowTest.calculateMedian(data, p, radius);
					//System.out.printf("Position %d, Radius %d : %g vs %g\n", p, radius, median2, median);
					Assert.assertEquals(String.format("Position %d, Radius %d", p, radius), median2, median, 1e-6);
				}
				for (int j = startData.length; j < data.length; j++, p++)
				{
					double median = mw.getMedian();
					mw.add(data[j]);
					double median2 = MedianWindowTest.calculateMedian(data, p, radius);
					//System.out.printf("Position %d, Radius %d : %g vs %g\n", p, radius, median2, median);
					Assert.assertEquals(String.format("Position %d, Radius %d", p, radius), median2, median, 1e-6);
				}
				for (int i = 2*radius+1; i-- > 0; p++)
				{
					double median = mw.getMedianYoungest(i + 1);
					double median2 = MedianWindowTest.calculateMedian(data, p, radius);
					//System.out.printf("Position %d, Radius %d : %g vs %g\n", p, radius, median2, median);
					Assert.assertEquals(String.format("Position %d, Radius %d", p, radius), median2, median, 1e-6);
				}
			}
		}
	}

	@Test
	public void canComputeMedianForDuplicateDataUsingDynamicLinkedList()
	{
		for (double value : values)
		{
			double[] data = mwt.createDuplicateData(dataSize, value);
			for (int radius : radii)
			{
				double[] startData = Arrays.copyOf(data, 2 * radius + 1);
				MedianWindowDLL mw = new MedianWindowDLL(startData);
				int p = 0;
				for (int i = 0; i < radius; i++, p++)
				{
					double median = mw.getMedianOldest(i + 1 + radius);
					double median2 = MedianWindowTest.calculateMedian(data, p, radius);
					//System.out.printf("Position %d, Radius %d : %g vs %g\n", p, radius, median2, median);
					Assert.assertEquals(String.format("Position %d, Radius %d", p, radius), median2, median, 1e-6);
				}
				for (int j = startData.length; j < data.length; j++, p++)
				{
					double median = mw.getMedian();
					mw.add(data[j]);
					double median2 = MedianWindowTest.calculateMedian(data, p, radius);
					//System.out.printf("Position %d, Radius %d : %g vs %g\n", p, radius, median2, median);
					Assert.assertEquals(String.format("Position %d, Radius %d", p, radius), median2, median, 1e-6);
				}
				for (int i = 2*radius+1; i-- > 0; p++)
				{
					double median = mw.getMedianYoungest(i + 1);
					double median2 = MedianWindowTest.calculateMedian(data, p, radius);
					//System.out.printf("Position %d, Radius %d : %g vs %g\n", p, radius, median2, median);
					Assert.assertEquals(String.format("Position %d, Radius %d", p, radius), median2, median, 1e-6);
				}
			}
		}
	}

	@Test
	public void isFasterThenMedianWindowUsingSortedCacheData()
	{
		for (int radius : speedRadii)
		{
			for (int increment : speedIncrement)
			{
				if (increment > radius)
					continue;
				isFaster(radius, increment);
			}
		}
	}

	private void isFaster(int radius, int increment)
	{
		int iterations = 20;
		double[][] data = new double[iterations][];
		for (int i = 0; i < iterations; i++)
			data[i] = mwt.createRandomData(dataSize);

		double[] m1 = new double[dataSize];
		// Initialise class
		int finalPosition = dataSize - radius;
		MedianWindow mw = new MedianWindow(data[0], radius);
		mw.setPosition(radius);
		long t1;
		if (increment == 1)
		{
			int j = 0;
			do
			{
				m1[j++] = mw.getMedian();
				mw.increment();
			} while (mw.getPosition() < finalPosition);

			long s1 = System.nanoTime();
			for (int iter = 0; iter < iterations; iter++)
			{
				mw = new MedianWindow(data[iter], radius);
				mw.setPosition(radius);
				do
				{
					mw.getMedian();
					mw.increment();
				} while (mw.getPosition() < finalPosition);
			}
			t1 = System.nanoTime() - s1;
		}
		else
		{
			int j = 0;
			do
			{
				m1[j++] = mw.getMedian();
				mw.increment(increment);
			} while (mw.getPosition() < finalPosition);

			long s1 = System.nanoTime();
			for (int iter = 0; iter < iterations; iter++)
			{
				mw = new MedianWindow(data[iter], radius);
				mw.setPosition(radius);
				do
				{
					mw.getMedian();
					mw.increment(increment);
				} while (mw.getPosition() < finalPosition);
			}
			t1 = System.nanoTime() - s1;
		}

		double[] m2 = new double[dataSize];
		double[] startData = Arrays.copyOf(data[0], 2 * radius + 1);
		MedianWindowDLL mw2 = new MedianWindowDLL(startData);
		long t2;
		if (increment == 1)
		{
			int k = 0;
			m2[k++] = mw2.getMedian();
			for (int j = startData.length; j < data[0].length; j++)
			{
				mw2.add(data[0][j]);
				m2[k++] = mw2.getMedian();
			}
			long s2 = System.nanoTime();
			for (int iter = 0; iter < iterations; iter++)
			{
				startData = Arrays.copyOf(data[iter], 2 * radius + 1);
				mw2 = new MedianWindowDLL(startData);
				mw2.getMedian();
				for (int j = startData.length; j < data[iter].length; j++)
				{
					mw2.add(data[iter][j]);
					mw2.getMedian();
				}
			}
			t2 = System.nanoTime() - s2;
		}
		else
		{
			final int limit = data[0].length - increment;
			int k = 0;
			m2[k++] = mw2.getMedian();
			for (int j = startData.length; j < limit; j += increment)
			{
				for (int i = 0; i < increment; i++)
					mw2.add(data[0][j + i]);
				m2[k++] = mw2.getMedian();
			}
			long s2 = System.nanoTime();
			for (int iter = 0; iter < iterations; iter++)
			{
				startData = Arrays.copyOf(data[iter], 2 * radius + 1);
				mw2 = new MedianWindowDLL(startData);
				mw2.getMedian();
				for (int j = startData.length; j < limit; j += increment)
				{
					for (int i = 0; i < increment; i++)
						mw2.add(data[iter][j + i]);
					mw2.getMedian();
				}
			}
			t2 = System.nanoTime() - s2;
		}

		Assert.assertArrayEquals(String.format("Radius %d, Increment %d", radius, increment), m1, m2, 1e-6);
		System.out.printf("Radius %d, Increment %d : Cached %d : DLL %d = %fx faster\n", radius, increment, t1, t2,
				(double) t1 / t2);

		// Only test the largest radii with a reasonable increment
		if (radius == speedRadii[speedRadii.length - 1] && increment <= radius / 2)
			Assert.assertTrue(String.format("Radius %d, Increment %d", radius, increment), t2 < t1);
	}
}
