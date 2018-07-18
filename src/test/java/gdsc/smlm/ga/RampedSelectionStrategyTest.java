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
package gdsc.smlm.ga;

import org.junit.Assert;
import org.junit.Test;

import gdsc.test.BaseTimingTask;
import gdsc.test.LogLevel;
import gdsc.test.TestLog;
import gdsc.test.TestSettings;
import gdsc.test.TimingResult;
import gdsc.test.TimingService;

@SuppressWarnings({ "javadoc" })
public class RampedSelectionStrategyTest
{
	@Test
	public void canSearchUsingActualKey()
	{
		final long[] sum = RampedSelectionStrategy.createSum(10);

		for (int i = 0; i < sum.length - 1; i++)
		{
			final long key = sum[i];
			final int j = RampedSelectionStrategy.search(sum, key);
			Assert.assertEquals(i + 1, j);
		}
	}

	@Test
	public void canBinarySearchUsingActualKey()
	{
		final long[] sum = RampedSelectionStrategy.createSum(10);

		for (int i = 0; i < sum.length - 1; i++)
		{
			final long key = sum[i];
			final int j = RampedSelectionStrategy.binarySearch(sum, key);
			Assert.assertEquals(i + 1, j);
		}
	}

	@Test
	public void canSearchUsingNotActualKey()
	{
		final long[] sum = RampedSelectionStrategy.createSum(10);

		for (int i = 0; i < sum.length; i++)
		{
			final long key = sum[i] - 1;
			final int j = RampedSelectionStrategy.search(sum, key);
			Assert.assertEquals(i, j);
		}
	}

	@Test
	public void canBinarySearchUsingNotActualKey()
	{
		final long[] sum = RampedSelectionStrategy.createSum(10);

		for (int i = 0; i < sum.length; i++)
		{
			final long key = sum[i] - 1;
			final int j = RampedSelectionStrategy.binarySearch(sum, key);
			Assert.assertEquals(i, j);
		}
	}

	@Test
	public void binarySearchEqualsSearch()
	{
		final long[] sum = RampedSelectionStrategy.createSum(100);
		for (int key = (int) sum[sum.length - 1]; key-- > 0;)
		{
			final int i = RampedSelectionStrategy.search(sum, key);
			final int j = RampedSelectionStrategy.binarySearch(sum, key);
			Assert.assertEquals(i, j);
		}
	}

	@Test
	public void speedTest50()
	{
		TestSettings.assumeLowComplexity();
		speedTest(50, false, 10);
	}

	@Test
	public void speedTest200()
	{
		TestSettings.assumeMediumComplexity();
		speedTest(200, true, 5);
	}

	@Test
	public void speedTest1000()
	{
		TestSettings.assumeMediumComplexity();
		speedTest(1000, true, 2);
	}

	// Too slow for common use
	@Test
	public void speedTest5000()
	{
		TestSettings.assumeHighComplexity();
		speedTest(5000, true, 1);
	}

	private static void speedTest(final int size, boolean faster, int runs)
	{
		final long[] sum = RampedSelectionStrategy.createSum(size);

		final TimingService ts = new TimingService(runs);

		ts.execute(new BaseTimingTask("search" + size)
		{
			@Override
			public Object getData(int i)
			{
				return sum;
			}

			@Override
			public Object run(Object data)
			{
				for (int key = (int) sum[sum.length - 1]; key-- > 0;)
					RampedSelectionStrategy.search(sum, key);
				return null;
			}
			@Override
			public int getSize()
			{
				return 1;
			}
		});

		ts.execute(new BaseTimingTask("binarySearch" + size)
		{
			@Override
			public Object getData(int i)
			{
				return sum[i];
			}

			@Override
			public Object run(Object data)
			{
				for (int key = (int) sum[sum.length - 1]; key-- > 0;)
					RampedSelectionStrategy.binarySearch(sum, key);
				return null;
			}

			@Override
			public int getSize()
			{
				return 1;
			}
		});

		final int n = ts.repeat();
		ts.repeat(n);
		if (TestSettings.allow(LogLevel.INFO))
			ts.report();

		final TimingResult slow = ts.get((faster) ? ts.getSize() - 2 : ts.getSize() - 1);
		final TimingResult fast = ts.get((faster) ? ts.getSize() - 1 : ts.getSize() - 2);
		TestLog.logSpeedTestResult(slow, fast);
	}
}
