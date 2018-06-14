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

import gdsc.test.TimingResult;
import gdsc.test.TimingService;
import gdsc.test.TimingTask;

public class RampedSelectionStrategyTest
{
	@Test
	public void canSearchUsingActualKey()
	{
		long[] sum = RampedSelectionStrategy.createSum(10);

		for (int i = 0; i < sum.length - 1; i++)
		{
			long key = sum[i];
			int j = RampedSelectionStrategy.search(sum, key);
			Assert.assertEquals(i + 1, j);
		}
	}

	@Test
	public void canBinarySearchUsingActualKey()
	{
		long[] sum = RampedSelectionStrategy.createSum(10);

		for (int i = 0; i < sum.length - 1; i++)
		{
			long key = sum[i];
			int j = RampedSelectionStrategy.binarySearch(sum, key);
			Assert.assertEquals(i + 1, j);
		}
	}

	@Test
	public void canSearchUsingNotActualKey()
	{
		long[] sum = RampedSelectionStrategy.createSum(10);

		for (int i = 0; i < sum.length; i++)
		{
			long key = sum[i] - 1;
			int j = RampedSelectionStrategy.search(sum, key);
			Assert.assertEquals(i, j);
		}
	}

	@Test
	public void canBinarySearchUsingNotActualKey()
	{
		long[] sum = RampedSelectionStrategy.createSum(10);

		for (int i = 0; i < sum.length; i++)
		{
			long key = sum[i] - 1;
			int j = RampedSelectionStrategy.binarySearch(sum, key);
			Assert.assertEquals(i, j);
		}
	}

	@Test
	public void binarySearchEqualsSearch()
	{
		long[] sum = RampedSelectionStrategy.createSum(100);
		for (int key = (int) sum[sum.length - 1]; key-- > 0;)
		{
			int i = RampedSelectionStrategy.search(sum, key);
			int j = RampedSelectionStrategy.binarySearch(sum, key);
			Assert.assertEquals(i, j);
		}
	}

	@Test
	public void speedTest50()
	{
		speedTest(50, false, 10);
	}

	@Test
	public void speedTest200()
	{
		speedTest(200, true, 5);
	}

	@Test
	public void speedTest1000()
	{
		speedTest(1000, true, 2);
	}

	// Too slow for common use
	//@Test
	public void speedTest5000()
	{
		speedTest(5000, true, 1);
	}

	private void speedTest(final int size, boolean faster, int runs)
	{
		final long[] sum = RampedSelectionStrategy.createSum(size);

		TimingService ts = new TimingService(runs);

		ts.execute(new TimingTask()
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
			public void check(int i, Object result)
			{
			}

			@Override
			public int getSize()
			{
				return 1;
			}

			@Override
			public String getName()
			{
				return "search" + size;
			}
		});

		ts.execute(new TimingTask()
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
			public void check(int i, Object result)
			{
			}

			@Override
			public int getSize()
			{
				return 1;
			}

			@Override
			public String getName()
			{
				return "binarySearch" + size;
			}
		});

		int n = ts.repeat();
		ts.repeat(n);

		ts.report();

		TimingResult slow = ts.get((faster) ? ts.getSize() - 2 : ts.getSize() - 1);
		TimingResult fast = ts.get((faster) ? ts.getSize() - 1 : ts.getSize() - 2);
		Assert.assertTrue(slow.getMin() > fast.getMin());
	}
}
