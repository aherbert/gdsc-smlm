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

import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Test;
import org.junit.internal.ArrayComparisonFailure;

import gdsc.test.TestSettings;

@SuppressWarnings({ "javadoc" })
public class IntBlockSumFilterTest extends AbstractFilterTest
{
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
	public static void sum(int[] data, int maxx, int maxy, int boxSize)
	{
		if (boxSize <= 0)
			return;

		final int n = (int) Math.ceil(boxSize);
		final int size = 2 * n + 1;

		final int[] out = new int[data.length];

		for (int y = 0; y < maxy; y++)
			for (int x = 0; x < maxx; x++)
			{
				int sum = 0;
				for (int yy = 0; yy < size; yy++)
				{
					final int yyy = y + yy - n;
					if (yyy < 0)
						//yyy = 0;
						continue;
					if (yyy >= maxy)
						//yyy = maxy - 1;
						continue;
					for (int xx = 0; xx < size; xx++)
					{
						final int xxx = x + xx - n;
						if (xxx < 0)
							//xxx = 0;
							continue;
						if (xxx >= maxx)
							//xxx = maxx - 1;
							continue;
						final int index = yyy * maxx + xxx;
						sum += data[index];
					}
				}
				out[y * maxx + x] = sum;
			}
		System.arraycopy(out, 0, data, 0, out.length);
	}

	/**
	 * Used to test the filter methods calculate the correct result
	 */
	private class BlockSumDataFilter extends IntDataFilter
	{
		public BlockSumDataFilter(String name, boolean isInterpolated)
		{
			super(name, isInterpolated);
		}

		IntBlockSumFilter f = new IntBlockSumFilter();

		@Override
		public void filter(int[] data, int width, int height, int boxSize)
		{
			f.rollingBlockFilter(data, width, height, boxSize);
		}

		@Override
		public void filterInternal(int[] data, int width, int height, int boxSize)
		{
			f.rollingBlockFilterInternal(data, width, height, boxSize);
		}

		@Override
		public void setWeights(float[] weights, int width, int height)
		{
			// Ignore weights
		}
	}

	private static void sumIsCorrect(int[] data, int width, int height, int boxSize, boolean internal,
			BlockSumDataFilter filter) throws ArrayComparisonFailure
	{
		final int[] data1 = data.clone();
		final int[] data2 = data.clone();

		sum(data1, width, height, boxSize);
		if (internal)
		{
			filter.filterInternal(data2, width, height, boxSize);
			intArrayEquals(data1, data2, width, height, boxSize, "Internal arrays do not match: [%dx%d] @ %d", width,
					height, boxSize);
		}
		else
		{
			filter.filter(data2, width, height, boxSize);
			intArrayEquals(data1, data2, width, height, 0, "Arrays do not match: [%dx%d] @ %d", width, height, boxSize);
		}
	}

	private static void checkIsCorrect(BlockSumDataFilter filter)
	{
		final RandomGenerator rg = TestSettings.getRandomGenerator();

		for (final int width : primes)
			for (final int height : primes)
			{
				final int[] data = createIntData(rg, width, height);

				for (final int boxSize : boxSizes)
					for (final boolean internal : checkInternal)
						sumIsCorrect(data, width, height, boxSize, internal, filter);
			}
	}

	@Test
	public void rollingBlockFilterIsCorrect()
	{
		final BlockSumDataFilter filter = new BlockSumDataFilter("rollingBlock", false)
		{
			@Override
			public void filter(int[] data, int width, int height, int boxSize)
			{
				f.rollingBlockFilter(data, width, height, boxSize);
			}

			@Override
			public void filterInternal(int[] data, int width, int height, int boxSize)
			{
				f.rollingBlockFilterInternal(data, width, height, boxSize);
			}
		};
		checkIsCorrect(filter);
	}

}
