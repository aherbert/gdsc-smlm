package gdsc.smlm.filters;

import org.junit.Assert;
import org.junit.Test;
import org.junit.internal.ArrayComparisonFailure;

public class IntBlockSumFilterTest
{
	private gdsc.core.utils.Random rand;

	int[] primes = new int[] { 113, 97, 53, 29 };
	//int[] primes = new int[] { 1024 };
	int[] boxSizes = new int[] { 15, 9, 5, 3, 2, 1 };
	boolean[] checkInternal = new boolean[] { true, false };

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

		int n = (int) Math.ceil(boxSize);
		int size = 2 * n + 1;

		int[] out = new int[data.length];

		for (int y = 0; y < maxy; y++)
		{
			for (int x = 0; x < maxx; x++)
			{
				int sum = 0;
				for (int yy = 0; yy < size; yy++)
				{
					int yyy = y + yy - n;
					if (yyy < 0)
						//yyy = 0;
						continue;
					if (yyy >= maxy)
						//yyy = maxy - 1;
						continue;
					for (int xx = 0; xx < size; xx++)
					{
						int xxx = x + xx - n;
						if (xxx < 0)
							//xxx = 0;
							continue;
						if (xxx >= maxx)
							//xxx = maxx - 1;
							continue;
						int index = yyy * maxx + xxx;
						sum += data[index];
					}
				}
				out[y * maxx + x] = sum;
			}
		}
		System.arraycopy(out, 0, data, 0, out.length);
	}

	private void intArrayEquals(String message, int[] data1, int[] data2, int maxx, int maxy, int boxSize)
	{
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
				if (data1[index] != data2[index])
				{
					Assert.fail(String.format("%s [%d,%d] %d != %d", message, x, y, data1[index], data2[index]));
				}
			}
		}
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

		}
	}

	private void sumIsCorrect(int[] data, int width, int height, int boxSize, boolean internal,
			BlockSumDataFilter filter) throws ArrayComparisonFailure
	{
		int[] data1 = data.clone();
		int[] data2 = data.clone();

		sum(data1, width, height, boxSize);
		if (internal)
		{
			filter.filterInternal(data2, width, height, boxSize);
			intArrayEquals(String.format("Internal arrays do not match: [%dx%d] @ %d", width, height, boxSize), data1,
					data2, width, height, boxSize);
		}
		else
		{
			filter.filter(data2, width, height, boxSize);
			intArrayEquals(String.format("Arrays do not match: [%dx%d] @ %d", width, height, boxSize), data1, data2,
					width, height, 0);
		}
	}

	private void checkIsCorrect(BlockSumDataFilter filter)
	{
		rand = new gdsc.core.utils.Random(-30051976);

		for (int width : primes)
			for (int height : primes)
			{
				int[] data = createData(width, height);

				for (int boxSize : boxSizes)
					for (boolean internal : checkInternal)
					{
						sumIsCorrect(data, width, height, boxSize, internal, filter);
					}
			}
	}

	@Test
	public void rollingBlockFilterIsCorrect()
	{
		BlockSumDataFilter filter = new BlockSumDataFilter("rollingBlock", false)
		{
			public void filter(int[] data, int width, int height, int boxSize)
			{
				f.rollingBlockFilter(data, width, height, (int) boxSize);
			}

			public void filterInternal(int[] data, int width, int height, int boxSize)
			{
				f.rollingBlockFilterInternal(data, width, height, (int) boxSize);
			}
		};
		checkIsCorrect(filter);
	}

	private int[] createData(int width, int height)
	{
		int[] data = new int[width * height];
		for (int i = data.length; i-- > 0;)
			data[i] = i;

		rand.shuffle(data);

		return data;
	}
}
