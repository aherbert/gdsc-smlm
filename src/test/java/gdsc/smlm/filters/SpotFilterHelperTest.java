package gdsc.smlm.filters;

import java.awt.Rectangle;

import org.junit.Assert;
import org.junit.Test;

public class SpotFilterHelperTest
{
	private gdsc.core.utils.Random rand;

	private Spot[] createData(int width, int height, int n)
	{
		if (n == 0)
			return new Spot[0];

		int[] data = new int[width * height];
		for (int i = n; i-- > 0;)
			data[i] = 1;

		rand.shuffle(data);

		Spot[] spots = new Spot[n];
		for (int i = 0, j = 0; i < data.length; i++)
		{
			if (data[i] == 1)
			{
				spots[j++] = new Spot(i % width, i / width, 1);
				if (j == n)
					break;
			}
		}

		return spots;
	}

	@Test
	public void canCountNeighbours()
	{
		rand = new gdsc.core.utils.Random(30051977);

		int width = 64, height = 64;
		int size = width * height;
		int n = 1; // Don't test simple case of no neighbours
		SpotFilterHelper h = new SpotFilterHelper();
		while (n < size / 4)
		{
			n *= 2;
			for (int loop = 0; loop < 5; loop++)
			{
				Spot[] spots = createData(width, height, n);
				for (int box : new int[] { 1, 2, 3, 4, 5 })
				{
					int[] e = countNeighbours(spots, width, height, box);
					int[] count = h.countNeighbours(spots, box);
					Assert.assertArrayEquals(e, count);
					int[] count2 = h.countNeighbours(spots, width, height, box);
					Assert.assertArrayEquals(e, count2);
				}
			}
		}
	}

	private int[] countNeighbours(Spot[] spots, int width, int height, int box)
	{
		short[] data = new short[width * height];
		for (int i = 0; i < spots.length; i++)
			data[spots[i].x + width * spots[i].y] = 1;
		Rectangle r = new Rectangle(width, height);
		Rectangle bounds = new Rectangle(2 * box + 1, 2 * box + 1);
		int[] count = new int[spots.length];
		for (int i = 0; i < spots.length; i++)
		{
			bounds.x = spots[i].x - box;
			bounds.y = spots[i].y - box;
			Rectangle limits = r.intersection(bounds);
			int sum = -1;
			for (int y = limits.y; y < limits.y + limits.height; y++)
				for (int x = limits.x, j = y * width + x; x < limits.x + limits.width; x++)
				{
					sum += data[j++];
				}
			count[i] = sum;
		}
		return count;
	}
}
