package gdsc.smlm.ij.utils;

import org.junit.Assert;
import org.junit.Test;

import ij.process.ImageProcessor;

public abstract class Image2DTest
{
	protected abstract Image2D createData(int w, int h);

	@Test
	public void canCrop()
	{
		canCrop(3, 4, 5, 6);
		canCrop(3, 4, 1, 1);
		canCrop(0, 0, 1, 2);
	}

	private void canCrop(int x, int y, int w, int h)
	{
		Image2D image = createData(x + w + 1, y + h + 1);
		ImageProcessor stack = image.getImageProcessor();

		Image2D croppedData = image.crop(x, y, w, h);
		Assert.assertEquals(croppedData.getWidth(), w);
		Assert.assertEquals(croppedData.getHeight(), h);

		Image2D croppedData2 = FloatImage2D.crop(stack, x, y, w, h, null);
		assertEquals(croppedData, croppedData2);

		ImageProcessor croppedStack = image.cropToProcessor(x, y, w, h);
		Assert.assertEquals(croppedStack.getWidth(), w);
		Assert.assertEquals(croppedStack.getHeight(), h);

		assertEquals(croppedData, new FloatImage2D(croppedStack));

		// Test it is the correct region
		ImageProcessor fp = image.getImageProcessor();

		fp.setRoi(x, y, w, h);
		fp = fp.crop();
		float[] e = (float[]) fp.getPixels();

		// Compare to the cropped stack 
		float[] o = (float[]) croppedStack.getPixels();
		Assert.assertArrayEquals(e, o, 0);

		// Compare to the cropped data
		for (int i = 0; i < e.length; i++)
			Assert.assertEquals(e[i], croppedData.get(i), 0);
	}

	@Test
	public void canInsert()
	{
		canInsert(3, 4, 5, 6);
		canInsert(3, 4, 1, 1);
		canInsert(0, 0, 1, 1);
	}

	private void canInsert(int x, int y, int w, int h)
	{
		// This test assumes that copy and crop work!
		for (int pad : new int[] { 0, 1 })
		{
			Image2D image = createData(x + w + pad, y + h + pad);
			Image2D image2 = image.copy();
			Image2D image3 = image.copy();

			Image2D blank = new FloatImage2D(w, h);

			image.insert(x, y, blank);

			Image2D croppedData = image.crop(x, y, w, h);

			assertEquals(croppedData, blank);

			ImageProcessor blankStack = blank.getImageProcessor();
			image2.insert(x, y, blankStack);

			croppedData = image2.crop(x, y, w, h);

			assertEquals(croppedData, blank);

			image3.insert(x, y, blankStack);

			croppedData = image3.crop(x, y, w, h);

			assertEquals(croppedData, blank);
		}
	}

	@Test
	public void canFindMin()
	{
		Image2D image = createData(2, 2);
		Assert.assertEquals(0, image.findMinIndex(0, 0, 2, 2));
		Assert.assertEquals(1, image.findMinIndex(1, 0, 2, 2));
		Assert.assertEquals(2, image.findMinIndex(0, 1, 2, 2));
		Assert.assertEquals(3, image.findMinIndex(1, 1, 2, 2));

		// Larger slices
		canFindMin(3, 4, 5, 6);
		canFindMin(3, 4, 1, 1);
		canFindMin(0, 0, 1, 1);
	}

	private void canFindMin(int x, int y, int w, int h)
	{
		// This test assumes that crop works!
		for (int pad : new int[] { 0, 1 })
		{
			Image2D image = createData(x + w + pad, y + h + pad);

			Image2D croppedData = image.crop(x, y, w, h);
			int i = findMinIndex(croppedData);
			int[] xy = croppedData.getXY(i);

			int j = image.findMinIndex(x, y, w, h);
			int[] xy2 = image.getXY(j);

			Assert.assertEquals(xy[0] + x, xy2[0]);
			Assert.assertEquals(xy[1] + y, xy2[1]);
		}
	}

	@Test
	public void canFindMax()
	{
		Image2D image = createData(2, 2);
		Assert.assertEquals(3, image.findMaxIndex(0, 0, 2, 2));
		Assert.assertEquals(2, image.findMaxIndex(0, 0, 1, 2));
		Assert.assertEquals(1, image.findMaxIndex(0, 0, 2, 1));
		Assert.assertEquals(0, image.findMaxIndex(0, 0, 1, 1));

		// Larger slices
		canFindMax(3, 4, 5, 6);
		canFindMax(3, 4, 1, 1);
		canFindMax(0, 0, 1, 1);
	}

	private void canFindMax(int x, int y, int w, int h)
	{
		// This test assumes that crop works!
		for (int pad : new int[] { 0, 1 })
		{
			Image2D image = createData(x + w + pad, y + h + pad);

			Image2D croppedData = image.crop(x, y, w, h);
			int i = findMaxIndex(croppedData);
			int[] xy = croppedData.getXY(i);

			int j = image.findMaxIndex(x, y, w, h);
			int[] xy2 = image.getXY(j);

			Assert.assertEquals(xy[0] + x, xy2[0]);
			Assert.assertEquals(xy[1] + y, xy2[1]);
		}
	}

	@Test
	public void canComputeSum()
	{
		// Bounds checks
		Image2D image = createData(2, 2);
		Assert.assertEquals(10, image.computeSum(0, 0, 2, 2), 0);
		Assert.assertEquals(0, image.computeSum(0, 0, 0, 0), 0);
		Assert.assertEquals(1, image.computeSum(0, 0, 1, 1), 0);
		Assert.assertEquals(2, image.computeSum(1, 0, 1, 1), 0);
		Assert.assertEquals(0, image.computeSum(-10, 0, 1, 1), 0);
		Assert.assertEquals(0, image.computeSum(10, 0, 1, 1), 0);
		Assert.assertEquals(0, image.computeSum(0, 10, 1, 1), 0);
		Assert.assertEquals(0, image.computeSum(0, -10, 1, 1), 0);

		// Larger slices
		canComputeSum(3, 4, 5, 6);
		canComputeSum(3, 4, 1, 1);
		canComputeSum(0, 0, 1, 1);
	}

	private void canComputeSum(int x, int y, int w, int h)
	{
		// This test assumes that crop works!
		for (int pad : new int[] { 0, 1 })
		{
			Image2D image = createData(x + w + pad, y + h + pad);

			Image2D croppedData = image.crop(x, y, w, h);
			double e = sum(croppedData);
			double o = image.computeSum(x, y, w, h);

			Assert.assertEquals(o, e, 0);
		}
	}

	@Test
	public void canComputeRollingSumTable()
	{
		int w = 2, h = 3;
		Image2D image = createData(w, h);
		double[] table = image.computeRollingSumTable(null);
		for (int hh = 1, i = 0; hh <= h; hh++)
			for (int ww = 1; ww <= w; ww++)
			{
				double e = image.computeSum(0, 0, ww, hh);
				double o = table[i++];
				Assert.assertEquals(e, o, 1e-3);
			}
	}

	@Test
	public void canComputeSumUsingTable()
	{
		// Bounds checks
		Image2D image = createData(2, 2);
		double[] table = image.computeRollingSumTable(null);
		//testComputeSum(36, image,table, 0, 0, 0, 2, 2, 2);
		testComputeSum(10, image, table, 0, 0, 5, 7);
		testComputeSum(0, image, table, 0, 0, 0, 0);
		testComputeSum(1, image, table, 0, 0, 1, 1);
		testComputeSum(2, image, table, 1, 0, 1, 1);
		testComputeSum(0, image, table, -10, 0, 1, 1);
		testComputeSum(0, image, table, 10, 0, 1, 1);
		testComputeSum(0, image, table, 0, 10, 1, 1);
		testComputeSum(0, image, table, 0, -10, 1, 1);

		// Larger slices
		canComputeSumUsingTable(3, 4, 5, 6);
		canComputeSumUsingTable(3, 4, 1, 1);
		canComputeSumUsingTable(0, 0, 1, 1);
	}

	private void canComputeSumUsingTable(int x, int y, int w, int h)
	{
		// This test assumes that crop works!
		for (int pad : new int[] { 0, 1 })
		{
			Image2D image = createData(x + w + pad, y + h + pad);
			double[] table = image.computeRollingSumTable(null);

			Image2D croppedData = image.crop(x, y, w, h);
			double e = sum(croppedData);

			testComputeSum(e, image, table, x, y, w, h);
		}
	}

	private void testComputeSum(double e, Image2D image, double[] table, int x, int y, int w, int h)
	{
		double o = image.computeSum(table, x, y, w, h);
		double o2 = image.computeSumFast(table, x, y, w, h);

		// This may be different due to floating point error
		// but we are adding integers so it should be OK
		Assert.assertEquals(e, o, 0);
		Assert.assertEquals(e, o2, 0);
	}

	private static void assertEquals(Image2D a, Image2D b)
	{
		for (int i = a.getDataLength(); i-- > 0;)
			Assert.assertEquals(a.get(i), b.get(i), 0);
	}

	private static int findMinIndex(Image2D image)
	{
		int min = 0;
		for (int i = image.getDataLength(); i-- > 0;)
			if (image.get(i) < image.get(min))
				min = i;
		return min;
	}

	private static int findMaxIndex(Image2D image)
	{
		int max = 0;
		for (int i = image.getDataLength(); i-- > 0;)
			if (image.get(i) > image.get(max))
				max = i;
		return max;
	}

	private static double sum(Image2D image)
	{
		double sum = 0;
		for (int i = image.getDataLength(); i-- > 0;)
			sum += image.get(i);
		return sum;
	}
}
