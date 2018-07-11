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
package gdsc.smlm.ij.utils;

import org.junit.Assert;
import org.junit.Test;

import ij.ImageStack;
import ij.process.ImageProcessor;

@SuppressWarnings({ "javadoc" })
public abstract class Image3DTest
{
	protected abstract Image3D createData(int w, int h, int d);

	protected abstract Image3D createEmptyData(int w, int h, int d);

	@Test
	public void canCrop()
	{
		canCrop(3, 4, 5, 6, 7, 8);
		canCrop(3, 4, 5, 1, 1, 1);
		canCrop(0, 0, 0, 1, 2, 3);
	}

	private void canCrop(int x, int y, int z, int w, int h, int d)
	{
		Image3D image = createData(x + w + 1, y + h + 1, z + d + 1);
		ImageStack stack = image.getImageStack();

		Image3D croppedData = image.crop(x, y, z, w, h, d);
		Assert.assertEquals(croppedData.getWidth(), w);
		Assert.assertEquals(croppedData.getHeight(), h);
		Assert.assertEquals(croppedData.getSize(), d);

		Image3D croppedData2 = FloatImage3D.crop(stack, x, y, z, w, h, d, null);
		assertEquals(croppedData, croppedData2);

		ImageStack croppedStack = image.cropToStack(x, y, z, w, h, d);
		Assert.assertEquals(croppedStack.getWidth(), w);
		Assert.assertEquals(croppedStack.getHeight(), h);
		Assert.assertEquals(croppedStack.getSize(), d);

		assertEquals(croppedData, new FloatImage3D(croppedStack));

		ImageStack croppedStack2 = Image3D.cropToStack(stack, x, y, z, w, h, d);
		assertEquals(croppedData, new FloatImage3D(croppedStack2));

		// Test it is the correct region
		ImageStack originalStack = image.getImageStack();
		for (int zz = 0; zz < d; zz++)
		{
			// Crop from the original data
			ImageProcessor fp = originalStack.getProcessor(1 + zz + z);
			fp.setRoi(x, y, w, h);
			fp = fp.crop();
			float[] e = (float[]) fp.getPixels();

			// Compare to the cropped stack
			float[] o = (float[]) croppedStack.getPixels(1 + zz);
			Assert.assertArrayEquals(e, o, 0);

			// Compare to the cropped data
			for (int i = 0, j = zz * e.length; i < e.length; i++, j++)
				Assert.assertEquals(e[i], croppedData.get(j), 0);
		}
	}

	@Test
	public void canInsert()
	{
		canInsert(3, 4, 5, 6, 7, 8);
		canInsert(3, 4, 5, 1, 1, 1);
		canInsert(0, 0, 0, 1, 2, 3);
	}

	private void canInsert(int x, int y, int z, int w, int h, int d)
	{
		// This test assumes that copy and crop work!
		for (int pad : new int[] { 0, 1 })
		{
			Image3D image = createData(x + w + pad, y + h + pad, z + d + pad);
			Image3D image2 = image.copy();
			Image3D image3 = image.copy();

			Image3D blank = createEmptyData(w, h, d);

			image.insert(x, y, z, blank);

			Image3D croppedData = image.crop(x, y, z, w, h, d);

			assertEquals(croppedData, blank);

			ImageStack blankStack = blank.getImageStack();
			image2.insert(x, y, z, blankStack);

			croppedData = image2.crop(x, y, z, w, h, d);

			assertEquals(croppedData, blank);

			for (int s = 0; s < blankStack.getSize(); s++)
				image3.insert(x, y, z + s, blankStack.getProcessor(1 + s));

			croppedData = image3.crop(x, y, z, w, h, d);

			assertEquals(croppedData, blank);
		}
	}

	@Test
	public void canFindMin()
	{
		Image3D image = createData(2, 2, 2);
		Assert.assertEquals(0, image.findMinIndex(0, 0, 0, 2, 2, 2));
		Assert.assertEquals(1, image.findMinIndex(1, 0, 0, 2, 2, 2));
		Assert.assertEquals(2, image.findMinIndex(0, 1, 0, 2, 2, 2));
		Assert.assertEquals(3, image.findMinIndex(1, 1, 0, 2, 2, 2));
		Assert.assertEquals(4, image.findMinIndex(0, 0, 1, 2, 2, 2));
		Assert.assertEquals(5, image.findMinIndex(1, 0, 1, 2, 2, 2));
		Assert.assertEquals(6, image.findMinIndex(0, 1, 1, 2, 2, 2));
		Assert.assertEquals(7, image.findMinIndex(1, 1, 1, 2, 2, 2));

		// Larger slices
		canFindMin(3, 4, 5, 6, 7, 8);
		canFindMin(3, 4, 5, 1, 1, 1);
		canFindMin(0, 0, 0, 1, 2, 3);
	}

	private void canFindMin(int x, int y, int z, int w, int h, int d)
	{
		// This test assumes that crop works!
		for (int pad : new int[] { 0, 1 })
		{
			Image3D image = createData(x + w + pad, y + h + pad, z + d + pad);

			Image3D croppedData = image.crop(x, y, z, w, h, d);
			int i = findMinIndex(croppedData);
			int[] xyz = croppedData.getXYZ(i);

			int j = image.findMinIndex(x, y, z, w, h, d);
			int[] xyz2 = image.getXYZ(j);

			Assert.assertEquals(xyz[0] + x, xyz2[0]);
			Assert.assertEquals(xyz[1] + y, xyz2[1]);
			Assert.assertEquals(xyz[2] + z, xyz2[2]);
		}
	}

	@Test
	public void canFindMax()
	{
		Image3D image = createData(2, 2, 2);
		Assert.assertEquals(7, image.findMaxIndex(0, 0, 0, 2, 2, 2));
		Assert.assertEquals(6, image.findMaxIndex(0, 0, 0, 1, 2, 2));
		Assert.assertEquals(5, image.findMaxIndex(0, 0, 0, 2, 1, 2));
		Assert.assertEquals(4, image.findMaxIndex(0, 0, 0, 1, 1, 2));
		Assert.assertEquals(3, image.findMaxIndex(0, 0, 0, 2, 2, 1));
		Assert.assertEquals(2, image.findMaxIndex(0, 0, 0, 1, 2, 1));
		Assert.assertEquals(1, image.findMaxIndex(0, 0, 0, 2, 1, 1));
		Assert.assertEquals(0, image.findMaxIndex(0, 0, 0, 1, 1, 1));

		// Larger slices
		canFindMax(3, 4, 5, 6, 7, 8);
		canFindMax(3, 4, 5, 1, 1, 1);
		canFindMax(0, 0, 0, 1, 2, 3);
	}

	private void canFindMax(int x, int y, int z, int w, int h, int d)
	{
		// This test assumes that crop works!
		for (int pad : new int[] { 0, 1 })
		{
			Image3D image = createData(x + w + pad, y + h + pad, z + d + pad);

			Image3D croppedData = image.crop(x, y, z, w, h, d);
			int i = findMaxIndex(croppedData);
			int[] xyz = croppedData.getXYZ(i);

			int j = image.findMaxIndex(x, y, z, w, h, d);
			int[] xyz2 = image.getXYZ(j);

			Assert.assertEquals(xyz[0] + x, xyz2[0]);
			Assert.assertEquals(xyz[1] + y, xyz2[1]);
			Assert.assertEquals(xyz[2] + z, xyz2[2]);
		}
	}

	@Test
	public void canComputeSum()
	{
		// Bounds checks
		Image3D image = createData(2, 2, 2);
		Assert.assertEquals(36, image.computeSum(0, 0, 0, 2, 2, 2), 0);
		Assert.assertEquals(0, image.computeSum(0, 0, 0, 0, 0, 0), 0);
		Assert.assertEquals(1, image.computeSum(0, 0, 0, 1, 1, 1), 0);
		Assert.assertEquals(2, image.computeSum(1, 0, 0, 1, 1, 1), 0);
		Assert.assertEquals(0, image.computeSum(-10, 0, 0, 1, 1, 1), 0);
		Assert.assertEquals(0, image.computeSum(10, 0, 0, 1, 1, 1), 0);
		Assert.assertEquals(0, image.computeSum(0, 10, 0, 1, 1, 1), 0);
		Assert.assertEquals(0, image.computeSum(0, -10, 0, 1, 1, 1), 0);
		Assert.assertEquals(0, image.computeSum(0, 0, 10, 1, 1, 1), 0);
		Assert.assertEquals(0, image.computeSum(0, 0, -10, 1, 1, 1), 0);

		// Larger slices
		canComputeSum(3, 4, 5, 6, 7, 8);
		canComputeSum(3, 4, 5, 1, 1, 1);
		canComputeSum(0, 0, 0, 1, 2, 3);
	}

	private void canComputeSum(int x, int y, int z, int w, int h, int d)
	{
		// This test assumes that crop works!
		for (int pad : new int[] { 0, 1 })
		{
			Image3D image = createData(x + w + pad, y + h + pad, z + d + pad);

			Image3D croppedData = image.crop(x, y, z, w, h, d);
			double e = sum(croppedData);
			double o = image.computeSum(x, y, z, w, h, d);

			Assert.assertEquals(o, e, 0);
		}
	}

	@Test
	public void canComputeRollingSumTable()
	{
		int w = 2, h = 3, d = 4;
		Image3D image = createData(w, h, d);
		double[] table = image.computeRollingSumTable(null);
		for (int dd = 1, i = 0; dd <= d; dd++)
			for (int hh = 1; hh <= h; hh++)
				for (int ww = 1; ww <= w; ww++)
				{
					double e = image.computeSum(0, 0, 0, ww, hh, dd);
					double o = table[i++];
					Assert.assertEquals(e, o, 1e-3);
				}
	}

	@Test
	public void canComputeSumUsingTable()
	{
		// Bounds checks
		Image3D image = createData(2, 2, 2);
		double[] table = image.computeRollingSumTable(null);
		//testComputeSum(36, image,table, 0, 0, 0, 2, 2, 2);
		testComputeSum(36, image, table, 0, 0, 0, 5, 7, 9);
		testComputeSum(0, image, table, 0, 0, 0, 0, 0, 0);
		testComputeSum(1, image, table, 0, 0, 0, 1, 1, 1);
		testComputeSum(2, image, table, 1, 0, 0, 1, 1, 1);
		testComputeSum(0, image, table, -10, 0, 0, 1, 1, 1);
		testComputeSum(0, image, table, 10, 0, 0, 1, 1, 1);
		testComputeSum(0, image, table, 0, 10, 0, 1, 1, 1);
		testComputeSum(0, image, table, 0, -10, 0, 1, 1, 1);
		testComputeSum(0, image, table, 0, 0, 10, 1, 1, 1);
		testComputeSum(0, image, table, 0, 0, -10, 1, 1, 1);

		// Larger slices
		canComputeSumUsingTable(3, 4, 5, 6, 7, 8);
		canComputeSumUsingTable(3, 4, 5, 1, 1, 1);
		canComputeSumUsingTable(0, 0, 0, 1, 2, 3);
	}

	private void canComputeSumUsingTable(int x, int y, int z, int w, int h, int d)
	{
		// This test assumes that crop works!
		for (int pad : new int[] { 0, 1 })
		{
			Image3D image = createData(x + w + pad, y + h + pad, z + d + pad);
			double[] table = image.computeRollingSumTable(null);

			Image3D croppedData = image.crop(x, y, z, w, h, d);
			double e = sum(croppedData);

			testComputeSum(e, image, table, x, y, z, w, h, d);
		}
	}

	private void testComputeSum(double e, Image3D image, double[] table, int x, int y, int z, int w, int h, int d)
	{
		double o = image.computeSum(table, x, y, z, w, h, d);
		double o2 = image.computeSumFast(table, x, y, z, w, h, d);

		// This may be different due to floating point error
		// but we are adding integers so it should be OK
		Assert.assertEquals(e, o, 0);
		Assert.assertEquals(e, o2, 0);
	}

	@Test
	public void canFill()
	{
		canFill(3, 4, 5, 6, 7, 8);
		canFill(3, 4, 5, 1, 1, 1);
		canFill(0, 0, 0, 1, 2, 3);
	}

	private void canFill(int x, int y, int z, int w, int h, int d)
	{
		// This test assumes that copy, crop and insert work!
		for (int pad : new int[] { 0, 1 })
		{
			Image3D image = createEmptyData(x + w + pad, y + h + pad, z + d + pad);
			image.fill(1);
			image.fill(x, y, z, w, h, d, 2);

			Image3D image2 = createEmptyData(x + w + pad, y + h + pad, z + d + pad);
			image2.fill(1);
			Image3D blank = createEmptyData(w, h, d);
			blank.fill(2);
			image2.insert(x, y, z, blank);

			Image3D croppedData = image.crop(x, y, z, w, h, d);

			assertEquals(croppedData, blank);

			assertEquals(image, image2);
		}
	}

	@Test
	public void canFillOutside()
	{
		canFillOutside(3, 4, 5, 6, 7, 8);
		canFillOutside(3, 4, 5, 1, 1, 1);
		canFillOutside(0, 0, 0, 1, 2, 3);
	}

	private void canFillOutside(int x, int y, int z, int w, int h, int d)
	{
		// This test assumes that copy, crop and insert work!
		for (int pad : new int[] { 0, 1 })
		{
			Image3D image = createEmptyData(x + w + pad, y + h + pad, z + d + pad);
			image.fill(1);
			image.fillOutside(x, y, z, w, h, d, 2);

			Image3D image2 = createEmptyData(x + w + pad, y + h + pad, z + d + pad);
			image2.fill(2);
			Image3D blank = createEmptyData(w, h, d);
			blank.fill(1);
			image2.insert(x, y, z, blank);

			Image3D croppedData = image.crop(x, y, z, w, h, d);

			assertEquals(croppedData, blank);

			assertEquals(image, image2);
		}
	}

	private static void assertEquals(Image3D a, Image3D b)
	{
		for (int i = a.getDataLength(); i-- > 0;)
			Assert.assertTrue(null, a.get(i) == b.get(i));
	}

	private static int findMinIndex(Image3D image)
	{
		int min = 0;
		for (int i = image.getDataLength(); i-- > 0;)
			if (image.get(i) < image.get(min))
				min = i;
		return min;
	}

	private static int findMaxIndex(Image3D image)
	{
		int max = 0;
		for (int i = image.getDataLength(); i-- > 0;)
			if (image.get(i) > image.get(max))
				max = i;
		return max;
	}

	private static double sum(Image3D image)
	{
		double sum = 0;
		for (int i = image.getDataLength(); i-- > 0;)
			sum += image.get(i);
		return sum;
	}

}
