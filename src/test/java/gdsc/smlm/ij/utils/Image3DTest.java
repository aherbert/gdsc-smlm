package gdsc.smlm.ij.utils;

import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.Maths;
import gdsc.core.utils.SimpleArrayUtils;
import ij.ImageStack;
import ij.process.ImageProcessor;

public class Image3DTest
{
	private Image3D createData(int w, int h, int d)
	{
		float[] data = SimpleArrayUtils.newArray(w * h * d, 0, 1f);
		return new Image3D(w, h, d, data);
	}

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

		Image3D croppedData = image.crop(x, y, z, w, h, d, null);
		Assert.assertEquals(croppedData.getWidth(), w);
		Assert.assertEquals(croppedData.getHeight(), h);
		Assert.assertEquals(croppedData.getSize(), d);

		Image3D croppedData2 = Image3D.crop(stack, x, y, z, w, h, d, null);
		Assert.assertArrayEquals(croppedData.getData(), croppedData2.getData(), 0);

		ImageStack croppedStack = image.cropToStack(x, y, z, w, h, d);
		Assert.assertEquals(croppedStack.getWidth(), w);
		Assert.assertEquals(croppedStack.getHeight(), h);
		Assert.assertEquals(croppedStack.getSize(), d);

		float[] croppedStackData = new Image3D(croppedStack).getData();
		Assert.assertArrayEquals(croppedData.getData(), croppedStackData, 0);

		ImageStack croppedStack2 = Image3D.cropToStack(stack, x, y, z, w, h, d);
		Assert.assertArrayEquals(croppedData.getData(), new Image3D(croppedStack2).getData(), 0);

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
			System.arraycopy(croppedData.getData(), zz * o.length, o, 0, o.length);
			Assert.assertArrayEquals(e, o, 0);
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

			Image3D blank = new Image3D(w, h, d);

			image.insert(x, y, z, blank);

			Image3D croppedData = image.crop(x, y, z, w, h, d, null);

			Assert.assertArrayEquals(croppedData.getData(), blank.getData(), 0);

			ImageStack blankStack = blank.getImageStack();
			image2.insert(x, y, z, blankStack);

			croppedData = image2.crop(x, y, z, w, h, d, null);

			Assert.assertArrayEquals(croppedData.getData(), blank.getData(), 0);

			for (int s = 0; s < blankStack.getSize(); s++)
				image3.insert(x, y, z + s, blankStack.getProcessor(1 + s));

			croppedData = image3.crop(x, y, z, w, h, d, null);

			Assert.assertArrayEquals(croppedData.getData(), blank.getData(), 0);
		}
	}
	
	@Test
	public void canComputeSum()
	{
		// Bounds checks
		Image3D image = createData(2, 2, 2);
		Assert.assertEquals(28, image.computeSum(0, 0, 0, 2, 2, 2), 0);
		Assert.assertEquals(0, image.computeSum(0, 0, 0, 0, 0, 0), 0);
		Assert.assertEquals(1, image.computeSum(1, 0, 0, 1, 1, 1), 0);
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
		// This test assumes that crop work!
		for (int pad : new int[] { 0, 1 })
		{
			Image3D image = createData(x + w + pad, y + h + pad, z + d + pad);

			Image3D croppedData = image.crop(x, y, z, w, h, d, null);
			double e = Maths.sum(croppedData.getData());
			double o = image.computeSum(x, y, z, w, h, d);
			
			Assert.assertEquals(o, e, 0);
		}
	}	
}
