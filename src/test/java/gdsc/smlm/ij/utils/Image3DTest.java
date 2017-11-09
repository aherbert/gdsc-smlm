package gdsc.smlm.ij.utils;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import ij.ImageStack;
import ij.process.ImageProcessor;

public class Image3DTest
{
	int size = 16;

	private Image3D createData()
	{
		RandomGenerator r = new Well19937c();
		float[] data = new float[size * size * size];
		for (int i = 0; i < size; i++)
		{
			data[i] = r.nextFloat();
		}
		return new Image3D(size, size, size, data);
	}

	@Test
	public void canCrop()
	{
		Image3D image = createData();
		ImageStack stack = image.getImageStack();
		
		int x = 3, y = 4, z = 5;
		int w = 6, h = 7, d = 8;
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
			croppedData.copySlice(zz, o, 0);
			System.arraycopy(croppedData.getData(), zz * o.length, o, 0, o.length);
			Assert.assertArrayEquals(e, o, 0);
		}
	}
}
