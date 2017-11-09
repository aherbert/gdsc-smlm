package gdsc.smlm.ij.utils;

import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.FloatEquality;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.smlm.function.StandardFloatValueProcedure;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.QuadraticAstigmatismZModel;
import gdsc.smlm.ij.utils.DHT3D;
import gdsc.smlm.ij.utils.Image3D;
import ij.ImageStack;
import ij.process.ImageProcessor;

public class DHT3DTest
{
	int size = 16;
	double centre = (size - 1) / 2.0;

	final static double gamma = 2;
	final static int zDepth = 5;
	protected QuadraticAstigmatismZModel zModel = new QuadraticAstigmatismZModel(gamma, zDepth);

	private DHT3D createData(double cx, double cy, double cz)
	{
		Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, size, GaussianFunctionFactory.FIT_ASTIGMATISM,
				zModel);
		int length = size * size;
		float[] data = new float[size * length];
		double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		a[Gaussian2DFunction.SIGNAL] = 1;
		a[Gaussian2DFunction.X_POSITION] = cx;
		a[Gaussian2DFunction.Y_POSITION] = cy;
		a[Gaussian2DFunction.X_SD] = 1;
		a[Gaussian2DFunction.Y_SD] = 1;
		StandardFloatValueProcedure p = new StandardFloatValueProcedure();
		for (int z = 0; z < size; z++)
		{
			a[Gaussian2DFunction.Z_POSITION] = z - cz;
			p.getValues(f, a, data, z * length);
		}
		return new DHT3D(size, size, size, data, false);
	}

	private DHT3D createData()
	{
		return createData(centre, centre, centre);
	}

	@Test
	public void canSwapOctants()
	{
		// This just tests that the swap of the DHT and the stack matches
		DHT3D dht = createData();
		ImageStack stack = dht.getImageStack();
		//gdsc.core.ij.Utils.display("Test", stack);
		dht.swapOctants();
		DHT3D.swapOctants(stack);

		float[] e = new DHT3D(stack).getData();
		float[] o = dht.getData();

		Assert.assertArrayEquals(e, o, 0);
	}

	@Test
	public void canConvolveAndDeconvolve()
	{
		DHT3D dht = createData();
		float[] pixels = dht.getData().clone();
		dht.transform();

		DHT3D convolved = dht.multiply(dht);
		DHT3D deconvolved = convolved.divide(dht);

		float[] e = dht.getData();
		float[] o = deconvolved.getData();
		for (int i = 0; i < e.length; i++)
			Assert.assertTrue(FloatEquality.almostEqualRelativeOrAbsolute(e[i], o[i], 1e-6f, 1e-6f));

		deconvolved.inverseTransform();

		// Test after reverse transform
		e = pixels;
		o = deconvolved.getData();

		for (int i = 0; i < e.length; i++)
			Assert.assertTrue(FloatEquality.almostEqualRelativeOrAbsolute(e[i], o[i], 1e-6f, 1e-6f));
	}

	@Test
	public void canCorrelate()
	{
		DHT3D dht = createData();
		dht.transform();

		// Centre of power spectrum
		int icentre = size / 2;

		for (int z = -1; z <= 1; z++)
			for (int y = -1; y <= 1; y++)
				for (int x = -1; x <= 1; x++)
				{
					DHT3D dht2 = createData(centre + x, centre + y, centre + z);
					dht2.transform();

					DHT3D correlation = dht2.conjugateMultiply(dht);
					correlation.inverseTransform();
					correlation.swapOctants();

					float[] pixels = correlation.getData();

					int i = SimpleArrayUtils.findMaxIndex(pixels);
					int[] xyz = correlation.getXYZ(i);

					// This is how far dht has to move to align with dht2.
					// To align dht2 with dht would be the opposite sign.
					int ox = xyz[0] - icentre;
					int oy = xyz[1] - icentre;
					int oz = xyz[2] - icentre;
					//System.out.printf("Shift [%d,%d,%d], centre [%d,%d,%d]\n", x, y, z, xyz[0], xyz[1], xyz[2]);
					Assert.assertEquals(x, ox);
					Assert.assertEquals(y, oy);
					Assert.assertEquals(z, oz);
				}
	}

	@Test
	public void canCrop()
	{
		DHT3D dht = createData();
		int x = 3, y = 4, z = 5;
		int w = 6, h = 7, d = 8;
		Image3D croppedData = dht.crop(x, y, z, w, h, d, null);
		Assert.assertEquals(croppedData.getWidth(), w);
		Assert.assertEquals(croppedData.getHeight(), h);
		Assert.assertEquals(croppedData.getSize(), d);
		ImageStack croppedStack = dht.cropToStack(x, y, z, w, h, d);
		Assert.assertEquals(croppedStack.getWidth(), w);
		Assert.assertEquals(croppedStack.getHeight(), h);
		Assert.assertEquals(croppedStack.getSize(), d);
		float[] croppedStackData = new Image3D(croppedStack).getData();

		Assert.assertArrayEquals(croppedData.getData(), croppedStackData, 0);

		// Test it is the correct region
		ImageStack originalStack = dht.getImageStack();
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
