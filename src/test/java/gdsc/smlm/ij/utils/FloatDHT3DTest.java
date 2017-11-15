package gdsc.smlm.ij.utils;

import java.util.Arrays;

import org.jtransforms.fft.FloatFFT_3D;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.FloatEquality;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.smlm.function.StandardFloatValueProcedure;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.smlm.function.gaussian.QuadraticAstigmatismZModel;
import ij.ImageStack;
import ij.process.FloatProcessor;

public class FloatDHT3DTest
{
	int size = 16;
	double centre = (size - 1) / 2.0;

	final static double gamma = 2;
	final static int zDepth = 5;
	protected QuadraticAstigmatismZModel zModel = new QuadraticAstigmatismZModel(gamma, zDepth);

	private FloatDHT3D createData(double cx, double cy, double cz)
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
		return new FloatDHT3D(size, size, size, data, false);
	}

	private FloatDHT3D createData()
	{
		return createData(centre, centre, centre);
	}

	private FloatDHT3D createOctants(int w, int h, int d)
	{
		return new FloatDHT3D(createOctantsStack(w, h, d));
	}

	static ImageStack createOctantsStack(int w, int h, int d)
	{
		int w_2 = w / 2;
		int h_2 = h / 2;
		int d_2 = d / 2;
		ImageStack stack = new ImageStack(w, h, d);
		FloatProcessor fp = new FloatProcessor(w, h);
		float[] pixels = (float[]) fp.getPixels();
		fill(fp, w_2, 0, w_2, h_2, 1);
		fill(fp, 0, 0, w_2, h_2, 2);
		fill(fp, 0, h_2, w_2, h_2, 3);
		fill(fp, w_2, h_2, w_2, h_2, 4);
		for (int z = 0; z < d_2; z++)
			stack.setPixels(pixels.clone(), 1 + z);
		fill(fp, w_2, 0, w_2, h_2, 5);
		fill(fp, 0, 0, w_2, h_2, 6);
		fill(fp, 0, h_2, w_2, h_2, 7);
		fill(fp, w_2, h_2, w_2, h_2, 8);
		for (int z = d_2; z < d; z++)
			stack.setPixels(pixels.clone(), 1 + z);
		return stack;
	}

	static void fill(FloatProcessor fp, int x, int y, int w, int h, double value)
	{
		fp.setRoi(x, y, w, h);
		fp.setValue(value);
		fp.fill();
	}

	@Test
	public void canSwapOctants()
	{
		FloatDHT3D dht;

		// Simple test
		float[] data = new float[] { 2, 1, 3, 4, 6, 5, 7, 8 };
		dht = new FloatDHT3D(2, 2, 2, data.clone(), false);
		dht.swapOctants();
		checkOctants(data, dht.getData());

		int[] test = new int[] { 2, 4, 6 };
		for (int w : test)
			for (int h : test)
				for (int d : test)
				{
					dht = createOctants(w, h, d);

					// This just tests that the swap of the DHT and the stack matches
					ImageStack stack = dht.getImageStack();
					//gdsc.core.ij.Utils.display("Test", stack);
					dht.swapOctants();
					FloatDHT3D.swapOctants(stack);

					float[] e = new FloatDHT3D(stack).getData();
					float[] o = dht.getData();

					Assert.assertArrayEquals(e, o, 0);
				}
	}

	private void checkOctants(float[] in, float[] out)
	{
		int[] swap = new int[9];
		swap[1] = 7;
		swap[2] = 8;
		swap[3] = 5;
		swap[4] = 6;
		swap[5] = 3;
		swap[6] = 4;
		swap[7] = 1;
		swap[8] = 2;
		for (int i = 0; i < in.length; i++)
			Assert.assertEquals(in[i], swap[(int) out[i]], 0);
	}

	@Test
	public void canConvolveAndDeconvolve()
	{
		FloatDHT3D dht = createData();
		float[] pixels = dht.getData().clone();
		dht.transform();

		FloatDHT3D convolved = dht.multiply(dht);
		FloatDHT3D deconvolved = convolved.divide(dht);

		float[] e = dht.getData();
		float[] o = deconvolved.getData();
		for (int i = 0; i < e.length; i++)
			Assert.assertTrue(FloatEquality.almostEqualRelativeOrAbsolute(e[i], o[i], 1e-6f, 1e-6f));

		deconvolved.inverseTransform();

		// Test after reverse transform
		e = pixels;
		o = deconvolved.getData();

		for (int i = 0; i < e.length; i++)
			Assert.assertTrue(FloatEquality.almostEqualRelativeOrAbsolute(e[i], o[i], 1e-7f, 1e-7f));
	}

	@Test
	public void canCorrelate()
	{
		FloatDHT3D dht = createData();
		dht.transform();

		// Centre of power spectrum
		int icentre = size / 2;

		for (int z = -1; z <= 1; z++)
			for (int y = -1; y <= 1; y++)
				for (int x = -1; x <= 1; x++)
				{
					FloatDHT3D dht2 = createData(centre + x, centre + y, centre + z);
					dht2.transform();

					FloatDHT3D correlation = dht2.conjugateMultiply(dht);
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
	public void canConvertToDFT()
	{
		FloatDHT3D dht = createData();
		float[] input = dht.getData().clone();
		dht.transform();

		FloatImage3D[] result = dht.toDFT(null, null);

		float rel = 1e-6f;
		float abs = 1e-6f;

		// Test reverse transform
		FloatDHT3D dht2 = FloatDHT3D.fromDFT(result[0], result[1], null);

		float[] e = dht.getData();
		float[] o = dht2.getData();
		for (int i = 0; i < e.length; i++)
			Assert.assertTrue(FloatEquality.almostEqualRelativeOrAbsolute(e[i], o[i], rel, abs));

		// Test verses full forward transform
		FloatFFT_3D fft = new FloatFFT_3D(dht.ns, dht.nr, dht.nc);
		float[] dft = Arrays.copyOf(input, 2 * e.length);
		fft.realForwardFull(dft);

		float[] or = result[0].getData();
		float[] oi = result[1].getData();
		for (int i = 0, j = 0; i < dft.length; i += 2, j++)
		{
			Assert.assertTrue(FloatEquality.almostEqualRelativeOrAbsolute(dft[i], or[j], rel, abs));
			Assert.assertTrue(FloatEquality.almostEqualRelativeOrAbsolute(dft[i + 1], oi[j], rel, abs));
		}
	}
}
