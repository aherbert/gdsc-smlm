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

import java.util.Arrays;

import org.jtransforms.fft.FloatFFT_2D;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.FloatEquality;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.smlm.function.StandardFloatValueProcedure;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import gdsc.test.TestAssert;
import ij.process.FHT2;
import ij.process.FloatProcessor;

public class FloatDHT2DTest
{
	int size = 16;
	double centre = (size - 1) / 2.0;

	private FloatDHT2D createData(double cx, double cy)
	{
		Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, size,
				GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);
		double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		a[Gaussian2DFunction.SIGNAL] = 1;
		a[Gaussian2DFunction.X_POSITION] = cx;
		a[Gaussian2DFunction.Y_POSITION] = cy;
		a[Gaussian2DFunction.X_SD] = 1.2;
		a[Gaussian2DFunction.Y_SD] = 1.1;
		StandardFloatValueProcedure p = new StandardFloatValueProcedure();
		p.getValues(f, a);
		return new FloatDHT2D(size, size, p.values, false);
	}

	private FloatDHT2D createData()
	{
		return createData(centre, centre);
	}

	private FloatDHT2D createQuadrants(int w, int h)
	{
		return new FloatDHT2D(createQuadrantsProcessor(w, h));
	}

	static FloatProcessor createQuadrantsProcessor(int w, int h)
	{
		int w_2 = w / 2;
		int h_2 = h / 2;
		FloatProcessor fp = new FloatProcessor(w, h);
		fill(fp, w_2, 0, w_2, h_2, 1);
		fill(fp, 0, 0, w_2, h_2, 2);
		fill(fp, 0, h_2, w_2, h_2, 3);
		fill(fp, w_2, h_2, w_2, h_2, 4);
		return fp;
	}

	static void fill(FloatProcessor fp, int x, int y, int w, int h, double value)
	{
		fp.setRoi(x, y, w, h);
		fp.setValue(value);
		fp.fill();
	}

	@Test
	public void canSwapQuadrants()
	{
		FloatDHT2D dht;

		// Simple test
		float[] data = new float[] { 2, 1, 3, 4 };
		dht = new FloatDHT2D(2, 2, data.clone(), false);
		dht.swapQuadrants();
		checkQuadrants(data, dht.getData());

		int[] test = new int[] { 2, 4, 6 };
		for (int w : test)
			for (int h : test)
			{
				dht = createQuadrants(w, h);

				float[] in = dht.getData().clone();

				dht.swapQuadrants();

				float[] o = dht.getData();

				checkQuadrants(in, o);
			}
	}

	private void checkQuadrants(float[] in, float[] out)
	{
		int[] swap = new int[9];
		swap[1] = 3;
		swap[2] = 4;
		swap[3] = 1;
		swap[4] = 2;
		for (int i = 0; i < in.length; i++)
			Assert.assertEquals(in[i], swap[(int) out[i]], 0);
	}

	@Test
	public void canConvolveAndDeconvolve()
	{
		FloatDHT2D dht = createData();

		float[] pixels = dht.getData().clone();
		dht.transform();

		FloatDHT2D copy = dht.copy();
		copy.initialiseFastMultiply();

		FloatDHT2D convolved = dht.multiply(dht);
		FloatDHT2D deconvolved = convolved.divide(dht);

		FloatDHT2D convolved2 = dht.multiply(copy);
		FloatDHT2D deconvolved2 = convolved.divide(copy);

		Assert.assertArrayEquals(convolved.getData(), convolved2.getData(), 0);
		Assert.assertArrayEquals(deconvolved.getData(), deconvolved2.getData(), 0);

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
		FloatDHT2D dht = createData();
		dht.transform();

		FloatDHT2D copy = dht.copy();
		copy.initialiseFastMultiply();

		// Centre of power spectrum
		int icentre = size / 2;

		for (int y = -1; y <= 1; y++)
			for (int x = -1; x <= 1; x++)
			{
				FloatDHT2D dht2 = createData(centre + x, centre + y);
				dht2.transform();

				FloatDHT2D correlation = dht2.conjugateMultiply(dht);
				FloatDHT2D correlation2 = dht2.conjugateMultiply(copy);
				Assert.assertArrayEquals(correlation.getData(), correlation2.getData(), 0);

				correlation.inverseTransform();
				correlation.swapQuadrants();

				float[] pixels = correlation.getData();

				int i = SimpleArrayUtils.findMaxIndex(pixels);
				int[] xy = correlation.getXY(i);

				// This is how far dht has to move to align with dht2.
				// To align dht2 with dht would be the opposite sign.
				int ox = xy[0] - icentre;
				int oy = xy[1] - icentre;
				//System.out.printf("Shift [%d,%d], centre [%d,%d]\n", x, y, xy[0], xy[1]);
				Assert.assertEquals(x, ox);
				Assert.assertEquals(y, oy);
			}
	}

	@Test
	public void canConvertToDFT()
	{
		FloatDHT2D dht = createData();
		float[] input = dht.getData().clone();
		dht.transform();

		FloatImage2D[] result = dht.toDFT(null, null);

		float rel = 1e-6f;
		float abs = 1e-6f;

		// Test reverse transform
		FloatDHT2D dht2 = FloatDHT2D.fromDFT(result[0], result[1], null);

		float[] e = dht.getData();
		float[] o = dht2.getData();
		for (int i = 0; i < e.length; i++)
			Assert.assertTrue(FloatEquality.almostEqualRelativeOrAbsolute(e[i], o[i], rel, abs));

		// Test verses full forward transform
		FloatFFT_2D fft = new FloatFFT_2D(dht.nr, dht.nc);
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

	@Test
	public void canMatchFHT2()
	{
		FloatDHT2D dht = createData();
		FloatDHT2D dht2 = createData(centre + 1, centre + 1);

		float[] pixels = dht.getData().clone();
		float[] pixels2 = dht2.getData().clone();

		FHT2 fht = new FHT2(pixels, size, false);
		FHT2 fht2 = new FHT2(pixels2, size, false);

		dht.transform();
		dht2.transform();
		fht.transform();
		fht2.transform();

		Assert.assertArrayEquals(fht.getData(), dht.getData(), 1e-6f);
		Assert.assertArrayEquals(fht2.getData(), dht2.getData(), 1e-6f);

		check("multiply", dht.multiply(dht2), fht.multiply(fht2), 1e-6f, 1e-6f);
		check("conjugateMultiply", dht.conjugateMultiply(dht2), fht.conjugateMultiply(fht2), 1e-6f, 1e-6f);
		check("divide", dht.divide(dht2), fht.divide(fht2), 1e-3f, 1e-3f);
	}

	private void check(String operation, FloatDHT2D dht, FHT2 fht, float rel, float abs)
	{
		float[] e = fht.getData();
		float[] o = dht.getData();
		for (int i = 0; i < e.length; i++)
			if (!FloatEquality.almostEqualRelativeOrAbsolute(e[i], o[i], rel, abs))
			{
				TestAssert.fail("%s [%d] %g vs %g = %g", operation, i, e[i], o[i],
						FloatEquality.relativeError(e[i], o[i]));
			}
	}
}
