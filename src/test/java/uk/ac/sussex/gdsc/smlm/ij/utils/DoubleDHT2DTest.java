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
package uk.ac.sussex.gdsc.smlm.ij.utils;

import java.util.Arrays;

import org.jtransforms.fft.DoubleFFT_2D;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;import uk.ac.sussex.gdsc.test.junit5.SeededTest;import uk.ac.sussex.gdsc.test.junit5.RandomSeed;import uk.ac.sussex.gdsc.test.junit5.SpeedTag;

import ij.process.FloatProcessor;
import uk.ac.sussex.gdsc.core.ij.process.FHT2;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.function.StandardValueProcedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.test.junit5.ExtraAssertions;

@SuppressWarnings({ "javadoc" })
public class DoubleDHT2DTest
{
	static int size = 16;
	static double centre = (size - 1) / 2.0;

	private static DoubleDHT2D createData(double cx, double cy)
	{
		final Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, size,
				GaussianFunctionFactory.FIT_ERF_FREE_CIRCLE, null);
		final double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
		a[Gaussian2DFunction.SIGNAL] = 1;
		a[Gaussian2DFunction.X_POSITION] = cx;
		a[Gaussian2DFunction.Y_POSITION] = cy;
		a[Gaussian2DFunction.X_SD] = 1.2;
		a[Gaussian2DFunction.Y_SD] = 1.1;
		final StandardValueProcedure p = new StandardValueProcedure();
		p.getValues(f, a);
		return new DoubleDHT2D(size, size, p.values, false);
	}

	private static DoubleDHT2D createData()
	{
		return createData(centre, centre);
	}

	private static DoubleDHT2D createQuadrants(int w, int h)
	{
		return new DoubleDHT2D(createQuadrantsProcessor(w, h));
	}

	static FloatProcessor createQuadrantsProcessor(int w, int h)
	{
		final int w_2 = w / 2;
		final int h_2 = h / 2;
		final FloatProcessor fp = new FloatProcessor(w, h);
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
		DoubleDHT2D dht;

		// Simple test
		final double[] data = new double[] { 2, 1, 3, 4 };
		dht = new DoubleDHT2D(2, 2, data.clone(), false);
		dht.swapQuadrants();
		checkQuadrants(data, dht.getData());

		final int[] test = new int[] { 2, 4, 6 };
		for (final int w : test)
			for (final int h : test)
			{
				dht = createQuadrants(w, h);

				final double[] in = dht.getData().clone();

				dht.swapQuadrants();

				final double[] o = dht.getData();

				checkQuadrants(in, o);
			}
	}

	private static void checkQuadrants(double[] in, double[] out)
	{
		final int[] swap = new int[9];
		swap[1] = 3;
		swap[2] = 4;
		swap[3] = 1;
		swap[4] = 2;
		for (int i = 0; i < in.length; i++)
			Assertions.assertEquals(in[i], swap[(int) out[i]]);
	}

	@Test
	public void canConvolveAndDeconvolve()
	{
		final DoubleDHT2D dht = createData();

		final double[] pixels = dht.getData().clone();
		dht.transform();

		final DoubleDHT2D copy = dht.copy();
		copy.initialiseFastMultiply();

		final DoubleDHT2D convolved = dht.multiply(dht);
		final DoubleDHT2D deconvolved = convolved.divide(dht);

		final DoubleDHT2D convolved2 = dht.multiply(copy);
		final DoubleDHT2D deconvolved2 = convolved.divide(copy);

		Assertions.assertArrayEquals(convolved.getData(), convolved2.getData());
		Assertions.assertArrayEquals(deconvolved.getData(), deconvolved2.getData());

		double[] e = dht.getData();
		double[] o = deconvolved.getData();
		for (int i = 0; i < e.length; i++)
			Assertions.assertTrue(DoubleEquality.almostEqualRelativeOrAbsolute(e[i], o[i], 1e-6, 1e-6));

		deconvolved.inverseTransform();

		// Test after reverse transform
		e = pixels;
		o = deconvolved.getData();

		for (int i = 0; i < e.length; i++)
			Assertions.assertTrue(DoubleEquality.almostEqualRelativeOrAbsolute(e[i], o[i], 1e-8, 1e-8));
	}

	@Test
	public void canCorrelate()
	{
		final DoubleDHT2D dht = createData();
		dht.transform();

		final DoubleDHT2D copy = dht.copy();
		copy.initialiseFastMultiply();

		// Centre of power spectrum
		final int icentre = size / 2;

		for (int y = -1; y <= 1; y++)
			for (int x = -1; x <= 1; x++)
			{
				final DoubleDHT2D dht2 = createData(centre + x, centre + y);
				dht2.transform();

				final DoubleDHT2D correlation = dht2.conjugateMultiply(dht);
				final DoubleDHT2D correlation2 = dht2.conjugateMultiply(copy);
				Assertions.assertArrayEquals(correlation.getData(), correlation2.getData());

				correlation.inverseTransform();
				correlation.swapQuadrants();

				final double[] pixels = correlation.getData();

				final int i = SimpleArrayUtils.findMaxIndex(pixels);
				final int[] xy = correlation.getXY(i);

				// This is how far dht has to move to align with dht2.
				// To align dht2 with dht would be the opposite sign.
				final int ox = xy[0] - icentre;
				final int oy = xy[1] - icentre;
				//System.out.printf("Shift [%d,%d], centre [%d,%d]\n", x, y, xy[0], xy[1]);
				Assertions.assertEquals(x, ox);
				Assertions.assertEquals(y, oy);
			}
	}

	@Test
	public void canConvertToDFT()
	{
		final DoubleDHT2D dht = createData();
		final double[] input = dht.getData().clone();
		dht.transform();

		final DoubleImage2D[] result = dht.toDFT(null, null);

		final double rel = 1e-6;
		final double abs = 1e-6;

		// Test reverse transform
		final DoubleDHT2D dht2 = DoubleDHT2D.fromDFT(result[0], result[1], null);

		final double[] e = dht.getData();
		final double[] o = dht2.getData();
		for (int i = 0; i < e.length; i++)
			Assertions.assertTrue(DoubleEquality.almostEqualRelativeOrAbsolute(e[i], o[i], rel, abs));

		// Test verses full forward transform
		final DoubleFFT_2D fft = new DoubleFFT_2D(dht.nr, dht.nc);
		final double[] dft = Arrays.copyOf(input, 2 * e.length);
		fft.realForwardFull(dft);

		final double[] or = result[0].getData();
		final double[] oi = result[1].getData();
		for (int i = 0, j = 0; i < dft.length; i += 2, j++)
		{
			Assertions.assertTrue(DoubleEquality.almostEqualRelativeOrAbsolute(dft[i], or[j], rel, abs));
			Assertions.assertTrue(DoubleEquality.almostEqualRelativeOrAbsolute(dft[i + 1], oi[j], rel, abs));
		}
	}

	@Test
	public void canMatchFHT2()
	{
		final DoubleDHT2D dht = createData();
		final DoubleDHT2D dht2 = createData(centre + 1, centre + 1);

		final float[] pixels = SimpleArrayUtils.toFloat(dht.getData());
		final float[] pixels2 = SimpleArrayUtils.toFloat(dht2.getData());

		final FHT2 fht = new FHT2(pixels, size, false);
		final FHT2 fht2 = new FHT2(pixels2, size, false);

		dht.transform();
		dht2.transform();
		fht.transform();
		fht2.transform();

		check("fht", fht.getData(), dht.getData(), 1e-6, 1e-6);
		check("fht2", fht2.getData(), dht2.getData(), 1e-6, 1e-6);

		check("multiply", dht.multiply(dht2), fht.multiply(fht2), 1e-6, 1e-6);
		check("conjugateMultiply", dht.conjugateMultiply(dht2), fht.conjugateMultiply(fht2), 1e-6, 1e-6);
		check("divide", dht.divide(dht2), fht.divide(fht2), 1e-3, 1e-3);
	}

	private static void check(String msg, float[] e, double[] o, double rel, double abs)
	{
		for (int i = 0; i < e.length; i++)
			if (!DoubleEquality.almostEqualRelativeOrAbsolute(e[i], o[i], rel, abs))
				ExtraAssertions.fail("%s [%d] %g vs %g = %g", msg, i, e[i], o[i], DoubleEquality.relativeError(e[i], o[i]));
	}

	private static void check(String operation, DoubleDHT2D dht, FHT2 fht, double rel, double abs)
	{
		final float[] e = fht.getData();
		final double[] o = dht.getData();
		check(operation, e, o, rel, abs);
	}
}
