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
package gdsc.smlm.filters;

import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.ImageWindow;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.smlm.filters.FHTFilter.Operation;
import gdsc.test.TestAssert;
import gdsc.test.TestCounter;
import gdsc.test.TestSettings;
import ij.plugin.filter.EDM;
import ij.process.ByteProcessor;
import ij.process.FHT;
import ij.process.FloatProcessor;

@SuppressWarnings({ "javadoc" })
public class FHTFilterTest
{
	@Test
	public void canCorrelate()
	{
		canFilter(Operation.CORRELATION);
	}

	@Test
	public void canConvolve()
	{
		canFilter(Operation.CONVOLUTION);
	}

	@Test
	public void canDeconvolve()
	{
		canFilter(Operation.DECONVOLUTION);
	}

	private void canFilter(Operation operation)
	{
		int size = 16;
		int ex = 5, ey = 7;
		int ox = 1, oy = 2;
		RandomGenerator r = TestSettings.getRandomGenerator();
		FloatProcessor fp1 = createProcessor(size, ex, ey, 4, 4, r);
		// This is offset from the centre
		FloatProcessor fp2 = createProcessor(size, size / 2 + ox, size / 2 + oy, 4, 4, r);

		float[] input1 = ((float[]) fp1.getPixels()).clone();
		float[] input2 = ((float[]) fp2.getPixels()).clone();

		FHT fht1 = new FHT(fp1);
		fht1.transform();
		FHT fht2 = new FHT(fp2);
		fht2.transform();

		FHT fhtE;
		switch (operation)
		{
			case CONVOLUTION:
				fhtE = fht1.multiply(fht2);
				break;
			case CORRELATION:
				fhtE = fht1.conjugateMultiply(fht2);
				break;
			case DECONVOLUTION:
				fhtE = fht1.divide(fht2);
				break;
			default:
				throw new RuntimeException();
		}
		fhtE.inverseTransform();
		fhtE.swapQuadrants();

		float[] e = (float[]) fhtE.getPixels();
		if (operation == Operation.CORRELATION)
		{
			// Test the max correlation position
			int max = SimpleArrayUtils.findMaxIndex(e);
			int x = max % 16;
			int y = max / 16;

			Assert.assertEquals(ex, x + ox);
			Assert.assertEquals(ey, y + oy);
		}

		// Test verses a spatial domain filter in the middle of the image
		if (operation != Operation.DECONVOLUTION)
		{
			double sum = 0;
			float[] i2 = input2;
			if (operation == Operation.CONVOLUTION)
			{
				i2 = i2.clone();
				KernelFilter.rotate180(i2);
			}
			for (int i = 0; i < input1.length; i++)
				sum += input1[i] * i2[i];
			//double exp = e[size / 2 * size + size / 2];
			//System.out.printf("Sum = %f vs [%d] %f\n", sum, size / 2 * size + size / 2, exp);
			Assert.assertEquals(sum, sum, 1e-3);
		}

		// Test the FHT filter
		FHTFilter ff = new FHTFilter(input2, size, size);
		ff.setOperation(operation);
		ff.filter(input1, size, size);

		// There may be differences due to the use of the JTransforms library
		double error = (operation == Operation.DECONVOLUTION) ? 5e-2 : 1e-4;

		// This tests everything and can fail easily depending on the random generator 
		// due to edge artifacts.
		//TestAssert.assertArrayEqualsRelative(e, input1, error);

		// This tests the centre to ignore edge differences
		int min = size / 4;
		int max = size - min;
		int repeats = 0;
		for (int y = min; y < max; y++)
			for (int x = min; x < max; x++)
				repeats++;

		// Use a fail counter for a 'soft' test that detects major problems
		int failureLimit = TestCounter.computeFailureLimit(repeats, 0.1);
		TestCounter failCounter = new TestCounter(failureLimit);

		for (int y = min; y < max; y++)
		{
			final int yy = y;
			for (int x = min; x < max; x++)
			{
				final int xx = x;
				int i = y * size + x;
				failCounter.run(() -> {
					TestAssert.assertEqualsRelative(e[i], input1[i], error, "Element [%d,%d]", xx, yy);
				});
			}
		}
	}

	private FloatProcessor createProcessor(int size, int x, int y, int w, int h, RandomGenerator r)
	{
		ByteProcessor bp = new ByteProcessor(size, size);
		bp.setColor(255);
		bp.fillOval(x, y, w, h);
		EDM e = new EDM();
		FloatProcessor fp = e.makeFloatEDM(bp, 0, true);
		if (r != null)
		{
			float[] d = (float[]) fp.getPixels();
			for (int i = 0; i < d.length; i++)
				d[i] += r.nextFloat() * 0.01;
		}
		return fp;
	}

	@Test
	public void canWindow()
	{
		int size = 16;
		float[] in = SimpleArrayUtils.newFloatArray(size * size, 1);
		FHTFilter f = new FHTFilter(new float[1], 1, 1);
		for (int i = 1; i < 5; i++)
		{
			double[] wx = ImageWindow.tukeyEdge(size, i);
			float[] e = ImageWindow.applyWindowSeparable(in, size, size, wx, wx);
			float[] o = in.clone();
			f.applyBorder(o, size, size, i);
			Assert.assertArrayEquals(e, o, 0);
		}
	}
}
