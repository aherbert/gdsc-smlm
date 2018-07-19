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
package uk.ac.sussex.gdsc.smlm.ij.frc;

import java.awt.Rectangle;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathArrays;
import org.junit.Assert;
import org.junit.Test;

import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.ij.results.IJImagePeakResults;
import uk.ac.sussex.gdsc.test.BaseTimingTask;
import uk.ac.sussex.gdsc.test.LogLevel;
import uk.ac.sussex.gdsc.test.TestSettings;
import uk.ac.sussex.gdsc.test.TimingService;
import uk.ac.sussex.gdsc.test.junit4.TestAssume;

@SuppressWarnings({ "javadoc" })
public class FRCTest
{
	@Test
	public void canComputeSine()
	{
		final int steps = 1000;
		final double delta = 2 * Math.PI / steps;
		for (int i = 0; i <= steps; i++)
		{
			final double a = i * delta;
			final double cosA = Math.cos(a);
			final double e = Math.sin(a);
			final double o = FRC.getSine(a, cosA);
			//System.out.printf("%f  %f ?= %f\n", a, e, o);
			Assert.assertTrue(DoubleEquality.almostEqualRelativeOrAbsolute(o, e, 1e-6, 1e-10));
		}
	}

	@Test
	public void canComputeMirrored()
	{
		// Sample lines through an image to create a structure.
		final int size = 1024;
		final double[][] data = new double[size * 2][];
		final RandomGenerator r = TestSettings.getRandomGenerator();
		for (int x = 0, y = 0, y2 = size, i = 0; x < size; x++, y++, y2--)
		{
			data[i++] = new double[] { x + r.nextGaussian() * 5, y + r.nextGaussian() * 5 };
			data[i++] = new double[] { x + r.nextGaussian() * 5, y2 + r.nextGaussian() * 5 };
		}
		// Create 2 images
		final Rectangle bounds = new Rectangle(0, 0, size, size);
		IJImagePeakResults i1 = createImage(bounds);
		IJImagePeakResults i2 = createImage(bounds);
		final int[] indices = SimpleArrayUtils.newArray(data.length, 0, 1);
		MathArrays.shuffle(indices, r);
		for (final int i : indices)
		{
			final IJImagePeakResults image = i1;
			i1 = i2;
			i2 = image;
			image.add((float) data[i][0], (float) data[i][1], 1);
		}
		i1.end();
		i2.end();
		final ImageProcessor ip1 = i1.getImagePlus().getProcessor();
		final ImageProcessor ip2 = i2.getImagePlus().getProcessor();
		// Test
		final FRC frc = new FRC();
		FloatProcessor[] fft1, fft2;
		fft1 = frc.getComplexFFT(ip1);
		fft2 = frc.getComplexFFT(ip2);

		final float[] dataA1 = (float[]) fft1[0].getPixels();
		final float[] dataB1 = (float[]) fft1[1].getPixels();
		final float[] dataA2 = (float[]) fft2[0].getPixels();
		final float[] dataB2 = (float[]) fft2[1].getPixels();

		final float[] numeratorE = new float[dataA1.length];
		final float[] absFFT1E = new float[dataA1.length];
		final float[] absFFT2E = new float[dataA1.length];

		FRC.compute(numeratorE, absFFT1E, absFFT2E, dataA1, dataB1, dataA2, dataB2);

		Assert.assertTrue("numeratorE", FRC.checkSymmetry(numeratorE, size));
		Assert.assertTrue("absFFT1E", FRC.checkSymmetry(absFFT1E, size));
		Assert.assertTrue("absFFT2E", FRC.checkSymmetry(absFFT2E, size));

		final float[] numeratorA = new float[dataA1.length];
		final float[] absFFT1A = new float[dataA1.length];
		final float[] absFFT2A = new float[dataA1.length];
		FRC.computeMirrored(size, numeratorA, absFFT1A, absFFT2A, dataA1, dataB1, dataA2, dataB2);

		//for (int y=0, i=0; y<size; y++)
		//	for (int x=0; x<size; x++, i++)
		//	{
		//		System.out.printf("[%d,%d = %d] %f ?= %f\n", x, y, i, numeratorE[i], numeratorA[i]);
		//	}

		Assert.assertArrayEquals("numerator", numeratorE, numeratorA, 0);
		Assert.assertArrayEquals("absFFT1", absFFT1E, absFFT1A, 0);
		Assert.assertArrayEquals("absFFT2", absFFT2E, absFFT2A, 0);

		FRC.computeMirroredFast(size, numeratorA, absFFT1A, absFFT2A, dataA1, dataB1, dataA2, dataB2);

		// Check this.
		for (int y = 1; y < size; y++)
			for (int x = 1, i = y * size + 1; x < size; x++, i++)
			{
				Assert.assertEquals("numerator", numeratorE[i], numeratorA[i], 0);
				Assert.assertEquals("absFFT1", absFFT1E[i], absFFT1A[i], 0);
				Assert.assertEquals("absFFT2", absFFT2E[i], absFFT2A[i], 0);
			}
	}

	private static IJImagePeakResults createImage(Rectangle bounds)
	{
		final IJImagePeakResults i1 = new IJImagePeakResults("1", bounds, 1);
		i1.setDisplayImage(false);
		i1.begin();
		return i1;
	}

	private abstract class MyTimingTask extends BaseTimingTask
	{
		public MyTimingTask(String name)
		{
			super(name);
		}

		@Override
		public int getSize()
		{
			return 1;
		}

		@Override
		public Object getData(int i)
		{
			return null;
		}
	}

	@Test
	public void computeSineIsFaster()
	{
		TestAssume.assumeHighComplexity();

		final int steps = 100000;
		final double delta = 2 * Math.PI / steps;
		final double[] a = new double[steps + 1];
		final double[] cosA = new double[steps + 1];
		for (int i = 0; i <= steps; i++)
		{
			a[i] = i * delta;
			cosA[i] = Math.cos(a[i]);
		}

		final TimingService ts = new TimingService(100);
		ts.execute(new MyTimingTask("sin")
		{
			@Override
			public Object run(Object data)
			{
				double d = 0;
				for (int i = 0; i < a.length; i++)
					d += Math.sin(a[i]);
				return d;
			}
		});
		ts.execute(new MyTimingTask("FastMath.sin")
		{
			@Override
			public Object run(Object data)
			{
				double d = 0;
				for (int i = 0; i < a.length; i++)
					d += FastMath.sin(a[i]);
				return d;
			}
		});
		ts.execute(new MyTimingTask("getSine")
		{
			@Override
			public Object run(Object data)
			{
				double d = 0;
				for (int i = 0; i < a.length; i++)
					d += FRC.getSine(a[i], cosA[i]);
				return d;
			}
		});

		final int size = ts.getSize();
		ts.repeat(size);
		if (TestSettings.allow(LogLevel.INFO))
			ts.report(size);

		Assert.assertTrue(ts.get(-1).getMean() < ts.get(-2).getMean());
		Assert.assertTrue(ts.get(-1).getMean() < ts.get(-3).getMean());
	}

	@Test
	public void computeMirroredIsFaster()
	{
		TestAssume.assumeMediumComplexity();

		// Sample lines through an image to create a structure.
		final int N = 2048;
		final double[][] data = new double[N * 2][];
		final RandomGenerator r = TestSettings.getRandomGenerator();
		for (int x = 0, y = 0, y2 = N, i = 0; x < N; x++, y++, y2--)
		{
			data[i++] = new double[] { x + r.nextGaussian() * 5, y + r.nextGaussian() * 5 };
			data[i++] = new double[] { x + r.nextGaussian() * 5, y2 + r.nextGaussian() * 5 };
		}
		// Create 2 images
		final Rectangle bounds = new Rectangle(0, 0, N, N);
		IJImagePeakResults i1 = createImage(bounds);
		IJImagePeakResults i2 = createImage(bounds);
		final int[] indices = SimpleArrayUtils.newArray(data.length, 0, 1);
		MathArrays.shuffle(indices, r);
		for (final int i : indices)
		{
			final IJImagePeakResults image = i1;
			i1 = i2;
			i2 = image;
			image.add((float) data[i][0], (float) data[i][1], 1);
		}
		i1.end();
		i2.end();
		final ImageProcessor ip1 = i1.getImagePlus().getProcessor();
		final ImageProcessor ip2 = i2.getImagePlus().getProcessor();
		// Test
		final FRC frc = new FRC();
		FloatProcessor[] fft1, fft2;
		fft1 = frc.getComplexFFT(ip1);
		fft2 = frc.getComplexFFT(ip2);

		final float[] dataA1 = (float[]) fft1[0].getPixels();
		final float[] dataB1 = (float[]) fft1[1].getPixels();
		final float[] dataA2 = (float[]) fft2[0].getPixels();
		final float[] dataB2 = (float[]) fft2[1].getPixels();

		final float[] numerator = new float[dataA1.length];
		final float[] absFFT1 = new float[dataA1.length];
		final float[] absFFT2 = new float[dataA1.length];

		final TimingService ts = new TimingService(10);
		ts.execute(new MyTimingTask("compute")
		{
			@Override
			public Object run(Object data)
			{
				FRC.compute(numerator, absFFT1, absFFT2, dataA1, dataB1, dataA2, dataB2);
				return null;
			}
		});
		ts.execute(new MyTimingTask("computeMirrored")
		{
			@Override
			public Object run(Object data)
			{
				FRC.computeMirrored(N, numerator, absFFT1, absFFT2, dataA1, dataB1, dataA2, dataB2);
				return null;
			}
		});
		ts.execute(new MyTimingTask("computeMirroredFast")
		{
			@Override
			public Object run(Object data)
			{
				FRC.computeMirroredFast(N, numerator, absFFT1, absFFT2, dataA1, dataB1, dataA2, dataB2);
				return null;
			}
		});

		final int size = ts.getSize();
		ts.repeat(size);
		if (TestSettings.allow(LogLevel.INFO))
			ts.report(size);

		Assert.assertTrue(ts.get(-1).getMean() < ts.get(-2).getMean());
		Assert.assertTrue(ts.get(-1).getMean() < ts.get(-3).getMean());
	}
}
