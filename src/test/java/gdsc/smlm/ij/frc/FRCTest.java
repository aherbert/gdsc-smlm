package gdsc.smlm.ij.frc;

import java.awt.Rectangle;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathArrays;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.ij.Utils;
import gdsc.core.test.BaseTimingTask;
import gdsc.core.test.TimingService;
import gdsc.core.utils.DoubleEquality;
import gdsc.smlm.ij.results.IJImagePeakResults;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class FRCTest
{
	@Test
	public void canComputeSine()
	{
		int steps = 1000;
		double delta = 2 * Math.PI / steps;
		for (int i = 0; i <= steps; i++)
		{
			double a = i * delta;
			double cosA = Math.cos(a);
			double e = Math.sin(a);
			double o = FRC.getSine(a, cosA);
			//System.out.printf("%f  %f ?= %f\n", a, e, o);
			Assert.assertTrue(DoubleEquality.almostEqualRelativeOrAbsolute(o, e, 1e-6, 1e-10));
		}
	}

	@Test
	public void canComputeMirrored()
	{
		// Sample lines through an image to create a structure.
		int size = 1024;
		double[][] data = new double[size * 2][];
		RandomGenerator r = new Well19937c(30051977);
		for (int x = 0, y = 0, y2 = size, i = 0; x < size; x++, y++, y2--)
		{
			data[i++] = new double[] { x + r.nextGaussian() * 5, y + r.nextGaussian() * 5 };
			data[i++] = new double[] { x + r.nextGaussian() * 5, y2 + r.nextGaussian() * 5 };
		}
		// Create 2 images
		Rectangle bounds = new Rectangle(0, 0, size, size);
		IJImagePeakResults i1 = createImage(bounds);
		IJImagePeakResults i2 = createImage(bounds);
		int[] indices = Utils.newArray(data.length, 0, 1);
		MathArrays.shuffle(indices, r);
		for (int i : indices)
		{
			IJImagePeakResults image = i1;
			i1 = i2;
			i2 = image;
			image.add((float) data[i][0], (float) data[i][1], 1);
		}
		i1.end();
		i2.end();
		ImageProcessor ip1 = i1.getImagePlus().getProcessor();
		ImageProcessor ip2 = i2.getImagePlus().getProcessor();
		// Test
		FRC frc = new FRC();
		FloatProcessor[] fft1, fft2;
		fft1 = frc.getComplexFFT(ip1);
		fft2 = frc.getComplexFFT(ip2);

		float[] dataA1 = (float[]) fft1[0].getPixels();
		float[] dataB1 = (float[]) fft1[1].getPixels();
		float[] dataA2 = (float[]) fft2[0].getPixels();
		float[] dataB2 = (float[]) fft2[1].getPixels();

		float[] numeratorE = new float[dataA1.length];
		float[] absFFT1E = new float[dataA1.length];
		float[] absFFT2E = new float[dataA1.length];

		FRC.compute(numeratorE, absFFT1E, absFFT2E, dataA1, dataB1, dataA2, dataB2);

		Assert.assertTrue("numeratorE", FRC.checkSymmetry(numeratorE, size));
		Assert.assertTrue("absFFT1E", FRC.checkSymmetry(absFFT1E, size));
		Assert.assertTrue("absFFT2E", FRC.checkSymmetry(absFFT2E, size));

		float[] numeratorA = new float[dataA1.length];
		float[] absFFT1A = new float[dataA1.length];
		float[] absFFT2A = new float[dataA1.length];
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
		for (int y=1; y<size; y++)
			for (int x=1, i=y*size+1; x<size; x++, i++)
			{
				Assert.assertEquals("numerator", numeratorE[i], numeratorA[i], 0);
				Assert.assertEquals("absFFT1", absFFT1E[i], absFFT1A[i], 0);
				Assert.assertEquals("absFFT2", absFFT2E[i], absFFT2A[i], 0);
			}
	}

	private IJImagePeakResults createImage(Rectangle bounds)
	{
		IJImagePeakResults i1 = new IJImagePeakResults("1", bounds, 1);
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

		public int getSize()
		{
			return 1;
		}

		public Object getData(int i)
		{
			return null;
		}
	}

	@Test
	public void computeSineIsFaster()
	{
		int steps = 100000;
		double delta = 2 * Math.PI / steps;
		final double[] a = new double[steps + 1];
		final double[] cosA = new double[steps + 1];
		for (int i = 0; i <= steps; i++)
		{
			a[i] = i * delta;
			cosA[i] = Math.cos(a[i]);
		}

		TimingService ts = new TimingService(100);
		ts.execute(new MyTimingTask("sin")
		{
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
			public Object run(Object data)
			{
				double d = 0;
				for (int i = 0; i < a.length; i++)
					d += FRC.getSine(a[i], cosA[i]);
				return d;
			}
		});

		ts.repeat(ts.getSize());
		ts.report();
	}
	
	@Test
	public void computeMirroredIsFaster()
	{
		// Sample lines through an image to create a structure.
		final int size = 2048;
		double[][] data = new double[size * 2][];
		RandomGenerator r = new Well19937c(30051977);
		for (int x = 0, y = 0, y2 = size, i = 0; x < size; x++, y++, y2--)
		{
			data[i++] = new double[] { x + r.nextGaussian() * 5, y + r.nextGaussian() * 5 };
			data[i++] = new double[] { x + r.nextGaussian() * 5, y2 + r.nextGaussian() * 5 };
		}
		// Create 2 images
		Rectangle bounds = new Rectangle(0, 0, size, size);
		IJImagePeakResults i1 = createImage(bounds);
		IJImagePeakResults i2 = createImage(bounds);
		int[] indices = Utils.newArray(data.length, 0, 1);
		MathArrays.shuffle(indices, r);
		for (int i : indices)
		{
			IJImagePeakResults image = i1;
			i1 = i2;
			i2 = image;
			image.add((float) data[i][0], (float) data[i][1], 1);
		}
		i1.end();
		i2.end();
		ImageProcessor ip1 = i1.getImagePlus().getProcessor();
		ImageProcessor ip2 = i2.getImagePlus().getProcessor();
		// Test
		FRC frc = new FRC();
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
		
		TimingService ts = new TimingService(10);
		ts.execute(new MyTimingTask("compute")
		{
			public Object run(Object data)
			{
				FRC.compute(numerator, absFFT1, absFFT2, dataA1, dataB1, dataA2, dataB2);
				return null;
			}
		});
		ts.execute(new MyTimingTask("computeMirrored")
		{
			public Object run(Object data)
			{
				FRC.computeMirrored(size, numerator, absFFT1, absFFT2, dataA1, dataB1, dataA2, dataB2);
				return null;
			}
		});
		ts.execute(new MyTimingTask("computeMirroredFast")
		{
			public Object run(Object data)
			{
				FRC.computeMirroredFast(size, numerator, absFFT1, absFFT2, dataA1, dataB1, dataA2, dataB2);
				return null;
			}
		});

		ts.repeat(ts.getSize());
		ts.report();
	}	
}
