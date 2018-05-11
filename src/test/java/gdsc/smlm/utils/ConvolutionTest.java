package gdsc.smlm.utils;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.util.FastMath;
import org.junit.Assert;
import org.junit.Assume;
import org.junit.Test;

import pl.edu.icm.jlargearrays.ConcurrencyUtils;

public class ConvolutionTest
{
	private static RandomGenerator random = new Well19937c();

	int sizeLoops = 10;
	int sLoops = 8;

	static
	{
		// Compare speeds when single threaded
		ConcurrencyUtils.setNumberOfThreads(1);
	}

	@Test
	public void canComputeConvolution()
	{
		int size = 10;
		for (int i = 0; i < sizeLoops; i++)
		{
			double s = 0.5;
			for (int j = 0; j < sLoops; j++)
			{
				double[] data = randomData(size);
				double[] kernel = createKernel(s);
				double[] r1 = Convolution.convolve(data, kernel);
				double[] r1b = Convolution.convolve(kernel, data);
				double[] r2 = Convolution.convolveFFT(data, kernel);
				double[] r2b = Convolution.convolveFFT(kernel, data);

				Assert.assertEquals(r1.length, r1b.length);
				Assert.assertEquals(r1.length, r2.length);
				Assert.assertEquals(r1.length, r2b.length);
				for (int k = 0; k < r1.length; k++)
				{
					double error = Math.abs(r1[k] * 1e-6);
					Assert.assertEquals(r1[k], r1b[k], error);
					Assert.assertEquals(r1[k], r2[k], error);
					Assert.assertEquals(r1[k], r2b[k], error);
				}

				s *= 2;
			}
			size *= 2;
		}
	}

	@Test
	public void canComputeDoubleConvolution()
	{
		int size = 10;
		for (int i = 0; i < sizeLoops; i++)
		{
			double s = 0.5;
			for (int j = 0; j < sLoops; j++)
			{
				double[] data1 = randomData(size);
				double[] data2 = randomData(size);
				double[] kernel = createKernel(s);
				double[] e1 = Convolution.convolve(kernel, data1);
				double[] e2 = Convolution.convolve(kernel, data2);
				double[][] r1 = Convolution.convolve(kernel, data1, data2);
				double[][] r2 = Convolution.convolveFFT(kernel, data1, data2);

				Assert.assertEquals(r1.length, 2);
				Assert.assertEquals(r2.length, 2);
				Assert.assertEquals(e1.length, r1[0].length);
				Assert.assertEquals(e2.length, r1[1].length);
				Assert.assertEquals(r1[0].length, r2[0].length);
				Assert.assertEquals(r1[1].length, r2[1].length);
				for (int k = 0; k < e1.length; k++)
				{
					double error = Math.abs(e1[k] * 1e-6);
					Assert.assertEquals(e1[k], r1[0][k], error);
					Assert.assertEquals(e2[k], r1[1][k], error);
					Assert.assertEquals(e1[k], r2[0][k], error);
					Assert.assertEquals(e2[k], r2[1][k], error);
				}

				s *= 2;
			}
			size *= 2;
		}
	}

	@Test
	public void doSpeedTest()
	{
		Assume.assumeTrue(false);
		
		int size = 10;
		for (int i = 0; i < sizeLoops; i++)
		{
			double s = 0.5;
			for (int j = 0; j < sLoops; j++)
			{
				speedTest(size, s);
				s *= 2;
			}
			size *= 2;
		}
	}

	private void speedTest(int size, double s)
	{
		final int RUNS = 1000;

		double[] data = randomData(size);
		double[] kernel = createKernel(s);

		// Warm up
		@SuppressWarnings("unused")
		double[] r1 = Convolution.convolve(kernel, data);
		@SuppressWarnings("unused")
		double[] r2 = Convolution.convolveFFT(kernel, data);

		long t1 = System.nanoTime();
		for (int i = 0; i < RUNS; i++)
			r1 = Convolution.convolve(kernel, data);
		t1 = System.nanoTime() - t1;

		long t2 = System.nanoTime();
		for (int i = 0; i < RUNS; i++)
			r2 = Convolution.convolveFFT(kernel, data);
		t2 = System.nanoTime() - t2;

		System.out.printf("Size=%d, s=%f (%d) [%d] : %d -> %d (%f)\n", size, s, kernel.length, size * kernel.length, t1,
				t2, t1 / (double) t2);
	}

	@Test
	public void doDoubleSpeedTest()
	{
		Assume.assumeTrue(false);
		
		int size = 10;
		for (int i = 0; i < sizeLoops; i++)
		{
			double s = 0.5;
			for (int j = 0; j < sLoops; j++)
			{
				doubleSpeedTest(size, s);
				s *= 2;
			}
			size *= 2;
		}
	}

	private void doubleSpeedTest(int size, double s)
	{
		final int RUNS = 1000;

		double[] data1 = randomData(size);
		double[] data2 = randomData(size);
		double[] kernel = createKernel(s);

		// Warm up
		@SuppressWarnings("unused")
		double[][] r1 = Convolution.convolve(kernel, data1, data2);
		@SuppressWarnings("unused")
		double[][] r2 = Convolution.convolveFFT(kernel, data1, data2);

		long t1 = System.nanoTime();
		for (int i = 0; i < RUNS; i++)
			r1 = Convolution.convolve(kernel, data1, data2);
		t1 = System.nanoTime() - t1;

		long t2 = System.nanoTime();
		for (int i = 0; i < RUNS; i++)
			r2 = Convolution.convolveFFT(kernel, data1, data2);
		t2 = System.nanoTime() - t2;

		System.out.printf("Size=%d, s=%f (%d) [%d] : %d -> %d (%f)\n", size, s, kernel.length, size * kernel.length, t1,
				t2, t1 / (double) t2);
	}
	
	@Test
	public void doSingleVsDoubleSpeedTest()
	{
		Assume.assumeTrue(false);
		
		int size = 10;
		for (int i = 0; i < sizeLoops/2; i++)
		{
			double s = 0.5;
			for (int j = 0; j < sLoops; j++)
			{
				singleVsDoubleSpeedTest(size, s);
				s *= 2;
			}
			size *= 2;
		}
	}

	private void singleVsDoubleSpeedTest(int size, double s)
	{
		final int RUNS = 1000;

		double[] data1 = randomData(size);
		double[] data2 = randomData(size);
		double[] kernel = createKernel(s);

		// Warm up
		@SuppressWarnings("unused")
		double[] r1;
		r1 = Convolution.convolve(kernel, data1);
		r1 = Convolution.convolve(kernel, data2);
		@SuppressWarnings("unused")
		double[][] r2 = Convolution.convolve(kernel, data1, data2);

		long t1 = System.nanoTime();
		for (int i = 0; i < RUNS; i++)
		{
			r1 = Convolution.convolve(kernel, data1);
			r1 = Convolution.convolve(kernel, data2);
		}
		t1 = System.nanoTime() - t1;

		long t2 = System.nanoTime();
		for (int i = 0; i < RUNS; i++)
			r2 = Convolution.convolve(kernel, data1, data2);
		t2 = System.nanoTime() - t2;

		System.out.printf("Size=%d, s=%f (%d) [%d] : %d -> %d (%f)\n", size, s, kernel.length, size * kernel.length, t1,
				t2, t1 / (double) t2);
	}
	
	@Test
	public void doSingleVsDoubleFFTSpeedTest()
	{
		Assume.assumeTrue(false);
		
		int size = 10;
		for (int i = 0; i < sizeLoops/2; i++)
		{
			double s = 0.5;
			for (int j = 0; j < sLoops; j++)
			{
				singleVsDoubleFFTSpeedTest(size, s);
				s *= 2;
			}
			size *= 2;
		}
	}

	private void singleVsDoubleFFTSpeedTest(int size, double s)
	{
		final int RUNS = 1000;

		double[] data1 = randomData(size);
		double[] data2 = randomData(size);
		double[] kernel = createKernel(s);

		// Warm up
		@SuppressWarnings("unused")
		double[] r1;
		r1 = Convolution.convolveFFT(kernel, data1);
		r1 = Convolution.convolveFFT(kernel, data2);
		@SuppressWarnings("unused")
		double[][] r2 = Convolution.convolveFFT(kernel, data1, data2);

		long t1 = System.nanoTime();
		for (int i = 0; i < RUNS; i++)
		{
			r1 = Convolution.convolveFFT(kernel, data1);
			r1 = Convolution.convolveFFT(kernel, data2);
		}
		t1 = System.nanoTime() - t1;

		long t2 = System.nanoTime();
		for (int i = 0; i < RUNS; i++)
			r2 = Convolution.convolveFFT(kernel, data1, data2);
		t2 = System.nanoTime() - t2;

		System.out.printf("Size=%d, s=%f (%d) [%d] : %d -> %d (%f)\n", size, s, kernel.length, size * kernel.length, t1,
				t2, t1 / (double) t2);
	}
	
	private double[] randomData(int size)
	{
		double[] data = new double[size];
		for (int i = 0; i < size; i++)
			data[i] = random.nextDouble();
		return data;
	}

	/**
	 * Create a Gaussian kernel of standard deviation s
	 * 
	 * @param s
	 * @return the kernel
	 */
	private double[] createKernel(double s)
	{
		final int radius = (int) Math.ceil(Math.abs(s) * 4) + 1;
		double[] kernel = new double[2 * radius + 1];
		final double norm = -0.5 / (s * s);
		for (int i = 0, j = radius, jj = radius; j < kernel.length; i++, j++, jj--)
			kernel[j] = kernel[jj] = FastMath.exp(norm * i * i);
		// Normalise
		double sum = 0;
		for (int j = 0; j < kernel.length; j++)
			sum += kernel[j];
		for (int j = 0; j < kernel.length; j++)
			kernel[j] /= sum;
		return kernel;
	}
}
