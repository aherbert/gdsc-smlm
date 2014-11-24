package gdsc.smlm.utils;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

import edu.emory.mathcs.utils.ConcurrencyUtils;

public class ConvolutionTest
{
	private static RandomGenerator random = new Well19937c();

	int sizeLoops = 8;
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
				double[] r2 = Convolution.convolveFFT(data, kernel);
				Assert.assertArrayEquals(
						String.format("Size=%d, s=%f (%d) [%d}", size, s, kernel.length, size * kernel.length), r1, r2,
						1e-10);

				s *= 2;
			}
			size *= 2;
		}
	}

	@Test
	public void doSpeedTest()
	{
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
		double[] r1 = Convolution.convolve(data, kernel);
		double[] r2 = Convolution.convolveFFT(data, kernel);

		long t1 = System.nanoTime();
		for (int i = 0; i < RUNS; i++)
			r1 = Convolution.convolve(data, kernel);
		t1 = System.nanoTime() - t1;

		long t2 = System.nanoTime();
		for (int i = 0; i < RUNS; i++)
			r2 = Convolution.convolveFFT(data, kernel);
		t2 = System.nanoTime() - t2;

		System.out.printf("Size=%d, s=%f (%d) [%d] : %d -> %d (%f)\n", size, s, kernel.length, size * kernel.length,
				t1, t2, t1 / (double) t2);
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
