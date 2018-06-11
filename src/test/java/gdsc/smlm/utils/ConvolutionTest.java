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
package gdsc.smlm.utils;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.util.FastMath;
import org.junit.Assert;
import org.junit.Assume;
import org.junit.Test;

import gdsc.smlm.utils.Convolution.ConvolutionValueProcedure;
import gdsc.smlm.utils.Convolution.DoubleConvolutionValueProcedure;
import gnu.trove.list.array.TDoubleArrayList;
import pl.edu.icm.jlargearrays.ConcurrencyUtils;

public class ConvolutionTest
{
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
		RandomGenerator random = new Well19937c(30051977);
		int size = 10;
		for (int i = 0; i < sizeLoops; i++)
		{
			double s = 0.5;
			for (int j = 0; j < sLoops; j++)
			{
				double[] data = randomData(random, size);
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
		RandomGenerator random = new Well19937c(30051977);
		int size = 10;
		for (int i = 0; i < sizeLoops; i++)
		{
			double s = 0.5;
			for (int j = 0; j < sLoops; j++)
			{
				double[] data1 = randomData(random, size);
				double[] data2 = randomData(random, size);
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
		RandomGenerator random = new Well19937c(30051977);
		final int RUNS = 1000;

		double[] data = randomData(random, size);
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
		RandomGenerator random = new Well19937c(30051977);
		final int RUNS = 1000;

		double[] data1 = randomData(random, size);
		double[] data2 = randomData(random, size);
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
		for (int i = 0; i < sizeLoops / 2; i++)
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
		RandomGenerator random = new Well19937c(30051977);
		final int RUNS = 1000;

		double[] data1 = randomData(random, size);
		double[] data2 = randomData(random, size);
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
		for (int i = 0; i < sizeLoops / 2; i++)
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
		RandomGenerator random = new Well19937c(30051977);
		final int RUNS = 1000;

		double[] data1 = randomData(random, size);
		double[] data2 = randomData(random, size);
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

	private double[] randomData(RandomGenerator random, int size)
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

	@Test
	public void canComputeScaledConvolution()
	{
		RandomGenerator random = new Well19937c(30051977);
		TDoubleArrayList list = new TDoubleArrayList();
		int size = 10;
		for (int i = 0; i < sizeLoops; i++)
		{
			double s = 0.5;
			for (int j = 0; j < sLoops; j++)
			{
				double[] data = randomData(random, size);
				double[] kernel = createKernel(s);

				for (int scale = 2; scale < 5; scale++)
				{
					double[] e = convolve(kernel, data, list, scale);
					double[] o = Convolution.convolve(kernel, data, scale);
					final double[] o2 = new double[o.length];
					Convolution.convolve(kernel, data, scale, new ConvolutionValueProcedure()
					{
						int i = 0;

						public boolean execute(double value)
						{
							o2[i++] = value;
							return true;
						}
					});

					Assert.assertArrayEquals(e, o, 0);
					Assert.assertArrayEquals(e, o2, 0);
				}

				s *= 2;
			}
			size *= 2;
		}
	}

	@Test
	public void canComputeDoubleScaledConvolution()
	{
		RandomGenerator random = new Well19937c(30051977);
		TDoubleArrayList list = new TDoubleArrayList();
		int size = 10;
		for (int i = 0; i < sizeLoops / 2; i++)
		{
			double s = 0.5;
			for (int j = 0; j < sLoops; j++)
			{
				double[] data1 = randomData(random, size);
				double[] data2 = randomData(random, size);
				double[] kernel = createKernel(s);

				for (int scale = 2; scale < 5; scale++)
				{
					double[] e1 = convolve(kernel, data1, list, scale);
					double[] e2 = convolve(kernel, data2, list, scale);
					double[][] o = Convolution.convolve(kernel, data1, data2, scale);
					final double[][] o2 = new double[2][o[0].length];
					Convolution.convolve(kernel, data1, data2, scale, new DoubleConvolutionValueProcedure()
					{
						int i = 0;

						public boolean execute(double value1, double value2)
						{
							o2[0][i] = value1;
							o2[1][i] = value2;
							i++;
							return true;
						}
					});

					Assert.assertArrayEquals(e1, o[0], 0);
					Assert.assertArrayEquals(e1, o2[0], 0);
					Assert.assertArrayEquals(e2, o[1], 0);
					Assert.assertArrayEquals(e2, o2[1], 0);
				}

				s *= 2;
			}
			size *= 2;
		}
	}

	private static double[] convolve(double[] kernel, double[] data, TDoubleArrayList list, int scale)
	{
		double[] data1 = scale(data, list, scale);
		return Convolution.convolve(kernel, data1);
	}

	private static double[] scale(double[] data, TDoubleArrayList list, int scale)
	{
		list.resetQuick();
		double[] fill = new double[scale - 1];
		list.add(data[0]);
		for (int i = 1; i < data.length; i++)
		{
			list.add(fill);
			list.add(data[i]);
		}
		return list.toArray();
	}

	@Test
	public void canComputeScaledConvolutionWithEarlyExit()
	{
		RandomGenerator random = new Well19937c(30051977);
		int size = 10;
		for (int i = 0; i < sizeLoops; i++)
		{
			double s = 0.5;
			for (int j = 0; j < sLoops; j++)
			{
				double[] data = randomData(random, size);
				double[] kernel = createKernel(s);

				for (int scale = 2; scale < 5; scale++)
				{
					double[] e = Convolution.convolve(kernel, data, scale);
					final double[] o = new double[e.length];
					final int limit = data.length;
					Convolution.convolve(kernel, data, scale, new ConvolutionValueProcedure()
					{
						int i = 0;

						public boolean execute(double value)
						{
							o[i++] = value;
							return i < limit;
						}
					});

					int k = 0;
					for (; k < limit; k++)
						Assert.assertEquals(e[k], o[k], 0);
					while (k < o.length)
						Assert.assertEquals(0, o[k++], 0);
				}

				s *= 2;
			}
			size *= 2;
		}
	}

	@Test
	public void canComputeDoubleScaledConvolutionWithEarlyExit()
	{
		RandomGenerator random = new Well19937c(30051977);
		int size = 10;
		for (int i = 0; i < sizeLoops; i++)
		{
			double s = 0.5;
			for (int j = 0; j < sLoops; j++)
			{
				double[] data1 = randomData(random, size);
				double[] data2 = randomData(random, size);
				double[] kernel = createKernel(s);

				for (int scale = 2; scale < 5; scale++)
				{
					double[][] e = Convolution.convolve(kernel, data1, data2, scale);
					final double[][] o = new double[2][e[0].length];
					final int limit = data1.length;
					Convolution.convolve(kernel, data1, data2, scale, new DoubleConvolutionValueProcedure()
					{
						int i = 0;

						public boolean execute(double value1, double value2)
						{
							o[0][i] = value1;
							o[1][i] = value1;
							i++;
							return i < limit;
						}
					});

					int k = 0;
					for (; k < limit; k++)
					{
						Assert.assertEquals(e[0][k], o[0][k], 0);
						Assert.assertEquals(e[0][k], o[1][k], 0);
					}
					while (k < o.length)
					{
						Assert.assertEquals(0, o[0][k], 0);
						Assert.assertEquals(0, o[1][k], 0);
						k++;
					}
				}

				s *= 2;
			}
			size *= 2;
		}
	}

	@Test
	public void doScaledSpeedTest()
	{
		Assume.assumeTrue(false);

		int size = 10;
		for (int scale = 4; scale <= 8; scale *= 2)
		{
			for (int i = 0; i < 4; i++)
			{
				double s = 0.5;
				for (int j = 0; j < 4; j++)
				{
					doScaledSpeedTest(size, s, scale);
					s *= 2;
				}
				size *= 2;
			}
		}
	}

	private void doScaledSpeedTest(int size, double s, int scale)
	{
		RandomGenerator random = new Well19937c(30051977);
		final int RUNS = 100;

		double[] data1 = randomData(random, size);
		double[] kernel = createKernel(s);
		TDoubleArrayList list = new TDoubleArrayList();

		// Warm up
		convolve(kernel, data1, list, scale);
		Convolution.convolve(kernel, data1, scale);

		long t1 = System.nanoTime();
		for (int i = 0; i < RUNS; i++)
		{
			convolve(kernel, data1, list, scale);
		}
		t1 = System.nanoTime() - t1;

		long t2 = System.nanoTime();
		for (int i = 0; i < RUNS; i++)
			Convolution.convolve(kernel, data1, scale);
		t2 = System.nanoTime() - t2;

		System.out.printf("Size=%d, s=%f, scale=%d (%d) [%d] : %d -> %d (%f)\n", size, s, scale, kernel.length,
				size * kernel.length, t1, t2, t1 / (double) t2);
	}
}
