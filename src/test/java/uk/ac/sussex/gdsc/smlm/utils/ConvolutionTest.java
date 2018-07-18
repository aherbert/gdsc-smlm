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
package uk.ac.sussex.gdsc.smlm.utils;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

import gnu.trove.list.array.TDoubleArrayList;
import pl.edu.icm.jlargearrays.ConcurrencyUtils;
import uk.ac.sussex.gdsc.smlm.utils.Convolution.ConvolutionValueProcedure;
import uk.ac.sussex.gdsc.smlm.utils.Convolution.DoubleConvolutionValueProcedure;
import uk.ac.sussex.gdsc.test.LogLevel;
import uk.ac.sussex.gdsc.test.TestAssert;
import uk.ac.sussex.gdsc.test.TestComplexity;
import uk.ac.sussex.gdsc.test.TestSettings;

@SuppressWarnings({ "javadoc" })
public class ConvolutionTest
{
	int sizeLoops = 8;
	int sLoops = 6;

	static
	{
		// Compare speeds when single threaded
		ConcurrencyUtils.setNumberOfThreads(1);
	}

	@Test
	public void canComputeConvolution()
	{
		final RandomGenerator random = TestSettings.getRandomGenerator();
		int size = 10;
		for (int i = 0; i < sizeLoops; i++)
		{
			double s = 0.5;
			for (int j = 0; j < sLoops; j++)
			{
				final double[] data = randomData(random, size);
				final double[] kernel = createKernel(s);
				final double[] r1 = Convolution.convolve(data, kernel);
				final double[] r1b = Convolution.convolve(kernel, data);
				final double[] r2 = Convolution.convolveFFT(data, kernel);
				final double[] r2b = Convolution.convolveFFT(kernel, data);

				Assert.assertEquals(r1.length, r1b.length);
				Assert.assertEquals(r1.length, r2.length);
				Assert.assertEquals(r1.length, r2b.length);

				TestAssert.assertArrayEqualsRelative("Spatial convolution doesn't match", r1, r1b, 1e-6);
				TestAssert.assertArrayEqualsRelative("FFT convolution doesn't match", r2, r2b, 1e-6);

				s *= 2;
			}
			size *= 2;
		}
	}

	@Test
	public void canComputeDoubleConvolution()
	{
		final RandomGenerator random = TestSettings.getRandomGenerator();
		int size = 10;
		for (int i = 0; i < sizeLoops; i++)
		{
			double s = 0.5;
			for (int j = 0; j < sLoops; j++)
			{
				final double[] data1 = randomData(random, size);
				final double[] data2 = randomData(random, size);
				final double[] kernel = createKernel(s);

				double[] e1, e2;
				double[][] r1;
				for (int fft = 0; fft < 2; fft++)
				{
					if (fft == 1)
					{
						e1 = Convolution.convolveFFT(kernel, data1);
						e2 = Convolution.convolveFFT(kernel, data2);
						r1 = Convolution.convolveFFT(kernel, data1, data2);
					}
					else
					{
						e1 = Convolution.convolve(kernel, data1);
						e2 = Convolution.convolve(kernel, data2);
						r1 = Convolution.convolve(kernel, data1, data2);
					}

					Assert.assertEquals(r1.length, 2);
					Assert.assertEquals(e1.length, r1[0].length);
					Assert.assertEquals(e2.length, r1[1].length);

					for (int k = 0; k < e1.length; k++)
					{
						// Exact match
						Assert.assertEquals(e1[k], r1[0][k], 0);
						Assert.assertEquals(e2[k], r1[1][k], 0);
					}
				}

				s *= 2;
			}
			size *= 2;
		}
	}

	@Test
	public void doSpeedTest()
	{
		TestSettings.assume(LogLevel.INFO, TestComplexity.MEDIUM);
		final RandomGenerator rg = TestSettings.getRandomGenerator();

		int size = 10;
		for (int i = 0; i < sizeLoops; i++)
		{
			double s = 0.5;
			for (int j = 0; j < sLoops; j++)
			{
				speedTest(rg, size, s);
				s *= 2;
			}
			size *= 2;
		}
	}

	private static void speedTest(RandomGenerator rg, int size, double s)
	{
		final int RUNS = 1000;

		final double[] data = randomData(rg, size);
		final double[] kernel = createKernel(s);

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
		TestSettings.assume(LogLevel.INFO, TestComplexity.MEDIUM);
		final RandomGenerator rg = TestSettings.getRandomGenerator();

		int size = 10;
		for (int i = 0; i < sizeLoops; i++)
		{
			double s = 0.5;
			for (int j = 0; j < sLoops; j++)
			{
				doubleSpeedTest(rg, size, s);
				s *= 2;
			}
			size *= 2;
		}
	}

	private static void doubleSpeedTest(RandomGenerator rg, int size, double s)
	{
		final int RUNS = 1000;

		final double[] data1 = randomData(rg, size);
		final double[] data2 = randomData(rg, size);
		final double[] kernel = createKernel(s);

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
		TestSettings.assume(LogLevel.INFO, TestComplexity.MEDIUM);

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

	private static void singleVsDoubleSpeedTest(int size, double s)
	{
		final RandomGenerator random = TestSettings.getRandomGenerator();
		final int RUNS = 1000;

		final double[] data1 = randomData(random, size);
		final double[] data2 = randomData(random, size);
		final double[] kernel = createKernel(s);

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
		TestSettings.assume(LogLevel.INFO, TestComplexity.MEDIUM);

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

	private static void singleVsDoubleFFTSpeedTest(int size, double s)
	{
		final RandomGenerator random = TestSettings.getRandomGenerator();
		final int RUNS = 1000;

		final double[] data1 = randomData(random, size);
		final double[] data2 = randomData(random, size);
		final double[] kernel = createKernel(s);

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

	private static double[] randomData(RandomGenerator random, int size)
	{
		final double[] data = new double[size];
		for (int i = 0; i < size; i++)
			data[i] = random.nextDouble();
		return data;
	}

	/**
	 * Create a Gaussian kernel of standard deviation s.
	 *
	 * @param s
	 *            the standard deviation
	 * @return the kernel
	 */
	private static double[] createKernel(double s)
	{
		final int radius = (int) Math.ceil(Math.abs(s) * 4) + 1;
		final double[] kernel = new double[2 * radius + 1];
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
		final RandomGenerator random = TestSettings.getRandomGenerator();
		final TDoubleArrayList list = new TDoubleArrayList();
		int size = 10;
		for (int i = 0; i < sizeLoops; i++)
		{
			double s = 0.5;
			for (int j = 0; j < sLoops; j++)
			{
				final double[] data = randomData(random, size);
				final double[] kernel = createKernel(s);

				for (int scale = 2; scale < 5; scale++)
				{
					final double[] e = convolve(kernel, data, list, scale);
					final double[] o = Convolution.convolve(kernel, data, scale);
					final double[] o2 = new double[o.length];
					Convolution.convolve(kernel, data, scale, new ConvolutionValueProcedure()
					{
						int i = 0;

						@Override
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
		final RandomGenerator random = TestSettings.getRandomGenerator();
		final TDoubleArrayList list = new TDoubleArrayList();
		int size = 10;
		for (int i = 0; i < sizeLoops / 2; i++)
		{
			double s = 0.5;
			for (int j = 0; j < sLoops; j++)
			{
				final double[] data1 = randomData(random, size);
				final double[] data2 = randomData(random, size);
				final double[] kernel = createKernel(s);

				for (int scale = 2; scale < 5; scale++)
				{
					final double[] e1 = convolve(kernel, data1, list, scale);
					final double[] e2 = convolve(kernel, data2, list, scale);
					final double[][] o = Convolution.convolve(kernel, data1, data2, scale);
					final double[][] o2 = new double[2][o[0].length];
					Convolution.convolve(kernel, data1, data2, scale, new DoubleConvolutionValueProcedure()
					{
						int i = 0;

						@Override
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
		final double[] data1 = scale(data, list, scale);
		return Convolution.convolve(kernel, data1);
	}

	private static double[] scale(double[] data, TDoubleArrayList list, int scale)
	{
		list.resetQuick();
		final double[] fill = new double[scale - 1];
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
		final RandomGenerator random = TestSettings.getRandomGenerator();
		int size = 10;
		final int sizeLoops = 4;
		final int sLoops = 2;
		for (int i = 0; i < sizeLoops; i++)
		{
			double s = 0.5;
			for (int j = 0; j < sLoops; j++)
			{
				final double[] data = randomData(random, size);
				final double[] kernel = createKernel(s);

				for (int scale = 2; scale < 5; scale++)
				{
					final double[] e = Convolution.convolve(kernel, data, scale);
					final double[] o = new double[e.length];
					final int limit = data.length;
					Convolution.convolve(kernel, data, scale, new ConvolutionValueProcedure()
					{
						int i = 0;

						@Override
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
		final RandomGenerator random = TestSettings.getRandomGenerator();
		int size = 10;
		final int sizeLoops = 4;
		final int sLoops = 2;
		for (int i = 0; i < sizeLoops; i++)
		{
			double s = 0.5;
			for (int j = 0; j < sLoops; j++)
			{
				final double[] data1 = randomData(random, size);
				final double[] data2 = randomData(random, size);
				final double[] kernel = createKernel(s);

				for (int scale = 2; scale < 5; scale++)
				{
					final double[][] e = Convolution.convolve(kernel, data1, data2, scale);
					final double[][] o = new double[2][e[0].length];
					final int limit = data1.length;
					Convolution.convolve(kernel, data1, data2, scale, new DoubleConvolutionValueProcedure()
					{
						int i = 0;

						@Override
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
		TestSettings.assume(LogLevel.INFO, TestComplexity.MEDIUM);

		int size = 10;
		for (int scale = 4; scale <= 8; scale *= 2)
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

	private static void doScaledSpeedTest(int size, double s, int scale)
	{
		final RandomGenerator random = TestSettings.getRandomGenerator();
		final int RUNS = 100;

		final double[] data1 = randomData(random, size);
		final double[] kernel = createKernel(s);
		final TDoubleArrayList list = new TDoubleArrayList();

		// Warm up
		convolve(kernel, data1, list, scale);
		Convolution.convolve(kernel, data1, scale);

		long t1 = System.nanoTime();
		for (int i = 0; i < RUNS; i++)
			convolve(kernel, data1, list, scale);
		t1 = System.nanoTime() - t1;

		long t2 = System.nanoTime();
		for (int i = 0; i < RUNS; i++)
			Convolution.convolve(kernel, data1, scale);
		t2 = System.nanoTime() - t2;

		System.out.printf("Size=%d, s=%f, scale=%d (%d) [%d] : %d -> %d (%f)\n", size, s, scale, kernel.length,
				size * kernel.length, t1, t2, t1 / (double) t2);
	}
}
