package gdsc.smlm.utils;

import java.util.Arrays;

import org.apache.commons.math3.util.FastMath;

import org.jtransforms.fft.DoubleFFT_1D;
import org.jtransforms.utils.CommonUtils;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Simple class to perform convolution
 */
public class Convolution
{
	/**
	 * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution">
	 * convolution</a> between two sequences.
	 * <p>
	 * The solution is obtained via straightforward computation of the convolution sum (and not via FFT). Whenever the
	 * computation needs an element that would be located at an index outside the input arrays, the value is assumed to
	 * be zero.
	 * <p>
	 * This has been taken from Apache Commons Math v3.3: org.apache.commons.math3.util.MathArrays
	 *
	 * @param x
	 *            First sequence.
	 *            Typically, this sequence will represent an input signal to a system.
	 * @param h
	 *            Second sequence.
	 *            Typically, this sequence will represent the impulse response of the system.
	 * @return the convolution of {@code x} and {@code h}.
	 *         This array's length will be {@code x.length + h.length - 1}.
	 * @throws IllegalArgumentException
	 *             if either {@code x} or {@code h} is {@code null} or either {@code x} or {@code h} is empty.
	 */
	public static double[] convolve(double[] x, double[] h) throws IllegalArgumentException
	{
		checkInput(x, h);

		final int xLen = x.length;
		final int hLen = h.length;

		// initialize the output array
		final int totalLength = xLen + hLen - 1;
		final double[] y = new double[totalLength];

		// straightforward implementation of the convolution sum
		for (int n = 0; n < totalLength; n++)
		{
			double yn = 0;
			int k = FastMath.max(0, n + 1 - xLen);
			int j = n - k;
			while (k < hLen && j >= 0)
			{
				yn += x[j--] * h[k++];
			}
			y[n] = yn;
		}

		return y;
	}

	/**
	 * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution">
	 * convolution</a> between two sequences.
	 * <p>
	 * The solution is obtained via multiplication in the frequency domain.
	 *
	 * @param x
	 *            First sequence.
	 *            Typically, this sequence will represent an input signal to a system.
	 * @param h
	 *            Second sequence.
	 *            Typically, this sequence will represent the impulse response of the system.
	 * @return the convolution of {@code x} and {@code h}.
	 *         This array's length will be {@code x.length + h.length - 1}.
	 * @throws IllegalArgumentException
	 *             if either {@code x} or {@code h} is {@code null} or either {@code x} or {@code h} is empty.
	 */
	public static double[] convolveFFT(double[] x, double[] h) throws IllegalArgumentException
	{
		checkInput(x, h);

		// This is not needed
		//if (x.length < h.length)
		//{
		//	// Swap so that the longest array is the signal
		//	final double[] tmp = x;
		//	x = h;
		//	h = tmp;
		//}

		final int xLen = x.length;
		final int hLen = h.length;
		final int totalLength = xLen + hLen - 1;

		// Get length to a power of 2
		int newL = CommonUtils.nextPow2(totalLength);

		// Double the new length for complex values in DoubleFFT_1D
		x = Arrays.copyOf(x, 2 * newL);
		h = Arrays.copyOf(h, x.length);

		DoubleFFT_1D fft = new DoubleFFT_1D(newL);

		// FFT
		fft.realForwardFull(x);
		fft.realForwardFull(h);

		// Complex multiply. Reuse data array
		for (int i = 0; i < x.length; i += 2)
		{
			int j = i + 1;
			double xi = x[i];
			double xj = x[j];
			double hi = h[i];
			double hj = h[j];
			h[i] = hi * xi - hj * xj;
			h[j] = hi * xj + hj * xi;
		}

		// Inverse FFT
		fft.complexInverse(h, true);

		// Fill result with real part
		final double[] y = new double[totalLength];
		for (int i = 0; i < totalLength; i++)
		{
			y[i] = h[2 * i];
		}
		return y;
	}

	/**
	 * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution">
	 * convolution</a> between two sequences.
	 * <p>
	 * The solution is obtained using either the spatial or frequency domain depending on the size. The switch is made
	 * when the min array length is above 127 and the product of the lengths is above 40000. Speed tests have
	 * been performed for single threaded FFT computation. The FFT library begins multi-threaded computation when the
	 * size of the array is above length 8192.
	 *
	 * @param x
	 *            First sequence.
	 *            Typically, this sequence will represent an input signal to a system.
	 * @param h
	 *            Second sequence.
	 *            Typically, this sequence will represent the impulse response of the system.
	 * @return the convolution of {@code x} and {@code h}.
	 *         This array's length will be {@code x.length + h.length - 1}.
	 * @throws IllegalArgumentException
	 *             if either {@code x} or {@code h} is {@code null} or either {@code x} or {@code h} is empty.
	 */
	public static double[] convolveFast(double[] x, double[] h) throws IllegalArgumentException
	{
		checkInput(x, h);
		// See Junit class ConvolveTest to determine when to switch to the FFT method.
		// This is not perfect for all length combinations but the switch will happen 
		// when the two methods are roughly the same speed.
		int min, max;
		if (x.length < h.length)
		{
			min = x.length;
			max = h.length;
		}
		else
		{
			min = h.length;
			max = x.length;
		}
		if (min >= 128 && (long) min * (long) max > 40000L)
			return convolveFFT(x, h);
		return convolve(x, h);
	}

	private static void checkInput(double[] x, double[] h)
	{
		if (x == null)
			throw new IllegalArgumentException("Input x is null");
		if (h == null)
			throw new IllegalArgumentException("Input h is null");

		if (x.length == 0 || h.length == 0)
		{
			throw new IllegalArgumentException("Input x or h have no length");
		}
	}

	/**
	 * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution">
	 * convolution</a> between one sequence and two other sequences.
	 * <p>
	 * The solution is obtained via straightforward computation of the convolution sum (and not via FFT). Whenever the
	 * computation needs an element that would be located at an index outside the input arrays, the value is assumed to
	 * be zero.
	 * <p>
	 * This has been adapted from Apache Commons Math v3.3: org.apache.commons.math3.util.MathArrays
	 *
	 * @param x
	 *            First sequence.
	 * @param h1
	 *            Second sequence 1.
	 * @param h2
	 *            Second sequence 2.
	 * @return the convolution of {@code x} and {@code h1} and {@code x} and {@code h2}.
	 *         This array's length will be [2][{@code x.length + h1.length - 1}].
	 * @throws IllegalArgumentException
	 *             If any input is null or empty. If h1 and h2 are different lengths.
	 */
	public static double[][] convolve(double[] x, double[] h1, double[] h2) throws IllegalArgumentException
	{
		checkInput(x, h1, h2);

		final int xLen = x.length;
		final int hLen = h1.length;

		// initialize the output array
		final int totalLength = xLen + hLen - 1;
		final double[][] y = new double[2][totalLength];

		// straightforward implementation of the convolution sum
		for (int n = 0; n < totalLength; n++)
		{
			double yn1 = 0, yn2 = 0;
			int k = FastMath.max(0, n + 1 - xLen);
			int j = n - k;
			while (k < hLen && j >= 0)
			{
				yn1 += x[j] * h1[k];
				yn2 += x[j] * h2[k];
				j--;
				k++;
			}
			y[0][n] = yn1;
			y[1][n] = yn2;
		}

		return y;
	}

	/**
	 * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution">
	 * convolution</a> between one sequence and two other sequences.
	 * <p>
	 * The solution is obtained via multiplication in the frequency domain.
	 *
	 * @param x
	 *            First sequence.
	 * @param h1
	 *            Second sequence 1.
	 * @param h2
	 *            Second sequence 2.
	 * @return the convolution of {@code x} and {@code h1} and {@code x} and {@code h2}.
	 *         This array's length will be [2][{@code x.length + h1.length - 1}].
	 * @throws IllegalArgumentException
	 *             If any input is null or empty. If h1 and h2 are different lengths.
	 */
	public static double[][] convolveFFT(double[] x, double[] h1, double[] h2) throws IllegalArgumentException
	{
		checkInput(x, h1, h2);

		final int xLen = x.length;
		final int hLen = h1.length;
		final int totalLength = xLen + hLen - 1;

		// Get length to a power of 2
		int newL = CommonUtils.nextPow2(totalLength);

		// Double the new length for complex values in DoubleFFT_1D
		x = Arrays.copyOf(x, 2 * newL);
		h1 = Arrays.copyOf(h1, x.length);
		h2 = Arrays.copyOf(h2, x.length);

		DoubleFFT_1D fft = new DoubleFFT_1D(newL);

		// FFT
		fft.realForwardFull(x);
		fft.realForwardFull(h1);
		fft.realForwardFull(h2);

		// Complex multiply. Reuse data array
		for (int i = 0; i < x.length; i += 2)
		{
			int j = i + 1;
			double xi = x[i];
			double xj = x[j];
			double hi = h1[i];
			double hj = h1[j];
			h1[i] = hi * xi - hj * xj;
			h1[j] = hi * xj + hj * xi;
			hi = h2[i];
			hj = h2[j];
			h2[i] = hi * xi - hj * xj;
			h2[j] = hi * xj + hj * xi;
		}

		// Inverse FFT
		fft.complexInverse(h1, true);
		fft.complexInverse(h2, true);

		// Fill result with real part
		final double[][] y = new double[2][totalLength];
		for (int i = 0; i < totalLength; i++)
		{
			y[0][i] = h1[2 * i];
			y[1][i] = h2[2 * i];
		}
		return y;
	}

	/**
	 * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution">
	 * convolution</a> between one sequence and two other sequences.
	 * <p>
	 * The solution is obtained using either the spatial or frequency domain depending on the size. The switch is made
	 * when the min array length is above 127 and the product of the lengths is above 40000. Speed tests have
	 * been performed for single threaded FFT computation. The FFT library begins multi-threaded computation when the
	 * size of the array is above length 8192.
	 *
	 * @param x
	 *            First sequence.
	 * @param h1
	 *            Second sequence 1.
	 * @param h2
	 *            Second sequence 2.
	 * @return the convolution of {@code x} and {@code h1} and {@code x} and {@code h2}.
	 *         This array's length will be [2][{@code x.length + h1.length - 1}].
	 * @throws IllegalArgumentException
	 *             If any input is null or empty. If h1 and h2 are different lengths.
	 */
	public static double[][] convolveFast(double[] x, double[] h1, double[] h2) throws IllegalArgumentException
	{
		checkInput(x, h1, h2);
		// See Junit class ConvolveTest to determine when to switch to the FFT method.
		// This is not perfect for all length combinations but the switch will happen 
		// when the two methods are roughly the same speed.
		int min, max;
		if (x.length < h1.length)
		{
			min = x.length;
			max = h1.length;
		}
		else
		{
			min = h1.length;
			max = x.length;
		}
		if (min >= 128 && (long) min * (long) max > 40000L)
			return convolveFFT(x, h1, h2);
		return convolve(x, h1, h2);
	}

	private static void checkInput(double[] x, double[] h1, double[] h2)
	{
		if (x == null)
			throw new IllegalArgumentException("Input x is null");
		if (h1 == null)
			throw new IllegalArgumentException("Input h1 is null");
		if (h2 == null)
			throw new IllegalArgumentException("Input h2 is null");

		if (x.length == 0 || h1.length == 0)
		{
			throw new IllegalArgumentException("Input x or h1 have no length");
		}
		if (h1.length != h2.length)
		{
			throw new IllegalArgumentException("Input h1 and h2 have different length");
		}
	}

	/**
	 * Create a 1-dimensional normalized Gaussian kernel with standard deviation sigma.
	 * To avoid a step due to the cutoff at a finite value, the near-edge values are
	 * replaced by a 2nd-order polynomial with its minimum=0 at the first out-of-kernel
	 * pixel. Thus, the kernel function has a smooth 1st derivative in spite of finite
	 * length.
	 *
	 * @param sigma
	 *            Standard deviation
	 * @param range
	 *            the range
	 * @param edgeCorrection
	 *            Set to true to perform the edge correction
	 * @return The kernel, decaying towards zero, which would be reached at the first out of kernel index
	 */
	public static double[] makeGaussianKernel(final double sigma, double range, boolean edgeCorrection)
	{
		// Limit range for the Gaussian
		if (range < 1)
			range = 1;
		else if (range > 38)
			range = 38;

		// Build half the kernel into the full kernel array. This is duplicated later.
		int kRadius = getGaussianHalfWidth(sigma, range) + 1;
		double[] kernel = new double[2 * kRadius - 1];

		kernel[0] = 1;
		final double s2 = sigma * sigma;
		for (int i = 1; i < kRadius; i++)
		{
			// Gaussian function
			kernel[i] = FastMath.exp(-0.5 * i * i / s2);
		}

		// Edge correction
		if (edgeCorrection && kRadius > 3)
		{
			double sqrtSlope = Double.MAX_VALUE;
			int r = kRadius;
			while (r > kRadius / 2)
			{
				r--;
				double a = Math.sqrt(kernel[r]) / (kRadius - r);
				if (a < sqrtSlope)
					sqrtSlope = a;
				else
					break;
			}
			//System.out.printf("Edge correction: s=%.3f, kRadius=%d, r=%d, sqrtSlope=%f\n", sigma, kRadius, r,
			//		sqrtSlope);
			for (int r1 = r + 2; r1 < kRadius; r1++)
				kernel[r1] = ((kRadius - r1) * (kRadius - r1) * sqrtSlope * sqrtSlope);
		}

		// Normalise
		double sum = kernel[0];
		for (int i = 1; i < kRadius; i++)
			sum += 2 * kernel[i];
		for (int i = 0; i < kRadius; i++)
			kernel[i] /= sum;

		// Create symmetrical
		System.arraycopy(kernel, 0, kernel, kRadius - 1, kRadius);
		for (int i = kRadius, j = i - 2; i < kernel.length; i++, j--)
		{
			kernel[j] = kernel[i];
		}
		return kernel;
	}

	/**
	 * Get half the width of the region smoothed by a Gaussian filter for the specified standard deviation. The full
	 * region size is 2N + 1
	 *
	 * @param sigma
	 *            the sigma
	 * @param range
	 *            the range
	 * @return The half width
	 */
	public static int getGaussianHalfWidth(double sigma, double range)
	{
		return (int) Math.ceil(sigma * range);
	}
}