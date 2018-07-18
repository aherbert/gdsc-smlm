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

import java.util.Arrays;

import org.apache.commons.math3.util.FastMath;
import org.jtransforms.fft.DoubleFFT_1D;
import org.jtransforms.utils.CommonUtils;

/**
 * Simple class to perform convolution
 */
public class Convolution
{
	/** The maximum size supported for scaled convolution. */
	public static final int MAX = 1 << 30;

	/**
	 * Interface to handle a convolution value
	 */
	public interface ConvolutionValueProcedure
	{
		/**
		 * Executes this procedure.
		 *
		 * @param value
		 *            the value of the convolution
		 * @return true, if further values should be computed
		 */
		public boolean execute(double value);
	}

	/**
	 * Interface to handle two convolution valuess
	 */
	public interface DoubleConvolutionValueProcedure
	{
		/**
		 * Executes this procedure.
		 *
		 * @param value1
		 *            the value of the convolution of the first input
		 * @param value2
		 *            the value of the convolution of the second input
		 * @return true, if further values should be computed
		 */
		public boolean execute(double value1, double value2);
	}

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

		// Straightforward implementation of the convolution sum
		for (int n = 0; n < totalLength; n++)
		{
			double yn = 0;
			int k = FastMath.max(0, n + 1 - xLen);
			int j = n - k;
			while (k < hLen && j >= 0)
				yn += x[j--] * h[k++];
			y[n] = yn;
		}

		return y;
	}

	/**
	 * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution">
	 * convolution</a> between two sequences.
	 * <p>
	 * The solution is obtained via straightforward computation of the convolution sum (and not via FFT). Whenever the
	 * computation needs an element that would be located at an index outside the input arrays, the value is assumed to
	 * be zero.
	 * <p>
	 * The convolution is computed dynamically and can be stopped.
	 *
	 * @param x
	 *            First sequence.
	 *            Typically, this sequence will represent an input signal to a system.
	 * @param h
	 *            Second sequence.
	 *            Typically, this sequence will represent the impulse response of the system.
	 * @param v
	 *            Output procedure for the convolution of {@code x} and {@code h}.
	 *            This total number of times this is called will be {@code x.length + h.length - 1}.
	 * @throws IllegalArgumentException
	 *             if either {@code x} or {@code h} is {@code null} or either {@code x} or {@code h} is empty.
	 */
	public static void convolve(double[] x, double[] h, ConvolutionValueProcedure v) throws IllegalArgumentException
	{
		checkInput(x, h);

		if (v == null)
			throw new IllegalArgumentException("No value procedure");

		final int xLen = x.length;
		final int hLen = h.length;

		// initialize the output array
		final int totalLength = xLen + hLen - 1;
		if (totalLength <= 0)
			throw new IllegalArgumentException("Unsupported size: " + ((long) xLen + hLen - 1));

		// Straightforward implementation of the convolution sum
		for (int n = 0; n < totalLength; n++)
		{
			double yn = 0;
			int k = FastMath.max(0, n + 1 - xLen);
			int j = n - k;
			while (k < hLen && j >= 0)
				yn += x[j--] * h[k++];
			if (!v.execute(yn))
				break;
		}
	}

	/**
	 * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution">
	 * convolution</a> between two sequences.
	 * <p>
	 * The solution is obtained via multiplication in the frequency domain. To reduce edge wrap artifacts the input
	 * signals should be windowed to zero at the ends.
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
		final int newL = CommonUtils.nextPow2(totalLength);

		// Double the new length for complex values in DoubleFFT_1D
		x = Arrays.copyOf(x, 2 * newL);
		h = Arrays.copyOf(h, x.length);

		final DoubleFFT_1D fft = new DoubleFFT_1D(newL);

		// FFT
		fft.realForwardFull(x);
		fft.realForwardFull(h);

		// Complex multiply. Reuse data array
		for (int i = 0; i < x.length; i += 2)
		{
			final int j = i + 1;
			final double xi = x[i];
			final double xj = x[j];
			final double hi = h[i];
			final double hj = h[j];
			h[i] = hi * xi - hj * xj;
			h[j] = hi * xj + hj * xi;
		}

		// Inverse FFT
		fft.complexInverse(h, true);

		// Fill result with real part
		final double[] y = new double[totalLength];
		for (int i = 0; i < totalLength; i++)
			y[i] = h[2 * i];
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

	/**
	 * Checks if convolution will use the FFT method.
	 *
	 * @param length1
	 *            the length 1
	 * @param length2
	 *            the length 2
	 * @return true, if using the FFT method
	 */
	public static boolean isFFT(int length1, int length2)
	{
		int min, max;
		if (length1 < length2)
		{
			min = length1;
			max = length2;
		}
		else
		{
			min = length2;
			max = length1;
		}
		if (min >= 128 && (long) min * (long) max > 40000L)
			return true;
		return false;
	}

	private static void checkInput(double[] x, double[] h)
	{
		if (x == null)
			throw new IllegalArgumentException("Input x is null");
		if (h == null)
			throw new IllegalArgumentException("Input h is null");

		if (x.length == 0 || h.length == 0)
			throw new IllegalArgumentException("Input x or h have no length");
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

		// Straightforward implementation of the convolution sum
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
	 * @param v
	 *            Output procedure for the convolution of {@code x} and {@code h1} or {@code h2}.
	 *            This total number of times this is called will be {@code x.length + h1.length - 1}.
	 * @throws IllegalArgumentException
	 *             If any input is null or empty. If h1 and h2 are different lengths.
	 */
	public static void convolve(double[] x, double[] h1, double[] h2, DoubleConvolutionValueProcedure v)
			throws IllegalArgumentException
	{
		checkInput(x, h1, h2);

		if (v == null)
			throw new IllegalArgumentException("No value procedure");

		final int xLen = x.length;
		final int hLen = h1.length;

		// initialize the output array
		final int totalLength = xLen + hLen - 1;
		if (totalLength <= 0)
			throw new IllegalArgumentException("Unsupported size: " + ((long) xLen + hLen - 1));

		// Straightforward implementation of the convolution sum
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
			if (!v.execute(yn1, yn2))
				break;
		}
	}

	/**
	 * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution">
	 * convolution</a> between one sequence and two other sequences.
	 * <p>
	 * The solution is obtained via multiplication in the frequency domain. To reduce edge wrap artifacts the input
	 * signals should be windowed to zero at the ends.
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
		final int newL = CommonUtils.nextPow2(totalLength);

		// Double the new length for complex values in DoubleFFT_1D
		x = Arrays.copyOf(x, 2 * newL);
		h1 = Arrays.copyOf(h1, x.length);
		h2 = Arrays.copyOf(h2, x.length);

		final DoubleFFT_1D fft = new DoubleFFT_1D(newL);

		// FFT
		fft.realForwardFull(x);
		fft.realForwardFull(h1);
		fft.realForwardFull(h2);

		// Complex multiply. Reuse data array
		for (int i = 0; i < x.length; i += 2)
		{
			final int j = i + 1;
			final double xi = x[i];
			final double xj = x[j];
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
			throw new IllegalArgumentException("Input x or h1 have no length");
		if (h1.length != h2.length)
			throw new IllegalArgumentException("Input h1 and h2 have different length");
	}

	/**
	 * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution">
	 * convolution</a> between two sequences.
	 * <p>
	 * The scale is used to increase the size of h dynamically to H with zero fill.
	 * The length of H is thus ((h.length-1) * scale + 1);
	 *
	 * @param x
	 *            First sequence.
	 * @param h
	 *            Second sequence.
	 * @param scale
	 *            the scale
	 * @return the convolution of {@code x} and {@code h}.
	 *         This array's length will be {@code x.length + H.length - 1}.
	 * @throws IllegalArgumentException
	 *             if either {@code x} or {@code h} is {@code null} or either {@code x} or {@code h} is empty.
	 * @throws IllegalArgumentException
	 *             if the scale is not strictly positive
	 */
	public static double[] convolve(double[] x, double[] h, int scale) throws IllegalArgumentException
	{
		checkInput(x, h);

		if (scale < 1)
			throw new IllegalArgumentException("Scale must be strictly positive");

		if (h.length == 1 || scale == 1)
			// No scaling
			return convolve(x, h);

		final int xLen = x.length;
		final int hLen = h.length;

		// initialize the output array
		final double HLen = (double) (h.length - 1) * scale + 1;
		final double totalLength = xLen + HLen - 1;
		if (totalLength > MAX)
			throw new IllegalArgumentException("Scale creates unsupported array size: " + totalLength);

		final int iTotalLength = (int) totalLength;
		final double[] y = new double[iTotalLength];

		// Convolution sum. x is reversed verses h.
		// h is scaled up with zeros.
		// This is equivalent to using x every interval of scale.
		for (int n = 0; n < iTotalLength; n++)
		{
			double yn = 0;
			// k is the index in the scaled up distribution H
			final int k = FastMath.max(0, n + 1 - xLen);
			// j is the index in the input distribution x
			int j = n - k;

			// k has to be scaled.
			// The modulus indicates how many values are zero
			// before the first non-zero value in H (in the descending direction).
			final int mod = k % scale;
			// kk is the index in input distribution h (in the descending direction).
			int kk = k / scale;
			// If there are non-zero value shift the indices
			if (mod != 0)
			{
				// Shift kk by one for the next non-zero value (in the ascending direction)
				kk++;
				// Shift j by the number of zero values (in the descending direction)
				j -= (scale - mod);
			}

			//int j = n - kk * scale;
			while (kk < hLen && j >= 0)
			{
				yn += x[j] * h[kk++];
				j -= scale;
			}
			y[n] = yn;
		}

		return y;
	}

	/**
	 * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution">
	 * convolution</a> between two sequences.
	 * <p>
	 * The scale is used to increase the size of h dynamically to H with zero fill.
	 * The length of H is thus ((h.length-1) * scale + 1);
	 * <p>
	 * The convolution is computed dynamically and can be stopped.
	 *
	 * @param x
	 *            First sequence.
	 * @param h
	 *            Second sequence.
	 * @param scale
	 *            the scale
	 * @param v
	 *            Output procedure for the convolution of {@code x} and {@code h}.
	 *            This total number of times this is called will be {@code x.length + H.length - 1}.
	 * @throws IllegalArgumentException
	 *             if either {@code x} or {@code h} is {@code null} or either {@code x} or {@code h} is empty.
	 * @throws IllegalArgumentException
	 *             if the scale is not strictly positive
	 */
	public static void convolve(double[] x, double[] h, int scale, ConvolutionValueProcedure v)
			throws IllegalArgumentException
	{
		// As above but dynamically output the result

		checkInput(x, h);

		if (scale < 1)
			throw new IllegalArgumentException("Scale must be strictly positive");
		if (v == null)
			throw new IllegalArgumentException("No value procedure");

		final int xLen = x.length;
		final int hLen = h.length;

		// For consistency just support up to the max for integers.
		// This could be changed to use long for the index.
		final double HLen = (double) (h.length - 1) * scale + 1;
		final double totalLength = xLen + HLen - 1;
		if (totalLength > MAX)
			throw new IllegalArgumentException("Scale creates unsupported size: " + totalLength);

		final int iTotalLength = (int) totalLength;

		for (int n = 0; n < iTotalLength; n++)
		{
			double yn = 0;
			final int k = FastMath.max(0, n + 1 - xLen);
			int j = n - k;
			final int mod = k % scale;
			int kk = k / scale;
			if (mod != 0)
			{
				kk++;
				j -= (scale - mod);
			}
			while (kk < hLen && j >= 0)
			{
				yn += x[j] * h[kk++];
				j -= scale;
			}
			if (!v.execute(yn))
				break;
		}
	}

	/**
	 * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution">
	 * convolution</a> between two sequences.
	 * <p>
	 * The scale is used to increase the size of h dynamically to H with zero fill.
	 * The length of H is thus ((h.length-1) * scale + 1);
	 *
	 * @param x
	 *            First sequence.
	 * @param h1
	 *            Second sequence 1.
	 * @param h2
	 *            Second sequence 2.
	 * @param scale
	 *            the scale
	 * @return the convolution of {@code x} and {@code h1} and {@code x} and {@code h2}.
	 *         This array's length will be [2][{@code x.length + h1.length - 1}].
	 * @throws IllegalArgumentException
	 *             If any input is null or empty. If h1 and h2 are different lengths.
	 * @throws IllegalArgumentException
	 *             if the scale is not strictly positive
	 */
	public static double[][] convolve(double[] x, double[] h1, double[] h2, int scale) throws IllegalArgumentException
	{
		checkInput(x, h1, h2);

		if (scale < 1)
			throw new IllegalArgumentException("Scale must be strictly positive");

		if (h1.length == 1 || scale == 1)
			// No scaling
			return convolve(x, h1, h2);

		final int xLen = x.length;
		final int hLen = h1.length;

		// initialize the output array
		final double HLen = (double) (h1.length - 1) * scale + 1;
		final double totalLength = xLen + HLen - 1;
		if (totalLength > MAX)
			throw new IllegalArgumentException("Scale creates unsupported array size: " + totalLength);

		final int iTotalLength = (int) totalLength;
		final double[][] y = new double[2][iTotalLength];

		// Convolution sum. x is reversed verses h.
		// h is scaled up with zeros.
		// This is equivalent to using x every interval of scale.
		for (int n = 0; n < iTotalLength; n++)
		{
			double yn1 = 0, yn2 = 0;
			// k is the index in the scaled up distribution H
			final int k = FastMath.max(0, n + 1 - xLen);
			// j is the index in the input distribution x
			int j = n - k;

			// k has to be scaled.
			// The modulus indicates how many values are zero
			// before the first non-zero value in H (in the descending direction).
			final int mod = k % scale;
			// kk is the index in input distribution h (in the descending direction).
			int kk = k / scale;
			// If there are non-zero value shift the indices
			if (mod != 0)
			{
				// Shift kk by one for the next non-zero value (in the ascending direction)
				kk++;
				// Shift j by the number of zero values (in the descending direction)
				j -= (scale - mod);
			}

			//int j = n - kk * scale;
			while (kk < hLen && j >= 0)
			{
				yn1 += x[j] * h1[kk];
				yn2 += x[j] * h2[kk];
				kk++;
				j -= scale;
			}
			y[0][n] = yn1;
			y[1][n] = yn2;
		}

		return y;
	}

	/**
	 * Calculates the <a href="http://en.wikipedia.org/wiki/Convolution">
	 * convolution</a> between two sequences.
	 * <p>
	 * The scale is used to increase the size of h dynamically to H with zero fill.
	 * The length of H is thus ((h1.length-1) * scale + 1);
	 * <p>
	 * The convolution is computed dynamically and can be stopped.
	 *
	 * @param x
	 *            First sequence.
	 * @param h1
	 *            Second sequence 1.
	 * @param h2
	 *            Second sequence 2.
	 * @param scale
	 *            the scale
	 * @param v
	 *            Output procedure for the convolution of {@code x} and {@code h1} or {@code h2}.
	 *            This total number of times this is called will be {@code x.length + h1.length - 1}.
	 * @throws IllegalArgumentException
	 *             If any input is null or empty. If h1 and h2 are different lengths.
	 * @throws IllegalArgumentException
	 *             if the scale is not strictly positive
	 * @throws IllegalArgumentException
	 *             if the output size is above the max size supported
	 */
	public static void convolve(double[] x, double[] h1, double[] h2, int scale, DoubleConvolutionValueProcedure v)
			throws IllegalArgumentException
	{
		checkInput(x, h1, h2);

		if (scale < 1)
			throw new IllegalArgumentException("Scale must be strictly positive");
		if (v == null)
			throw new IllegalArgumentException("No value procedure");

		final int xLen = x.length;
		final int hLen = h1.length;

		// For consistency just support up to the max for integers.
		// This could be changed to use long for the index.
		final double HLen = (double) (h1.length - 1) * scale + 1;
		final double totalLength = xLen + HLen - 1;
		if (totalLength > MAX)
			throw new IllegalArgumentException("Scale creates unsupported size: " + totalLength);

		final int iTotalLength = (int) totalLength;

		for (int n = 0; n < iTotalLength; n++)
		{
			double yn1 = 0, yn2 = 0;
			final int k = FastMath.max(0, n + 1 - xLen);
			int j = n - k;
			final int mod = k % scale;
			int kk = k / scale;
			if (mod != 0)
			{
				kk++;
				j -= (scale - mod);
			}
			while (kk < hLen && j >= 0)
			{
				yn1 += x[j] * h1[kk];
				yn2 += x[j] * h2[kk];
				kk++;
				j -= scale;
			}
			if (!v.execute(yn1, yn2))
				break;
		}
	}
}
