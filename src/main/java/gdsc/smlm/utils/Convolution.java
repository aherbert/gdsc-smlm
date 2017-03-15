package gdsc.smlm.utils;

import java.util.Arrays;

import org.apache.commons.math3.util.FastMath;

import org.jtransforms.fft.DoubleFFT_1D;

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

		if (x.length < h.length)
		{
			// Swap so that the longest array is the signal
			final double[] tmp = x;
			x = h;
			h = tmp;
		}

		final int xLen = x.length;
		final int hLen = h.length;
		final int totalLength = xLen + hLen - 1;

		// Get length to a power of 2
		int newL = 2;
		while (newL < totalLength) //xLen)
			newL *= 2;

		// Double the new length for complex values in DoubleFFT_1D
		x = Arrays.copyOf(x, 2 * newL);
		h = Arrays.copyOf(h, x.length);
		double[] tmp = new double[x.length];

		DoubleFFT_1D fft = new DoubleFFT_1D(newL);

		// FFT
		fft.realForwardFull(x);
		fft.realForwardFull(h);

		// Complex multiply
		for (int i = 0; i < newL; i++)
		{
			final int ii = 2 * i;
			tmp[ii] = x[ii] * h[ii] - x[ii + 1] * h[ii + 1];
			tmp[ii + 1] = x[ii] * h[ii + 1] + x[ii + 1] * h[ii];
		}

		// Inverse FFT
		fft.complexInverse(tmp, true);

		// Fill result with real part
		final double[] y = new double[totalLength];
		for (int i = 0; i < totalLength; i++)
		{
			y[i] = tmp[2 * i];
		}
		return y;
	}

	private static void checkInput(double[] x, double[] h)
	{
		if (x == null)
			throw new IllegalArgumentException("Input x is null");
		if (h == null)
			throw new IllegalArgumentException("Input g is null");

		if (x.length == 0 || h.length == 0)
		{
			throw new IllegalArgumentException("Input x or h have no length");
		}
	}
}