package gdsc.smlm.utils;

import java.util.Arrays;

import org.apache.commons.math3.util.FastMath;

import gdsc.core.utils.Maths;
import gnu.trove.list.array.TDoubleArrayList;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Store a Gaussian kernel for use in convolution
 */
public class GaussianKernel implements Cloneable
{
	/** The maximum size of half the Gaussian kernel. */
	public static final int HALF_WIDTH_LIMIT = 1 << 28;

	/** The standard deviation. */
	final public double s;

	final private double var2;

	private int currentScale;

	private TDoubleArrayList halfKernel;

	/**
	 * Instantiates a new gaussian kernel.
	 *
	 * @param s
	 *            The standard deviation
	 */
	public GaussianKernel(double s)
	{
		// Check not infinite or NaN
		if (!(s > 0 && s < Double.MAX_VALUE))
			throw new IllegalArgumentException("Standard deviation must be positive");
		this.s = s;
		currentScale = 1;
		halfKernel = new TDoubleArrayList();
		// Initialise the first value exp(-0)
		halfKernel.add(1);
		// Precompute exponential scaling factor
		var2 = -(s * s * 2);
	}

	/**
	 * Create a 1-dimensional normalized Gaussian kernel with standard deviation scale*s.
	 * To avoid a step due to the cutoff at a finite value, the near-edge values are
	 * replaced by a 2nd-order polynomial with its minimum=0 at the first out-of-kernel
	 * pixel. Thus, the kernel function has a smooth 1st derivative in spite of finite
	 * length.
	 *
	 * @param scale
	 *            the scale (must be a power of 2)
	 * @param range
	 *            the range (in units of standard deviation)
	 * @param edgeCorrection
	 *            Set to true to perform the edge correction
	 * @return The kernel, decaying towards zero, which would be reached at the first out of kernel index
	 */
	public double[] getGaussianKernel(int scale, double range, boolean edgeCorrection)
	{
		if (!Maths.isPow2(scale))
			throw new IllegalArgumentException("Scale must be a power of 2: " + scale);

		// Limit range for the Gaussian
		if (range < 1)
			range = 1;
		else if (range > 38)
			range = 38;

		int kRadius = getGaussianHalfWidth(s * scale, range) + 1;
		increaseScale(scale);
		increaseKernel(kRadius);
		//increaseKernel(scale, kRadius);

		// Create kernel
		// Note: The stored values in the halfKernel are always non-zero.
		double[] kernel = new double[2 * kRadius - 1];
		kernel[0] = 1;
		if (currentScale == scale)
		{
			for (int i = 1, size = Math.min(kRadius, halfKernel.size()); i < size; i++)
			{
				kernel[i] = halfKernel.getQuick(i);
			}
		}
		else
		{
			final int upsample = currentScale / scale;
			final double step = 1.0 / scale;
			for (int i = 1, j = upsample; i < kRadius; i++, j += upsample)
			{
				// Just in case down-scaling requires a different end point in the kernel
				// check the size
				if (j < halfKernel.size())
				{
					kernel[i] = halfKernel.getQuick(j);
				}
				else
				{
					kernel[i] = FastMath.exp(Maths.pow2(i * step) / var2);
					// Check if zero
					if (kernel[i] == 0)
						break;
				}
			}
		}

		return buildKernel(kernel, kRadius, edgeCorrection);
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
		double limit = Math.ceil(sigma * range);
		// Ensure the kernel is clipped to the size of an array
		return (limit < HALF_WIDTH_LIMIT) ? (int) limit : HALF_WIDTH_LIMIT;
	}

	private void increaseScale(int scale)
	{
		if (currentScale < scale)
		{
			// Up sample the current kernel
			final int upsample = scale / currentScale;
			currentScale = scale;

			final double[] g = halfKernel.toArray();
			final int kRadius = g.length;
			final double step = 1.0 / currentScale;
			halfKernel.resetQuick();
			halfKernel.add(1);

			for (int i = 1, j = 0; i < kRadius; i++, j += upsample)
			{
				for (int k = 1; k < upsample; k++)
					halfKernel.add(FastMath.exp(Maths.pow2((j + k) * step) / var2));
				halfKernel.add(g[i]);
			}
		}
	}

	private void increaseKernel(int kRadius)
	{
		if (halfKernel.size() < kRadius)
		{
			final double step = 1.0 / currentScale;
			for (int i = halfKernel.size(); i < kRadius; i++)
			{
				double v = FastMath.exp(Maths.pow2(i * step) / var2);
				if (v == 0)
					break;
				halfKernel.add(v);
			}
		}
	}

	@SuppressWarnings("unused")
	private void increaseKernel(int scale, int kRadius)
	{
		if (currentScale < scale)
		{
			currentScale = scale;
			halfKernel.resetQuick();
			halfKernel.add(1);
		}

		if (halfKernel.size() < kRadius)
		{
			final double step = 1.0 / currentScale;
			for (int i = halfKernel.size(); i < kRadius; i++)
			{
				double v = FastMath.exp(Maths.pow2(i * step) / var2);
				if (v == 0)
					break;
				halfKernel.add(v);
			}
		}
	}

	private static double[] buildKernel(double[] kernel, int kRadius, boolean edgeCorrection)
	{
		// Clip in the event that zeros occurred during computation
		if (kernel[kRadius - 1] == 0)
		{
			while (kernel[--kRadius] == 0)
				;
			if (kRadius == 1)
				return new double[] { 1 };
			kernel = Arrays.copyOf(kernel, 2 * kRadius - 1);
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
			if (kernel[i] == 0)
				break;
		}

		return buildKernel(kernel, kRadius, edgeCorrection);
	}

	/**
	 * Create a 1-dimensional normalized Gaussian kernel with standard deviation sigma.
	 * The kernel is constructed using the Error function (Erf) to compute the sum of
	 * the Gaussian from x-0.5 to x+0.5 for each x sample point.
	 *
	 * @param sigma
	 *            Standard deviation
	 * @param range
	 *            the range
	 * @return The Erf kernel
	 */
	public static double[] makeErfGaussianKernel(double sigma, double range)
	{
		// Limit range for the Gaussian
		if (range < 1)
			range = 1;
		else if (range > 38)
			range = 38;

		// Build half the kernel into the full kernel array. This is duplicated later.
		int kRadius = getGaussianHalfWidth(sigma, range) + 1;
		double[] kernel = new double[2 * kRadius - 1];

		if (kRadius == 1)
		{
			kernel[0] = 1;
			return kernel;
		}

		// Use the error function to get the integral of the Gaussian.
		final double sqrt_var_by_2 = Math.sqrt(sigma * sigma * 2);

		double upper = org.apache.commons.math3.special.Erf.erf(-0.5 / sqrt_var_by_2);
		for (int i = 0; i < kRadius; i++)
		{
			double lower = upper;
			upper = org.apache.commons.math3.special.Erf.erf((i + 0.5) / sqrt_var_by_2);
			kernel[i] = (upper - lower) * 0.5;
			if (kernel[i] == 0)
				break;
		}

		return buildKernel(kernel, kRadius, false);
	}

	@Override
	public GaussianKernel clone()
	{
		try
		{
			GaussianKernel k = (GaussianKernel) super.clone();
			k.halfKernel = new TDoubleArrayList(this.halfKernel);
			return k;
		}
		catch (CloneNotSupportedException e)
		{
			return new GaussianKernel(s);
		}
	}
}