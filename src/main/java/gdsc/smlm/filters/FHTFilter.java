package gdsc.smlm.filters;

import gdsc.core.utils.ImageWindow;
import gdsc.core.utils.Maths;
import ij.process.FHT2;
import ij.process.FloatProcessor;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Computes a convolution/correlation in the frequency domain using a Fast Hartley Tranform. An option edge window
 * function can be applied.
 */
public class FHTFilter extends BaseFilter
{
	private final float[] kernel;
	private final int kw;
	private final int kh;
	private final int kN; // Next power of 2 for the kernel

	private FHT2 fht = null;
	private float[] tmp;
	private boolean convolution = false;

	// Cache the window function
	private double[] w = null;
	private int edge = -1;

	/**
	 * Instantiates a new FHT filter. It is assumed that the kernel has an appropriate window function applied to the
	 * edges.
	 *
	 * @param kernel
	 *            the kernel
	 * @param kw
	 *            the kernel width
	 * @param kh
	 *            the kernel height
	 * @throws IllegalArgumentException
	 *             if the kernel width or height does not match the kernel size
	 */
	public FHTFilter(float[] kernel, int kw, int kh) throws IllegalArgumentException
	{
		checkKernel(kernel, kw, kh);
		this.kernel = kernel.clone();
		this.kw = kw;
		this.kh = kh;
		kN = Maths.nextPow2(Math.max(kw, kh));
		double scale = getScale(kernel);
		// Scale
		if (scale != 1)
		{
			for (int i = 0; i < kernel.length; i++)
				kernel[i] *= scale;
		}
	}

	private static void checkKernel(float[] kernel, int kw, int kh)
	{
		if (kw < 1 || kh < 1)
			throw new IllegalArgumentException("Kernel width & height must be positive");
		if (kw * kh != kernel.length)
			throw new IllegalArgumentException("Kernel width x height != kernel length");
	}

	/**
	 * Gets the scale of the kernel (i.e. 1/sum). This is used to scale the kernel to sum to 1.
	 *
	 * @param kernel
	 *            the kernel
	 * @return the scale
	 */
	public static double getScale(float[] kernel)
	{
		double sum = 0.0;
		for (int i = 0; i < kernel.length; i++)
			sum += kernel[i];
		if (sum != 0.0)
			return 1.0 / sum;
		return 1.0;
	}

	/**
	 * Compute the filter.
	 * <p>
	 * Note: the input data is destructively modified
	 *
	 * @param data
	 *            The input/output data (packed in YX order)
	 * @param maxx
	 *            The width of the data
	 * @param maxy
	 *            The height of the data
	 */
	public void filter(float[] data, final int maxx, final int maxy)
	{
		filterInternal(data, maxx, maxy, 0);
	}

	/**
	 * Compute the filter.
	 * <p>
	 * Note: the input data is destructively modified
	 *
	 * @param data
	 *            the data
	 * @param maxx
	 *            the maxx
	 * @param maxy
	 *            the maxy
	 * @param border
	 *            the border to use for the Tukey edge window
	 */
	public void filter(float[] data, final int maxx, final int maxy, int border)
	{
		if (border < 0)
		{
			border = 0;
		}
		else
		{
			border = Maths.min(border, maxx / 2, maxy / 2);
		}
		filterInternal(data, maxx, maxy, border);
	}

	/**
	 * Compute the filter.
	 * <p>
	 * Note: the input data is destructively modified
	 *
	 * @param data
	 *            the data
	 * @param maxx
	 *            the maxx
	 * @param maxy
	 *            the maxy
	 * @param border
	 *            the border
	 */
	private void filterInternal(float[] data, final int maxx, final int maxy, int border)
	{
		initialiseKernel(maxx, maxy);

		FHT2 fht2 = createFHT(data, maxx, maxy, border);
		int maxN = fht.getWidth();

		FHT2 result = (convolution) ? fht2.multiply(fht, tmp) : fht2.conjugateMultiply(fht, tmp);
		// Transform using the kernel FHT with precomputed tables
		this.fht.rc2DFHT(tmp, true, maxN);
		// result.inverseTransform();
		result.swapQuadrants();
		if (maxx < maxN || maxy < maxN)
		{
			int x = (maxN - maxx) / 2;
			int y = (maxN - maxy) / 2;
			for (int to = 0, from = y * maxN + x; to < data.length; to += maxx, from += maxN)
			{
				System.arraycopy(tmp, from, data, to, maxx);
			}
		}
		else
		{
			System.arraycopy(tmp, 0, data, 0, tmp.length);
		}
	}

	/**
	 * Initialise kernel. It is recreated only if the FHT size has changed.
	 *
	 * @param maxx
	 *            the width of the target image
	 * @param maxy
	 *            the height of the target image
	 */
	private void initialiseKernel(int maxx, int maxy)
	{
		int maxN = Maths.nextPow2(Maths.max(maxx, maxy, kN));
		if (tmp == null || tmp.length != maxN * maxN)
			tmp = new float[maxN * maxN];
		if (fht != null && maxN == fht.getWidth())
			// Already initialised
			return;
		int size = maxN * maxN;
		// No window function for the kernel so just create a new FHT
		if (kw < maxN || kh < maxN)
		{
			// Too small so insert in the middle
			fht = new FHT2(new float[size], maxN, false);
			int x = (maxN - maxx) / 2;
			int y = (maxN - maxy) / 2;
			fht.insert(new FloatProcessor(kw, kh, kernel), x, y);
		}
		else
		{
			// Clone to avoid destroying data
			fht = new FHT2(kernel.clone(), maxN, false);
		}
		fht.transform();
		// This is used for the output complex multiple of the two FHTs
		tmp = new float[size];
	}

	/**
	 * Creates the FHT.
	 *
	 * @param data
	 *            the image data
	 * @param maxx
	 *            the image width
	 * @param maxy
	 *            the image height
	 * @param border
	 *            the border
	 * @return the fht2
	 */
	private FHT2 createFHT(float[] data, int maxx, int maxy, int border)
	{
		if (border != 0)
			applyBorderInternal(data, maxx, maxy, border);

		int maxN = fht.getWidth();
		if (maxx < maxN || maxy < maxN)
		{
			// Too small so insert in the middle of a new processor
			int x = (maxN - maxx) / 2;
			int y = (maxN - maxy) / 2;
			FloatProcessor fp = new FloatProcessor(maxN, maxN);
			fp.insert(new FloatProcessor(maxx, maxy, data), x, y);
			data = (float[]) fp.getPixels();
		}

		// Tranform using the kernel FHT using the precomputed tables
		this.fht.rc2DFHT(data, false, maxN);
		return new FHT2(data, maxN, true);
	}

	/**
	 * Apply a border to the kernel image using a Tukey window.
	 *
	 * @param kernel
	 *            the kernel
	 * @param kw
	 *            the kernel width
	 * @param kh
	 *            the kernel height
	 * @param border
	 *            the border
	 * @throws IllegalArgumentException
	 *             if the kernel width or height does not match the kernel size
	 */
	public void applyBorder(float[] kernel, int kw, int kh, int border)
	{
		if (border <= 0)
			return;
		checkKernel(kernel, kw, kh);
		applyBorderInternal(kernel, kw, kh, border);
	}

	/**
	 * Apply a border to the kernel image using a Tukey window.
	 *
	 * @param kernel
	 *            the kernel
	 * @param kw
	 *            the kernel width
	 * @param kh
	 *            the kernel height
	 * @param border
	 *            the border (must be positive)
	 */
	private void applyBorderInternal(float[] kernel, int kw, int kh, int border)
	{
		if (edge != border)
		{
			edge = border;
			// Get the whole window but we only use part of it.
			w = ImageWindow.tukeyEdge(Math.min(kw, kh), border);
		}

		// Assume that the border will only be a fraction of the image and perform 
		// selective weighting

		for (int b = 0; b < border; b++)
		{
			// index = y * kw + x

			// Rows: y = b or kh-b; x = 0
			double weight = w[b];
			for (int i = 0, lx = b * kw, ux = (kh - b) * kw; i < kw; i++, lx++, ux++)
			{
				kernel[lx] *= weight;
				kernel[ux] *= weight;
			}
			// Columns: y = 0; x = b or kw-b
			for (int i = 0, ly = b, uy = kw - b; i < kh; i++, ly += kw, uy += kw)
			{
				kernel[ly] *= weight;
				kernel[uy] *= weight;
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public FHTFilter clone()
	{
		FHTFilter fht = (FHTFilter) clone();
		fht.tmp = null; // This cannot be shared across instances
		return fht;
	}

	/**
	 * Gets the kernal.
	 *
	 * @return the kernal
	 */
	public float[] getKernal()
	{
		return kernel.clone();
	}

	/**
	 * Gets the kernal width.
	 *
	 * @return the kernal width
	 */
	public int getKernalWidth()
	{
		return kw;
	}

	/**
	 * Gets the kernal height.
	 *
	 * @return the kernal height
	 */
	public int getKernalHeight()
	{
		return kh;
	}

	public boolean isConvolution()
	{
		return convolution;
	}

	public void setConvolution(boolean convolution)
	{
		this.convolution = convolution;
	}
}