package gdsc.smlm.filters;

import org.jtransforms.dht.FloatDHT_2D;
import org.jtransforms.utils.CommonUtils;

import gdsc.core.utils.ImageWindow;
import gdsc.core.utils.Maths;
import ij.process.FHT2;

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

	private FloatDHT_2D dht = null;
	private FHT2 kernelFht = null;
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

		FHT2 dataFht = createFHT(data, maxx, maxy, border);
		int maxN = kernelFht.getWidth();

		FHT2 result = (convolution) ? dataFht.multiply(kernelFht, tmp) : dataFht.conjugateMultiply(kernelFht, tmp);

		// Do the transform using JTransforms as it is faster
		dht.inverse(result.getData(), true);

		//result.inverseTransform();

		result.swapQuadrants();
		if (maxx < maxN || maxy < maxN)
		{
			int x = getInsert(maxN, maxx);
			int y = getInsert(maxN, maxy);
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
	 * Initialise the kernel FHT. It is recreated only if the target size has changed.
	 *
	 * @param maxx
	 *            the width of the target image
	 * @param maxy
	 *            the height of the target image
	 */
	public void initialiseKernel(int maxx, int maxy)
	{
		int maxN = Maths.nextPow2(Maths.max(maxx, maxy, kN));
		if (tmp == null || tmp.length != maxN * maxN)
			tmp = new float[maxN * maxN];
		if (kernelFht != null && maxN == kernelFht.getWidth())
			// Already initialised
			return;
		int size = maxN * maxN;
		// No window function for the kernel so just create a new FHT
		float[] data;
		if (kw < maxN || kh < maxN)
		{
			// Too small so insert in the middle
			data = new float[size];
			int x = getInsert(maxN, kw);
			int y = getInsert(maxN, kh);
			insert(kernel, kw, kh, data, maxN, x, y);
		}
		else
		{
			// Clone to avoid destroying data
			data = kernel.clone();
		}

		// Do the transform using JTransforms as it is faster. Do not allow multi-threading.
		CommonUtils.setThreadsBeginN_2D(Long.MAX_VALUE);
		dht = new FloatDHT_2D(maxN, maxN);
		CommonUtils.resetThreadsBeginN();
		dht.forward(data);
		kernelFht = new FHT2(data, maxN, true);

		//kernelFht = new FHT2(data, maxN, false);
		//kernelFht.transform();		

		kernelFht.initialiseFastMultiply();

		// This is used for the output complex multiple of the two FHTs
		tmp = new float[size];
	}

	private static int getInsert(int maxN, int width)
	{
		// Note the FHT power spectrum centre is at n/2 of an even sized image.
		// So we must insert the centre at that point. To do this we check for odd/even
		// and offset if necessary. This allows the FHTFilter to match the correlation
		// result from a filter in the spatial domain performed using the KernelFilter class. 
		int diff = maxN - width;
		return ((diff & 1) == 1) ? (diff + 1) / 2 : diff / 2;
	}

	private static void insert(float[] source, int maxx, int maxy, float[] dest, int maxN, int x, int y)
	{
		for (int from = 0, to = y * maxN + x; from < source.length; from += maxx, to += maxN)
		{
			System.arraycopy(source, from, dest, to, maxx);
		}
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

		int maxN = kernelFht.getWidth();
		if (maxx < maxN || maxy < maxN)
		{
			// Too small so insert in the middle of a new processor
			float[] data2 = new float[maxN * maxN];
			int x = getInsert(maxN, maxx);
			int y = getInsert(maxN, maxy);
			insert(data, maxx, maxy, data2, maxN, x, y);
			data = data2;
		}

		// Do the transform using JTransforms as it is faster
		dht.forward(data);
		FHT2 result = new FHT2(data, maxN, true);

		//FHT2 result = new FHT2(data, maxN, false);
		// Copy the initialised tables
		//result.copyTables(kernelFht);
		//result.transform();

		return result;
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
		// Border will be positive, ensure no overflow of kernel size
		border = Maths.min(border, kw / 2, kh / 2);

		if (edge != border)
		{
			edge = border;
			// Cache the whole window but we only use part of it.
			w = ImageWindow.tukeyEdge(Math.min(kw, kh), border);
		}

		// Assume that the border will only be a fraction of the image and perform 
		// selective weighting

		for (int b = 0, ri = -1, rj = kernel.length; b < border; b++)
		{
			double weight = w[b];
			// index = y * kw + x
			// Rows (ri|rj): y = b or kh-b-1; x = 0
			// Columns (ri|rj): y = 0; x = b or kw-b-1
			for (int i = 0, ci = b, cj = kw - b - 1; i < kh; i++, ci += kw, cj += kw)
			{
				kernel[++ri] *= weight;
				kernel[--rj] *= weight;
				kernel[ci] *= weight;
				kernel[cj] *= weight;
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
		FHTFilter fht = (FHTFilter) super.clone();
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

	/**
	 * Checks if is a convolution filter. The default is correlation.
	 *
	 * @return true, if is a convolution filter
	 */
	public boolean isConvolution()
	{
		return convolution;
	}

	/**
	 * Sets the convolution flag. The default is correlation.
	 *
	 * @param convolution
	 *            the new convolution flag
	 */
	public void setConvolution(boolean convolution)
	{
		this.convolution = convolution;
	}
}