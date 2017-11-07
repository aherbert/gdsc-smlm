package gdsc.smlm.ij.utils;

import org.jtransforms.dht.FloatDHT_3D;

import ij.ImageStack;
import ij.process.ImageProcessor;
import pl.edu.icm.jlargearrays.LargeArray;

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
 * Wrapper to compute the discrete Hartley transform on 3D data. This uses the JTransforms library.
 */
public class DHT3D
{
	public final int ns, nr, nc;

	private boolean isFrequencyDomain;
	private final FloatDHT_3D dht;
	private final float[] data;

	/**
	 * Instantiates a new 3D discrete Hartley transform
	 *
	 * @param stack
	 *            the stack
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if the combined dimensions is too large for an array
	 */
	public DHT3D(ImageStack stack) throws IllegalArgumentException
	{
		ns = stack.getSize();
		nr = stack.getHeight();
		nc = stack.getWidth();

		long size = (long) ns * nr * nc;
		if (size > LargeArray.getMaxSizeOf32bitArray())
			throw new IllegalArgumentException("3D data too large");

		dht = new FloatDHT_3D(ns, nr, nc);

		data = new float[(int) size];

		int nr_nc = nr * nc;
		if (stack.getBitDepth() == 32)
		{
			for (int s = 0; s < ns; s++)
			{
				System.arraycopy(stack.getPixels(s + 1), 0, data, s * nr_nc, nr_nc);
			}
		}
		else
		{
			for (int s = 1, i = 0; s <= ns; s++)
			{
				ImageProcessor ip = stack.getProcessor(s);
				for (int j = 0; i < nr_nc; j++)
					data[i++] = ip.getf(j);
			}
		}
	}

	private DHT3D(float[] data, int ns, int nr, int nc, boolean isFrequencyDomain, FloatDHT_3D dht)
	{
		this.ns = ns;
		this.nr = nr;
		this.nc = nc;
		this.dht = dht; // This can be reused across objects
		this.data = data;
		this.isFrequencyDomain = isFrequencyDomain;
	}

	/**
	 * Convert to an image stack.
	 *
	 * @return the image stack
	 */
	public ImageStack getImageStack()
	{
		ImageStack stack = new ImageStack(nc, nr);
		int nr_nc = nr * nc;
		for (int s = 0; s < ns; s++)
		{
			float[] pixels = new float[nr_nc];
			System.arraycopy(data, s * nr_nc, pixels, 0, nr_nc);
			stack.addSlice(null, pixels);
		}
		return stack;
	}

	/**
	 * Performs a forward transform, converting this image into the frequency domain.
	 * 
	 * @throws IllegalArgumentException
	 *             If already in the frequency domain
	 */
	public void transform() throws IllegalArgumentException
	{
		if (isFrequencyDomain)
			throw new IllegalArgumentException("Already frequency domain DHT");
		dht.forward(data);
		isFrequencyDomain = true;
	}

	/**
	 * Performs an inverse transform, converting this image into the space domain.
	 *
	 * @throws IllegalArgumentException
	 *             If already in the space domain
	 */
	public void inverseTransform() throws IllegalArgumentException
	{
		if (!isFrequencyDomain)
			throw new IllegalArgumentException("Already space domain DHT");
		dht.inverse(data, true);
		isFrequencyDomain = false;
	}

	/**
	 * Checks if is frequency domain.
	 *
	 * @return true, if is frequency domain
	 */
	public boolean isFrequencyDomain()
	{
		return isFrequencyDomain;
	}

	/**
	 * Returns the image resulting from the point by point Hartley multiplication
	 * of this image and the specified image. Both images are assumed to be in
	 * the frequency domain. Multiplication in the frequency domain is equivalent
	 * to convolution in the space domain.
	 *
	 * @param dht
	 *            the dht
	 * @return the result
	 * @throws IllegalArgumentException
	 *             if the dht is not the same dimensions
	 */
	public DHT3D multiply(DHT3D dht) throws IllegalArgumentException
	{
		checkDHT(dht);

		float[] tmp = new float[data.length];
		float[] h1 = this.data;
		float[] h2 = dht.data;
		int nr_nc = nr * nc;

		for (int s = 0, ns_s = 0, i = 0; s < ns; s++, ns_s = ns - s)
		{
			for (int r = 0, nr_r = 0; r < nr; r++, nr_r = nr - r)
			{
				for (int c = 0, nc_c = 0; c < nc; c++, nc_c = nc - c, i++)
				{
					// This is actually doing for 3D data stored as x[slices][rows][columns]
					// https://en.wikipedia.org/wiki/Discrete_Hartley_transform
					//h2e = (h2[s][r][c] + h2[Ns-s][Nr-r][Nr-c]) / 2;
					//h2o = (h2[s][r][c] - h2[Ns-s][Nr-r][Nr-c]) / 2;
					//tmp[s][r][c] = (float) (h1[s][r][c] * h2e + h1[Ns-s][Nr-r][Nc-c] * h2o);
					int j = ns_s * nr_nc + nr_r * nc + nc_c;
					double h2e = (h2[i] + h2[j]) / 2;
					double h2o = (h2[i] - h2[j]) / 2;
					tmp[i] = (float) (h1[i] * h2e + h1[j] * h2o);
				}
			}
		}

		return new DHT3D(tmp, ns, nr, nc, true, this.dht);
	}

	/**
	 * Returns the image resulting from the point by point Hartley conjugate
	 * multiplication of this image and the specified image. Both images are
	 * assumed to be in the frequency domain. Conjugate multiplication in
	 * the frequency domain is equivalent to correlation in the space domain.
	 *
	 * @param dht
	 *            the dht
	 * @return the result
	 * @throws IllegalArgumentException
	 *             if the dht is not the same dimensions
	 */
	public DHT3D conjugateMultiply(DHT3D dht) throws IllegalArgumentException
	{
		checkDHT(dht);

		float[] tmp = new float[data.length];
		float[] h1 = this.data;
		float[] h2 = dht.data;
		int nr_nc = nr * nc;

		for (int s = 0, ns_s = 0, i = 0; s < ns; s++, ns_s = ns - s)
		{
			for (int r = 0, nr_r = 0; r < nr; r++, nr_r = nr - r)
			{
				for (int c = 0, nc_c = 0; c < nc; c++, nc_c = nc - c, i++)
				{
					int j = ns_s * nr_nc + nr_r * nc + nc_c;
					double h2e = (h2[i] + h2[j]) / 2;
					double h2o = (h2[i] - h2[j]) / 2;
					// As per multiply but reverse the addition sign for the conjugate  
					tmp[i] = (float) (h1[i] * h2e - h1[j] * h2o);
				}
			}
		}

		return new DHT3D(tmp, ns, nr, nc, true, this.dht);
	}

	/**
	 * Returns the image resulting from the point by point Hartley division
	 * of this image by the specified image. Both images are assumed to be in
	 * the frequency domain. Division in the frequency domain is equivalent
	 * to deconvolution in the space domain.
	 *
	 * @param dht
	 *            the dht
	 * @return the result
	 * @throws IllegalArgumentException
	 *             if the dht is not the same dimensions or in the frequency domain
	 */
	public DHT3D divide(DHT3D dht) throws IllegalArgumentException
	{
		checkDHT(dht);

		float[] tmp = new float[data.length];
		float[] h1 = this.data;
		float[] h2 = dht.data;
		int nr_nc = nr * nc;

		for (int s = 0, ns_s = 0, i = 0; s < ns; s++, ns_s = ns - s)
		{
			for (int r = 0, nr_r = 0; r < nr; r++, nr_r = nr - r)
			{
				for (int c = 0, nc_c = 0; c < nc; c++, nc_c = nc - c, i++)
				{
					// This is a copy of the divide operation in ij.process.FHT
					int j = ns_s * nr_nc + nr_r * nc + nc_c;
					double mag = h2[i] * h2[i] + h2[j] * h2[j];
					if (mag < 1e-20)
						mag = 1e-20;
					double h2e = (h2[i] + h2[j]);
					double h2o = (h2[i] - h2[j]);
					tmp[i] = (float) ((h1[i] * h2e - h1[j] * h2o) / mag);
				}
			}
		}

		return new DHT3D(tmp, ns, nr, nc, true, this.dht);
	}

	private void checkDHT(DHT3D dht)
	{
		if (dht.ns != ns || dht.nr != nr || dht.nc != nc)
			throw new IllegalArgumentException("Dimension mismatch");
		if (!dht.isFrequencyDomain || !isFrequencyDomain)
			throw new IllegalArgumentException("Require frequency domain DHT");
	}

	/**
	 * Swap quadrants 1+7, 2+8, 4+6, 3+5 of the specified image stack
	 * so the power spectrum origin is at the center of the image.
	 * 
	 * <pre>
	 *        2----1
	 *       /    /
	 *      3----4
	 * 
	 *        6----5
	 *       /    / 
	 *      7----8
	 * </pre>
	 * 
	 * Requires even dimensions.
	 *
	 * @param stack
	 *            the stack
	 */
	public static void swapQuadrants(ImageStack stack)
	{
		if (stack.getBitDepth() != 32)
			throw new IllegalArgumentException("Require float stack");
		int ns = stack.getSize();
		int nr = stack.getHeight();
		int nc = stack.getWidth();
		if ((ns & 1) == 1 || (nr & 1) == 1 || (nc & 1) == 1)
			throw new IllegalArgumentException("Require even dimensions");

		// TODO - Finish this

		//		int width = ip.getWidth();
		//		float[] pixels = (float[]) ip.getPixels();
		//		int size = width / 2;
		//		float[] a = new float[size * size];
		//		float[] b = new float[size * size];

		//		crop(pixels, width, a, size, 0, size);
		//		crop(pixels, width, b, 0, size, size);
		//		insert(pixels, width, b, size, 0, size);
		//		insert(pixels, width, a, 0, size, size);
		//		crop(pixels, width, a, 0, 0, size);
		//		crop(pixels, width, b, size, size, size);
		//		insert(pixels, width, b, 0, 0, size);
		//		insert(pixels, width, a, size, size, size);
	}

	/**
	 * Crop from the source.
	 *
	 * @param source
	 *            the source pixels
	 * @param width
	 *            the width of the source pixels
	 * @param x
	 *            the source x location
	 * @param y
	 *            the source y location
	 * @param w
	 *            the source width
	 * @param h
	 *            the source height
	 * @param buffer
	 *            the buffer pixels
	 */
	private static void crop(float[] source, int width, int x, int y, int w, int h, float[] buffer)
	{
		for (int ys = y + h; ys-- > y;)
		{
			int si = ys * width + x;
			int bi = (ys - y) * w;
			System.arraycopy(source, si, buffer, bi, w);
		}
	}

	/**
	 * Insert into the source.
	 *
	 * @param source
	 *            the source pixels
	 * @param width
	 *            the width of the source pixels
	 * @param x
	 *            the source x location
	 * @param y
	 *            the source y location
	 * @param w
	 *            the source width
	 * @param h
	 *            the source height
	 * @param buffer
	 *            the buffer pixels
	 */
	private static void insert(float[] source, int width, int x, int y, int w, int h, float[] buffer)
	{
		for (int ys = y + h; ys-- > y;)
		{
			int si = ys * width + x;
			int bi = (ys - y) * w;
			System.arraycopy(buffer, bi, source, si, w);
		}
	}
}