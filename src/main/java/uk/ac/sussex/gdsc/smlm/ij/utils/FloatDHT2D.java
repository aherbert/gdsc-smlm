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
package uk.ac.sussex.gdsc.smlm.ij.utils;

import org.jtransforms.dht.FloatDHT_2D;

import ij.process.ImageProcessor;
import pl.edu.icm.jlargearrays.LargeArray;

/**
 * Wrapper to compute the discrete Hartley transform on 2D data. This uses the JTransforms library.
 */
public class FloatDHT2D extends FloatImage2D
{
	private boolean isFrequencyDomain;
	private final FloatDHT_2D dht;
	// Used for fast multiply operations
	private double[] h2e, h2o, mag;
	private int[] jj;

	/**
	 * Instantiates a new 2D discrete Hartley transform
	 *
	 * @param image
	 *            the stack
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if the combined dimensions is too large for an array
	 */
	public FloatDHT2D(ImageProcessor image) throws IllegalArgumentException
	{
		super(image);
		LargeArray.setMaxSizeOf32bitArray(MAX_SIZE_OF_32_BIT_ARRAY);
		dht = new FloatDHT_2D(nr, nc);
	}

	/**
	 * Instantiates a new 2D discrete Hartley transform.
	 *
	 * @param nc
	 *            the number of columns
	 * @param nr
	 *            the number of rows
	 * @param data
	 *            the data
	 * @param isFrequencyDomain
	 *            Set to true if in the frequency domain
	 * @throws IllegalArgumentException
	 *             If any dimension is less than 2, or if the data is not the correct length
	 */
	public FloatDHT2D(int nc, int nr, float[] data, boolean isFrequencyDomain) throws IllegalArgumentException
	{
		super(nc, nr, data);
		LargeArray.setMaxSizeOf32bitArray(MAX_SIZE_OF_32_BIT_ARRAY);
		dht = new FloatDHT_2D(nr, nc);
		this.isFrequencyDomain = isFrequencyDomain;
	}

	/**
	 * Instantiates a new 2D discrete Hartley transform.
	 *
	 * @param nc
	 *            the number of columns
	 * @param nr
	 *            the number of rows
	 * @param data
	 *            the data
	 * @param isFrequencyDomain
	 *            the is frequency domain
	 * @param dht
	 *            the dht
	 */
	private FloatDHT2D(int nc, int nr, float[] data, boolean isFrequencyDomain, FloatDHT_2D dht)
	{
		super(nc, nr, data, false);
		this.isFrequencyDomain = isFrequencyDomain;
		this.dht = dht; // This can be reused across objects
	}

	/**
	 * Return a copy of the 2D discrete Hartley transform.
	 *
	 * @return the copy
	 */
	@Override
	public FloatDHT2D copy()
	{
		final FloatDHT2D copy = new FloatDHT2D(nc, nr, data.clone(), isFrequencyDomain, dht);
		copy.h2e = h2e;
		copy.h2o = h2o;
		copy.jj = jj;
		copy.mag = mag;
		return copy;
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
		resetFastOperations();
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
		resetFastOperations();
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
	 * Initialise fast operations for {@link #multiply(FloatDHT2D)} and {@link #conjugateMultiply(FloatDHT2D)}. This
	 * pre-computes
	 * the values needed for the operations.
	 * <p>
	 * Note: This initialises the FHT object for use as the argument to the operation, for example if a convolution
	 * kernel is to be applied to many FHT objects.
	 */
	public void initialiseFastMultiply()
	{
		if (h2e == null)
		{
			// Do this on new arrays for thread safety (i.e. concurrent initialisation)
			final float[] h2 = getData();
			final double[] h2e = new double[h2.length];
			final double[] h2o = new double[h2e.length];
			final int[] jj = new int[h2e.length];
			for (int r = 0, nr_m_r = 0, i = 0; r < nr; r++, nr_m_r = nr - r)
				for (int c = 0, nc_m_c = 0; c < nc; c++, nc_m_c = nc - c, i++)
				{
					final int j = nr_m_r * nc + nc_m_c;
					h2e[i] = ((double) h2[i] + (double) h2[j]) / 2.0;
					h2o[i] = ((double) h2[i] - (double) h2[j]) / 2.0;
					jj[i] = j;
				}
			this.h2o = h2o;
			this.jj = jj;
			// Assign at the end for thread safety (i.e. concurrent initialisation)
			this.h2e = h2e;
		}
	}

	/**
	 * Initialise fast operations for {@link #multiply(FloatDHT2D)}, {@link #conjugateMultiply(FloatDHT2D)} and
	 * {@link #divide(FloatDHT2D)}. This pre-computes the values needed for the operations.
	 * <p>
	 * Note: This initialises the FHT object for use as the argument to the operation, for example if a deconvolution
	 * kernel is to be applied to many FHT objects.
	 */
	public void initialiseFastOperations()
	{
		initialiseFastMultiply();
		if (mag == null)
		{
			// Do this on new arrays for thread safety (i.e. concurrent initialisation)
			final double[] mag = new double[h2e.length];
			final float[] h2 = getData();
			for (int i = 0; i < h2.length; i++)
				// Note that pre-computed h2e and h2o are divided by 2 so we also
				// divide the magnitude by 2 to allow reuse of the pre-computed values
				// in the divide operation (which does not require h2e/2 and h2o/2)
				mag[i] = Math.max(1e-20, h2[i] * h2[i] + h2[jj[i]] * h2[jj[i]]) / 2;
			this.mag = mag;
		}
	}

	/**
	 * Checks if is initialised for fast multiply.
	 *
	 * @return true, if is fast multiply
	 */
	public boolean isFastMultiply()
	{
		return h2e != null;
	}

	/**
	 * Checks if is initialised for fast operations.
	 *
	 * @return true, if is fast operations
	 */
	public boolean isFastOperations()
	{
		return mag != null;
	}

	private void resetFastOperations()
	{
		h2e = null;
		h2o = null;
		jj = null;
		mag = null;
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
	public FloatDHT2D multiply(FloatDHT2D dht) throws IllegalArgumentException
	{
		return multiply(dht, null);
	}

	/**
	 * Returns the image resulting from the point by point Hartley multiplication
	 * of this image and the specified image. Both images are assumed to be in
	 * the frequency domain. Multiplication in the frequency domain is equivalent
	 * to convolution in the space domain.
	 *
	 * @param dht
	 *            the dht
	 * @param tmp
	 *            the tmp buffer to use for the result
	 * @return the result
	 * @throws IllegalArgumentException
	 *             if the dht is not the same dimensions
	 */
	public FloatDHT2D multiply(FloatDHT2D dht, float[] tmp) throws IllegalArgumentException
	{
		checkDHT(dht);
		return (dht.isFastMultiply()) ? multiply(dht.h2e, dht.h2o, dht.jj, tmp) : multiply(dht.getData(), tmp);
	}

	/**
	 * Returns the image resulting from the point by point Hartley multiplication
	 * of this image and the specified image. Both images are assumed to be in
	 * the frequency domain. Multiplication in the frequency domain is equivalent
	 * to convolution in the space domain.
	 *
	 * @param h2
	 *            the second FHT
	 * @param tmp
	 *            the tmp buffer to use for the result
	 * @return the result
	 * @throws IllegalArgumentException
	 *             if the dht is not the same dimensions
	 */
	private FloatDHT2D multiply(float[] h2, float[] tmp) throws IllegalArgumentException
	{
		final float[] h1 = this.data;
		if (tmp == null || tmp.length != h1.length)
			tmp = new float[h1.length];

		for (int r = 0, nr_m_r = 0, i = 0; r < nr; r++, nr_m_r = nr - r)
			for (int c = 0, nc_m_c = 0; c < nc; c++, nc_m_c = nc - c, i++)
			{
				// This is actually doing for 2D data stored as x[rows][columns]
				// https://en.wikipedia.org/wiki/Discrete_Hartley_transform
				//h2e = (h2[r][c] + h2[Nr-r][Nr-c]) / 2;
				//h2o = (h2[r][c] - h2[Nr-r][Nr-c]) / 2;
				//tmp[r][c] = (float) (h1[r][c] * h2e + h1[Nr-r][Nc-c] * h2o);
				final int j = nr_m_r * nc + nc_m_c;
				final double h2e = ((double) h2[i] + (double) h2[j]) / 2.0;
				final double h2o = ((double) h2[i] - (double) h2[j]) / 2.0;
				tmp[i] = (float) (h1[i] * h2e + h1[j] * h2o);
			}

		return new FloatDHT2D(nc, nr, tmp, true, this.dht);
	}

	/**
	 * Returns the image resulting from the point by point Hartley multiplication
	 * of this image and the specified image. Both images are assumed to be in
	 * the frequency domain. Multiplication in the frequency domain is equivalent
	 * to convolution in the space domain.
	 *
	 * @param h2e
	 *            the pre-initialised h2e value
	 * @param h2o
	 *            the pre-initialised h2o value
	 * @param jj
	 *            the pre-initialised j index
	 * @param tmp
	 *            the buffer for the result (can be null)
	 * @return the result
	 */
	private FloatDHT2D multiply(double[] h2e, double[] h2o, int[] jj, float[] tmp)
	{
		final float[] h1 = getData();
		if (tmp == null || tmp.length != h1.length)
			tmp = new float[h1.length];
		for (int i = 0; i < h1.length; i++)
			tmp[i] = (float) (h1[i] * h2e[i] + h1[jj[i]] * h2o[i]);
		return new FloatDHT2D(nc, nr, tmp, true, this.dht);
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
	public FloatDHT2D conjugateMultiply(FloatDHT2D dht) throws IllegalArgumentException
	{
		return conjugateMultiply(dht, null);
	}

	/**
	 * Returns the image resulting from the point by point Hartley conjugate
	 * multiplication of this image and the specified image. Both images are
	 * assumed to be in the frequency domain. Conjugate multiplication in
	 * the frequency domain is equivalent to correlation in the space domain.
	 *
	 * @param dht
	 *            the dht
	 * @param tmp
	 *            the tmp buffer to use for the result
	 * @return the result
	 * @throws IllegalArgumentException
	 *             if the dht is not the same dimensions
	 */
	public FloatDHT2D conjugateMultiply(FloatDHT2D dht, float[] tmp) throws IllegalArgumentException
	{
		checkDHT(dht);
		return (dht.isFastMultiply()) ? conjugateMultiply(dht.h2e, dht.h2o, dht.jj, tmp)
				: conjugateMultiply(dht.getData(), tmp);
	}

	/**
	 * Returns the image resulting from the point by point Hartley conjugate
	 * multiplication of this image and the specified image. Both images are
	 * assumed to be in the frequency domain. Conjugate multiplication in
	 * the frequency domain is equivalent to correlation in the space domain.
	 *
	 * @param h2
	 *            the second DHT
	 * @param tmp
	 *            the tmp buffer to use for the result
	 * @return the result
	 * @throws IllegalArgumentException
	 *             if the dht is not the same dimensions
	 */
	private FloatDHT2D conjugateMultiply(float[] h2, float[] tmp) throws IllegalArgumentException
	{
		final float[] h1 = this.data;
		if (tmp == null || tmp.length != h1.length)
			tmp = new float[h1.length];

		for (int r = 0, nr_m_r = 0, i = 0; r < nr; r++, nr_m_r = nr - r)
			for (int c = 0, nc_m_c = 0; c < nc; c++, nc_m_c = nc - c, i++)
			{
				final int j = nr_m_r * nc + nc_m_c;
				final double h2e = ((double) h2[i] + (double) h2[j]) / 2.0;
				final double h2o = ((double) h2[i] - (double) h2[j]) / 2.0;
				// As per multiply but reverse the addition sign for the conjugate
				tmp[i] = (float) (h1[i] * h2e - h1[j] * h2o);
			}

		return new FloatDHT2D(nc, nr, tmp, true, this.dht);
	}

	/**
	 * Returns the image resulting from the point by point Hartley conjugate
	 * multiplication of this image and the specified image. Both images are
	 * assumed to be in the frequency domain. Conjugate multiplication in
	 * the frequency domain is equivalent to correlation in the space domain.
	 *
	 * @param h2e
	 *            the pre-initialised h2e value
	 * @param h2o
	 *            the pre-initialised h2o value
	 * @param jj
	 *            the pre-initialised j index
	 * @param tmp
	 *            the buffer for the result (can be null)
	 * @return the fht2
	 */
	private FloatDHT2D conjugateMultiply(double[] h2e, double[] h2o, int[] jj, float[] tmp)
	{
		final float[] h1 = getData();
		if (tmp == null || tmp.length != h1.length)
			tmp = new float[h1.length];
		for (int i = 0; i < h1.length; i++)
			tmp[i] = (float) (h1[i] * h2e[i] - h1[jj[i]] * h2o[i]);
		return new FloatDHT2D(nc, nr, tmp, true, this.dht);
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
	public FloatDHT2D divide(FloatDHT2D dht) throws IllegalArgumentException
	{
		return divide(dht, null);
	}

	/**
	 * Returns the image resulting from the point by point Hartley division
	 * of this image by the specified image. Both images are assumed to be in
	 * the frequency domain. Division in the frequency domain is equivalent
	 * to deconvolution in the space domain.
	 *
	 * @param dht
	 *            the dht
	 * @param tmp
	 *            the tmp buffer to use for the result
	 * @return the result
	 * @throws IllegalArgumentException
	 *             if the dht is not the same dimensions or in the frequency domain
	 */
	public FloatDHT2D divide(FloatDHT2D dht, float[] tmp) throws IllegalArgumentException
	{
		checkDHT(dht);
		return (dht.isFastOperations()) ? divide(dht.h2e, dht.h2o, dht.jj, dht.mag, tmp) : divide(dht.getData(), tmp);
	}

	/**
	 * Returns the image resulting from the point by point Hartley division
	 * of this image by the specified image. Both images are assumed to be in
	 * the frequency domain. Division in the frequency domain is equivalent
	 * to deconvolution in the space domain.
	 *
	 * @param h2
	 *            the second DHT
	 * @param tmp
	 *            the tmp buffer to use for the result
	 * @return the result
	 * @throws IllegalArgumentException
	 *             if the dht is not the same dimensions or in the frequency domain
	 */
	private FloatDHT2D divide(float[] h2, float[] tmp) throws IllegalArgumentException
	{
		final float[] h1 = this.data;
		if (tmp == null || tmp.length != h1.length)
			tmp = new float[h1.length];

		for (int r = 0, nr_m_r = 0, i = 0; r < nr; r++, nr_m_r = nr - r)
			for (int c = 0, nc_m_c = 0; c < nc; c++, nc_m_c = nc - c, i++)
			{
				// This is a copy of the divide operation in ij.process.FHT
				final int j = nr_m_r * nc + nc_m_c;
				final double h2i = h2[i];
				final double h2j = h2[j];
				double mag = h2i * h2i + h2j * h2j;
				if (mag < 1e-20)
					mag = 1e-20;
				final double h2e = (h2i + h2j);
				final double h2o = (h2i - h2j);
				tmp[i] = (float) ((h1[i] * h2e - h1[j] * h2o) / mag);
			}

		return new FloatDHT2D(nc, nr, tmp, true, this.dht);
	}

	/**
	 * Returns the image resulting from the point by point Hartley division
	 * of this image by the specified image. Both images are assumed to be in
	 * the frequency domain. Division in the frequency domain is equivalent
	 * to deconvolution in the space domain.
	 *
	 * @param h2e
	 *            the pre-initialised h2e value
	 * @param h2o
	 *            the pre-initialised h2o value
	 * @param jj
	 *            the pre-initialised j index
	 * @param mag
	 *            the pre-initialised magnitude value
	 * @param tmp
	 *            the buffer for the result (can be null)
	 */
	private FloatDHT2D divide(double[] h2e, double[] h2o, int[] jj, double[] mag, float[] tmp)
	{
		final float[] h1 = getData();
		if (tmp == null || tmp.length != h1.length)
			tmp = new float[h1.length];
		for (int i = 0; i < h1.length; i++)
			tmp[i] = (float) ((h1[i] * h2e[i] - h1[jj[i]] * h2o[i]) / mag[i]);
		return new FloatDHT2D(nc, nr, tmp, true, this.dht);
	}

	/**
	 * Check the DHT matches the dimensions of this DHT. Check both are in the frequency domain.
	 *
	 * @param dht
	 *            the dht
	 * @throws IllegalArgumentException
	 *             If multiplication is not possible
	 */
	private void checkDHT(FloatDHT2D dht) throws IllegalArgumentException
	{
		if (dht.nr != nr || dht.nc != nc)
			throw new IllegalArgumentException("Dimension mismatch");
		if (!dht.isFrequencyDomain || !isFrequencyDomain)
			throw new IllegalArgumentException("Require frequency domain DHT");
	}

	/**
	 * Converts this DHT to a discrete Fourier transform (DFT) and returns it as a two image pair. The image is assumed
	 * to be in the frequency domain.
	 *
	 * @param real
	 *            the buffer to use for the real component (can be null)
	 * @param imaginary
	 *            the buffer to use for the imaginary component (can be null)
	 * @return [real, imaginary]
	 * @throws IllegalArgumentException
	 *             if not in the frequency domain
	 * @see <A href=
	 *      "https://en.wikipedia.org/wiki/Hartley_transform#Relation_to_Fourier_transform">https://en.wikipedia.org/
	 *      wiki/Hartley_transform#Relation_to_Fourier_transform</a>
	 */
	public FloatImage2D[] toDFT(float[] real, float[] imaginary) throws IllegalArgumentException
	{
		if (!isFrequencyDomain)
			throw new IllegalArgumentException("Require frequency domain DHT");

		final float[] h1 = this.data;
		if (real == null || real.length != h1.length)
			real = new float[h1.length];
		if (imaginary == null || imaginary.length != h1.length)
			imaginary = new float[h1.length];

		for (int r = 0, nr_m_r = 0, i = 0; r < nr; r++, nr_m_r = nr - r)
			for (int c = 0, nc_m_c = 0; c < nc; c++, nc_m_c = nc - c, i++)
			{
				// This is a copy of the getComplexTransform operation in ij.process.FHT
				final int j = nr_m_r * nc + nc_m_c;
				real[i] = (h1[i] + h1[j]) * 0.5f;
				imaginary[i] = (-h1[i] + h1[j]) * 0.5f;
			}

		return new FloatImage2D[] { new FloatImage2D(nc, nr, real, false), new FloatImage2D(nc, nr, imaginary, false) };
	}

	/**
	 * Convert a discrete Fourier transform (DFT) to a DHT
	 *
	 * @param real
	 *            the real component
	 * @param imaginary
	 *            the imaginary component
	 * @param tmp
	 *            the tmp buffer to use for the result
	 * @return the DHT
	 * @throws IllegalArgumentException
	 *             If there is a dimension mismatch
	 */
	public static FloatDHT2D fromDFT(FloatImage2D real, FloatImage2D imaginary, float[] tmp)
			throws IllegalArgumentException
	{
		if (real.nr != imaginary.nr || real.nc != imaginary.nc)
			throw new IllegalArgumentException("Dimension mismatch");

		final float[] re = real.getData();
		final float[] im = imaginary.getData();
		if (tmp == null || tmp.length != re.length)
			tmp = new float[re.length];

		final int nc = real.nc;
		final int nr = real.nr;

		for (int r = 0, nr_m_r = 0, i = 0; r < nr; r++, nr_m_r = nr - r)
			for (int c = 0, nc_m_c = 0; c < nc; c++, nc_m_c = nc - c, i++)
			{
				final int j = nr_m_r * nc + nc_m_c;
				// Reverse the toDFT() method
				// re = (a+b)/2
				// im = (-a+b)/2
				// b = re + im
				// a = 2*re - b
				tmp[j] = re[i] + im[i];
				tmp[i] = 2 * re[i] - tmp[j];
			}

		return new FloatDHT2D(nc, nr, tmp, true, new FloatDHT_2D(nr, nc));
	}

	/**
	 * Returns the absolute value (amplitude) of the Hartley transform. The image is assumed to be in
	 * the frequency domain.
	 *
	 * @param tmp
	 *            the tmp buffer to use for the result
	 * @return the result
	 * @throws IllegalArgumentException
	 *             if not in the frequency domain
	 */
	public FloatImage2D getAbsoluteValue(float[] tmp) throws IllegalArgumentException
	{
		if (!isFrequencyDomain)
			throw new IllegalArgumentException("Require frequency domain DHT");

		final float[] h1 = this.data;
		if (tmp == null || tmp.length != h1.length)
			tmp = new float[h1.length];

		for (int r = 0, nr_m_r = 0, i = 0; r < nr; r++, nr_m_r = nr - r)
			for (int c = 0, nc_m_c = 0; c < nc; c++, nc_m_c = nc - c, i++)
			{
				// This is a copy of the amplitude operation in ij.process.FHT
				final int j = nr_m_r * nc + nc_m_c;
				tmp[i] = (float) Math.sqrt(h1[i] * h1[i] + h1[j] * h1[j]);
			}

		return new FloatImage2D(nc, nr, tmp, false);
	}

	/**
	 * Swap quadrants 1 and 3 and 2 and 4 of image
	 * so the power spectrum origin is at the center of the image.
	 *
	 * <pre>
	    2 1
	    3 4
	 * </pre>
	 */
	public void swapQuadrants()
	{
		swapQuadrants(this);
		resetFastOperations();
	}

	/**
	 * Swap quadrants 1 and 3 and 2 and 4 of the specified ImageProcessor
	 * so the power spectrum origin is at the center of the image.
	 *
	 * <pre>
	    2 1
	    3 4
	 * </pre>
	 *
	 * @param image
	 *            The image (must be even dimensions)
	 * @throws IllegalArgumentException
	 *             If not even dimensions
	 */
	public static void swapQuadrants(FloatImage2D image) throws IllegalArgumentException
	{
		// This is a specialised version to allow using a float buffer and
		// optimised for even sized images

		final int ny = image.getHeight();
		final int nx = image.getWidth();
		if ((ny & 1) == 1 || (nx & 1) == 1)
			throw new IllegalArgumentException("Require even dimensions");

		final int ny_2 = ny / 2;
		final int nx_2 = nx / 2;

		final float[] tmp = new float[nx];
		final float[] a = image.getData();

		//@formatter:off
		// We swap: 0 <=> nx_2, 0 <=> ny_2
		// 1 <=> 3
		swap(a, a, nx, nx_2,    0,    0, ny_2, nx_2, ny_2, tmp);
		// 2 <=> 4
		swap(a, a, nx,    0,    0, nx_2, ny_2, nx_2, ny_2, tmp);
		//@formatter:on
	}

	/**
	 * Swap the rectangle pixel values from a with b.
	 * <p>
	 * No bounds checks are performed so use with care!
	 *
	 * @param a
	 *            the a pixels
	 * @param b
	 *            the b pixels (must match a.length)
	 * @param width
	 *            the width of each set of pixels
	 * @param ax
	 *            the x origin from a
	 * @param ay
	 *            the y origin from a
	 * @param bx
	 *            the x origin from b
	 * @param by
	 *            the b origin from b
	 * @param w
	 *            the width of the rectangle to swap
	 * @param h
	 *            the height of the rectangle to swap
	 * @param tmp
	 *            the tmp buffer (must be at least width in length)
	 */
	public static void swap(float[] a, float[] b, int width, int ax, int ay, int bx, int by, int w, int h, float[] tmp)
	{
		for (int ayy = ay + h, byy = by + h - 1; ayy-- > ay; byy--)
		{
			final int ai = ayy * width + ax;
			final int bi = byy * width + bx;
			System.arraycopy(a, ai, tmp, 0, w);
			System.arraycopy(b, bi, a, ai, w);
			System.arraycopy(tmp, 0, b, bi, w);
		}
	}
}