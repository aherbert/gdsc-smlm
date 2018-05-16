package gdsc.smlm.function;

import gdsc.core.utils.Maths;
import gdsc.smlm.utils.Convolution;

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
 * Calculate the Fisher information for a continuous Poisson-Gaussian distribution.
 * <p>
 * Uses the equation of Chao, et al (2013) Nature Methods, 10, 335-338, SI Eq S7.
 * <p>
 * Performs a convolution with a finite Gaussian kernel. The Gaussian is constructed using a range of the standard
 * deviation (s) and sampled at least every s/2.
 * <p>
 * An optimisation is used to avoid computation on tiny Gaussian kernels (i.e. too small to be computed)
 * This will occur when the Gaussian standard deviation is less than 0.02. The result is no convolution and the result
 * computes Poisson Fisher information.
 * <p>
 * An optimisation is used when the mean of the Poisson is above a threshold. In this case the Poisson can be
 * approximated as a Gaussian and the Fisher information is returned for the Gaussian-Gaussian convolution.
 */
public class RealPoissonGaussianFisherInformation extends PoissonGaussianFisherInformation
{
	/** The default scale for the kernel. */
	private final int defaultScale;

	/**
	 * The Gaussian convolution kernels for different scaling. The scale is 2^index, e.g. 1, 2, 4, 8, 16, 32, 64, 128.
	 */
	private final double[][] kernel;

	/**
	 * Instantiates a new real poisson gaussian fisher information.
	 *
	 * @param s
	 *            the standard deviation of the Gaussian
	 * @throws IllegalArgumentException
	 *             If the standard deviation is not strictly positive
	 */
	public RealPoissonGaussianFisherInformation(double s) throws IllegalArgumentException
	{
		this(s, 7);
	}

	/**
	 * Instantiates a new real poisson gaussian fisher information.
	 *
	 * @param s
	 *            the standard deviation of the Gaussian
	 * @param range
	 *            the range of the Gaussian kernel (in SD units). This is clipped to the range
	 *            1-38 to provide a meaningful convolution.
	 * @throws IllegalArgumentException
	 *             If the standard deviation is not strictly positive
	 */
	public RealPoissonGaussianFisherInformation(double s, double range) throws IllegalArgumentException
	{
		super(s, range);

		// Check if the Gaussian standard deviation is above the threshold for computation.
		// Also check if the gaussian filter will touch more than one Poisson value.
		// Otherwise convolution is not possible.
		// The limit s==0.04 is based on the scale being 4/s = 100. Do not support scaling 
		// greater than this. It is unlikely anyway.
		if (s >= 0.04 && s * range >= 1)
		{
			// Determine how much to up-sample so that the convolution with the Gaussian
			// uses multiple values of the Gaussian.
			defaultScale = getScale(s);

			// Store the Gaussian kernels for convolution:
			// 1, 2, 4, 8, 16, 32, 64, 128
			kernel = new double[8][];
		}
		else
		{
			defaultScale = 0;
			kernel = null;
		}
	}

	private static int getScale(double s)
	{
		double scale = Math.ceil(4 / s);
		if (scale > 128)
			return 128;
		return Maths.nextPow2((int) scale);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.PoissonGaussianFisherInformation#getGaussianKernel(double)
	 */
	@Override
	protected int getKernelScale(double t)
	{
		if (defaultScale == 0)
			return 0;
		// Choose the kernel. A small mean requires more Gaussian samples.
		// Note the default scale is the minimum required to sample at 0.5 SD units.
		// Find the same for the Poisson using its variance.
		// mean 4 => scale = 1
		// mean <4 => scale = 2
		// mean <1 => scale >= 4
		// This may have to be changed.
		return Math.max(defaultScale, getScale(t / 2));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.PoissonGaussianFisherInformation#getGaussianKernel(int)
	 */
	@Override
	protected double[] getGaussianKernel(int scale)
	{
		// Get the Gaussian kernel
		int index = log2(scale);
		if (kernel[index] == null)
		{
			//kernel[index] = Convolution.makeGaussianKernel(s * scale, range, true);
			kernel[index] = Convolution.makeGaussianKernel(s * scale, range, false);
			// This does not work. The Poisson is discrete and so must be convolved
			// with a continuous Gaussian. Using the Error function is equivalent
			// to integrating the continuous Gaussian with the same value of the 
			// Poisson over the range x-0.5 to x+0.5. However the Poisson is zero
			// in all that range except for x.
			//kernel[index] = Convolution.makeErfGaussianKernel(s * scale, range);
		}
		return kernel[index];
	}

	private int log2(int scale)
	{
		int bits = 30;
		while ((scale & (1 << bits)) == 0)
			bits--;
		return bits;
	}
}