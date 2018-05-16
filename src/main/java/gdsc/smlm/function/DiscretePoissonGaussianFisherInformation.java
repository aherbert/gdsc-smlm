package gdsc.smlm.function;

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
 * Calculate the Fisher information for a discrete Poisson-Gaussian distribution.
 * <p>
 * Uses the equation of Chao, et al (2013) Nature Methods, 10, 335-338, SI Eq S7. This has been modified so the
 * convolution with the Gaussian uses the error function.
 * <p>
 * Performs a convolution with a finite Gaussian kernel. The Gaussian is constructed using a range of the standard
 * deviation (s) and sampled using the error function in the interval x-0.5 to x+0.5 to generate a discrete Gaussian
 * PMF. This computed Fisher information is thus for a discrete Poisson-Gaussian PMF.
 * <p>
 * An optimisation is used when the mean of the Poisson is above a threshold. In this case the Poisson can be
 * approximated as a Gaussian and the Fisher information is returned for the Gaussian-Gaussian convolution.
 */
public class DiscretePoissonGaussianFisherInformation extends PoissonGaussianFisherInformation
{
	/**
	 * The Gaussian convolution kernel.
	 */
	private final double[] kernel;

	/**
	 * Instantiates a new real poisson gaussian fisher information.
	 *
	 * @param s
	 *            the standard deviation of the Gaussian
	 * @throws IllegalArgumentException
	 *             If the standard deviation is not strictly positive
	 */
	public DiscretePoissonGaussianFisherInformation(double s) throws IllegalArgumentException
	{
		// With no upscaling the convolution can handle a large range.
		this(s, 8);
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
	public DiscretePoissonGaussianFisherInformation(double s, double range) throws IllegalArgumentException
	{
		super(s, range);
		kernel = createKernel(s, range);
	}

	private static double[] createKernel(double s, double range)
	{
		// Get the range of the kernel
		if (range * s < 0.5)
		{
			// Convolution will not cover more than 1 sample of the Poisson PMF
			return null;
		}

		return Convolution.makeErfGaussianKernel(s, range);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.PoissonGaussianFisherInformation#getGaussianKernel(double)
	 */
	@Override
	protected int getKernelScale(double t)
	{
		return (kernel == null) ? 0 : 1;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.PoissonGaussianFisherInformation#getGaussianKernel(int)
	 */
	@Override
	protected double[] getGaussianKernel(int scale)
	{
		return kernel;
	}
}