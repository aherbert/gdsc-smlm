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
	/**
	 * The Gaussian convolution kernels for different ranges.
	 */
	private final double[][] kernel;

	/**
	 * Instantiates a new real poisson gaussian fisher information.
	 *
	 * @param s
	 *            the standard deviation of the Gaussian
	 * @throws IllegalArgumentException
	 *             If the standard deviation is not strictly positive
	 * @throws IllegalArgumentException
	 *             If the sampling is below 1
	 * @throws IllegalArgumentException
	 *             If the maximum kernel size after scaling is too large
	 */
	public RealPoissonGaussianFisherInformation(double s) throws IllegalArgumentException
	{
		this(s, PoissonGaussianFisherInformation.DEFAULT_SAMPLING);
	}

	/**
	 * Instantiates a new real poisson gaussian fisher information.
	 *
	 * @param s
	 *            the standard deviation of the Gaussian
	 * @param sampling
	 *            The number of Gaussian samples to take per standard deviation
	 * @throws IllegalArgumentException
	 *             If the standard deviation is not strictly positive
	 * @throws IllegalArgumentException
	 *             If the sampling is below 1
	 * @throws IllegalArgumentException
	 *             If the maximum kernel size after scaling is too large
	 */
	public RealPoissonGaussianFisherInformation(double s, double sampling) throws IllegalArgumentException
	{
		super(s, sampling);
		kernel = new double[39][];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.PoissonGaussianFisherInformation#getGaussianKernel(int)
	 */
	@Override
	protected double[] getGaussianKernel(int scale, int range)
	{
		// Get the Gaussian kernel
		if (kernel[range] == null)
		{
			kernel[range] = Convolution.makeGaussianKernel(s * scale, range, true);
			//kernel[range] = Convolution.makeGaussianKernel(s * scale, range, false);
			
			// This does not work. The Poisson is discrete and so must be convolved
			// with a continuous Gaussian. Using the Error function is equivalent
			// to integrating the continuous Gaussian with the same value of the 
			// Poisson over the range x-0.5 to x+0.5.
			
			// Lack of granularity in the kernel makes the A^2/P function incorrect. 
			
			//kernel[range] = Convolution.makeErfGaussianKernel(s * scale, range);
		}
		return kernel[range];
	}
}