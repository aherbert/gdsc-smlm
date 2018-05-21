package gdsc.smlm.function;

import gdsc.smlm.utils.Convolution;
import gdsc.smlm.utils.GaussianKernel;

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
 * Calculate the Fisher information for a continuous Poisson-Gamma-Gaussian distribution.
 * <p>
 * Uses a modified form of the equation of Chao, et al (2013) Nature Methods, 10, 335-338, SI Eq S8.
 * <p>
 * Performs a convolution with a finite Gaussian kernel.
 */
public class RealPoissonGammaGaussianFisherInformation extends PoissonGammaGaussianFisherInformation
{
	/**
	 * The Gaussian convolution kernels for different ranges.
	 */
	private final double[][][] kernel;

	/**
	 * Instantiates a new real poisson gaussian fisher information.
	 *
	 * @param m
	 *            the gain multiplication factor
	 * @param s
	 *            the standard deviation of the Gaussian
	 * @throws IllegalArgumentException
	 *             If the standard deviation is not strictly positive
	 * @throws IllegalArgumentException
	 *             If the sampling is below 1
	 * @throws IllegalArgumentException
	 *             If the maximum kernel size after scaling is too large
	 */
	public RealPoissonGammaGaussianFisherInformation(double m, double s) throws IllegalArgumentException
	{
		this(m, s, PoissonGammaGaussianFisherInformation.DEFAULT_SAMPLING);
	}

	/**
	 * Instantiates a new real poisson gaussian fisher information.
	 *
	 * @param m
	 *            the gain multiplication factor
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
	public RealPoissonGammaGaussianFisherInformation(double m, double s, double sampling)
			throws IllegalArgumentException
	{
		super(m, s, sampling);
		kernel = new double[LOG_2_MAX_SCALE + 1][39][];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.PoissonGammaGaussianFisherInformation#getUnitGaussianKernel(int)
	 */
	@Override
	protected double[] getUnitGaussianKernel(int scale, int range)
	{
		int index = getIndex(scale);
		if (kernel[index][range] == null)
		{
			kernel[index][range] = GaussianKernel.makeGaussianKernel(scale, range, true);

			// This does not work as the lack of granularity in the 
			// kernel makes the A^2/P function incorrect.

			//kernel[range] = GaussianKernel.makeErfGaussianKernel(scale, range);
		}
		return kernel[index][range];
	}
}