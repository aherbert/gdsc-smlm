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
 * Calculate the Fisher information for a continuous Poisson-Gamma-Gaussian distribution.
 * <p>
 * Uses a modified form of the equation of Chao, et al (2013) Nature Methods, 10, 335-338, SI Eq S8.
 * <p>
 * Performs a convolution with a finite Gaussian kernel.
 */
public class RealPoissonGammaGaussianFisherInformation extends PoissonGammaGaussianFisherInformation
{
	/**
	 * The Gaussian convolution kernels for different scaling. The scale is 2^index, e.g. 1, 2, 4, 8, 16, 32, 64, 128.
	 */
	private final double[][] kernel;

	/**
	 * The Erf Gaussian convolution kernels for different scaling. The scale is 2^index, e.g. 1, 2, 4, 8, 16, 32, 64, 128.
	 */
	private final double[][] erfKernel;
	
	/**
	 * Instantiates a new real poisson gaussian fisher information.
	 *
	 * @param m
	 *            the gain multiplication factor
	 * @param s
	 *            the standard deviation of the Gaussian
	 * @throws IllegalArgumentException
	 *             If the standard deviation is not strictly positive
	 */
	public RealPoissonGammaGaussianFisherInformation(double m, double s) throws IllegalArgumentException
	{
		this(m, s, 6);
	}

	/**
	 * Instantiates a new real poisson gaussian fisher information.
	 *
	 * @param m
	 *            the gain multiplication factor
	 * @param s
	 *            the standard deviation of the Gaussian
	 * @param range
	 *            the range of the Gaussian kernel (in SD units). This is clipped to the range
	 *            1-38 to provide a meaningful convolution.
	 * @throws IllegalArgumentException
	 *             If the standard deviation is not strictly positive
	 */
	public RealPoissonGammaGaussianFisherInformation(double m, double s, double range) throws IllegalArgumentException
	{
		super(m, s, range);

		// Store the Gaussian kernels for convolution:
		// 1, 2, 4, 8, 16, 32, 64, 128
		kernel = new double[8][];
		erfKernel = new double[kernel.length][];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.PoissonGammaGaussianFisherInformation#getUnitGaussianKernel(int)
	 */
	@Override
	protected double[] getUnitGaussianKernel(int scale)
	{
		// Get the Gaussian kernel
		int index = log2(scale);
		if (kernel[index] == null)
			kernel[index] = Convolution.makeGaussianKernel(scale, range, false);
		return kernel[index];
	}

	/* (non-Javadoc)
	 * @see gdsc.smlm.function.PoissonGammaGaussianFisherInformation#getUnitErfGaussianKernel(int)
	 */
	@Override
	protected double[] getUnitErfGaussianKernel(int scale)
	{
		// Get the Gaussian kernel
		int index = log2(scale);
		if (erfKernel[index] == null)
			erfKernel[index] = Convolution.makeErfGaussianKernel(scale, range);
		return erfKernel[index];
	}
	
	private int log2(int scale)
	{
		int bits = 8;
		while ((scale & (1 << bits)) == 0)
			bits--;
		return bits;
	}
}