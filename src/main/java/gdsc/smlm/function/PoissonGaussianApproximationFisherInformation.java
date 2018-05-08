package gdsc.smlm.function;

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
 * Calculate the Fisher information for a Poisson-Gaussian distribution.
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
public class PoissonGaussianApproximationFisherInformation implements FisherInformation
{
	/** The variance of the Gaussian. */
	public final double variance;

	/**
	 * Instantiates a new poisson gaussian fisher information.
	 *
	 * @param s
	 *            the standard deviation of the Gaussian
	 * @throws IllegalArgumentException
	 *             If the standard deviation is not strictly positive
	 */
	public PoissonGaussianApproximationFisherInformation(double s) throws IllegalArgumentException
	{
		if (!(s > 0 && s <= Double.MAX_VALUE))
			throw new IllegalArgumentException("Gaussian variance must be strictly positive");
		this.variance = s * s;
	}

	/**
	 * {@inheritDoc}
	 * <p>
	 * Gets the approximate Poisson-Gaussian Fisher information.
	 * Approximate the Poisson as a Gaussian (u=t, var=t) and convolve with a Gaussian (u=0,var=s*s).
	 * Gaussian-Gaussian convolution: var1 * var2 => var = var1+var2.
	 * The Fisher information of Gaussian mean is 1/variance.
	 * The Poisson-Gaussian Fisher information is therefore 1 / (t + s*s).
	 * 
	 * @see gdsc.smlm.function.FisherInformation#getFisherInformation(double)
	 */
	public double getFisherInformation(double t) throws IllegalArgumentException
	{
		if (t <= 0)
			throw new IllegalArgumentException("Poisson mean must be positive");
		return 1.0 / (t + variance);
	}
}