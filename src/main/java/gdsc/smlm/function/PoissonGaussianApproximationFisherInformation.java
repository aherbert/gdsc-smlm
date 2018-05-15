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
 * Calculate the Fisher information for a Poisson-Gaussian distribution using an approximation of the Poisson (mean=t)
 * as a Gaussian (u=t, var=t).
 */
public class PoissonGaussianApproximationFisherInformation extends BasePoissonFisherInformation
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
	
	@Override
	protected void postClone()
	{
		// Nothing to do
	}
}