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
 * Calculate the Fisher information for a Poisson-distributed random variable using an interpolation
 * of the alpha scale parameter. The alpha scale parameter is the ratio between the Fisher information
 * for a Poisson distribution and the Fisher information of another Poisson-based distribution,
 * e.g. a Poisson-Gaussian convolution.
 */
public class InterpolatedPoissonFisherInformation extends BasePoissonFisherInformation
{
	// TODO - store a polynomial spline function to interpolate the alpha parameter.
	// This must use a log x scale (log(Poisson mean)).

	// At the minimum the alpha is fixed to the min.

	// At the maximum the Fisher information can use a fixed alpha or a approximation
	// function.

	// Within the range the poisson mean is converted to a log scale and alpha is
	// interpolated.

	// The Fisher information is returned using the alpha multiplied by the 
	// Poisson Fisher information.

	/**
	 * Instantiates a new poisson gaussian fisher information.
	 */
	public InterpolatedPoissonFisherInformation() throws IllegalArgumentException
	{
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
		// Poisson fisher information
		double I = 1.0 / t;
		if (I != Double.POSITIVE_INFINITY)
			I *= getAlpha(t);
		return I;
	}

	private double getAlpha(double t)
	{
		// TODO
		return 1;
	}

	@Override
	protected void postClone()
	{
		// Nothing to do
	}
}