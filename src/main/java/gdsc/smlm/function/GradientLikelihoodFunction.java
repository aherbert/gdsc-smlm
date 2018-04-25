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
 * Provides functions to compute the gradient of the likelihood for a univariate distribution.
 */
public interface GradientLikelihoodFunction extends LikelihoodFunction
{
	/**
	 * Compute the likelihood of an observation x given a parameter value theta.
	 * <p>
	 * This is the probability mass function P(X=x|θ) or the probability density function f(x|θ) depending on parameter
	 * θ.
	 *
	 * @param o
	 *            The observed value (x)
	 * @param t
	 *            The parameter value (θ)
	 * @param dp_dt
	 *            the gradient d P(X=x|θ) dθ
	 * @return The likelihood
	 */
	public double likelihood(final double o, final double t, final double[] dp_dt);
}