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
 * Provides functions to compute the log likelihood for a distribution.
 */
public interface LogLikelihoodFunction
{
	/**
	 * Compute the log-likelihood of an observation given an expected value
	 * 
	 * @param o
	 *            The observed count
	 * @param e
	 *            The expected count
	 * @return The log-likelihood
	 */
	public double logLikelihood(final double o, final double e);
}