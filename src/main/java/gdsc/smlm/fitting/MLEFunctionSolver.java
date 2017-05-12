package gdsc.smlm.fitting;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Defines methods to fit a function with coefficients (a) using maximum likelihood estimation.
 */
public interface MLEFunctionSolver extends FunctionSolver
{
	/**
	 * Gets the log likelihood.
	 *
	 * @return the log likelihood
	 */
	public double getLogLikelihood();

	/**
	 * Gets the log likelihood ratio.
	 *
	 * @return the log likelihood ratio
	 */
	public double getLogLikelihoodRatio();

	/**
	 * Gets the probability Q that a value of the log likelihood ratio as poor as the value should occur by chance.
	 *
	 * @return the q-value
	 */
	public double getQ();
}
