package gdsc.smlm.fitting;

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
 * Calculator for the Fisher information, a symmetric positive definite matrix containing the amount of information that
 * an observable random variable X carries about an unknown parameter θ of a distribution that models X.
 * <p>
 * The calculator will compute the Fisher Information Matrix (I) using numerical gradients:
 * 
 * <pre>
 * Iij = E [ (d log f(X;θ) / dθi) (d log f(X;θ) / dθj) | θ ]
 * E = Expected value
 * f(X;θ) = Likelihood function for data X given parameters θ
 * </pre>
 */
public interface FisherInformationCalculator
{
	/**
	 * Compute the Fisher information, a symmetric positive definite matrix containing the amount of information that
	 * an observable random variable X carries about an unknown parameter θ of a distribution that models X.
	 *
	 * @param parameters
	 *            the parameters (θ)
	 * @return the fisher information matrix
	 */
	public FisherInformationMatrix compute(double[] parameters);
}