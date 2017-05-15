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
 * Defines methods to fit a function with coefficients (a) using weighted least-squares estimation.
 */
public interface WLSEFunctionSolver extends FunctionSolver
{
	/**
	 * Gets the chi squared.
	 *
	 * @return the chi squared
	 */
	public double getChiSquared();

	/**
	 * Gets the probability Q that a value of chi-square as poor as the value should occur by chance.
	 * <p>
	 * A low value indicates greater statistical significance, i.e. greater confidence that the observed deviation from
	 * the null hypothesis is significant, with the null model being that the fit is good. The confidence in rejecting
	 * the null hypothesis is 100 * (1 - q) percent.
	 *
	 * @return the q-value
	 */
	public double getQ();
}
