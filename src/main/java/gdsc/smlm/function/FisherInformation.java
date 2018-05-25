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
 * Calculate the Fisher information: the amount of information that an observable random variable X carries about an
 * unknown parameter θ of a distribution that models X.
 * 
 * <pre>
 * I = E [ (d log f(X;θ) / dθi) (d log f(X;θ) / dθj) | θ ]
 * E = Expected value
 * </pre>
 */
public interface FisherInformation
{
	/**
	 * Gets the fisher information: the amount of information that
	 * an observable random variable X carries about an unknown parameter θ of a distribution that models X.
	 *
	 * @param t
	 *            parameter θ of a distribution that models X
	 * @return the fisher information
	 * @throws IllegalArgumentException
	 *             if the parameter is not in the valid range
	 */
	public double getFisherInformation(double t) throws IllegalArgumentException;

	/**
	 * Checks if the parameter θ is in a valid range to compute a representable value.
	 * <p>
	 * If not true then it would be expected that
	 * {@link #getFisherInformation(double)} will: throw an exception; compute zero; or compute infinity.
	 *
	 * @param t
	 *            parameter θ of a distribution that models X
	 * @return true, if a representable value can be computed
	 */
	public boolean isValid(double t);
}