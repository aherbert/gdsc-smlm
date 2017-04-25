package gdsc.smlm.function;

import org.apache.commons.math3.util.Pair;

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
 * Defines the non-linear fitting function
 */
public interface ExtendedNonLinearFunction extends NonLinearFunction
{
	/**
	 * Compute the values of all the data points
	 * 
	 * @param variables
	 * @return The values
	 */
	public double[] computeValues(double[] variables);

	/**
	 * Compute the Jacobian of the gradients of all the data points
	 * 
	 * @param variables
	 * @return The Jacobian
	 */
	public double[][] computeJacobian(double[] variables);

	/**
	 * Return true if the function can compute the values and the Jacobian.
	 *
	 * @return true, if the function can compute the values and the Jacobian.
	 */
	public boolean canComputeValuesAndJacobian();

	/**
	 * Compute the values and the Jacobian of the gradients of all the data points
	 * 
	 * @param variables
	 * @return The values and the Jacobian
	 */
	public Pair<double[], double[][]> computeValuesAndJacobian(double[] variables);
}
