/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package uk.ac.sussex.gdsc.smlm.function;

import uk.ac.sussex.gdsc.smlm.utils.Pair;

/**
 * Defines the non-linear fitting function
 */
public interface ExtendedNonLinearFunction extends NonLinearFunction
{
	/**
	 * Compute the values of all the data points.
	 *
	 * @param variables
	 *            the variables
	 * @return The values
	 */
	public double[] computeValues(double[] variables);

	/**
	 * Compute the Jacobian of the gradients of all the data points.
	 *
	 * @param variables
	 *            the variables
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
	 * Compute the values and the Jacobian of the gradients of all the data points.
	 *
	 * @param variables
	 *            the variables
	 * @return The values and the Jacobian
	 */
	public Pair<double[], double[][]> computeValuesAndJacobian(double[] variables);
}
