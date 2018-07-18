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

/**
 * Interface for functions to produce a value, first and second partial derivatives
 */
public interface ExtendedGradient2Procedure
{
	/**
	 * Executes this procedure.
	 *
	 * @param value
	 *            the value of the function
	 * @param dy_da
	 *            Partial first derivative of function with respect to each coefficient (a)
	 * @param d2y_dadb
	 *            Partial second derivative of function with respect to each coefficient pair (a,b). Packed linearly
	 *            with size n*n with n the number of coefficients.
	 */
	public void executeExtended(double value, double[] dy_da, double[] d2y_dadb);
}
