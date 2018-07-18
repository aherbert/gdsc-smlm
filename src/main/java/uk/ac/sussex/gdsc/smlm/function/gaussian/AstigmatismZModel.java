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
package uk.ac.sussex.gdsc.smlm.function.gaussian;

/**
 * Implements a astigmatism model of a 2D Gaussian function, where z-depth determines the x and y width.
 */
public interface AstigmatismZModel extends Cloneable
{
	/**
	 * Gets the standard deviation in dimension x for the given z-depth.
	 *
	 * @param z
	 *            the z
	 * @return
	 * 		the standard deviation in dimension x
	 */
	public double getSx(double z);

	/**
	 * Gets the standard deviation and partial derivatives in dimension x for the given z-depth.
	 *
	 * @param z
	 *            the z
	 * @param ds_dz
	 *            the first derivative of s given z
	 * @return
	 * 		the standard deviation in dimension x
	 */
	public double getSx(double z, double[] ds_dz);

	/**
	 * Gets the standard deviation and partial derivatives in dimension x for the given z-depth.
	 *
	 * @param z
	 *            the z
	 * @param ds_dz
	 *            the first and second derivative of s given z
	 * @return
	 * 		the standard deviation in dimension x
	 */
	public double getSx2(double z, double[] ds_dz);

	/**
	 * Gets the standard deviation in dimension y for the given z-depth.
	 *
	 * @param z
	 *            the z
	 * @return
	 * 		the standard deviation in dimension y
	 */
	public double getSy(double z);

	/**
	 * Gets the standard deviation and partial derivatives in dimension y for the given z-depth.
	 *
	 * @param z
	 *            the z
	 * @param ds_dz
	 *            the first derivative of s given z
	 * @return
	 * 		the standard deviation in dimension y
	 */
	public double getSy(double z, double[] ds_dz);

	/**
	 * Gets the standard deviation and partial derivatives in dimension y for the given z-depth.
	 *
	 * @param z
	 *            the z
	 * @param ds_dz
	 *            the first and second derivative of s given z
	 * @return
	 * 		the standard deviation in dimension y
	 */
	public double getSy2(double z, double[] ds_dz);
}
