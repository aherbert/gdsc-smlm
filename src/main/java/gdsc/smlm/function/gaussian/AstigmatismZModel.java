package gdsc.smlm.function.gaussian;

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
