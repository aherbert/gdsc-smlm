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
 * Implements a no astigmatism model of a 2D Gaussian function. The width is constant.
 */
public class NullAstigmatismZModel implements AstigmatismZModel
{
	/** The width in x. */
	public final double sx;
	/** The width in y. */
	public final double sy;

	/**
	 * Instantiates a new null astigmatism Z model.
	 *
	 * @param sx
	 *            the sx
	 * @param sy
	 *            the sy
	 * @throws IllegalArgumentException
	 *             if the widths are not positive
	 */
	public NullAstigmatismZModel(double sx, double sy) throws IllegalArgumentException
	{
		if (!(sx > 0 && sy > 0))
			throw new IllegalArgumentException("Width must be positive");
		this.sx = sx;
		this.sy = sy;
	}

	public double getSx(double z)
	{
		return sx;
	}

	public double getSx(double z, double[] ds_dz)
	{
		ds_dz[0] = 0;
		return sx;
	}

	public double getSx2(double z, double[] ds_dz)
	{
		ds_dz[0] = 0;
		ds_dz[1] = 0;
		return sx;
	}

	public double getSy(double z)
	{
		return sy;
	}

	public double getSy(double z, double[] ds_dz)
	{
		ds_dz[0] = 0;
		return sy;
	}

	public double getSy2(double z, double[] ds_dz)
	{
		ds_dz[0] = 0;
		ds_dz[1] = 0;
		return sy;
	}
}
