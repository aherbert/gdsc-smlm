package gdsc.smlm.model;

import xal.tools.math.BesselFunction;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Simple class to calculate AiryPattern
 */
public class AiryPattern
{
	/**
	 * Convert Airy radius to Gaussian approximation
	 * <p>
	 * See Abraham, et al (2009) Opt Express 17: 23352.
	 */
	private static final double FACTOR = 1.323;

	/**
	 * Calculate the intensity of the AiryPattern at distance x from the centre
	 * 
	 * @param x
	 * @return The intensity
	 */
	public static double intensity(final double x)
	{
		if (x == 0)
			return 1;
		final double y = BesselFunction.J1(x) / x;
		return 4.0 * y * y;
	}

	/**
	 * Calculate the intensity of the AiryPattern at distance x from the centre using a Gaussian approximation
	 * 
	 * @param x
	 * @return The intensity
	 */
	public static double intensityGaussian(double x)
	{
		if (x == 0)
			return 1;
		x /= FACTOR;
		return Math.exp(-0.5 * (x * x));
	}

	/**
	 * Calculate the total power of the AiryPattern at distance x from the centre
	 * 
	 * @param x
	 * @return The total power
	 */
	public static double power(final double x)
	{
		if (x == 0)
			return 0;
		final double j0 = BesselFunction.J0(x);
		final double j1 = BesselFunction.J1(x);
		return 1.0 - j0 * j0 - j1 * j1;
	}
}