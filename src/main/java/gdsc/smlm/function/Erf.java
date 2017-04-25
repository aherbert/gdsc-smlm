package gdsc.smlm.function;

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

import org.apache.commons.math3.util.FastMath;

/**
 * Class for computing the error function
 * <p>
 * The implementation is based upon an approximation:
 * Winitzki, Sergei (6 February 2008). "A handy approximation for the error function and its inverse".
 * http://sites.google.com/site/winitzki/sergei-winitzkis-files/erf-approx.pdf
 */
public class Erf
{
	private final static double a = 0.147;
	private final static double four_over_pi = 4.0 / Math.PI;

	/**
	 * Returns the error function.
	 *
	 * <p>
	 * erf(x) = 2/&radic;&pi; <sub>0</sub>&int;<sup>x</sup> e<sup>-t<sup>2</sup></sup>dt
	 * </p>
	 *
	 * <p>
	 * This implementation computes erf(x) using the approximation by Sergei Winitzki. The maximum error is about
	 * 0.00012 for all x.
	 * </p>
	 *
	 * <p>
	 * The value returned is always between -1 and 1 (inclusive).
	 * If {@code abs(x) > 40}, then {@code erf(x)} is indistinguishable from
	 * either 1 or -1 as a double, so the appropriate extreme value is returned.
	 * </p>
	 *
	 * @param x
	 *            the value.
	 * @return the error function erf(x)
	 */
	public static double erf(double x)
	{
		if (FastMath.abs(x) > 40)
		{
			return x > 0 ? 1 : -1;
		}
		final double x2 = x * x;
		final double ax2 = a * x2;
		final double ret = Math.sqrt(1 - FastMath.exp(-x2 * (four_over_pi + ax2) / (1 + ax2)));
		return x < 0 ? -ret : ret;
	}

	/**
	 * Returns the difference between erf(x1) and erf(x2).
	 *
	 * @param x1
	 *            the first value
	 * @param x2
	 *            the second value
	 * @return erf(x2) - erf(x1)
	 */
	public static double erf(double x1, double x2)
	{
		return erf(x2) - erf(x1);
	}
}
