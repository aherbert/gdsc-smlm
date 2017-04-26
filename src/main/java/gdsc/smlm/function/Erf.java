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
 * Class for computing the error function using approximations.
 * <p>
 * Methods in this class are ordered by speed, not accuracy. The default call to {@link #erf(double)} is a good
 * compromise between both.
 * 
 * @see https://en.wikipedia.org/wiki/Error_function#Approximation_with_elementary_functions
 * @see Winitzki, Sergei (6 February 2008). "A handy approximation for the error function and its inverse".
 *      http://sites.google.com/site/winitzki/sergei-winitzkis-files/erf-approx.pdf
 */
public class Erf
{
	/**
	 * Returns the error function.
	 *
	 * <p>
	 * erf(x) = 2/&radic;&pi; <sub>0</sub>&int;<sup>x</sup> e<sup>-t<sup>2</sup></sup>dt
	 * </p>
	 *
	 * <p>
	 * This implementation computes erf(x) using the approximation by Abramowitz and Stegun. The maximum absolute error
	 * is about 5e-4 for all x.
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
	public static double erf0(double x)
	{
		final boolean negative = (x < 0);
		if (negative)
			x = -x;
		if (x > 40)
			return negative ? -1 : 1;

		final double x2 = x * x;
		final double ret = 1 - 1 / power4(1.0 + 0.278393 * x + 0.230389 * x2 + 0.000972 * x2 * x + 0.078108 * x2 * x2);

		return (negative) ? -ret : ret;
	}

	private static double power4(double d)
	{
		d = d * d; // power 2
		return d * d;
	}

	/**
	 * Returns the difference between erf(x1) and erf(x2). Uses the fast approximation {@link #erf0(double)}.
	 *
	 * @param x1
	 *            the first value
	 * @param x2
	 *            the second value
	 * @return erf(x2) - erf(x1)
	 */
	public static double erf0(double x1, double x2)
	{
		return erf0(x2) - erf0(x1);
	}

	/**
	 * Returns the error function.
	 *
	 * <p>
	 * erf(x) = 2/&radic;&pi; <sub>0</sub>&int;<sup>x</sup> e<sup>-t<sup>2</sup></sup>dt
	 * </p>
	 *
	 * <p>
	 * This implementation computes erf(x) using the approximation by Abramowitz and Stegun. The maximum absolute error
	 * is about 3e-7 for all x.
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
		final boolean negative = (x < 0);
		if (negative)
			x = -x;
		if (x > 40)
			return negative ? -1 : 1;

		final double x2 = x * x;
		final double x3 = x2 * x;
		final double ret = 1 - 1 / power16(1.0 + 0.0705230784 * x + 0.0422820123 * x2 + 0.0092705272 * x3 +
				0.0001520143 * x2 * x2 + 0.0002765672 * x2 * x3 + 0.0000430638 * x3 * x3);

		return (negative) ? -ret : ret;
	}

	private static double power16(double d)
	{
		d = d * d; // power2
		d = d * d; // power4
		d = d * d; // power8
		return d * d;
	}

	/**
	 * Returns the difference between erf(x1) and erf(x2). Uses the fast approximation {@link #erf(double)}.
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

	private final static double four_over_pi = 4.0 / Math.PI;

	/**
	 * Returns the error function.
	 *
	 * <p>
	 * erf(x) = 2/&radic;&pi; <sub>0</sub>&int;<sup>x</sup> e<sup>-t<sup>2</sup></sup>dt
	 * </p>
	 *
	 * <p>
	 * This implementation computes erf(x) using the approximation by Sergei Winitzki. This involves a sqrt() and exp()
	 * function call and so is slower than {@link #erf(double)} and {@link #erf0(double)}. The maximum absolute error is
	 * about 0.00012 for all x.
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
	public static double erf2(double x)
	{
		final boolean negative = (x < 0);
		if (negative)
			x = -x;
		if (x > 40)
			return negative ? -1 : 1;

		final double x2 = x * x;
		final double ax2 = 0.147 * x2;
		final double ret = Math.sqrt(1 - FastMath.exp(-x2 * (four_over_pi + ax2) / (1 + ax2)));

		return negative ? -ret : ret;
	}

	/**
	 * Returns the difference between erf(x1) and erf(x2). Uses the fast approximation {@link #erf2(double)}.
	 *
	 * @param x1
	 *            the first value
	 * @param x2
	 *            the second value
	 * @return erf(x2) - erf(x1)
	 */
	public static double erf2(double x1, double x2)
	{
		return erf2(x2) - erf2(x1);
	}
}
