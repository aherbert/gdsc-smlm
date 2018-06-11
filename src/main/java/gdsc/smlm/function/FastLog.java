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
package gdsc.smlm.function;

/**
 * Base class for an algorithm computing a fast approximation to the log function
 */
public abstract class FastLog
{
	/**
	 * Natural logarithm of 2
	 */
	public static final double LN2 = Math.log(2);
	/**
	 * Natural logarithm of 2
	 */
	public static final float LN2F = (float) LN2;

	/** The default value for n (the number of bits to keep from the mantissa) */
	public static final int N = 13;

	/**
	 * Gets the scale to convert from log2 to logB (using the given base) by multiplication.
	 * 
	 * <pre>
	 * scale = Math.log(2) / Math.log(base)
	 * </pre>
	 *
	 * @return the scale
	 */
	public static double getScale(double base)
	{
		if (!(base > 0 && base != Double.POSITIVE_INFINITY))
			throw new IllegalArgumentException("Base must be a real positive number");
		return LN2 / Math.log(base);
	}

	/**
	 * Calculate the logarithm with base 2.
	 *
	 * <pre>
	 * return Math.log(val) / Math.log(2)
	 * </pre>
	 * 
	 * Math.log(2) is pre-computed.
	 * 
	 * @param val
	 *            the input value
	 * @return the log2 of the value
	 */
	public static double exactLog2(double val)
	{
		return Math.log(val) / LN2;
	}

	/**
	 * Gets the base.
	 *
	 * @return the base
	 */
	public abstract double getBase();

	/**
	 * Gets the scale to convert from log2 to logB (the base given by {@link #getBase()}) by multiplication.
	 *
	 * @return the scale
	 */
	public abstract double getScale();

	/**
	 * Gets the number of most significant bits to keep from the mantissa, i.e. the binary precision of the floating
	 * point number before computing the log.
	 *
	 * @return the number of most significant bits to keep
	 */
	public abstract int getN();

	/**
	 * Calculate the logarithm to base 2, handling special cases.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN or less than zero, then the result is NaN.
	 * <li>If the argument is positive infinity, then the result is positive infinity.
	 * <li>If the argument is positive zero or negative zero, then the result is negative infinity.
	 * </ul>
	 *
	 * @param x
	 *            the argument
	 * @return log(x)
	 */
	public abstract float log2(float x);

	/**
	 * Calculate the logarithm to base 2. Requires the argument be finite and strictly positive.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN, then the result is incorrect.
	 * <li>If the argument is negative, then the result is incorrect.
	 * <li>If the argument is positive infinity, then the result is incorrect.
	 * <li>If the argument is positive zero or negative zero, then the result is incorrect.
	 * </ul>
	 * <p>
	 * Sub-classes may handle some of the special cases.
	 * 
	 * @param x
	 *            the argument (must be strictly positive)
	 * @return log(x)
	 */
	public abstract float fastLog2(float x);

	/**
	 * Calculate the logarithm to the configured base, handling special cases.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN or less than zero, then the result is NaN.
	 * <li>If the argument is positive infinity, then the result is positive infinity.
	 * <li>If the argument is positive zero or negative zero, then the result is negative infinity.
	 * </ul>
	 *
	 * @param x
	 *            the argument
	 * @return log(x)
	 * @see #getBase()
	 */
	public abstract float log(float x);

	/**
	 * Calculate the logarithm to the configured base. Requires the argument be finite and strictly positive.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN, then the result is incorrect.
	 * <li>If the argument is negative, then the result is incorrect.
	 * <li>If the argument is positive infinity, then the result is incorrect.
	 * <li>If the argument is positive zero or negative zero, then the result is incorrect.
	 * </ul>
	 * <p>
	 * Sub-classes may handle some of the special cases.
	 * 
	 * @param x
	 *            the argument (must be strictly positive)
	 * @return log(x)
	 * @see #getBase()
	 */
	public abstract float fastLog(float x);

	/**
	 * Calculate the logarithm to base 2, handling special cases.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN or less than zero, then the result is NaN.
	 * <li>If the argument is positive infinity, then the result is positive infinity.
	 * <li>If the argument is positive zero or negative zero, then the result is negative infinity.
	 * </ul>
	 *
	 * @param x
	 *            the argument
	 * @return log(x)
	 */
	public abstract float log2(double x);

	/**
	 * Calculate the logarithm to base 2. Requires the argument be finite and strictly positive.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN, then the result is incorrect.
	 * <li>If the argument is negative, then the result is incorrect.
	 * <li>If the argument is positive infinity, then the result is incorrect.
	 * <li>If the argument is positive zero or negative zero, then the result is incorrect.
	 * </ul>
	 * <p>
	 * Sub-classes may handle some of the special cases.
	 * 
	 * @param x
	 *            the argument (must be strictly positive)
	 * @return log(x)
	 */
	public abstract float fastLog2(double x);

	/**
	 * Calculate the logarithm to base 2, handling special cases.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN or less than zero, then the result is NaN.
	 * <li>If the argument is positive infinity, then the result is positive infinity.
	 * <li>If the argument is positive zero or negative zero, then the result is negative infinity.
	 * </ul>
	 * <p>
	 * This can be re-implemented if {@link #fastLog(double)} can be returned in double precision.
	 *
	 * @param x
	 *            the argument
	 * @return log(x)
	 */
	public double log2D(double x)
	{
		return log2(x);
	}

	/**
	 * Calculate the logarithm to base 2. Requires the argument be finite and strictly positive.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN, then the result is incorrect.
	 * <li>If the argument is negative, then the result is incorrect.
	 * <li>If the argument is positive infinity, then the result is incorrect.
	 * <li>If the argument is positive zero or negative zero, then the result is incorrect.
	 * </ul>
	 * <p>
	 * Sub-classes may handle some of the special cases.
	 * <p>
	 * This can be re-implemented if {@link #fastLog(double)} can be returned in double precision.
	 * 
	 * @param x
	 *            the argument (must be strictly positive)
	 * @return log(x)
	 */
	public double fastLog2D(double x)
	{
		return fastLog2(x);
	}

	/**
	 * Calculate the logarithm to the configured base, handling special cases.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN or less than zero, then the result is NaN.
	 * <li>If the argument is positive infinity, then the result is positive infinity.
	 * <li>If the argument is positive zero or negative zero, then the result is negative infinity.
	 * </ul>
	 *
	 * @param x
	 *            the argument
	 * @return log(x)
	 * @see #getBase()
	 */
	public abstract float log(double x);

	/**
	 * Calculate the logarithm to the configured base. Requires the argument be finite and strictly positive.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN, then the result is incorrect.
	 * <li>If the argument is negative, then the result is incorrect.
	 * <li>If the argument is positive infinity, then the result is incorrect.
	 * <li>If the argument is positive zero or negative zero, then the result is incorrect.
	 * </ul>
	 * <p>
	 * Sub-classes may handle some of the special cases.
	 * 
	 * @param x
	 *            the argument (must be strictly positive)
	 * @return log(x)
	 * @see #getBase()
	 */
	public abstract float fastLog(double x);

	/**
	 * Calculate the logarithm to the configured base, handling special cases.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN or less than zero, then the result is NaN.
	 * <li>If the argument is positive infinity, then the result is positive infinity.
	 * <li>If the argument is positive zero or negative zero, then the result is negative infinity.
	 * </ul>
	 * <p>
	 * This can be re-implemented if {@link #log(double)} can be returned in double precision.
	 *
	 * @param x
	 *            the argument
	 * @return log(x)
	 * @see #getBase()
	 */
	public double logD(double x)
	{
		return log(x);
	}

	/**
	 * Calculate the logarithm to the configured base. Requires the argument be finite and strictly positive.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN, then the result is incorrect.
	 * <li>If the argument is negative, then the result is incorrect.
	 * <li>If the argument is positive infinity, then the result is incorrect.
	 * <li>If the argument is positive zero or negative zero, then the result is incorrect.
	 * </ul>
	 * <p>
	 * Sub-classes may handle some of the special cases.
	 * <p>
	 * This can be re-implemented if {@link #fastLog(double)} can be returned in double precision.
	 * 
	 * @param x
	 *            the argument (must be strictly positive)
	 * @return log(x)
	 * @see #getBase()
	 */
	public double fastLogD(double x)
	{
		return fastLog(x);
	}
}
