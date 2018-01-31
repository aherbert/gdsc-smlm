package gdsc.smlm.function;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (c) 2004-2016 Hanns Holger Rutz. All rights reserved.
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Implementation of the ICSILog algorithm
 * as described in O. Vinyals, G. Friedland, N. Mirghafori
 * "Revisiting a basic function on current CPUs: A fast logarithm implementation
 * with adjustable accuracy" (2007).
 * <p>
 * This class is based on the original algorithm description and a Java implementation by Hanns Holger Rutz.
 *
 * @see <a href=
 *      "https://www.javatips.net/api/Eisenkraut-master/src/main/java/de/sciss/eisenkraut/math/FastLog.java">https://www.javatips.net/api/Eisenkraut-master/src/main/java/de/sciss/eisenkraut/math/FastLog.java</a>
 * @see <a href=
 *      "http://www.icsi.berkeley.edu/pubs/techreports/TR-07-002.pdf">http://www.icsi.berkeley.edu/pubs/techreports/TR-07-002.pdf</a>
 */
public class FFastLog
{
	/**
	 * Natural logarithm of 2
	 */
	public static final double LN2 = Math.log(2);

	/** The default value for q (the number of bits to remove from the mantissa) */
	public static final int Q = 11;

	private final int q, qM1;
	private final float[] data;
	private float scale;

	/**
	 * Create a new natural logarithm calculation instance. This will
	 * hold the pre-calculated log values for base E
	 * and a table size of 11.
	 */
	public FFastLog()
	{
		this(Math.E, Q);
	}

	/**
	 * Create a new logarithm calculation instance. This will
	 * hold the pre-calculated log values for base E
	 * and a table size depending on a given mantissa quantization.
	 * 
	 * @param q
	 *            the quantization, the number of bits to remove
	 *            from the mantissa. for q = 11, the table storage
	 *            requires 32 KB.
	 */
	public FFastLog(int q)
	{
		this(Math.E, q);
	}

	/**
	 * Create a new logarithm calculation instance. This will
	 * hold the pre-calculated log values for a given base
	 * and a table size depending on a given mantissa quantization.
	 * 
	 * @param base
	 *            the logarithm base (e.g. 2 for log duals, 10 for
	 *            decibels calculations, Math.E for natural log)
	 * @param q
	 *            the quantization, the number of bits to remove
	 *            from the mantissa. for q = 11, the table storage
	 *            requires 32 KB.
	 */
	public FFastLog(double base, int q)
	{
		if (!(base > 0 && base != Double.POSITIVE_INFINITY))
			throw new IllegalArgumentException("Base must be a real positive number");
		if (q < 0 || q > 23)
			throw new IllegalArgumentException("Q must be in the range 0<=q<=23");

		final int p = 1 << (24 - q);

		this.q = q;
		qM1 = q - 1;
		scale = (float) (LN2 / Math.log(base));
		data = new float[p];

		for (int i = 0; i < p; i++)
		{
			// Store log2 value of a range of floating point numbers using a limited
			// precision mantissa (m).

			// See Float.intBitsToFloat(int):
			// int s = ((bits >> 31) == 0) ? 1 : -1;
			// int e = ((bits >> 23) & 0xff);
			// int m = (e == 0) ?
			//                 (bits & 0x7fffff) << 1 :
			//                 (bits & 0x7fffff) | 0x800000;
			//
			// Then the floating-point result equals the value of the mathematical
			// expression s x m x 2^(e-150)

			// int m = i << q
			// int e = 0
			// log2(m x 2^(e-150)) == log2(m) + log2(2^-150) = log2(m) -150 

			data[i] = (float) (log2(i << q) - 150);
		}
	}

	/**
	 * Calculate the logarithm with base 2.
	 *
	 * @param val
	 *            the input value
	 * @return the log2 of the value
	 */
	public static double log2(double val)
	{
		return Math.log(val) / LN2;
	}

	/**
	 * Calculate the logarithm to the base given in the constructor. Requires the argument be finite and positive.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN, then the result is incorrect.
	 * <li>If the argument is negative, then the result is incorrect (log(-x)).
	 * <li>If the argument is positive infinity, then the result is incorrect (Math.log(Float.MAX_VALUE)).
	 * <li>If the argument is positive zero or negative zero, then the result is negative infinity.
	 * </ul>
	 * 
	 * @param x
	 *            the argument. must be positive!
	 * @return log( x )
	 */
	public float fastLog(float x)
	{
		final int bits = Float.floatToRawIntBits(x);

		// Note the documentation from Float.intBitsToFloat(int):
		// int e = ((bits >> 23) & 0xff);
		// int m = (e == 0) ?
		//                 (bits & 0x7fffff) << 1 :
		//                 (bits & 0x7fffff) | 0x800000;

		final int e = (bits >> 23) & 0xff;
		// raw mantissa, conversion is done with the bit shift to reduce precision
		final int m = (bits & 0x7fffff);

		return (e == 0 ? data[m >> qM1] : e + data[((m | 0x00800000) >> q)]) * scale;
	}

	/**
	 * Calculate the logarithm to the base given in the constructor, handling special cases.
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
	 * @return log( x )
	 */
	public float log(float x)
	{
		// Basic implementation
		//if (x > 0)
		//	return (x == Float.POSITIVE_INFINITY) ? Float.POSITIVE_INFINITY : fastLog(x);
		//return (x == 0) ? Float.NEGATIVE_INFINITY : Float.NaN;

		// Re-implement to avoid float comparisons (which will be slower than int comparisons) 
		final int bits = Float.floatToRawIntBits(x);

		// Note the documentation from Float.intBitsToFloat(int):
		// int e = ((bits >> 23) & 0xff);
		// int m = (e == 0) ?
		//                 (bits & 0x7fffff) << 1 :
		//                 (bits & 0x7fffff) | 0x800000;

		final boolean negative = (bits >> 31) != 0;
		final int e = (bits >> 23) & 0xff;
		// raw mantissa, conversion is done with the bit shift to reduce precision
		final int m = (bits & 0x7fffff);

		if (e == 255)
		{
			// All bits set is a special case
			if (m != 0)
				return Float.NaN;
			// +/- Infinity
			return (negative) ? Float.NaN : Float.POSITIVE_INFINITY;
		}

		if (negative)
			// Only -0 is allowed
			return (e == 0 && m == 0) ? Float.NEGATIVE_INFINITY : Float.NaN;

		return (e == 0 ? data[m >> qM1] : e + data[((m | 0x00800000) >> q)]) * scale;
	}

	/**
	 * Calculate the logarithm to the base given in the constructor. Requires the argument be finite and positive.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN, then the result is incorrect.
	 * <li>If the argument is negative, then the result is incorrect (log(-x)).
	 * <li>If the argument is positive infinity, then the result is incorrect (Math.log(Float.MAX_VALUE)).
	 * <li>If the argument is positive zero or negative zero, then the result is negative infinity.
	 * </ul>
	 * 
	 * @param x
	 *            the argument. must be positive!
	 * @return log( x )
	 */
	public double fastLog(double x)
	{
		return fastLog((float) x);
	}

	/**
	 * Calculate the logarithm to the base given in the constructor, handling special cases.
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
	 * @return log( x )
	 */
	public double log(double x)
	{
		if (x > 0)
			return (x == Double.POSITIVE_INFINITY) ? Double.POSITIVE_INFINITY : fastLog(x);
		return (x == 0) ? Double.NEGATIVE_INFINITY : Double.NaN;
	}
}
