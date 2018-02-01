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
 * This class is based on the original algorithm description and a Java implementation by Hanns Holger Rutz. It has been
 * adapted for use with double-precision data.
 * <p>
 * Note: log(float) is not provided as it is dynamically cast up to a double and float values can be represented as a
 * double. If your computation generates a float value to be logged then use {@link FFastLog }.
 *
 * @see <a href=
 *      "https://www.javatips.net/api/Eisenkraut-master/src/main/java/de/sciss/eisenkraut/math/FastLog.java">https://www.javatips.net/api/Eisenkraut-master/src/main/java/de/sciss/eisenkraut/math/FastLog.java</a>
 * @see <a href=
 *      "http://www.icsi.berkeley.edu/pubs/techreports/TR-07-002.pdf">http://www.icsi.berkeley.edu/pubs/techreports/TR-07-002.pdf</a>
 */
public class DFastLog extends BaseFastLog
{
	private final double base;
	private final int q, q_minus_1;
	private final float[] data;
	private final double scale;

	/**
	 * Create a new natural logarithm calculation instance.
	 */
	public DFastLog()
	{
		this(Math.E, N);
	}

	/**
	 * Create a new logarithm calculation instance. This will
	 * hold the pre-calculated log values for base E
	 * and a table size depending on a given mantissa quantization.
	 * 
	 * @param n
	 *            The number of bits to keep from the mantissa.
	 *            Table storage = 2^(n+1) * 4 bytes, e.g. 64Kb for n=13
	 */
	public DFastLog(int n)
	{
		this(Math.E, n);
	}

	/**
	 * Create a new logarithm calculation instance. This will
	 * hold the pre-calculated log values for a given base
	 * and a table size depending on a given mantissa quantization.
	 * 
	 * @param base
	 *            the logarithm base (e.g. 2 for log duals, 10 for
	 *            decibels calculations, Math.E for natural log)
	 * @param n
	 *            The number of bits to keep from the mantissa.
	 *            Table storage = 2^(n+1) * 4 bytes, e.g. 64Kb for n=13
	 */
	public DFastLog(double base, int n)
	{
		// Note: Only support the max precision we can store in a table.
		// Note: A large table would be very costly to compute!
		if (n < 0 || n > 30)
			throw new IllegalArgumentException("N must be in the range 0<=n<=30");
		scale = getScale(base);
		this.base = base;

		final int size = 1 << (n + 1);

		q = 52 - n;
		q_minus_1 = q - 1;
		data = new float[size];

		for (int i = 0; i < size; i++)
		{
			// Store log2 value of a range of floating point numbers using a limited
			// precision mantissa (m). The purpose of this code is to enumerate all 
			// possible mantissas of a double with limited precision (52-q). Note the
			// 53 for a 52-bit mantissa comes from the fact that the mantissa 
			// represents the digits of a binary number after the binary-point: .10101010101
			// It is assumed that the digit before the point is a 1 if the exponent
			// is non-zero. Otherwise the binary point is moved to the right of the first
			// digit (i.e. a bit shift left).

			// See Double.longBitsToDouble(int):
			// int s = ((bits >>> 63) == 0) ? 1 : -1;
			// int e = (int)((bits >>> 52) & 0x7ffL);
			// long m = (e == 0) ?
			//                 (bits & 0xfffffffffffffL) << 1 :
			//                 (bits & 0xfffffffffffffL) | 0x10000000000000L;
			//
			// Then the floating-point result equals the value of the mathematical
			// expression s x m x 2^(e-1075)

			// For a precision of n=(52-q)=6
			// We enumerate:
			// ( .000000 to  .111111) * 2^q (i.e. q additional zeros) 
			// (1.000000 to 1.111111) * 2^q

			// The bit shift is performed on integer data to construct the desired mantissa
			// which is then converted to a double for the call to exactLog2(double).

			// int m = i << q
			// log2(m x 2^(e-1075)) == log2(m) + log2(2^e-1075) = log2(m) +(e-1075)
			// We subtract the -1075 here so that the log2(float) can be reconstructed
			// from the table of log2(m) + e.

			data[i] = (float) (exactLog2(((long) i) << q) - 1075);
		}

		// We need the complete table to do this.
		// Comment out for production code since the tolerance is variable.
		//for (int i = 1; i < size; i++)
		//{
		//	double value = ((long) i) << q;
		//	double log2 = data[i] + 1075;
		//	assert Math.abs((log2 - fastLog2(value)) / log2) < 1e-5 : String.format("[%d] log2(%g)  %g != %g  %g", i,
		//			value, log2, fastLog2(value), gdsc.core.utils.DoubleEquality.relativeError(log2, fastLog2(value)));
		//}
	}


	/**
	 * Gets the base.
	 *
	 * @return the base
	 */
	public double getBase()
	{
		return base;
	}

	/**
	 * Gets the scale to convert from log2 to logB (the base given in the constructor) by multiplication.
	 *
	 * @return the scale
	 */
	public double getScale()
	{
		return scale;
	}

	/**
	 * Gets the number of most significant bits to keep from the mantissa, i.e. the binary precision of the floating point number before
	 * computing the log.
	 *
	 * @return the number of most significant bits to keep
	 */
	public int getN()
	{
		return 52 - q;
	}
	
	/**
	 * Gets the number of least signification bits to ignore from the mantissa of a double (52-bits).
	 *
	 * @return the number of least signification bits to ignore
	 */
	public int getQ()
	{
		return q;
	}
	
	/**
	 * Calculate the logarithm to the base given in the constructor. Requires the argument be finite and positive.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN, then the result is incorrect.
	 * <li>If the argument is negative, then the result is incorrect (log(-x)).
	 * <li>If the argument is positive infinity, then the result is incorrect (Math.log(Double.MAX_VALUE)).
	 * <li>If the argument is positive zero or negative zero, then the result is negative infinity.
	 * </ul>
	 * 
	 * @param x
	 *            the argument. must be positive!
	 * @return log( x )
	 */
	public double fastLog2(double x)
	{
		final long bits = Double.doubleToRawLongBits(x);

		// Note the documentation from Double.longBitsToDouble(int):
		// int e = (int)((bits >>> 52) & 0x7ffL);
		// long m = (e == 0) ?
		//                 (bits & 0xfffffffffffffL) << 1 :
		//                 (bits & 0xfffffffffffffL) | 0x10000000000000L;

		final int e = (int) ((bits >>> 52) & 0x7ffL);
		// raw mantissa, conversion is done with the bit shift to reduce precision
		final long m = (bits & 0xfffffffffffffL);

		return (e == 0 ? data[(int) (m >>> q_minus_1)] : e + data[(int) ((m | 0x10000000000000L) >>> q)]);
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
	public double log2(double x)
	{
		// Basic implementation
		//if (x > 0)
		//	return (x == Double.POSITIVE_INFINITY) ? Double.POSITIVE_INFINITY : fastLog(x);
		//return (x == 0) ? Double.NEGATIVE_INFINITY : Double.NaN;

		// Re-implement to avoid double comparisons (which will be slower than int comparisons) 
		final long bits = Double.doubleToRawLongBits(x);

		// Note the documentation from Double.longBitsToDouble(int):
		// int e = (int)((bits >>> 52) & 0x7ffL);
		// long m = (e == 0) ?
		//                 (bits & 0xfffffffffffffL) << 1 :
		//                 (bits & 0xfffffffffffffL) | 0x10000000000000L;

		final boolean negative = (bits >>> 63) != 0;
		final int e = (int) ((bits >>> 52) & 0x7ffL);
		// raw mantissa, conversion is done with the bit shift to reduce precision
		final long m = (bits & 0xfffffffffffffL);

		if (e == 2047)
		{
			// All bits set is a special case
			if (m != 0)
				return Double.NaN;
			// +/- Infinity
			return (negative) ? Double.NaN : Double.POSITIVE_INFINITY;
		}

		if (negative)
			// Only -0 is allowed
			return (e == 0 && m == 0) ? Double.NEGATIVE_INFINITY : Double.NaN;

		return (e == 0 ? data[(int) (m >>> q_minus_1)] : e + data[(int) ((m | 0x10000000000000L) >>> q)]);
	}

	/**
	 * Calculate the logarithm to the base given in the constructor. Requires the argument be finite and positive.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN, then the result is incorrect.
	 * <li>If the argument is negative, then the result is incorrect (log(-x)).
	 * <li>If the argument is positive infinity, then the result is incorrect (Math.log(Double.MAX_VALUE)).
	 * <li>If the argument is positive zero or negative zero, then the result is negative infinity.
	 * </ul>
	 * 
	 * @param x
	 *            the argument. must be positive!
	 * @return log( x )
	 */
	public double fastLog(double x)
	{
		return fastLog2(x) * scale;
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
		return log2(x) * scale;
	}
}
