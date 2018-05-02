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
public class FFastLog extends FastLog
{
	/** The base. */
	private final double base;
	/** The number of bits to remove from a float mantissa. */
	private final int q;
	/** (q-1). */
	private final int q_minus_1;
	/** The number of bits to remove from a double mantissa. */
	private final int qd;
	/** (qd-1). */
	private final int qd_minus_1;
	/**
	 * The table of the log2 value of binary number 1.0000... to 1.1111...., depending on the precision.
	 * The table has had the float bias (127) and mantissa size (23) pre-subtracted (i.e. -150 in total).
	 */
	private final float[] data;
	/** The scale used to convert the log2 to the logB (using the base). */
	private final float scale;

	/**
	 * Create a new natural logarithm calculation instance.
	 */
	public FFastLog()
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
	public FFastLog(int n)
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
	public FFastLog(double base, int n)
	{
		if (n < 0 || n > 23)
			throw new IllegalArgumentException("N must be in the range 0<=n<=23");
		scale = (float) getScale(base);
		this.base = base;

		final int size = 1 << (n + 1);

		q = 23 - n;
		q_minus_1 = q - 1;
		qd = 52 - n;
		qd_minus_1 = qd - 1;
		data = new float[size];

		for (int i = 0; i < size; i++)
		{
			// Store log2 value of a range of floating point numbers using a limited
			// precision mantissa (m). The purpose of this code is to enumerate all 
			// possible mantissas of a float with limited precision (23-q). Note the
			// 24 for a 23-bit mantissa comes from the fact that the mantissa 
			// represents the digits of a binary number after the binary-point: .10101010101
			// It is assumed that the digit before the point is a 1 if the exponent
			// is non-zero. Otherwise the binary point is moved to the right of the first
			// digit (i.e. a bit shift left).

			// See Float.intBitsToFloat(int):
			// int s = ((bits >> 31) == 0) ? 1 : -1;
			// int e = ((bits >> 23) & 0xff); // Unsigned exponent
			// int m = (e == 0) ?
			//                 (bits & 0x7fffff) << 1 :
			//                 (bits & 0x7fffff) | 0x800000;
			//
			// Then the floating-point result equals the value of the mathematical
			// expression s x m x 2^(e-150).

			// For a precision of n=(23-q)=6
			// We enumerate:
			// ( .000000 to  .111111) * 2^q (i.e. q additional zeros) 
			// (1.000000 to 1.111111) * 2^q

			// The bit shift is performed on integer data to construct the desired mantissa
			// which is then converted to a double for the call to exactLog2(double).

			// int m = i << q
			// log2(m x 2^(e-150)) == log2(m) + log2(2^e-150) = log2(m) +(e-150)
			// We subtract the -150 here so that the log2(float) can be reconstructed
			// from the table of log2(m) + e.

			data[i] = (float) (exactLog2(i << q) - 150);
		}

		// We need the complete table to do this
		// Comment out for production code since the tolerance is variable.
		//for (int i = 1; i < size; i++)
		//{
		//	float value = i << q;
		//	float log2 = data[i] + 150f;
		//	assert Math.abs((log2 - fastLog2(value)) / log2) < 1e-6f : String.format("[%d] log2(%g)  %g != %g  %g", i,
		//			value, log2, fastLog2(value), gdsc.core.utils.FloatEquality.relativeError(log2, fastLog2(value)));
		//}
	}

	@Override
	public double getBase()
	{
		return base;
	}

	@Override
	public double getScale()
	{
		return scale;
	}

	@Override
	public int getN()
	{
		return 23 - q;
	}

	@Override
	public float log2(float x)
	{
		final int bits = Float.floatToRawIntBits(x);

		// Note the documentation from Float.intBitsToFloat(int):
		// int s = ((bits >> 31) == 0) ? 1 : -1;
		// int e = ((bits >> 23) & 0xff);
		// int m = (e == 0) ?
		//                 (bits & 0x7fffff) << 1 :
		//                 (bits & 0x7fffff) | 0x800000;
		// Then the floating-point result equals the value of the mathematical
		// expression s x m x 2^(e-150):
		// e-127 is the unbiased exponent. 23 is the mantissa precision
		// = s x m x 2^(e-127-23) 

		final int e = (bits >> 23) & 0xff;
		// raw mantissa, conversion is done with the bit shift to reduce precision
		final int m = (bits & 0x7fffff);

		if (e == 255)
		{
			// All bits set is a special case
			if (m != 0)
				return Float.NaN;
			// +/- Infinity
			return ((bits >> 31) != 0) ? Float.NaN : Float.POSITIVE_INFINITY;
		}

		if ((bits >> 31) != 0)
		{
			// Only -0 is allowed
			return (e == 0 && m == 0) ? Float.NEGATIVE_INFINITY : Float.NaN;
		}

		return (e == 0 ? data[m >>> q_minus_1] : e + data[((m | 0x00800000) >>> q)]);
	}

	/**
	 * Calculate the logarithm using base 2. Requires the argument be finite and positive.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN, then the result is incorrect (>fastLog2(Float.MAX_VALUE)).
	 * <li>If the argument is negative, then the result is incorrect (fastLog2(-x)).
	 * <li>If the argument is positive infinity, then the result is incorrect (fastLog2(Float.MAX_VALUE)).
	 * <li>If the argument is positive zero or negative zero, then the result is negative infinity.
	 * </ul>
	 * 
	 * @param x
	 *            the argument (must be strictly positive)
	 * @return log2( x )
	 */
	@Override
	public float fastLog2(float x)
	{
		final int bits = Float.floatToRawIntBits(x);
		final int e = (bits >> 23) & 0xff;
		final int m = (bits & 0x7fffff);
		return (e == 0 ? data[m >>> q_minus_1] : e + data[((m | 0x00800000) >>> q)]);
	}

	@Override
	public float log(float x)
	{
		// Re-implement to avoid float comparisons (which will be slower than int comparisons) 
		final int bits = Float.floatToRawIntBits(x);
		final int e = (bits >> 23) & 0xff;
		final int m = (bits & 0x7fffff);
		if (e == 255)
		{
			if (m != 0)
				return Float.NaN;
			return ((bits >> 31) != 0) ? Float.NaN : Float.POSITIVE_INFINITY;
		}
		if ((bits >> 31) != 0)
		{
			return (e == 0 && m == 0) ? Float.NEGATIVE_INFINITY : Float.NaN;
		}
		return (e == 0 ? data[m >>> q_minus_1] : e + data[((m | 0x00800000) >>> q)]) * scale;
	}

	/**
	 * Calculate the logarithm to the base given in the constructor. Requires the argument be finite and positive.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN, then the result is incorrect (>fastLog(Float.MAX_VALUE)).
	 * <li>If the argument is negative, then the result is incorrect (fastLog(-x)).
	 * <li>If the argument is positive infinity, then the result is incorrect (fastLog(Float.MAX_VALUE)).
	 * <li>If the argument is positive zero or negative zero, then the result is negative infinity.
	 * </ul>
	 * 
	 * @param x
	 *            the argument (must be strictly positive)
	 * @return log( x )
	 */
	@Override
	public float fastLog(float x)
	{
		final int bits = Float.floatToRawIntBits(x);
		final int e = (bits >> 23) & 0xff;
		final int m = (bits & 0x7fffff);
		return (e == 0 ? data[m >>> q_minus_1] : e + data[((m | 0x00800000) >>> q)]) * scale;
	}

	@Override
	public float log2(double x)
	{
		final long bits = Double.doubleToRawLongBits(x);

		// Note the documentation from Double.longBitsToDouble(int):
		// int s = ((bits >> 63) == 0) ? 1 : -1;
		// int e = (int)((bits >>> 52) & 0x7ffL);
		// long m = (e == 0) ?
		//                 (bits & 0xfffffffffffffL) << 1 :
		//                 (bits & 0xfffffffffffffL) | 0x10000000000000L;
		// Then the floating-point result equals the value of the mathematical
		// expression s x m x 2^(e-1075):
		// e-1023 is the unbiased exponent. 52 is the mantissa precision
		// = s x m x 2^(e-1023-52) 

		// Get the biased exponent
		final int e = (int) ((bits >>> 52) & 0x7ffL);
		// Mantissa
		final long m = (bits & 0xfffffffffffffL);

		if (e == 2047)
		{
			// All bits set is a special case
			if (m != 0)
				return Float.NaN;
			// +/- Infinity
			return ((bits >> 63) != 0L) ? Float.NaN : Float.POSITIVE_INFINITY;
		}

		// Check for negatives
		if ((bits >> 63) != 0L)
		{
			// Only -0 is allowed
			return (e == 0 && m == 0L) ? Float.NEGATIVE_INFINITY : Float.NaN;
		}

		// We must subtract -1075 from the exponent.
		// However -150 has been pre-subtracted in the table.
		// and the mantissa has 29 more digits of significance.
		// So take away 1075-150-29 = 896.
		return (e == 0 ? data[(int) (m >>> qd_minus_1)] - 896
				: e - 896 + data[(int) ((m | 0x10000000000000L) >>> qd)]);
	}

	/**
	 * Calculate the logarithm using base 2. Requires the argument be finite and positive.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN, then the result is incorrect (>fastLog2(Float.MAX_VALUE)).
	 * <li>If the argument is negative, then the result is incorrect (fastLog2(-x)).
	 * <li>If the argument is positive infinity, then the result is incorrect (fastLog2(Float.MAX_VALUE)).
	 * <li>If the argument is positive zero or negative zero, then the result is negative infinity.
	 * </ul>
	 * 
	 * @param x
	 *            the argument (must be strictly positive)
	 * @return log( x )
	 */
	@Override
	public float fastLog2(double x)
	{
		final long bits = Double.doubleToRawLongBits(x);
		final int e = (int) ((bits >>> 52) & 0x7ffL);
		final long m = (bits & 0xfffffffffffffL);
		return (e == 0 ? data[(int) (m >>> qd_minus_1)] - 896 : e - 896 + data[(int) ((m | 0x10000000000000L) >>> qd)]);
	}

	@Override
	public float log(double x)
	{
		final long bits = Double.doubleToRawLongBits(x);

		// Note the documentation from Double.longBitsToDouble(int):
		// int s = ((bits >> 63) == 0) ? 1 : -1;
		// int e = (int)((bits >>> 52) & 0x7ffL);
		// long m = (e == 0) ?
		//                 (bits & 0xfffffffffffffL) << 1 :
		//                 (bits & 0xfffffffffffffL) | 0x10000000000000L;
		// Then the floating-point result equals the value of the mathematical
		// expression s x m x 2^(e-1075):
		// e-1023 is the unbiased exponent. 52 is the mantissa precision
		// = s x m x 2^(e-1023-52) 

		// Get the biased exponent
		final int e = (int) ((bits >>> 52) & 0x7ffL);
		// Mantissa
		final long m = (bits & 0xfffffffffffffL);

		if (e == 2047)
		{
			// All bits set is a special case
			if (m != 0)
				return Float.NaN;
			// +/- Infinity
			return ((bits >> 63) != 0L) ? Float.NaN : Float.POSITIVE_INFINITY;
		}

		// Check for negatives
		if ((bits >> 63) != 0L)
		{
			// Only -0 is allowed
			return (e == 0 && m == 0L) ? Float.NEGATIVE_INFINITY : Float.NaN;
		}

		// We must subtract -1075 from the exponent.
		// However -150 has been pre-subtracted in the table.
		// and the mantissa has 29 more digits of significance.
		// So take away 1075-150-29 = 896.
		return (e == 0 ? data[(int) (m >>> qd_minus_1)] - 896
				: e - 896 + data[(int) ((m | 0x10000000000000L) >>> qd)]) * scale;
	}

	/**
	 * Calculate the logarithm to the base given in the constructor. Requires the argument be finite and positive.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN, then the result is incorrect (>fastLog(Float.MAX_VALUE)).
	 * <li>If the argument is negative, then the result is incorrect (fastLog(-x)).
	 * <li>If the argument is positive infinity, then the result is incorrect (fastLog(Float.MAX_VALUE)).
	 * <li>If the argument is positive zero or negative zero, then the result is negative infinity.
	 * </ul>
	 * 
	 * @param x
	 *            the argument (must be strictly positive)
	 * @return log( x )
	 */
	@Override
	public float fastLog(double x)
	{
		final long bits = Double.doubleToRawLongBits(x);
		final int e = (int) ((bits >>> 52) & 0x7ffL);
		final long m = (bits & 0xfffffffffffffL);
		return (e == 0 ? data[(int) (m >>> qd_minus_1)] - 896
				: e - 896 + data[(int) ((m | 0x10000000000000L) >>> qd)]) * scale;
	}
}
