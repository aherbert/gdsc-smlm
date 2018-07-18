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
package uk.ac.sussex.gdsc.smlm.function;

/**
 * Implementation of the ICSILog algorithm
 * as described in O. Vinyals, G. Friedland, N. Mirghafori
 * "Revisiting a basic function on current CPUs: A fast logarithm implementation
 * with adjustable accuracy" (2007).
 * <p>
 * This class is based on the original algorithm description.
 * <p>
 * The algorithm has been changed to detects when the unbiased exponent is zero and maintains the precision.
 * <p>
 * Using look-up table the relative error ((fastLog(x)-Math.log(x))/Math.log(x)) is large (e>>1) when the input value x
 * is close to 1. So the algorithm detects values close to 1 and uses Math.log instead.
 * <p>
 * This is a copy of TurboLog but implements rounding on the mantissa. This allows this class to achieve the same error
 * as TurboLog(n) using (n-1), i.e. half the table size.
 *
 * @see <a href=
 *      "http://www.icsi.berkeley.edu/pubs/techreports/TR-07-002.pdf">http://www.icsi.berkeley.edu/pubs/techreports/TR-
 *      07-002.pdf</a>
 */
public class TurboLog2 extends TurboLog
{
	/** The number of bits to remove from a float mantissa. */
	private final int q;
	/** The number of bits to remove from a double mantissa. */
	private final int qd;
	/**
	 * The table of the log value of the floating point mantissa (a binary number 1.0000... to 1.1111....), depending on
	 * the precision.
	 */
	private final float[] logMantissa;

	/** The number to add to a float mantissa for rounding. */
	private final int roundF;
	/** The number to add to a double mantissa for rounding. */
	private final long roundD;

	/**
	 * Create a new natural logarithm calculation instance. This will
	 * hold the pre-calculated log values for base E
	 * and the default table size.
	 */
	public TurboLog2()
	{
		// Since rounding doubles the precision we can reduce the default table size
		this(N - 1);
	}

	/**
	 * Create a new natural logarithm calculation instance. This will
	 * hold the pre-calculated log values for base E
	 * and a table size depending on a given mantissa precision.
	 *
	 * @param n
	 *            The number of bits to keep from the mantissa.
	 *            Table storage = 2^n * 4 bytes, e.g. 32Kb for n=13.
	 */
	public TurboLog2(int n)
	{
		// Store log value of a range of floating point numbers using a limited
		// precision mantissa (m). The purpose of this code is to enumerate all
		// possible mantissas of a float with limited precision (23-q). Note the
		// mantissa represents the digits of a binary number after the binary-point: .10101010101.
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
		// expression s x m x 2^(e-150):
		// e-127 is the unbiased exponent. 23 is the mantissa precision
		// = s x m x 2^(e-127-23)

		// E.g. For a precision of n=(23-q)=6
		// We enumerate:
		// (1.000000 to 1.111111)

		// The mantissa is incremented using an integer representation to allow
		// exact enumeration. This is then converted to a float for the call to
		// log(double).

		q = 23 - n;
		qd = 52 - n;
		int x = 0x3F800000; // Set the exponent to 0 so the float value=1.0
		//assert Float.intBitsToFloat(x) == 1.0f : "value is not 1.0f";
		final int inc = 1 << q; // Amount to increase the mantissa

		final int size = 1 << n;
		// Add an extra value in case the final mantissa is rounded up
		logMantissa = new float[size + 1];
		for (int i = 0; i < size; i++)
		{
			final float value = Float.intBitsToFloat(x);
			final float logv = (float) Math.log(value);
			logMantissa[i] = logv;
			x += inc;

			//assert logv == fastLog(value) : String.format("[%d] data[i](%g)  %g != %g  %g", i, value, logv,
			//		fastLog2(value), uk.ac.sussex.gdsc.core.utils.FloatEquality.relativeError(logv, fastLog2(value)));
		}
		// For rounding the final mantissa up
		logMantissa[size] = logMantissa[size - 1];

		// To round a mantissa add a number corresponding to the first insignificant digit.
		// E.g. for n=10, q=13, number to add is 1 << (q-1)
		//   1101010101.............
		// +           1000000000000
		if (q != 0)
		{
			roundF = 1 << (q - 1);
			roundD = 1L << (qd - 1);
		}
		else
		{
			roundF = 0;
			roundD = 0;
		}
	}

	@Override
	public int getN()
	{
		return 23 - q;
	}

	@Override
	public float log(float x)
	{
		final int bits = Float.floatToRawIntBits(x);
		final int e = (bits >>> 23) & 0xff;
		final int m = (bits & 0x7fffff);

		// Edge case for NaN and +/- Infinity
		if (e == 255)
		{
			if (m != 0)
				return Float.NaN;
			return ((bits & 0x80000000) != 0) ? Float.NaN : Float.POSITIVE_INFINITY;
		}

		// Edge case for negatives
		if ((bits & 0x80000000) != 0)
			// Only allow -0
			return (e == 0 && m == 0) ? Float.NEGATIVE_INFINITY : Float.NaN;

		// Note the documentation from Float.intBitsToFloat(int):
		// int s = ((bits >> 31) == 0) ? 1 : -1;
		// int e = ((bits >> 23) & 0xff); // Unsigned exponent
		// int m = (e == 0) ?
		//                 (bits & 0x7fffff) << 1 :
		//                 (bits & 0x7fffff) | 0x800000;
		//
		// Then the floating-point result equals the value of the mathematical
		// expression s x m x 2^(e-150):
		// e-127 is the unbiased exponent. 23 is the mantissa precision
		// = s x m x 2^(e-127-23)
		//
		// Here we have m as an index to the log of the mantissa including
		// the binary point. So we just need to compute
		// log(m x 2^(e-127))
		// = log(m) + log(2^(e-127))
		// = log(m) + (e-127) * log(2)

		// Check the exponent
		if (e == 0)
			return (m == 0) ? Float.NEGATIVE_INFINITY : computeSubnormal(m << 1);

		// When the value is close to 1 then the relative error can be very large
		if ((e == 126 && m >= lowerBoundMantissaF) || (e == 127 && m <= upperBoundMantissaF))
			return (float) Math.log(x);

		// Round the mantissa
		return logMantissa[(m + roundF) >>> q] + logExpF[e];
	}

	/**
	 * Calculate the natural logarithm. Requires the argument be finite and positive.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN, then the result is incorrect (>fastLog(Float.MAX_VALUE)).
	 * <li>If the argument is negative, then the result is incorrect (fastLog(-x)).
	 * <li>If the argument is positive infinity, then the result is incorrect (>fastLog(Float.MAX_VALUE)).
	 * <li>If the argument is positive zero or negative zero, then the result is negative infinity.
	 * </ul>
	 *
	 * @param x
	 *            the argument (must be strictly positive)
	 * @return log(x)
	 */
	@Override
	public float fastLog(float x)
	{
		// As above but no checks for NaN or infinity
		final int bits = Float.floatToRawIntBits(x);
		final int e = ((bits >>> 23) & 0xff);
		final int m = (bits & 0x7fffff);
		if (e == 0)
			return (m == 0) ? Float.NEGATIVE_INFINITY : computeSubnormal(m << 1);
		if ((e == 126 && m >= lowerBoundMantissaF) || (e == 127 && m <= upperBoundMantissaF))
			return (float) Math.log(x);
		return logMantissa[(m + roundF) >>> q] + logExpF[e];
	}

	@Override
	public float log(double x)
	{
		final long bits = Double.doubleToRawLongBits(x);
		final int e = (int) ((bits >>> 52) & 0x7ffL);
		final long m = (bits & 0xfffffffffffffL);

		// Edge case for NaN and +/- Infinity
		if (e == 2047)
		{
			if (m != 0L)
				return Float.NaN;
			return ((bits & 0x8000000000000000L) != 0L) ? Float.NaN : Float.POSITIVE_INFINITY;
		}

		// Edge case for negatives
		if ((bits & 0x8000000000000000L) != 0L)
			// Only allow -0
			return (e == 0 && m == 0L) ? Float.NEGATIVE_INFINITY : Float.NaN;

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
		//
		// Here we have m as an index to the log of the mantissa including
		// the binary point. So we just need to compute
		// log(m x 2^(e-1023))
		// = log(m) + log(2^(e-1023))
		// = log(m) + (e-1023) * log(2)

		// Check the exponent
		if (e == 0)
			return (m == 0L) ? Float.NEGATIVE_INFINITY : computeSubnormalF(m << 1);

		// When the value is close to 1 then the relative error can be very large
		if ((e == 1022 && m >= lowerBoundMantissa) || (e == 1023 && m <= upperBoundMantissa))
			return (float) Math.log(x);

		// Round the mantissa
		return logMantissa[(int) ((m + roundD) >>> qd)] + logExpD[e];
	}

	/**
	 * Calculate the natural logarithm. Requires the argument be finite and positive.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN, then the result is incorrect (>fastLog(Double.MAX_VALUE)).
	 * <li>If the argument is negative, then the result is incorrect (fastLog(-x)).
	 * <li>If the argument is positive infinity, then the result is incorrect (>fastLog(Double.MAX_VALUE)).
	 * <li>If the argument is positive zero or negative zero, then the result is negative infinity.
	 * </ul>
	 *
	 * @param x
	 *            the argument (must be strictly positive)
	 * @return log(x)
	 */
	@Override
	public float fastLog(double x)
	{
		// As above but no checks for NaN or infinity
		final long bits = Double.doubleToRawLongBits(x);
		final int e = (int) ((bits >>> 52) & 0x7ffL);
		final long m = (bits & 0xfffffffffffffL);
		if (e == 0)
			return (m == 0L) ? Float.NEGATIVE_INFINITY : computeSubnormalF(m << 1);
		if ((e == 1022 && m >= lowerBoundMantissa) || (e == 1023 && m <= upperBoundMantissa))
			return (float) Math.log(x);
		// Round the mantissa
		return logMantissa[(int) ((m + roundD) >>> qd)] + logExpD[e];
	}

	@Override
	public double logD(double x)
	{
		final long bits = Double.doubleToRawLongBits(x);
		final int e = (int) ((bits >>> 52) & 0x7ffL);
		final long m = (bits & 0xfffffffffffffL);

		// Edge case for NaN and +/- Infinity
		if (e == 2047)
		{
			if (m != 0L)
				return Double.NaN;
			return ((bits & 0x8000000000000000L) != 0L) ? Double.NaN : Double.POSITIVE_INFINITY;
		}

		// Edge case for negatives
		if ((bits & 0x8000000000000000L) != 0L)
			// Only allow -0
			return (e == 0 && m == 0L) ? Double.NEGATIVE_INFINITY : Double.NaN;

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
		//
		// Here we have m as an index to the log of the mantissa including
		// the binary point. So we just need to compute
		// log(m x 2^(e-1023))
		// = log(m) + log(2^(e-1023))
		// = log(m) + (e-1023) * log(2)

		// Check the exponent
		if (e == 0)
			return (m == 0L) ? Double.NEGATIVE_INFINITY : computeSubnormal(m << 1);

		// When the value is close to 1 then the relative error can be very large
		if ((e == 1022 && m >= lowerBoundMantissa) || (e == 1023 && m <= upperBoundMantissa))
			return Math.log(x);

		// Round the mantissa
		//return logMantissa[(int) ((m+roundD) >>> qd)] + logExpD[e];
		return logMantissa[(int) ((m + roundD) >>> qd)] + (e - 1023) * LN2;
	}

	/**
	 * Calculate the natural logarithm. Requires the argument be finite and positive.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN, then the result is incorrect (>fastLog(Double.MAX_VALUE)).
	 * <li>If the argument is negative, then the result is incorrect (fastLog(-x)).
	 * <li>If the argument is positive infinity, then the result is incorrect (>fastLog(Double.MAX_VALUE)).
	 * <li>If the argument is positive zero or negative zero, then the result is negative infinity.
	 * </ul>
	 *
	 * @param x
	 *            the argument (must be strictly positive)
	 * @return log(x)
	 */
	@Override
	public double fastLogD(double x)
	{
		// As above but no checks for NaN or infinity
		final long bits = Double.doubleToRawLongBits(x);
		final int e = (int) ((bits >>> 52) & 0x7ffL);
		final long m = (bits & 0xfffffffffffffL);
		if (e == 0)
			return (m == 0L) ? Double.NEGATIVE_INFINITY : computeSubnormal(m << 1);
		if ((e == 1022 && m >= lowerBoundMantissa) || (e == 1023 && m <= upperBoundMantissa))
			return Math.log(x);
		//return logMantissa[(int) ((m+roundD) >>> qd)] + logExpD[e];
		return logMantissa[(int) ((m + roundD) >>> qd)] + (e - 1023) * LN2;
	}
}
