package gdsc.smlm.function;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
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
 * This class is based on the original algorithm description.
 * <p>
 * When the unbiased exponent is zero a conversion from float to double is made to preserve the precision. If already a
 * double then the full Math.log function is used.
 * <p>
 * The relative error ((fastLog(x)-Math.log(x))/Math.log(x)) is large (e~0.76) when the input value x is close to 1.
 *
 * @see <a href=
 *      "http://www.icsi.berkeley.edu/pubs/techreports/TR-07-002.pdf">http://www.icsi.berkeley.edu/pubs/techreports/TR-
 *      07-002.pdf</a>
 */
public class TurboLog extends FastLog
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
	/**
	 * The table of the log value of the unbaised float exponent (an integer from -127 to 128).
	 */
	private static final float[] logExpF;
	/**
	 * The table of the log value of the unbaised double exponent (an integer from -1023 to 1024).
	 */
	private static final float[] logExpD;

	static
	{
		// Note: the exponent is already in base 2. Just multiply by ln(2) to convert to base E
		logExpF = new float[256]; // 8-bit exponent
		for (int i = 0; i < logExpF.length; i++)
			logExpF[i] = (float) ((i - 127) * LN2);
		logExpD = new float[2048]; // 11-bit exponent
		for (int i = 0; i < logExpD.length; i++)
			logExpD[i] = (float) ((i - 1023) * LN2);
	}

	/**
	 * Create a new natural logarithm calculation instance. This will
	 * hold the pre-calculated log values for base E
	 * and the default table size.
	 *
	 * @param dataType
	 *            the data type
	 */
	public TurboLog()
	{
		this(N);
	}

	/**
	 * Create a new natural logarithm calculation instance. This will
	 * hold the pre-calculated log values for base E
	 * and a table size depending on a given mantissa precision.
	 *
	 * @param n
	 *            The number of bits to keep from the mantissa.
	 *            Table storage = 2^n * 4 bytes, e.g. 32Kb for n=13.
	 * @param dataType
	 *            the data type
	 */
	public TurboLog(int n)
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
		int inc = 1 << q; // Amount to increase the mantissa

		final int size = 1 << n;
		logMantissa = new float[size];
		for (int i = 0; i < size; i++)
		{
			float value = Float.intBitsToFloat(x);
			float logv = (float) Math.log(value);
			logMantissa[i] = logv;
			x += inc;

			//assert gdsc.core.utils.FloatEquality.almostEqualRelativeOrAbsolute(logv, fastLog2(value), 1e-6f, 1e-16f)
			assert logv == fastLog(value) : String.format("[%d] data[i](%g)  %g != %g  %g", i, value, logv,
					fastLog2(value), gdsc.core.utils.FloatEquality.relativeError(logv, fastLog2(value)));
		}
	}

	@Override
	public int getN()
	{
		return 23 - q;
	}

	@Override
	public double getScale()
	{
		return LN2;
	}

	@Override
	public double getBase()
	{
		return Math.E;
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
			return ((bits >>> 31) != 0) ? Float.NaN : Float.POSITIVE_INFINITY;
		}

		// Edge case for negatives
		if ((bits >>> 31) != 0)
		{
			// Only allow -0
			return (e == 0 && m == 0) ? Float.NEGATIVE_INFINITY : Float.NaN;
		}

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
		{
			// When the exponent is zero there in no assumed leading 1.
			// So the look-up table is invalid. In this case a cast to a double
			// should restore the leading 1 as the exponent has more precision.
			// The only case when this is not true is if m==0 and thus the 
			// value is zero.
			if (m == 0)
				return Float.NEGATIVE_INFINITY;
			//return fastLog((double) x);
			
			// TODO
			// See FastMath.log for how to normalise the sub-normal number			

			// Re-implement double version here as we assume that e will not be zero again.
			final long lbits = Double.doubleToLongBits(x);
			final int le = (int) ((lbits >>> 52) & 0x7ffL);
			final long lm = (lbits & 0xfffffffffffffL);
			return logMantissa[(int) (lm >>> qd)] + logExpD[le];
		}

		// TODO - fix this for the double version too
		
		// When the value is close to 1 then the relative error can be very large
		if ((e == 126 || e == 127) && x < 1.01f && x > 0.99f)
		{
			return (float) Math.log(x);
		}

		return logMantissa[m >>> q] + logExpF[e];
	}

	/**
	 * Calculate the natural logarithm. Requires the argument be finite and strictly positive.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN, then the result is incorrect (>fastLog(Float.MAX_VALUE)).
	 * <li>If the argument is negative, then the result is incorrect (fastLog(-x)).
	 * <li>If the argument is positive infinity, then the result is incorrect (>fastLog(Float.MAX_VALUE)).
	 * <li>If the argument is positive zero or negative zero, then the result is incorrect (fastLog((double)x)).
	 * </ul>
	 * 
	 * @param x
	 *            the argument (must be strictly positive)
	 * @return log(x)
	 */
	public float fastLog(float x)
	{
		final int bits = Float.floatToRawIntBits(x);
		final int e = ((bits >>> 23) & 0xff);
		if (e == 0)
		{
			// Maintain precision when there is no leading 1.
			// Assume that x is strictly positive and do not check for zero.
			//return fastLog((double) x);

			// Re-implement double version here as we assume that e will not be zero again.
			final long lbits = Double.doubleToLongBits(x);
			final int le = (int) ((lbits >>> 52) & 0x7ffL);
			final long lm = (lbits & 0xfffffffffffffL);
			return logMantissa[(int) (lm >>> qd)] + logExpD[le];
		}

		// When the value is close to 1 then the relative error can be very large
		if ((e == 126 || e == 127) && x < 1.01f && x > 0.99f)
		{
			return (float) Math.log(x);
		}

		final int m = (bits & 0x7fffff);
		return logMantissa[m >>> q] + logExpF[e];
	}

	@Override
	public float log(double x)
	{
		final long bits = Double.doubleToLongBits(x);
		final int e = (int) ((bits >>> 52) & 0x7ffL);
		final long m = (bits & 0xfffffffffffffL);

		// Edge case for NaN and +/- Infinity
		if (e == 2047)
		{
			if (m != 0L)
				return Float.NaN;
			return ((bits >>> 63) != 0L) ? Float.NaN : Float.POSITIVE_INFINITY;
		}

		// Edge case for negatives
		if ((bits >>> 63) != 0L)
		{
			// Only allow -0
			return (e == 0 && m == 0L) ? Float.NEGATIVE_INFINITY : Float.NaN;
		}

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
		{
			// When the exponent is zero there in no assumed leading 1.
			// So the look-up table is invalid. In this case resort to full precision.
			return (float) Math.log(x);
		}
		else
		{
			return logMantissa[(int) (m >>> qd)] + logExpD[e];
		}
	}

	/**
	 * Calculate the natural logarithm. Requires the argument be finite and positive.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN, then the result is incorrect (>fastLog(Double.MAX_VALUE)).
	 * <li>If the argument is negative, then the result is incorrect (fastLog(-x)).
	 * <li>If the argument is positive infinity, then the result is incorrect (>fastLog(Double.MAX_VALUE)).
	 * <li>If the argument is positive zero or negative zero, then the result is Double.NEGATIVE_INFINITY.
	 * </ul>
	 * 
	 * @param x
	 *            the argument (must be strictly positive)
	 * @return log(x)
	 */
	public float fastLog(double x)
	{
		final long bits = Double.doubleToLongBits(x);
		final int e = (int) ((bits >>> 52) & 0x7ffL);
		if (e == 0)
		{
			// Maintain precision when there is no leading 1.
			return (float) Math.log(x);
		}
		final long m = (bits & 0xfffffffffffffL);
		return logMantissa[(int) (m >>> qd)] + logExpD[e];
	}

	@Override
	public double logD(double x)
	{
		final long bits = Double.doubleToLongBits(x);
		final int e = (int) ((bits >>> 52) & 0x7ffL);
		final long m = (bits & 0xfffffffffffffL);

		// Edge case for NaN and +/- Infinity
		if (e == 2047)
		{
			if (m != 0L)
				return Double.NaN;
			return ((bits >>> 63) != 0L) ? Double.NaN : Double.POSITIVE_INFINITY;
		}

		// Edge case for negatives
		if ((bits >>> 63) != 0L)
		{
			// Only allow -0
			return (e == 0 && m == 0L) ? Double.NEGATIVE_INFINITY : Double.NaN;
		}

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
		{
			// When the exponent is zero there in no assumed leading 1.
			// So the look-up table is invalid. In this case resort to full precision.
			return Math.log(x);
		}
		else
		{
			return logMantissa[(int) (m >>> qd)] + logExpD[e];
		}
	}

	/**
	 * Calculate the natural logarithm. Requires the argument be finite and positive.
	 * <p>
	 * Special cases:
	 * <ul>
	 * <li>If the argument is NaN, then the result is incorrect (>fastLog(Double.MAX_VALUE)).
	 * <li>If the argument is negative, then the result is incorrect (fastLog(-x)).
	 * <li>If the argument is positive infinity, then the result is incorrect (>fastLog(Double.MAX_VALUE)).
	 * <li>If the argument is positive zero or negative zero, then the result is Double.NEGATIVE_INFINITY.
	 * </ul>
	 * 
	 * @param x
	 *            the argument (must be strictly positive)
	 * @return log(x)
	 */
	public double fastLogD(double x)
	{
		final long bits = Double.doubleToLongBits(x);
		final int e = (int) ((bits >>> 52) & 0x7ffL);
		if (e == 0)
		{
			// Maintain precision when there is no leading 1.
			return Math.log(x);
		}
		final long m = (bits & 0xfffffffffffffL);
		return logMantissa[(int) (m >>> qd)] + logExpD[e];
	}

	// We don't support other bases so do a simple conversion for log2 for the super-class method

	@Override
	public float log2(float x)
	{
		return log(x) / LN2F;
	}

	@Override
	public float fastLog2(float x)
	{
		return fastLog(x) / LN2F;
	}

	@Override
	public double log2D(double x)
	{
		return log(x) / LN2;
	}

	@Override
	public double fastLog2D(double x)
	{
		return fastLog(x) / LN2;
	}

	@Override
	public float log2(double x)
	{
		return log(x) / LN2F;
	}

	@Override
	public float fastLog2(double x)
	{
		return fastLog(x) / LN2F;
	}
}
