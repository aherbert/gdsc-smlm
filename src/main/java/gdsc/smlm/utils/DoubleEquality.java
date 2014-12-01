package gdsc.smlm.utils;

import org.apache.commons.math3.util.FastMath;

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
 * Provides equality functions for double doubleing polong numbers
 * <p>
 * Adapted from http://www.cygnus-software.com/papers/comparingdoubles/comparingdoubles.htm
 */
public class DoubleEquality
{
	private static final double RELATIVE_ERROR = 1e-2f;
	private static final double ABSOLUTE_ERROR = 1e-10f;
	private static final long SIGNIFICANT_DIGITS = 3;
	
	private double maxRelativeError;
	private double maxAbsoluteError;
	private long maxUlps;

	/**
	 * Default constructor
	 */
	public DoubleEquality()
	{
		init(RELATIVE_ERROR, ABSOLUTE_ERROR, SIGNIFICANT_DIGITS);
	}

	/**
	 * Override constructor
	 * 
	 * @param maxRelativeError
	 * @param maxAbsoluteError
	 */
	public DoubleEquality(double maxRelativeError, double maxAbsoluteError)
	{
		init(maxRelativeError, maxAbsoluteError, SIGNIFICANT_DIGITS);
	}

	/**
	 * Override constructor
	 * 
	 * @param maxRelativeError
	 * @param maxAbsoluteError
	 * @param significantDigits
	 */
	public DoubleEquality(double maxRelativeError, double maxAbsoluteError, long significantDigits)
	{
		init(maxRelativeError, maxAbsoluteError, significantDigits);
	}

	/**
	 * Override constructor
	 * 
	 * @param significantDigits
	 * @param maxAbsoluteError
	 */
	public DoubleEquality(long significantDigits, double maxAbsoluteError)
	{
		init(RELATIVE_ERROR, maxAbsoluteError, significantDigits);
	}
	
	private void init(double maxRelativeError, double maxAbsoluteError, long significantDigits)
	{
		this.maxRelativeError = maxRelativeError;
		this.maxAbsoluteError = maxAbsoluteError;
		setSignificantDigits(significantDigits);
	}

	/**
	 * Compares two doubles are within the configured errors.
	 * 
	 * @param A
	 * @param B
	 * @return True if equal
	 */
	public boolean almostEqualRelativeOrAbsolute(double A, double B)
	{
		return almostEqualRelativeOrAbsolute(A, B, maxRelativeError, maxAbsoluteError);
	}

	/**
	 * Compares two doubles are within the configured number of bits variation using long comparisons.
	 * 
	 * @param A
	 * @param B
	 * @return True if equal
	 */
	public boolean almostEqualComplement(double A, double B)
	{
		return almostEqualComplement(A, B, maxUlps, maxAbsoluteError);
	}

	/**
	 * Compares two doubles within the configured number of bits variation using long comparisons.
	 * 
	 * @param A
	 * @param B
	 * @return -1, 0 or 1
	 */
	public long compareComplement(double A, double B)
	{
		return compareComplement(A, B, maxUlps);
	}

	/**
	 * Compares two double arrays are within the configured errors.
	 * 
	 * @param A
	 * @param B
	 * @return True if equal
	 */
	public boolean almostEqualRelativeOrAbsolute(double[] A, double[] B)
	{
		for (int i = 0; i < A.length; i++)
			if (!almostEqualRelativeOrAbsolute(A[i], B[i], maxRelativeError, maxAbsoluteError))
				return false;
		return true;
	}

	/**
	 * Compares two double arrays are within the configured number of bits variation using long comparisons.
	 * 
	 * @param A
	 * @param B
	 * @return True if equal
	 */
	public boolean almostEqualComplement(double[] A, double[] B)
	{
		for (int i = 0; i < A.length; i++)
			if (!almostEqualComplement(A[i], B[i], maxUlps, maxAbsoluteError))
				return false;
		return true;
	}
	
	/**
	 * Compares two doubles are within the specified errors.
	 * 
	 * @param A
	 * @param B
	 * @param maxRelativeError
	 *            The relative error allowed between the numbers
	 * @param maxAbsoluteError
	 *            The absolute error allowed between the numbers. Should be a small number (e.g. 1e-10)
	 * @return True if equal
	 */
	public static boolean almostEqualRelativeOrAbsolute(double A, double B, double maxRelativeError, double maxAbsoluteError)
	{
		// Check the two numbers are within an absolute distance.
		final double difference = Math.abs(A - B);
		if (difference <= maxAbsoluteError)
			return true;
		final double size = FastMath.max(Math.abs(A), Math.abs(B));
		if (difference <= size * maxRelativeError)
			return true;
		return false;
	}

	/**
	 * Compute the relative error between two doubles.
	 * 
	 * @param A
	 * @param B
	 * @return The relative error
	 */
	public static double relativeError(double A, double B)
	{
		if (Math.abs(B) > Math.abs(A))
			return Math.abs((A - B) / B);
		else
			return Math.abs((A - B) / A);
	}

	/**
	 * Compares two doubles are within the specified number of bits variation using long comparisons.
	 * 
	 * @param A
	 * @param B
	 * @param maxUlps
	 *            How many representable doubles we are willing to accept between A and B
	 * @param maxAbsoluteError
	 *            The absolute error allowed between the numbers. Should be a small number (e.g. 1e-10)
	 * @return True if equal
	 */
	public static boolean almostEqualComplement(double A, double B, long maxUlps, double maxAbsoluteError)
	{
		// Make sure maxUlps is non-negative and small enough that the
		// default NAN won't compare as equal to anything.
		//assert (maxUlps > 0 && maxUlps < 4 * 1024 * 1024);
		
		if (Math.abs(A - B) < maxAbsoluteError)
			return true;
		if (complement(A, B) <= maxUlps)
			return true;
		return false;
	}

	/**
	 * Compares two doubles within the specified number of bits variation using long comparisons.
	 * 
	 * @param A
	 * @param B
	 * @param maxUlps
	 *            How many representable doubles we are willing to accept between A and B
	 * @return -1, 0 or 1
	 */
	public static int compareComplement(double A, double B, long maxUlps)
	{
		long c = signedComplement(A, B);
		if (c < -maxUlps)
			return -1;
		if (c > maxUlps)
			return 1;
		return 0;
	}

	/**
	 * @param maxRelativeError
	 *            the maxRelativeError to set
	 */
	public void setMaxRelativeError(double maxRelativeError)
	{
		this.maxRelativeError = maxRelativeError;
	}

	/**
	 * @return the maxRelativeError
	 */
	public double getMaxRelativeError()
	{
		return maxRelativeError;
	}

	/**
	 * @param maxAbsoluteError
	 *            the maxAbsoluteError to set
	 */
	public void setMaxAbsoluteError(double maxAbsoluteError)
	{
		this.maxAbsoluteError = maxAbsoluteError;
	}

	/**
	 * @return the maxAbsoluteError
	 */
	public double getMaxAbsoluteError()
	{
		return maxAbsoluteError;
	}

	/**
	 * @param maxUlps
	 *            the maximum error in terms of Units in the Last Place
	 */
	public void setMaxUlps(long maxUlps)
	{
		this.maxUlps = maxUlps;
	}

	/**
	 * @return the maximum error in terms of Units in the Last Place
	 */
	public long getMaxUlps()
	{
		return maxUlps;
	}
	
	/**
	 * Set the maximum error in terms of Units in the Last Place using the number of decimal significant digits
	 * 
	 * @param significantDigits The number of significant digits for comparisons
	 */
	public void setSignificantDigits(long significantDigits)
	{
		this.maxUlps = getUlps(significantDigits);
	}
	
	// The following methods are different between the FloatEquality and DoubleEquality class
	
	/**
	 * Compute the number of representable doubles until a difference in significant digits
	 * <p>
	 * The number of doubles are computed between Math.power(10, sig) and 1 + Math.power(10, sig)
	 * 
	 * @param significantDigits
	 *            The significant digits
	 * @return The number of representable doubles (Units in the Last Place)
	 */
	public static long getUlps(long significantDigits)
	{
		long value1 = (long)Math.pow(10.0, significantDigits-1);
		long value2 = value1 + 1;
		long ulps = Double.doubleToRawLongBits((double)value2) - Double.doubleToRawLongBits((double)value1);
		return (ulps < 0) ? 0 : ulps;
	}

	/**
	 * Compute the number of bits variation using long comparisons.
	 * 
	 * @param A
	 * @param B
	 *            
	 * @return How many representable doubles we are between A and B
	 */
	public static long complement(double A, double B)
	{
		long aInt = Double.doubleToRawLongBits(A);
		// Make aInt lexicographically ordered as a twos-complement long
		if (aInt < 0)
			aInt = 0x8000000000000000L - aInt;
		// Make bInt lexicographically ordered as a twos-complement long
		long bInt = Double.doubleToRawLongBits(B);
		if (bInt < 0)
			bInt = 0x8000000000000000L - bInt;
		return Math.abs(aInt - bInt);
	}

	/**
	 * Compute the number of bits variation using long comparisons.
	 * 
	 * @param A
	 * @param B
	 *            
	 * @return How many representable doubles we are between A and B
	 */
	public static long signedComplement(double A, double B)
	{
		long aInt = Double.doubleToRawLongBits(A);
		// Make aInt lexicographically ordered as a twos-complement long
		if (aInt < 0)
			aInt = 0x8000000000000000L - aInt;
		// Make bInt lexicographically ordered as a twos-complement long
		long bInt = Double.doubleToRawLongBits(B);
		if (bInt < 0)
			bInt = 0x8000000000000000L - bInt;
		return aInt - bInt;
	}
}
