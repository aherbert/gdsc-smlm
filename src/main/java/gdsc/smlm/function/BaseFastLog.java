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
 * Base class for the ICSILog algorithm
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
abstract class BaseFastLog
{
	/**
	 * Natural logarithm of 2
	 */
	public static final double LN2 = Math.log(2);

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
}
