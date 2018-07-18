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

import java.util.Arrays;

import org.apache.commons.math3.special.Gamma;

/**
 * Compute the log of n!
 */
public class LogFactorial
{
	private static double[] table;

	/** All long-representable factorials */
	static final long[] FACTORIALS = new long[] { 1l, 1l, 2l, 6l, 24l, 120l, 720l, 5040l, 40320l, 362880l, 3628800l,
			39916800l, 479001600l, 6227020800l, 87178291200l, 1307674368000l, 20922789888000l, 355687428096000l,
			6402373705728000l, 121645100408832000l, 2432902008176640000l };

	static
	{
		table = new double[FACTORIALS.length];
		for (int k = 0; k < FACTORIALS.length; k++)
			table[k] = Math.log(FACTORIALS[k]);
	}

	private static final Object lock = new Object();

	/**
	 * Gets the maximum N that is tabulated.
	 *
	 * @return the max N
	 */
	public static int getTableMaxN()
	{
		return table.length - 1;
	}

	/**
	 * Increase the tabulated values up to a max n. Does nothing if already above the given n.
	 *
	 * @param n
	 *            the n
	 */
	public static void increaseTableMaxN(int n)
	{
		if (getTableMaxN() < n)
			synchronized (lock)
			{
				table = increaseSize(table, n);
			}
	}

	private static double[] increaseSize(double[] table, int n)
	{
		final double[] newTable = Arrays.copyOf(table, n + 1);

		// Using gamma for consistency with non-tabulated values
		int k = table.length - 1;
		while (k < n)
		{
			k++;
			newTable[k] = Gamma.logGamma(k + 1);
		}

		return newTable;
	}

	/**
	 * Reduces the tabulated values down to a max n. Does nothing if already below the given n.
	 *
	 * @param n
	 *            the new table max N
	 */
	public static void reduceTableMaxN(int n)
	{
		if (getTableMaxN() > n)
		{
			// Keep the representable factorials
			n = Math.max(n, FACTORIALS.length);

			synchronized (lock)
			{
				table = Arrays.copyOf(table, n + 1);
			}
		}
	}

	/**
	 * Compute the log of n!. Uses tabulated values or the gamma function if n is large.
	 *
	 * @param n
	 *            the n (must be positive)
	 * @return log(n!)
	 * @throws ArrayIndexOutOfBoundsException
	 *             if n is negative
	 */
	public static double logF(int n) throws ArrayIndexOutOfBoundsException
	{
		final double[] logF = table;
		if (n < logF.length)
			return logF[n];
		return Gamma.logGamma(n + 1);
	}

	/**
	 * Compute the log of k!. Uses the gamma function
	 *
	 * @param k
	 *            the k
	 * @return log(k!)
	 */
	public static double logF(double k)
	{
		if (k <= 1)
			return 0;
		return Gamma.logGamma(k + 1);
	}

	/////////////////////////////////////
	// Instance with pre-computed values
	/////////////////////////////////////

	private double[] objectTable;
	private final Object objectLock = new Object();

	/**
	 * Instantiates a new log factorial using the current static table.
	 */
	public LogFactorial()
	{
		// Copy the static values already present
		objectTable = table;
	}

	/**
	 * Instantiates a new log factorial using the given table size.
	 *
	 * @param n
	 *            the n
	 */
	public LogFactorial(int n)
	{
		if (n < 0)
			throw new IllegalArgumentException("N must be positive");

		// Copy the static values already present
		final double[] masterTable = LogFactorial.table;
		objectTable = Arrays.copyOf(masterTable, n + 1);

		// Fill in the rest (if required)
		int k = masterTable.length - 1;
		while (k < n)
		{
			k++;
			objectTable[k] = Gamma.logGamma(k + 1);
		}
	}

	/**
	 * Gets the maximum N that is tabulated.
	 *
	 * @return the max N
	 */
	public int getMaxN()
	{
		return objectTable.length - 1;
	}

	/**
	 * Increase the tabulated values up to a max n. Does nothing if already above the given n.
	 *
	 * @param n
	 *            the n
	 */
	public void increaseMaxN(int n)
	{
		if (getMaxN() < n)
			synchronized (objectLock)
			{
				objectTable = increaseSize(objectTable, n);
			}
	}

	/**
	 * Ensure the table contains values for the specified range of N. This can be called before using the object with a
	 * known range of n.
	 *
	 * @param minN
	 *            the min N
	 * @param maxN
	 *            the max N
	 */
	public void ensureRange(int minN, int maxN)
	{
		if (getMaxN() < maxN)
			synchronized (objectLock)
			{
				// Resize but do not compute.
				// Use the master table if it is bigger as that has all values pre-computed.
				final double[] masterTable = LogFactorial.table;
				objectTable = Arrays.copyOf((masterTable.length > objectTable.length) ? masterTable : objectTable,
						maxN + 1);
			}

		// Check range has pre-computed values
		for (int n = Math.max(2, minN); n <= maxN; n++)
			if (objectTable[n] == 0)
				objectTable[n] = Gamma.logGamma(n + 1);
	}

	/**
	 * Reduces the tabulated values down to a max n. Does nothing if already below the given n.
	 *
	 * @param n
	 *            the new table max N
	 */
	public void reduceMaxN(int n)
	{
		if (getMaxN() > n)
		{
			// Keep the representable factorials
			n = Math.max(n, FACTORIALS.length);

			synchronized (lock)
			{
				objectTable = Arrays.copyOf(objectTable, n + 1);
			}
		}
	}

	/**
	 * Get the log of n! using tabulated values.
	 *
	 * @param n
	 *            the n (must be positive)
	 * @return log(n!)
	 * @throws ArrayIndexOutOfBoundsException
	 *             if n is outside the table bounds
	 * @see #getMaxN()
	 */
	public double getLogF(int n) throws ArrayIndexOutOfBoundsException
	{
		return objectTable[n];
	}
}
