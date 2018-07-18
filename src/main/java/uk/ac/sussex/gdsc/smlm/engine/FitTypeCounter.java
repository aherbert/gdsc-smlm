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
package uk.ac.sussex.gdsc.smlm.engine;

import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.math3.util.FastMath;

/**
 * Count the type of fit that was performed
 */
public class FitTypeCounter
{
	private final AtomicInteger[] count;

	/**
	 * Instantiates a new fit type counter.
	 */
	public FitTypeCounter()
	{
		this.count = new AtomicInteger[(int) FastMath.pow(2, FitType.NO_OF_FLAGS)];
		for (int i = 0; i < count.length; i++)
			count[i] = new AtomicInteger();
	}

	/**
	 * @return the total number of flags for the type of fit
	 */
	public int size()
	{
		return count.length;
	}

	/**
	 * Add a single count of the given type.
	 *
	 * @param fitType
	 *            the fit type
	 */
	public void add(FitType fitType)
	{
		count[fitType.getFlags()].incrementAndGet();
	}

	/**
	 * Add a value of the given type to the count.
	 *
	 * @param fitType
	 *            the fit type
	 * @param value
	 *            the value
	 */
	public void add(FitType fitType, int value)
	{
		count[fitType.getFlags()].addAndGet(value);
	}

	/**
	 * Get the count of the given type.
	 *
	 * @param fitType
	 *            the fit type
	 * @return The count
	 */
	public int get(FitType fitType)
	{
		return count[fitType.getFlags()].get();
	}

	/**
	 * Get the count of the given type
	 *
	 * @param flags
	 *            The flags that must be set
	 * @return The count
	 */
	public int getSet(int flags)
	{
		if (flags == 0)
			return count[0].get();

		int total = 0;
		for (int i = 1; i < count.length; i++)
			if ((i & flags) == flags)
				total += count[i].get();
		return total;
	}

	/**
	 * Get the count of the given type
	 *
	 * @param flags
	 *            The flags that must be unset
	 * @return The count
	 */
	public int getUnset(int flags)
	{
		int total = 0;
		if (flags == 0)
			// Count all but zero
			for (int i = 1; i < count.length; i++)
				total += count[i].get();
		else
			for (int i = 0; i < count.length; i++)
				if ((i & flags) == 0)
					total += count[i].get();
		return total;
	}

	/**
	 * Get the count of the given type
	 *
	 * @param setFlags
	 *            The flags that must be set
	 * @param unsetFlags
	 *            The flags that must be unset
	 * @return The count
	 */
	public int get(int setFlags, int unsetFlags)
	{
		// Check set flags and unset flags do not clash
		// If they do then return 0
		// TODO - check this works...
		if ((setFlags & unsetFlags) != 0)
			return 0;

		int total = 0;
		for (int i = 0; i < count.length; i++)
			if ((i & setFlags) == setFlags && (i & unsetFlags) == 0)
				total += count[i].get();
		return total;
	}

	/**
	 * @return The total count
	 */
	public int getTotal()
	{
		int total = 0;
		for (int i = 0; i < count.length; i++)
			total += count[i].get();
		return total;
	}
}
