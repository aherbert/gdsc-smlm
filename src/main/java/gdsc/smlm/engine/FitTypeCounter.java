package gdsc.smlm.engine;

import java.util.concurrent.atomic.AtomicInteger;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2015 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Count the type of fit that was performed
 */
public class FitTypeCounter
{
	private AtomicInteger[] count;

	public FitTypeCounter()
	{
		this.count = new AtomicInteger[(int) Math.pow(2, FitType.NO_OF_FLAGS)];
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
	 * Add a single count of the given type
	 * 
	 * @param fitType
	 */
	public void add(FitType fitType)
	{
		count[fitType.getFlags()].incrementAndGet();
	}

	/**
	 * Add a value of the given type to the count
	 * 
	 * @param fitType
	 * @param value
	 */
	public void add(FitType fitType, int value)
	{
		count[fitType.getFlags()].addAndGet(value);
	}

	/**
	 * Get the count of the given type
	 * 
	 * @param fitType
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
		{
			// Count all but zero
			for (int i = 1; i < count.length; i++)
				total += count[i].get();
		}
		else
		{
			for (int i = 0; i < count.length; i++)
				if ((i & flags) == 0)
					total += count[i].get();
		}
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