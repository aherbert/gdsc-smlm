package gdsc.smlm.results.count;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Base class for fail counters
 */
public abstract class BaseFailCounter implements FailCounter
{
	private String description;

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#getDescription()
	 */
	public String getDescription()
	{
		if (description == null)
			description = generateDescription();
		return description;
	}

	/**
	 * Generate the description.
	 *
	 * @return the description
	 */
	protected abstract String generateDescription();

	/**
	 * Check the number is positive.
	 *
	 * @param c
	 *            the number
	 * @throws IllegalStateException
	 *             if not positive
	 */
	protected static void checkPositive(int c) throws IllegalStateException
	{
		if (c < 0)
			throw new IllegalStateException("Negative count: " + c);
	}
}