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
 * Combine the result of two fail counters using an AND operator
 */
public class CombinedAndFailCounter extends CombinedFailCounter
{
	/**
	 * Instantiates a new combined or fail counter.
	 *
	 * @param c1
	 *            the first counter
	 * @param c2
	 *            the second counter
	 */
	public CombinedAndFailCounter(FailCounter c1, FailCounter c2)
	{
		super(c1, c2);
	}

	@Override
	protected String getOperator()
	{
		return "&&";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#isOK()
	 */
	public boolean isOK()
	{
		return c1.isOK() && c2.isOK();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#newCounter()
	 */
	public FailCounter newCounter()
	{
		return new CombinedAndFailCounter(c1.newCounter(), c2.newCounter());
	}

	/**
	 * Join the fail counters. 
	 * <p>
	 * If both are not null then return a combined fail counter. 
	 * <p>
	 * If either are null then a single counter will be returned.
	 * <p>
	 * If both are null then null will be returned.
	 *
	 * @param c1
	 *            the first counter
	 * @param c2
	 *            the second counter
	 * @return the fail counter
	 */
	public static FailCounter join(FailCounter c1, FailCounter c2)
	{
		if (c1 != null)
			return (c2 != null) ? new CombinedAndFailCounter(c1, c2) : c1;
		return c2;
	}
}