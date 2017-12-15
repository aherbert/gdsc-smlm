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
 * Combine the result of two fail counters using an OR operator
 */
public class CombinedOrFailCounter extends CombinedFailCounter
{
	/**
	 * Instantiates a new combined or fail counter.
	 *
	 * @param c1
	 *            the first counter
	 * @param c2
	 *            the second counter
	 */
	public CombinedOrFailCounter(FailCounter c1, FailCounter c2)
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
		return c1.isOK() || c2.isOK();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#newCounter()
	 */
	public FailCounter newCounter()
	{
		return new CombinedOrFailCounter(c1.newCounter(), c2.newCounter());
	}
}