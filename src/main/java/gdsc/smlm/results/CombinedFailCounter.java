package gdsc.smlm.results;

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
 * Combine the result of two fail counters
 */
public abstract class CombinedFailCounter extends BaseFailCounter
{
	protected final FailCounter c1, c2;

	/**
	 * Instantiates a new combined fail counter.
	 *
	 * @param c1
	 *            the first counter
	 * @param c2
	 *            the second counter
	 */
	public CombinedFailCounter(FailCounter c1, FailCounter c2)
	{
		if (c1 == null || c2 == null)
			throw new NullPointerException();
		this.c1 = c1;
		this.c2 = c2;
	}

	@Override
	protected String generateDescription()
	{
		StringBuilder sb = new StringBuilder();
		add(sb, c1);
		sb.append(" ").append(getOperator()).append(" ");
		add(sb, c2);
		return sb.toString();
	}

	private void add(StringBuilder sb, FailCounter c)
	{
		if (c instanceof CombinedFailCounter)
		{
			sb.append("(");
			sb.append(c.getDescription());
			sb.append(")");
		}
		else
		{
			sb.append(c.getDescription());
		}
	}

	/**
	 * Get the string representation of the operator used to combine the two fail counters. This is used in the filter
	 * name.
	 * 
	 * @return The operator
	 */
	protected abstract String getOperator();

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#pass()
	 */
	public void pass()
	{
		c1.pass();
		c2.pass();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#pass(int)
	 */
	public void pass(int n)
	{
		c1.pass(n);
		c2.pass(n);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#fail()
	 */
	public void fail()
	{
		c1.fail();
		c2.fail();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#fail(int)
	 */
	public void fail(int n)
	{
		c1.fail(n);
		c2.fail(n);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#reset()
	 */
	public void reset()
	{
		c1.reset();
		c2.reset();
	}
}