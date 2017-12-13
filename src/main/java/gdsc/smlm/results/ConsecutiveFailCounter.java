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
 * Stop evaluating when a number of consecutive failures occurs.
 */
public class ConsecutiveFailCounter extends BaseFailCounter
{
	/** The fail count. */
	private int failCount = 0;

	/** The number of allowed failures. */
	private final int allowedFailures;

	/**
	 * Instantiates a new consecutive fail counter.
	 *
	 * @param allowedFailures
	 *            the number of allowed failures
	 */
	private ConsecutiveFailCounter(int allowedFailures)
	{
		this.allowedFailures = allowedFailures;
	}

	@Override
	protected String generateDescription()
	{
		return "consecutiveFailures=" + allowedFailures;
	}

	/**
	 * Instantiates a new consecutive fail counter.
	 *
	 * @param allowedFailures
	 *            the number of allowed failures
	 */
	public static ConsecutiveFailCounter create(int allowedFailures)
	{
		return new ConsecutiveFailCounter(Math.max(0, allowedFailures));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#pass()
	 */
	public void pass()
	{
		failCount = 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#pass(int)
	 */
	public void pass(int n)
	{
		failCount = 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#fail()
	 */
	public void fail()
	{
		failCount++;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#fail(int)
	 */
	public void fail(int n)
	{
		if (n < 0)
			throw new IllegalArgumentException("Number of fails must be positive");
		failCount += n;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#isOK()
	 */
	public boolean isOK()
	{
		return failCount <= allowedFailures;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#newCounter()
	 */
	public FailCounter newCounter()
	{
		return new ConsecutiveFailCounter(allowedFailures);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#reset()
	 */
	public void reset()
	{
		failCount = 0;
	}

	/**
	 * Gets the fail count.
	 *
	 * @return the fail count
	 */
	public int getFailCount()
	{
		return failCount;
	}

	/**
	 * Gets the number of allowed failures.
	 *
	 * @return the number of allowed failures.
	 */
	public int getAllowedFailures()
	{
		return allowedFailures;
	}
}
