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
 * Stop evaluating when a number of cumulative failures occurs. The counting is weighted so that fails increment and
 * passes decrement different amounts.
 */
public class WeightedFailCounter extends BaseFailCounter
{
	/** The fail count. Use a long to avoid overflow errors. */
	private long failCount = 0L;

	/** The number of allowed failures. */
	private final int allowedFailures;

	/** The amount to increment for a fail. */
	private final int failIncrement;

	/** The amount to decrement for a pass. */
	private final int passDecrement;

	/**
	 * Instantiates a new weighted fail counter.
	 *
	 * @param allowedFailures
	 *            the number of allowed failures
	 * @param failIncrement
	 *            the fail increment
	 * @param passDecrement
	 *            the pass decrement
	 */
	private WeightedFailCounter(int allowedFailures, int failIncrement, int passDecrement)
	{
		this.allowedFailures = allowedFailures;
		this.failIncrement = failIncrement;
		this.passDecrement = passDecrement;
	}

	@Override
	protected String generateDescription()
	{
		return String.format("cumulativeFailures=%d;fail+%d;pass-%d", allowedFailures, failIncrement, passDecrement);
	}

	/**
	 * Instantiates a new weighted fail counter.
	 *
	 * @param allowedFailures
	 *            the number of allowed failures
	 * @param failIncrement
	 *            the amount to increment for a fail (set to 1 if below 1)
	 * @param passDecrement
	 *            the amount to decrement for a pass (set to 0 if below 0)
	 * @return the weighted fail counter
	 */
	public static WeightedFailCounter create(int allowedFailures, int failIncrement, int passDecrement)
	{
		return new WeightedFailCounter(Math.max(0, allowedFailures), Math.max(1, failIncrement),
				Math.max(0, passDecrement));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#pass()
	 */
	public void pass()
	{
		failCount -= passDecrement;
		if (failCount < 0L)
			failCount = 0L;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#pass(int)
	 */
	public void pass(int n)
	{
		if (n < 0)
			throw new IllegalArgumentException("Number of passes must be positive");
		failCount -= n * passDecrement;
		if (failCount < 0L)
			failCount = 0L;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#fail()
	 */
	public void fail()
	{
		failCount += failIncrement;
		if (failCount < 0L)
			throw new IllegalStateException("Unable to increment");
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
		failCount += n * failIncrement;
		if (failCount < 0L)
			throw new IllegalStateException("Unable to increment");
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
		return new WeightedFailCounter(allowedFailures, failIncrement, passDecrement);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#reset()
	 */
	public void reset()
	{
		failCount = 0L;
	}

	/**
	 * Gets the fail count.
	 *
	 * @return the fail count
	 */
	public long getFailCount()
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

	/**
	 * Gets the amount to increment for a fail.
	 *
	 * @return the fail increment
	 */
	public int getFailIncrement()
	{
		return failIncrement;
	}

	/**
	 * Gets the amount to decrement for a pass
	 *
	 * @return the pass decrement
	 */
	public int getPassDecrement()
	{
		return passDecrement;
	}
}
