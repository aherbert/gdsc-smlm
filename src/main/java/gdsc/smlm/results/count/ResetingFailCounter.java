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
 * Stop evaluating when a number of cumulative failures occurs. The failures count is reset to a fraction of the current
 * value for each pass.
 */
public class ResetingFailCounter extends BaseFailCounter
{
	/** The fail count. */
	private int failCount = 0;

	/** The number of allowed failures. */
	private final int allowedFailures;

	/** The fraction of the current failures count to reset to for a pass. */
	private final double resetFraction;

	/**
	 * Instantiates a new weighted fail counter.
	 *
	 * @param allowedFailures
	 *            the number of allowed failures
	 * @param resetFraction
	 *            the reset fraction
	 */
	private ResetingFailCounter(int allowedFailures, double resetFraction)
	{
		this.allowedFailures = allowedFailures;
		this.resetFraction = resetFraction;
	}

	@Override
	protected String generateDescription()
	{
		return String.format("allowedFailures=%d;resetFraction=%f", allowedFailures, resetFraction);
	}

	/**
	 * Instantiates a new weighted fail counter.
	 *
	 * @param allowedFailures
	 *            the number of allowed failures
	 * @param resetFraction
	 *            The fraction of the current failures count to reset to for a pass.
	 * @return the weighted fail counter
	 */
	public static ResetingFailCounter create(int allowedFailures, double resetFraction)
	{
		if (!(resetFraction >= 0 && resetFraction <= 1))
			throw new IllegalArgumentException("Reset must be in the range 0-1");
		return new ResetingFailCounter(Math.max(0, allowedFailures), resetFraction);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#pass()
	 */
	public void pass()
	{
		failCount = (int) (failCount * resetFraction);
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
		while (n-- > 0)
		{
			pass();
			if (failCount == 0)
				break;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#fail()
	 */
	public void fail()
	{
		if (failCount == Integer.MAX_VALUE)
			throw new IllegalStateException("Unable to increment");
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
		if (Integer.MAX_VALUE - n < failCount)
			throw new IllegalStateException("Unable to increment");
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
		return new ResetingFailCounter(allowedFailures, resetFraction);
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
	 * Gets the fraction of the current failures count to reset to for a pass.
	 *
	 * @return the reset fraction
	 */
	public double getResetFraction()
	{
		return resetFraction;
	}

}
