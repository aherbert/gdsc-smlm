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
package gdsc.smlm.results.count;

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
	@Override
	public void pass()
	{
		failCount = 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#pass(int)
	 */
	@Override
	public void pass(int n)
	{
		failCount = 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#fail()
	 */
	@Override
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
	@Override
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
	@Override
	public boolean isOK()
	{
		return failCount <= allowedFailures;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#newCounter()
	 */
	@Override
	public FailCounter newCounter()
	{
		return new ConsecutiveFailCounter(allowedFailures);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.FailCounter#reset()
	 */
	@Override
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
