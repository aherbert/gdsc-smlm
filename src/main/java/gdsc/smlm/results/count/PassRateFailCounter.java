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
 * Stop evaluating when the pass rate falls below a set fraction. A minimum number of pass/fail counts can be specified.
 */
public class PassRateFailCounter extends BaseFailCounter
{
	/** The pass count. */
	private int passCount = 0;

	/** The fail count. */
	private int failCount = 0;

	/** The number of allowed counts. */
	private final int allowedCounts;

	/** The pass rate. */
	private final double passRate;

	/**
	 * Instantiates a new pass rate fail counter.
	 *
	 * @param allowedCounts
	 *            the number of allowed counts before the pass rate is evaluated
	 * @param passRate
	 *            the reset fraction
	 */
	private PassRateFailCounter(int allowedCounts, double passRate)
	{
		this.allowedCounts = allowedCounts;
		this.passRate = passRate;
	}

	@Override
	protected String generateDescription()
	{
		return String.format("allowedCounts=%d;passRate=%f", allowedCounts, passRate);
	}

	/**
	 * Instantiates a new pass rate fail counter.
	 *
	 * @param allowedCounts
	 *            the number of allowed counts before the pass rate is evaluated
	 * @param passRate
	 *            The fraction of the current failures count to reset to for a pass.
	 * @return the pass rate fail counter
	 */
	public static PassRateFailCounter create(int allowedCounts, double passRate)
	{
		if (!(passRate >= 0 && passRate <= 1))
			throw new IllegalArgumentException("Reset must be in the range 0-1");
		return new PassRateFailCounter(Math.max(0, allowedCounts), passRate);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.FailCounter#pass()
	 */
	@Override
	public void pass()
	{
		passCount++;
		if (passCount < 0)
			throw new IllegalStateException("Unable to increment");
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.FailCounter#pass(int)
	 */
	@Override
	public void pass(int n)
	{
		if (n < 0)
			throw new IllegalArgumentException("Number of passes must be positive");
		passCount += n;
		if (passCount < 0)
			throw new IllegalStateException("Unable to increment");
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.FailCounter#fail()
	 */
	@Override
	public void fail()
	{
		failCount++;
		if (failCount < 0)
			throw new IllegalStateException("Unable to increment");
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
		failCount += n;
		if (failCount < 0)
			throw new IllegalStateException("Unable to increment");
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.FailCounter#isOK()
	 */
	@Override
	public boolean isOK()
	{
		double total = failCount + passCount;
		return total <= allowedCounts || passCount / total >= passRate;
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.FailCounter#newCounter()
	 */
	@Override
	public FailCounter newCounter()
	{
		return new PassRateFailCounter(allowedCounts, passRate);
	}

	/*
	 * (non-Javadoc)
	 *
	 * @see gdsc.smlm.results.FailCounter#reset()
	 */
	@Override
	public void reset()
	{
		passCount = 0;
		failCount = 0;
	}

	/**
	 * Gets the pass count.
	 *
	 * @return the pass count
	 */
	public int getPassCount()
	{
		return passCount;
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
	 * Gets the number of allowed counts before the pass rate is evaluated.
	 *
	 * @return the allowed counts
	 */
	public int getAllowedCounts()
	{
		return allowedCounts;
	}

	/**
	 * Gets the pass rate.
	 *
	 * @return the pass rate
	 */
	public double getPassRate()
	{
		return passRate;
	}
}
