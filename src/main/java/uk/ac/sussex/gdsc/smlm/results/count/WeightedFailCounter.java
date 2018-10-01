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
package uk.ac.sussex.gdsc.smlm.results.count;

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
        return String.format("weightedFailures=%d;fail+%d;pass-%d", allowedFailures, failIncrement, passDecrement);
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

    /** {@inheritDoc} */
    @Override
    public void pass()
    {
        failCount -= passDecrement;
        if (failCount < 0L)
            failCount = 0L;
    }

    /** {@inheritDoc} */
    @Override
    public void pass(int n)
    {
        if (n < 0)
            throw new IllegalArgumentException("Number of passes must be positive");
        failCount -= n * passDecrement;
        if (failCount < 0L)
            failCount = 0L;
    }

    /** {@inheritDoc} */
    @Override
    public void fail()
    {
        failCount += failIncrement;
        if (failCount < 0L)
            throw new IllegalStateException("Unable to increment");
    }

    /** {@inheritDoc} */
    @Override
    public void fail(int n)
    {
        if (n < 0)
            throw new IllegalArgumentException("Number of fails must be positive");
        failCount += n * failIncrement;
        if (failCount < 0L)
            throw new IllegalStateException("Unable to increment");
    }

    /** {@inheritDoc} */
    @Override
    public boolean isOK()
    {
        return failCount <= allowedFailures;
    }

    /** {@inheritDoc} */
    @Override
    public FailCounter newCounter()
    {
        return new WeightedFailCounter(allowedFailures, failIncrement, passDecrement);
    }

    /** {@inheritDoc} */
    @Override
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
