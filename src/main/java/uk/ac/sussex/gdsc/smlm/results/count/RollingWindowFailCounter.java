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

import uk.ac.sussex.gdsc.core.utils.BooleanRollingArray;

/**
 * Stop evaluating when a number of failures occurs within a window.
 */
public class RollingWindowFailCounter extends BaseFailCounter
{
    /** The rolling array stores true for a failure. */
    private final BooleanRollingArray rollingArray;

    /** The number of allowed failures. */
    private final int allowedFailures;

    /**
     * Instantiates a new rolling window fail counter.
     *
     * @param allowedFailures
     *            the number of allowed failures
     * @param window
     *            the window size
     */
    private RollingWindowFailCounter(int allowedFailures, int window)
    {
        this.allowedFailures = allowedFailures;
        rollingArray = new BooleanRollingArray(window);
    }

    @Override
    protected String generateDescription()
    {
        return "rollingFailures=" + allowedFailures + "/" + getWindow();
    }

    /**
     * Instantiates a new rolling window fail counter.
     *
     * @param allowedFailures
     *            the number of allowed failures
     * @param window
     *            the window size
     * @return the rolling window fail counter
     * @throws IllegalArgumentException
     *             If the window is not strictly positive, or the window is smaller that the allowed failures
     */
    public static RollingWindowFailCounter create(int allowedFailures, int window) throws IllegalArgumentException
    {
        if (window < 1)
            throw new IllegalArgumentException("Window must be strictly positive");
        if (window < allowedFailures)
            throw new IllegalArgumentException("Window must be larger than the allowed failures");
        return new RollingWindowFailCounter(Math.max(0, allowedFailures), window);
    }

    /** {@inheritDoc} */
    @Override
    public void pass()
    {
        rollingArray.add(false);
    }

    /** {@inheritDoc} */
    @Override
    public void pass(int n)
    {
        if (n < 0)
            throw new IllegalArgumentException("Number of passes must be positive");
        rollingArray.add(false, n);
    }

    /** {@inheritDoc} */
    @Override
    public void fail()
    {
        rollingArray.add(true);
    }

    /** {@inheritDoc} */
    @Override
    public void fail(int n)
    {
        if (n < 0)
            throw new IllegalArgumentException("Number of fails must be positive");
        rollingArray.add(true, n);
    }

    /** {@inheritDoc} */
    @Override
    public boolean isOK()
    {
        return (rollingArray.isFull()) ? getFailCount() <= allowedFailures : true;
    }

    /** {@inheritDoc} */
    @Override
    public FailCounter newCounter()
    {
        return new RollingWindowFailCounter(allowedFailures, getWindow());
    }

    /** {@inheritDoc} */
    @Override
    public void reset()
    {
        rollingArray.clear();
    }

    /**
     * Gets the window size.
     *
     * @return the window size
     */
    public int getWindow()
    {
        return rollingArray.getCapacity();
    }

    /**
     * Gets the fail count within the current window.
     *
     * @return the fail count
     */
    public int getFailCount()
    {
        return rollingArray.getTrueCount();
    }

    /**
     * Gets the current window size. This may be smaller than the window size if not enough pass/fail events have been
     * registered.
     *
     * @return the current window size
     */
    public int getCurrentWindowSize()
    {
        return rollingArray.getCount();
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
