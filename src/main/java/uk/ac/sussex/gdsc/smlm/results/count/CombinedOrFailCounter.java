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
     * @see uk.ac.sussex.gdsc.smlm.results.FailCounter#isOK()
     */
    @Override
    public boolean isOK()
    {
        return c1.isOK() || c2.isOK();
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.results.FailCounter#newCounter()
     */
    @Override
    public FailCounter newCounter()
    {
        return new CombinedOrFailCounter(c1.newCounter(), c2.newCounter());
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
            return (c2 != null) ? new CombinedOrFailCounter(c1, c2) : c1;
        return c2;
    }
}
