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
package uk.ac.sussex.gdsc.smlm.results.filter;

import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Filter results using the combination of two filters. Results must pass both filters
 */
public class AndFilter extends CombinedFilter
{
    /**
     * Instantiates a new and filter.
     *
     * @param filter1
     *            the filter 1
     * @param filter2
     *            the filter 2
     */
    public AndFilter(Filter filter1, Filter filter2)
    {
        super(filter1, filter2);
    }

    @Override
    protected String getOperator()
    {
        return "&&";
    }

    @Override
    public boolean accept(PeakResult peak)
    {
        return accept1(peak) && accept2(peak);
    }

    @Override
    public int validate(final PreprocessedPeakResult peak)
    {
        if (accept1(peak) && accept2(peak))
            return 0;
        // We get here if filter 1 failed; or filter 1 passed but filter 2 failed.
        return (result1 == 0) ? result2 : result1;
    }

    @Override
    protected Filter createFilter(Filter f1, Filter f2)
    {
        return new AndFilter(f1, f2);
    }

    @Override
    public Filter clone()
    {
        return new AndFilter(filter1.clone(), filter2.clone());
    }
}
