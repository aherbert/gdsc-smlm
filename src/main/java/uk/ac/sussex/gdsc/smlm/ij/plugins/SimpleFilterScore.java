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
package uk.ac.sussex.gdsc.smlm.ij.plugins;

import uk.ac.sussex.gdsc.smlm.results.filter.DirectFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.FilterScore;

/**
 * Store the filter score used in benchmarking.
 */
public class SimpleFilterScore extends FilterScore
{
    /** The result. */
    final FilterScoreResult r;

    /**
     * Instantiates a new simple filter score.
     *
     * @param r
     *            the result
     * @param allSameType
     *            the all same type
     * @param criteriaPassed
     *            the criteria passed
     */
    public SimpleFilterScore(FilterScoreResult r, boolean allSameType, boolean criteriaPassed)
    {
        super(r.filter, r.score, r.criteria, allSameType, criteriaPassed);
        this.r = r;
    }

    @Override
    protected int compareParameters(FilterScore that)
    {
        // We only use this class with DirectFilter
        return ((DirectFilter) that.filter).weakestUnsafe((DirectFilter) this.filter);
    }
}
