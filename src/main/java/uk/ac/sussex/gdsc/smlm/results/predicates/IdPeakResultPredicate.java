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
package uk.ac.sussex.gdsc.smlm.results.predicates;

import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Test a result using the id.
 */
public class IdPeakResultPredicate implements PeakResultPredicate
{
    /** The id. */
    private final int id;

    /**
     * Instantiates a new id peak result predicate.
     *
     * @param id
     *            the id
     */
    public IdPeakResultPredicate(int id)
    {
        this.id = id;
    }

    /** {@inheritDoc} */
    @Override
    public boolean test(PeakResult t)
    {
        return t.getId() == id;
    }
}
