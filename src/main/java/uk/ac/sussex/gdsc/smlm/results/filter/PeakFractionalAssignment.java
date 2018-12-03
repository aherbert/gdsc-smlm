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

import uk.ac.sussex.gdsc.core.match.ImmutableFractionalAssignment;

/**
 * Extends the fractional assignment to add a reference to the peak result.
 */
public class PeakFractionalAssignment extends ImmutableFractionalAssignment
{
    /** The peak result. */
    public final PreprocessedPeakResult peakResult;

    /**
     * Instantiates a new custom fractional assignment.
     *
     * @param targetId
     *            the target id
     * @param predictedId
     *            the predicted id
     * @param distance
     *            the distance
     * @param score
     *            the score
     * @param peakResult
     *            the peak result
     */
    public PeakFractionalAssignment(int targetId, int predictedId, double distance, double score,
            PreprocessedPeakResult peakResult)
    {
        super(targetId, predictedId, distance, score);
        this.peakResult = peakResult;
    }
}
