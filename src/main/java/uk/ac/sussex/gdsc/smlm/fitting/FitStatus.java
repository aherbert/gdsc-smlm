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
package uk.ac.sussex.gdsc.smlm.fitting;

/**
 * Define the status of a fit result.
 */
public enum FitStatus
{
    //@formatter:off
    /** OK. */
    OK{ @Override public String getName() { return "OK"; }},
    /** Singular non-linear model */
    SINGULAR_NON_LINEAR_MODEL{ @Override public String getName() { return "Singular non-linear model"; }},
    /** Singular non-linear solution */
    SINGULAR_NON_LINEAR_SOLUTION{ @Override public String getName() { return "Singular non-linear solution"; }},
    /** Invalid gradients. */
    INVALID_GRADIENTS{ @Override public String getName() { return "Invalid gradients"; }},
    /** Failed to converge. */
    FAILED_TO_CONVERGE{ @Override public String getName() { return "Failed to converge"; }},
    /** Too many iterations. */
    TOO_MANY_ITERATIONS{ @Override public String getName() { return "Too many iterations"; }},
    /** Too many evaluations. */
    TOO_MANY_EVALUATIONS{ @Override public String getName() { return "Too many evaluations"; }},
    /** Invalid likelihood. */
    INVALID_LIKELIHOOD{ @Override public String getName() { return "Invalid likelihood"; }},
    /** Bad parameters. */
    BAD_PARAMETERS{ @Override public String getName() { return "Bad parameters"; }},
    /** Failed to estimate width. */
    FAILED_TO_ESTIMATE_WIDTH{ @Override public String getName() { return "Failed to estimate width"; }},
    /** Coordinates moved. */
    COORDINATES_MOVED{ @Override public String getName() { return "Coordinates moved"; }},
    /** Outside fit region. */
    OUTSIDE_FIT_REGION{ @Override public String getName() { return "Outside fit region"; }},
    /** Insufficient signal. */
    INSUFFICIENT_SIGNAL{ @Override public String getName() { return "Insufficient signal"; }},
    /** Insufficient Signal-to-Noise Ratio (SNR) */
    INSUFFICIENT_SNR{ @Override public String getName() { return "Insufficient SNR"; }},
    /** Width diverged. */
    WIDTH_DIVERGED{ @Override public String getName() { return "Width diverged"; }},
    /** Z-coordinate moved */
    Z_MOVED{ @Override public String getName() { return "Z-coordinate moved"; }},
    /** Insufficient precision. */
    INSUFFICIENT_PRECISION{ @Override public String getName() { return "Insufficient precision"; }},
    /** Neighbour overlap. */
    NEIGHBOUR_OVERLAP{ @Override public String getName() { return "Neighbour overlap"; }},
    /** Failed smart filter. */
    FAILED_SMART_FILTER{ @Override public String getName() { return "Failed smart filter"; }},
    /** Drift to another result. */
    DRIFT_TO_ANOTHER_RESULT{ @Override public String getName() { return "Drift to another result"; }},
    /** Failed validation. */
    FAILED_VALIDATION{ @Override public String getName() { return "Failed validation"; }},
    /** No model improvement. */
    NO_MODEL_IMPROVEMENT{ @Override public String getName() { return "No model improvement"; }},
    /** Line search error. */
    LINE_SEARCH_ERROR{ @Override public String getName() { return "Line search error"; }},
    /** Unknown. */
    UNKNOWN{ @Override public String getName() { return "Unknown"; }};
	//@formatter:on

    @Override
    public String toString()
    {
        return getName();
    }

    /**
     * Gets the name.
     *
     * @return the name
     */
    abstract public String getName();
}
