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
package uk.ac.sussex.gdsc.smlm.engine.filter;

import java.util.List;

import uk.ac.sussex.gdsc.smlm.fitting.FitResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;

/**
 * Provides a system for filtering the fitted results using distance to a set of target coordinates.
 * <p>
 * Note that the target coordinates should be relative to the fitting region bounds, not the bounds of the data frame.
 *
 * @deprecated Filtering of the results is no longer supported
 */
@Deprecated
public abstract class ResultFilter
{
    /** The filter. */
    protected List<float[]> filter;

    /** The distance squared. */
    protected float d2;

    /** The number of maxima. */
    protected int nMaxima;

    /** The filtered count. */
    protected int filteredCount = 0;

    /** The filtered fit results. */
    protected FitResult[] filteredFitResults;

    /** The filtered indices. */
    protected int[] filteredIndices;

    /** The peak results. */
    protected List<PeakResult> peakResults;

    /**
     * Instantiates a new result filter.
     *
     * @param filter
     *            The list of target coordinates (relative to the fitting region bounds)
     * @param d
     *            The distance
     * @param nMaxima
     *            The potential number of maxima that will be fitted
     */
    public ResultFilter(List<float[]> filter, float d, int nMaxima)
    {
        if (filter == null)
            throw new IllegalArgumentException("null filter list");
        for (final float[] f : filter)
            if (f == null)
                throw new IllegalArgumentException("null array used for filter element");
        this.filter = filter;
        d2 = d * d;
        this.nMaxima = nMaxima;
    }

    /**
     * Pass in a list of fitted peaks to be filtered. Called when fitting was successful
     *
     * @param fitResult
     *            The output from the fitting routine
     * @param maxIndex
     *            The source index that was fitted
     * @param results
     *            The fitted peaks
     */
    public abstract void filter(FitResult fitResult, int maxIndex, PeakResult... results);

    /**
     * Pass in a starting coordinate to be filtered. Called when fitting was unsuccessful but the starting point can
     * still be filtered.
     *
     * @param fitResult
     *            The output from the fitting routine
     * @param maxIndex
     *            The source index that was fitted
     * @param x
     *            The x position of the source index
     * @param y
     *            The y position of the source index
     */
    public abstract void filter(FitResult fitResult, int maxIndex, float x, float y);

    /**
     * Called when all results have been input and the final filtered results are required.
     */
    public abstract void finalise();

    /**
     * @return The number of results that pass the filter.
     */
    public int getFilteredCount()
    {
        return filteredCount;
    }

    /**
     * @return The indices that were fitted that pass the filter.
     */
    public int[] getMaxIndices()
    {
        return filteredIndices;
    }

    /**
     * @return The output from the fitting routine of any positions that pass the filter. This can includes
     *         positions that were within distance of the target coordinates but that were not fitted (e.g. due to
     *         failure to converge, etc)
     */
    public FitResult[] getFitResults()
    {
        return filteredFitResults;
    }

    /**
     * @return The fitted peaks that pass the filter.
     */
    public List<PeakResult> getResults()
    {
        return peakResults;
    }

}
