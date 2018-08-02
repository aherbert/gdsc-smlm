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
package uk.ac.sussex.gdsc.smlm.filters;

import java.util.Arrays;
import java.util.List;

/**
 * Identifies candidate spots (local maxima) in an image.
 */
public abstract class SpotFilter implements Cloneable
{
    /**
     * Checks if the filter is weighted, i.e. supports {@link #setWeights(float[], int, int)}.
     *
     * @return true, if is weighted
     */
    public abstract boolean isWeighted();

    /**
     * Sets the weights of the data. This should be called before {@link #list(float[], int, int)} or
     * {@link #rank(float[], int, int)} is called with data samples.
     * <p>
     * Calling this in advance allows efficient caching of pre-computed weightings.
     *
     * @param weights
     *            the weights of the data (can be null)
     * @param width
     *            The width of the data
     * @param height
     *            The height of the data
     */
    public abstract void setWeights(final float[] weights, final int width, final int height);

    /**
     * Checks for weights. Weights are set using {@link #setWeights(float[], int, int)}.
     *
     * @return true, if successful
     */
    public abstract boolean hasWeights();

    /**
     * Find the candidate spots in the data.
     *
     * @param data
     *            The data
     * @param width
     *            The width of the data
     * @param height
     *            The height of the data
     * @return The candidate spots
     */
    protected abstract Spot[] find(final float[] data, final int width, final int height);

    /**
     * Get the pre-processed data produced by the find method
     *
     * @return The pre-processed data produced by the {@link #find(float[], int, int)} method
     */
    public abstract float[] getPreprocessedData();

    /**
     * List the candidate spots in the data. The list will be in the order the candidates are found.
     *
     * @param data
     *            The data
     * @param width
     *            The width of the data
     * @param height
     *            The height of the data
     * @return The candidate spots (may be an empty array but will not be null)
     */
    public Spot[] list(float[] data, int width, int height)
    {
        final Spot[] spots = find(data, width, height);
        return (spots == null) ? new Spot[0] : spots;
    }

    /**
     * List and then rank the candidate spots in the data. The list will be in the order defined by sorting the
     * candidates.
     *
     * @param data
     *            The data
     * @param width
     *            The width of the data
     * @param height
     *            The height of the data
     * @return The candidate spots (may be an empty array but will not be null)
     */
    public Spot[] rank(float[] data, int width, int height)
    {
        final Spot[] spots = find(data, width, height);
        if (spots == null)
            return new Spot[0];
        Arrays.sort(spots);
        return spots;
    }

    /**
     * Return true if the intensity value of the candidate spots is absolute, i.e. is the height of the candidate
     * maximum using the original data scale. Otherwise the intensity is relative, for example this could be relative to
     * the local background.
     *
     * @return True if the intensity value of the candidate spots is absolute
     */
    public abstract boolean isAbsoluteIntensity();

    /*
     * (non-Javadoc)
     *
     * @see java.lang.Object#clone()
     */
    @Override
    public SpotFilter clone()
    {
        try
        {
            return (SpotFilter) super.clone();
        }
        catch (final CloneNotSupportedException e)
        {
            return null;
        }
    }

    /**
     * @return A description of the filter and parameters
     */
    public String getDescription()
    {
        return getName() + ": " + Arrays.toString(getParameters().toArray());
    }

    /**
     * @return The name of the filter
     */
    public abstract String getName();

    /**
     * @return The parameters of the filter
     */
    public abstract List<String> getParameters();

    /**
     * Get the width spread of data used to process each position
     *
     * @return The spread
     */
    public abstract double getSpread();
}
