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

/**
 * Defines no normalisation of an area around a point
 */
public class NonNormaliser implements Normaliser
{
    /** An instance for applying no normalisation */
    public static final NonNormaliser INSTANCE = new NonNormaliser();

    /**
     * Normalise the sum.
     *
     * @param sum
     *            the sum
     * @param index
     *            the index
     * @return the normalised value
     */
    @Override
    public float normalise(double sum, int index)
    {
        return (float) sum;
    }

    /**
     * Normalise the sum.
     *
     * @param sum
     *            the sum
     * @param index
     *            the index
     * @return the normalised value
     */
    @Override
    public float normalise(float sum, int index)
    {
        return sum;
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.filters.Normaliser#normalise(float[], int)
     */
    @Override
    public void normalise(float[] data, int size)
    {
        // Do nothing
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.filters.Normaliser#normalise(float[], float[], int)
     */
    @Override
    public void normalise(float[] data, float[] out, int size)
    {
        System.arraycopy(data, 0, out, 0, size);
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.filters.Normaliser#normalise(float[], int, int, int)
     */
    @Override
    public void normalise(float[] data, int maxx, int maxy, int border)
    {
        // Do nothing
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.filters.Normaliser#normalise(float[], float[], int, int, int)
     */
    @Override
    public void normalise(float[] data, float[] out, int maxx, int maxy, int border)
    {
        final int width = maxx - 2 * border;
        if (width <= 0)
            return;
        final int max = (maxy - border) * maxx;
        for (int i = border * maxx + border; i < max; i += maxx)
            System.arraycopy(data, i, out, i, width);
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.filters.Normaliser#normalise(double[], float[], int)
     */
    @Override
    public void normalise(double[] data, float[] out, int size)
    {
        for (int i = 0; i < size; i++)
            out[i] = (float) (data[i]);
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.filters.Normaliser#normalise(double[], float[], int, int, int)
     */
    @Override
    public void normalise(double[] data, float[] out, int maxx, int maxy, int border)
    {
        final int xlimit = maxx - border;
        final int ylimit = maxy - border;
        for (int y = border; y < ylimit; y++)
            for (int x = border, i = y * maxx + border; x < xlimit; x++, i++)
                out[i] = (float) (data[i]);
    }
}
