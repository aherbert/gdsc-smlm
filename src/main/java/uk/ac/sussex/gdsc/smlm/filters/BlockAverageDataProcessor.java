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

import org.apache.commons.math3.util.FastMath;

import uk.ac.sussex.gdsc.core.ij.Utils;

/**
 * Identifies candidate spots (local maxima) in an image. The image is smoothed with an average box filter.
 */
public class BlockAverageDataProcessor extends DataProcessor
{
    private final double smooth;
    private BlockMeanFilter filter = new BlockMeanFilter();

    /**
     * Constructor
     *
     * @param border
     *            The border to ignore for maxima
     * @param smooth
     *            The smoothing width to apply to the data
     */
    public BlockAverageDataProcessor(int border, double smooth)
    {
        super(border);
        this.smooth = convert(smooth);
    }

    /**
     * Convert the smoothing parameter to the value which is used for the BlockMeanFilter.
     * We only use int smoothing. Values below zero are set to zero.
     *
     * @param smooth
     *            the smoothing parameter
     * @return The adjusted value
     * @see BlockMeanFilter
     */
    public static double convert(double smooth)
    {
        if (smooth < 0)
            return 0;
        return (int) smooth;
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.filters.DataProcessor#isWeighted()
     */
    @Override
    public boolean isWeighted()
    {
        return true;
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.filters.DataProcessor#setWeights(float[], int, int)
     */
    @Override
    public void setWeights(float[] weights, int width, int height)
    {
        if (smooth > 0)
            filter.setWeights(weights, width, height);
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.filters.DataProcessor#hasWeights()
     */
    @Override
    public boolean hasWeights()
    {
        return filter.hasWeights();
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.filters.DataProcessor#process(float[], int, int)
     */
    @Override
    public float[] process(float[] data, int width, int height)
    {
        float[] smoothData = data;
        if (smooth > 0)
        {
            // Smoothing destructively modifies the data so create a copy
            smoothData = Arrays.copyOf(data, width * height);

            // Check upper limits are safe
            final int tmpSmooth = FastMath.min((int) smooth, FastMath.min(width, height) / 2);

            if (tmpSmooth <= getBorder())
                filter.rollingBlockFilterInternal(smoothData, width, height, tmpSmooth);
            else
                filter.rollingBlockFilter(smoothData, width, height, tmpSmooth);
        }
        return smoothData;
    }

    /**
     * Gets the smooth.
     *
     * @return the smoothing width
     */
    public double getSmooth()
    {
        return smooth;
    }

    /*
     * (non-Javadoc)
     *
     * @see java.lang.Object#clone()
     */
    @Override
    public BlockAverageDataProcessor clone()
    {
        final BlockAverageDataProcessor f = (BlockAverageDataProcessor) super.clone();
        // Ensure the object is duplicated and not passed by reference.
        f.filter = filter.clone();
        return f;
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.filters.DataProcessor#getName()
     */
    @Override
    public String getName()
    {
        return "Block Average";
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.filters.DataProcessor#getParameters()
     */
    @Override
    public List<String> getParameters()
    {
        final List<String> list = super.getParameters();
        list.add("smooth = " + Utils.rounded(smooth));
        return list;
    }

    /*
     * (non-Javadoc)
     *
     * @see uk.ac.sussex.gdsc.smlm.filters.DataProcessor#getSpread()
     */
    @Override
    public double getSpread()
    {
        return 2 * smooth + 1;
    }
}
