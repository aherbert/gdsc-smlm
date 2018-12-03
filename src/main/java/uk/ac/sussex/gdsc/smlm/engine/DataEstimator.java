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
package uk.ac.sussex.gdsc.smlm.engine;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.stat.ranking.NaNStrategy;

import uk.ac.sussex.gdsc.core.threshold.AutoThreshold;
import uk.ac.sussex.gdsc.core.threshold.AutoThreshold.Method;
import uk.ac.sussex.gdsc.core.threshold.FloatHistogram;
import uk.ac.sussex.gdsc.core.threshold.Histogram;
import uk.ac.sussex.gdsc.core.utils.NoiseEstimator;
import uk.ac.sussex.gdsc.core.utils.Statistics;

/**
 * Provide methods to estimate parameters of the data. Data is partitioned into foreground and background using
 * thresholding. The background region must be a minimum fraction of the total data. If this is not achieved then the
 * estimates are made using all the data.
 */
public class DataEstimator
{
    private final int ESTIMATE_LARGE_ENOUGH = 0;
    private final int ESTIMATE_BACKGROUND = 1;
    private final int ESTIMATE_NOISE = 2;
    private final int ESTIMATE_THRESHOLD = 3;
    private final int ESTIMATE_BACKGROUND_SIZE = 4;

    private final float[] data;
    private Histogram h;
    private final int width, height;
    private float fraction = 0.25f;
    private int histogramSize = 2048;
    private AutoThreshold.Method thresholdMethod = Method.DEFAULT;
    private float[] estimate = null;

    /**
     * Create a new DataEstimator.
     *
     * @param data
     *            The data
     * @param width
     *            The width of the data
     * @param height
     *            The height of the data
     */
    public DataEstimator(float[] data, int width, int height)
    {
        if (data == null)
            throw new IllegalArgumentException("Input data must not be null");
        if (data.length < width * height)
            throw new IllegalArgumentException("Input data must not be smaller than width * height");
        this.data = data;
        this.width = width;
        this.height = height;
    }

    /**
     * Gets a clone of the data.
     *
     * @return the data (or null)
     */
    public float[] getData()
    {
        return data.clone();
    }

    /**
     * Checks if the background region is large enough to produce estimates.
     *
     * @return true, if there is a background
     */
    public boolean isBackgroundRegion()
    {
        getEstimate();
        return estimate[ESTIMATE_LARGE_ENOUGH] == 1;
    }

    /**
     * Gets the background. This is the mean of the data in the background region.
     *
     * @return the background
     */
    public float getBackground()
    {
        getEstimate();
        return estimate[ESTIMATE_BACKGROUND];
    }

    /**
     * Gets the noise. This is the standard deviation of the data in the background region.
     *
     * @return the noise
     */
    public float getNoise()
    {
        getEstimate();
        return estimate[ESTIMATE_NOISE];
    }

    /**
     * Gets the threshold for the background region.
     *
     * @return the noise
     */
    public float getThreshold()
    {
        getEstimate();
        return estimate[ESTIMATE_THRESHOLD];
    }

    /**
     * Gets the size of the background region.
     *
     * @return the size
     */
    public float getBackgroundSize()
    {
        getEstimate();
        return estimate[ESTIMATE_BACKGROUND_SIZE];
    }

    private void getEstimate()
    {
        if (estimate == null)
        {
            estimate = new float[5];

            if (h == null)
            {
                h = FloatHistogram.buildHistogram(data.clone(), true);
                h = h.compact(histogramSize);
            }

            // Threshold the data
            final float t = estimate[ESTIMATE_THRESHOLD] = h.getAutoThreshold(thresholdMethod);

            // Get stats below the threshold
            Statistics stats = new Statistics();
            for (int i = h.minBin; i <= h.maxBin; i++)
            {
                if (h.getValue(i) >= t)
                    break;
                stats.add(h.histogramCounts[i], h.getValue(i));
            }

            // Check if background region is large enough
            estimate[ESTIMATE_BACKGROUND_SIZE] = stats.getN();
            if (stats.getN() > fraction * data.length)
                // Background region is large enough
                estimate[ESTIMATE_LARGE_ENOUGH] = 1;
            else
                // Recompute with all the data
                stats = Statistics.create(data);

            estimate[ESTIMATE_BACKGROUND] = (float) stats.getMean();
            estimate[ESTIMATE_NOISE] = (float) stats.getStandardDeviation();
        }
    }

    /**
     * Estimate the noise in the all the data.
     *
     * @param method
     *            the method
     * @return the noise
     */
    public float getNoise(NoiseEstimator.Method method)
    {
        final NoiseEstimator ne = NoiseEstimator.wrap(data, width, height);
        return (float) ne.getNoise(method);
    }

    /**
     * Get the percentile value of the data.
     *
     * @param percentile
     *            The percentile
     * @return the percentile value
     */
    public float getPercentile(double percentile)
    {
        // Check the input
        if (percentile <= 0)
            percentile = Double.MIN_NORMAL;
        if (percentile > 100)
            percentile = 100;

        // The data should not have NaN so we ignore them for speed.
        final Percentile p = new Percentile(percentile).withNaNStrategy(NaNStrategy.FIXED);
        final int size = width * height;
        final double[] values = new double[size];
        for (int i = 0; i < size; i++)
            values[i] = data[i];
        return (float) p.evaluate(values);
    }

    /**
     * Gets the fraction of the data size that the background region must achieve to be used.
     *
     * @return the fraction
     */
    public float getFraction()
    {
        return fraction;
    }

    /**
     * Sets the fraction of the data size that the background region must achieve to be used.
     *
     * @param fraction
     *            the new fraction
     */
    public void setFraction(float fraction)
    {
        this.fraction = fraction;
        this.estimate = null;
    }

    /**
     * Gets the threshold method.
     *
     * @return the threshold method
     */
    public AutoThreshold.Method getThresholdMethod()
    {
        return thresholdMethod;
    }

    /**
     * Sets the threshold method.
     *
     * @param thresholdMethod
     *            the new threshold method
     */
    public void setThresholdMethod(AutoThreshold.Method thresholdMethod)
    {
        this.thresholdMethod = thresholdMethod;
        this.estimate = null;
    }

    /**
     * SGet the size of the histogram used to compute the threshold.
     *
     * @return the histogram size
     */
    public int getHistogramSize()
    {
        return histogramSize;
    }

    /**
     * Set the size of the histogram used to compute the threshold.
     *
     * @param histogramSize
     *            the histogram size
     */
    public void setHistogramSize(int histogramSize)
    {
        this.histogramSize = histogramSize;
        this.estimate = null;
        this.h = null;
    }
}
