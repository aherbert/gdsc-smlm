package gdsc.smlm.engine;

import gdsc.core.threshold.AutoThreshold;
import gdsc.core.threshold.AutoThreshold.Method;
import gdsc.core.threshold.FloatHistogram;
import gdsc.core.threshold.Histogram;
import gdsc.core.utils.NoiseEstimator;
import gdsc.core.utils.Statistics;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Provide methods to estimate parameters of the data. Data is partitioned into foreground and background using
 * thresholding. The background region must be a minimum fraction of the total data. If this is not achieved then the
 * estimates are made using all the data.
 */
public class DataEstimator
{
	final private float[] data;
	private Histogram h;
	final int width, height;
	private float fraction = 0.25f;
	private int histogramSize = 2048;
	private AutoThreshold.Method thresholdMethod = Method.DEFAULT;
	private float[] estimate = null;

	/**
	 * Create a new DataEstimator
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
	 * Gets a clone of the data
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
		return estimate[0] == 1;
	}

	/**
	 * Gets the background. This is the mean of the data in the background region.
	 *
	 * @return the background
	 */
	public float getBackground()
	{
		getEstimate();
		return estimate[1];
	}

	/**
	 * Gets the noise. This is the standard deviation of the data in the background region.
	 *
	 * @return the noise
	 */
	public float getNoise()
	{
		getEstimate();
		return estimate[2];
	}

	/**
	 * Gets the threshold for the background region.
	 *
	 * @return the noise
	 */
	public float getThreshold()
	{
		getEstimate();
		return estimate[3];
	}

	private void getEstimate()
	{
		if (estimate == null)
		{
			estimate = new float[4];

			if (h == null)
			{
				h = FloatHistogram.buildHistogram(data.clone(), true);
				h = h.compact(histogramSize);
			}

			// Threshold the data
			final float t = estimate[3] = h.getAutoThreshold(thresholdMethod);

			// Get stats below the threshold
			Statistics stats = new Statistics();
			for (int i = h.minBin; i <= h.maxBin; i++)
			{
				if (h.getValue(i) >= t)
					break;
				stats.add(h.h[i], h.getValue(i));
			}

			// Check if background region is large enough
			if (stats.getN() > fraction * data.length)
			{
				// Background region is large enough
				estimate[0] = 1;
			}
			else
			{
				// Recompute with all the data
				stats = new Statistics(data);
			}

			// Background
			estimate[1] = (float) stats.getMean();
			// Noise
			estimate[2] = (float) stats.getStandardDeviation();
		}
	}

	/**
	 * Estimate the noise in the all the data
	 *
	 * @param method
	 *            the method
	 * @return the noise
	 */
	public float getNoise(NoiseEstimator.Method method)
	{
		NoiseEstimator ne = new NoiseEstimator(data, width, height);
		return (float) ne.getNoise(method);
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
	 * SGet the size of the histogram used to compute the threshold
	 * 
	 * @return the histogram size
	 */
	public int getHistogramSize()
	{
		return histogramSize;
	}

	/**
	 * Set the size of the histogram used to compute the threshold
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