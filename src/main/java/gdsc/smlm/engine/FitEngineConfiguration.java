package gdsc.smlm.engine;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import gdsc.smlm.fitting.FitConfiguration;
import gdsc.smlm.utils.NoiseEstimator;
import gdsc.smlm.utils.NoiseEstimator.Method;

/**
 * Specifies the configuration for the fit engine
 */
public class FitEngineConfiguration implements Cloneable
{
	private FitConfiguration fitConfiguration;

	private double smooth = 0.5;
	private double smooth2 = 3;
	private double search = 1;
	private double fitting = 3;
	private int failuresLimit = 3;
	private boolean includeNeighbours = true;
	private double neighbourHeightThreshold = 0.3;
	private double residualsThreshold = 1;
	private NoiseEstimator.Method noiseMethod = Method.QUICK_RESIDUALS_LEAST_MEAN_OF_SQUARES;
	private DataFilter dataFilter = DataFilter.MEAN;

	/**
	 * Constructor
	 * 
	 * @param fitConfiguration
	 */
	public FitEngineConfiguration(FitConfiguration fitConfiguration)
	{
		this.fitConfiguration = fitConfiguration;
	}

	/**
	 * @return the smoothing window size
	 */
	public double getSmooth()
	{
		return smooth;
	}

	/**
	 * @param smooth
	 *            the size of the smoothing window. The actual window is calculated dynamically in conjunction with the
	 *            peak widths.
	 */
	public void setSmooth(double smooth)
	{
		this.smooth = smooth;
	}

	/**
	 * @return the second smoothing window size
	 */
	public double getSmooth2()
	{
		return smooth2;
	}

	/**
	 * Use this parameter to perform a difference-of-smoothing (Top-hat) filter to identify peaks. The second smoothing
	 * window should be larger than the first. The second smoothed image is subtracted from the first to create a
	 * difference-of-smoothing such as a difference-of-Gaussians band-pass filter.
	 * 
	 * @param smooth2
	 *            the size of the second smoothing window. The actual window is calculated dynamically in conjunction
	 *            with the peak widths.
	 */
	public void setSmooth2(double smooth2)
	{
		this.smooth2 = smooth2;
	}

	/**
	 * @return the size of the region to search for local maxima
	 */
	public double getSearch()
	{
		return search;
	}

	/**
	 * @param search
	 *            the size of the region to search for local maxima. The actual window is calculated dynamically
	 *            in conjunction with the peak widths.
	 */
	public void setSearch(double search)
	{
		this.search = search;
	}

	/**
	 * @return the fitting window size
	 */
	public double getFitting()
	{
		return fitting;
	}

	/**
	 * @param fitting
	 *            the size of the window used for fitting around a maxima. The actual window is calculated dynamically
	 *            in conjunction with the peak widths.
	 */
	public void setFitting(double fitting)
	{
		this.fitting = fitting;
	}

	/**
	 * @return the failuresLimit
	 */
	public int getFailuresLimit()
	{
		return failuresLimit;
	}

	/**
	 * @param failuresLimit
	 *            the number of consecutive failures that stops the fitting process on the frame
	 */
	public void setFailuresLimit(int failuresLimit)
	{
		this.failuresLimit = failuresLimit;
	}

	/**
	 * @return the fitConfiguration
	 */
	public FitConfiguration getFitConfiguration()
	{
		return fitConfiguration;
	}

	/**
	 * @return the includeNeighbours
	 */
	public boolean isIncludeNeighbours()
	{
		return includeNeighbours;
	}

	/**
	 * Include neighbour maxima in the fitting.
	 * <p>
	 * Use this option when the fitting search region is large relative the the smoothing, thus other peaks may be
	 * within the region used for fitting.
	 * 
	 * @param includeNeighbours
	 */
	public void setIncludeNeighbours(boolean includeNeighbours)
	{
		this.includeNeighbours = includeNeighbours;
	}

	/**
	 * @return the neighbourHeightThreshold
	 */
	public double getNeighbourHeightThreshold()
	{
		return neighbourHeightThreshold;
	}

	/**
	 * @param neighbourHeightThreshold
	 *            Set the height threshold that determines if a neighbour peak should be fitted (specified as a fraction
	 *            of the central peak relative to the background)
	 */
	public void setNeighbourHeightThreshold(double neighbourHeightThreshold)
	{
		this.neighbourHeightThreshold = neighbourHeightThreshold;
	}

	/**
	 * @return the residuals threshold
	 */
	public double getResidualsThreshold()
	{
		return residualsThreshold;
	}

	/**
	 * @param residualsThreshold
	 *            Set the threshold for the residuals analysis that determines if a two-kernel model should be fitted
	 */
	public void setResidualsThreshold(double residualsThreshold)
	{
		this.residualsThreshold = residualsThreshold;
	}

	/**
	 * @return the method used to estimate the image noise
	 */
	public NoiseEstimator.Method getNoiseMethod()
	{
		return noiseMethod;
	}

	/**
	 * @param noiseMethod
	 *            Set the method used to estimate the image noise
	 */
	public void setNoiseMethod(NoiseEstimator.Method noiseMethod)
	{
		this.noiseMethod = noiseMethod;
	}

	/**
	 * @param noiseMethod
	 *            Set the method used to estimate the image noise
	 */
	public void setNoiseMethod(int noiseMethod)
	{
		if (noiseMethod >= 0 && noiseMethod < NoiseEstimator.Method.values().length)
		{
			setNoiseMethod(NoiseEstimator.Method.values()[noiseMethod]);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	public Object clone()
	{
		try
		{
			FitEngineConfiguration f = (FitEngineConfiguration) super.clone();
			// Ensure the object is duplicated and not passed by reference.
			f.fitConfiguration = (FitConfiguration) fitConfiguration.clone();
			return f;
		}
		catch (CloneNotSupportedException e)
		{
			// Ignore
		}
		return null;
	}

	/**
	 * Ensure that the internal state of the object is initialised. This is used after deserialisation since some state
	 * is not saved but restored from other property values.
	 */
	public void initialiseState()
	{
		if (fitConfiguration == null)
			fitConfiguration = new FitConfiguration();
		fitConfiguration.initialiseState();
		if (noiseMethod == null)
			noiseMethod = Method.QUICK_RESIDUALS_LEAST_MEAN_OF_SQUARES;
		if (dataFilter == null)
			dataFilter = DataFilter.MEAN;
	}

	/**
	 * @return the filter to apply to the data before identifying local maxima
	 */
	public DataFilter getDataFilter()
	{
		return dataFilter;
	}

	/**
	 * @param DataFilter
	 *            the filter to apply to the data before identifying local maxima
	 */
	public void setDataFilter(DataFilter dataFilter)
	{
		this.dataFilter = dataFilter;
	}

	/**
	 * @param DataFilter
	 *            the filter to apply to the data before identifying local maxima
	 */
	public void setDataFilter(int dataFilter)
	{
		if (dataFilter >= 0 && dataFilter < DataFilter.values().length)
		{
			setDataFilter(DataFilter.values()[dataFilter]);
		}
	}
}