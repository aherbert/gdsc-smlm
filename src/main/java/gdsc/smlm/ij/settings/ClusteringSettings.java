package gdsc.smlm.ij.settings;

import gdsc.core.clustering.ClusteringAlgorithm;
import gdsc.smlm.data.config.UnitConfig.TimeUnit;
import gdsc.smlm.ij.settings.ClusteringSettingsHelper.OptimiserPlot;
import gdsc.smlm.results.TraceManager;
import gdsc.smlm.results.TraceManager.TraceMode;

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

/**
 * Contain the settings for the clustering algorithm
 */
public class ClusteringSettings
{
	public double distanceThreshold = 50;
	public double distanceExclusion = 0;
	private double timeThreshold = 1;
	private TimeUnit timeUnit = TimeUnit.SECOND;
	public int traceMode = TraceMode.LATEST_FORERUNNER.ordinal();
	public int clusteringAlgorithm = ClusteringAlgorithm.PAIRWISE.ordinal();
	public int pulseInterval = 0;
	public int pulseWindow = 0;
	public boolean splitPulses = false;
	public double blinkingRate = 1;
	public boolean optimise = false;
	public double minDistanceThreshold = 0;
	public double maxDistanceThreshold = 500;
	/** The min time threshold for optimisation (time is in frames). */
	public int minTimeThreshold = 0;
	/** The max time threshold for optimisation (time is in frames). */
	public int maxTimeThreshold = 20;
	public int optimiserSteps = 10;
	public int optimiserPlot = OptimiserPlot.BILINEAR.ordinal();
	public boolean saveTraces = false;
	public boolean showHistograms = false;
	public boolean saveTraceData = false;
	public String traceDataDirectory = "";
	public int histogramBins = 100;
	public boolean removeOutliers = false;
	public boolean refitOption = false;
	// Options for tracing diffusion
	public int minimumTraceLength = 6;
	public boolean truncate = false;
	public boolean internalDistances = true;
	public boolean subSampledDistances = false;
	public boolean ignoreEnds = true;
	public boolean precisionCorrection = true;
	public boolean msdCorrection = true;
	public boolean mle = true;
	public int fitLength = 6;
	public int fitRestarts = 3;
	public int jumpDistance = 1;

	public int getOptimiserPlot()
	{
		return optimiserPlot;
	}

	public void setOptimiserPlot(int optimiserPlot)
	{
		this.optimiserPlot = optimiserPlot;
	}

	public int getTraceMode()
	{
		return traceMode;
	}

	public void setTraceMode(int traceMode)
	{
		this.traceMode = traceMode;
	}

	public TimeUnit getTimeUnit()
	{
		if (timeUnit == null)
			timeUnit = TimeUnit.SECOND;
		return timeUnit;
	}

	public void setTimeUnit(TimeUnit timeUnit)
	{
		this.timeUnit = timeUnit;
	}

	public void setTimeUnit(int timeUnit)
	{
		if (timeUnit < 0 || timeUnit >= TimeUnit.values().length)
			this.timeUnit = TimeUnit.SECOND;
		else
			this.timeUnit = TimeUnit.values()[timeUnit];
	}

	public int getClusteringAlgorithm()
	{
		return clusteringAlgorithm;
	}

	public void setClusteringAlgorithm(int clusteringAlgorithm)
	{
		this.clusteringAlgorithm = clusteringAlgorithm;
	}

	/**
	 * Gets the time threshold.
	 * <p>
	 * The units are specified in {@link #getTimeUnit()}
	 *
	 * @return the time threshold
	 */
	public double getTimeThreshold()
	{
		return timeThreshold;
	}

	/**
	 * Sets the time threshold.
	 * <p>
	 * The units are specified in {@link #getTimeUnit()}
	 *
	 * @param timeThreshold
	 *            the new time threshold
	 */
	public void setTimeThreshold(double timeThreshold)
	{
		this.timeThreshold = timeThreshold;
	}
}
