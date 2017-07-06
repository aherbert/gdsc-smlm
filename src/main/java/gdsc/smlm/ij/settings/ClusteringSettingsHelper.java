package gdsc.smlm.ij.settings;

import gdsc.core.clustering.ClusteringAlgorithm;
import gdsc.smlm.data.config.UnitConfig.TimeUnit;
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
public class ClusteringSettingsHelper
{
	public enum OptimiserPlot
	{
		//@formatter:off
		NONE{ public String getName() { return "None"; }}, 
		NEAREST_NEIGHBOUR{ public String getName() { return "Nearest neighbour"; }}, 
		BILINEAR{ public String getName() { return "Bi-linear"; }};
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

	public double distanceThreshold = 50;
	public double distanceExclusion = 0;
	public double timeThreshold = 1;
	public TimeUnit timeUnit = TimeUnit.SECOND;
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

	public static OptimiserPlot getOptimiserPlot(int optimiserPlot)
	{
		if (optimiserPlot < 0 || optimiserPlot >= OptimiserPlot.values().length)
			return OptimiserPlot.NONE;
		return OptimiserPlot.values()[optimiserPlot];
	}

	public static TraceMode getTraceMode(int traceMode)
	{
		if (traceMode < 0 || traceMode >= TraceMode.values().length)
			return TraceMode.LATEST_FORERUNNER;
		return TraceMode.values()[traceMode];
	}

	public static ClusteringAlgorithm getClusteringAlgorithm(int clusteringAlgorithm)
	{
		if (clusteringAlgorithm < 0 || clusteringAlgorithm >= ClusteringAlgorithm.values().length)
			return ClusteringAlgorithm.PAIRWISE;
		return ClusteringAlgorithm.values()[clusteringAlgorithm];
	}
}
