package gdsc.smlm.ij.settings;

import gdsc.smlm.results.TraceManager;
import gdsc.smlm.results.TraceManager.TraceMode;
import gdsc.smlm.results.clustering.ClusteringAlgorithm;

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
	public enum OptimiserPlot
	{
		NONE("None"), NEAREST_NEIGHBOUR("Nearest neighbour"), BILINEAR("Bi-linear");

		private String name;

		private OptimiserPlot(String name)
		{
			this.name = name;
		}

		@Override
		public String toString()
		{
			return name;
		}
	}
	
	public double distanceThreshold = 50;
	public double distanceExclusion = 0;
	public double timeThreshold = 5;
	private TraceManager.TraceMode traceMode = TraceMode.LATEST_FORERUNNER; 
	private ClusteringAlgorithm clusteringAlgorithm = ClusteringAlgorithm.PAIRWISE; 
	public int pulseInterval = 0;
	public int pulseWindow = 0;
	public boolean splitPulses = false;
	public double blinkingRate = 1;
	public boolean optimise = false;
	public double minDistanceThreshold = 0;
	public double maxDistanceThreshold = 500;
	public double minTimeThreshold = 0;
	public double maxTimeThreshold = 20;
	public int optimiserSteps = 10;
	private OptimiserPlot optimiserPlot = OptimiserPlot.BILINEAR;
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
	
	public OptimiserPlot getOptimiserPlot()
	{
		if (optimiserPlot == null)
			optimiserPlot = OptimiserPlot.NONE;
		return optimiserPlot;
	}
	
	public void setOptimiserPlot(OptimiserPlot optimiserPlot)
	{
		this.optimiserPlot = optimiserPlot;
	}
	
	public void setOptimiserPlot(int optimiserPlot)
	{
		if (optimiserPlot < 0 || optimiserPlot >= OptimiserPlot.values().length)
			this.optimiserPlot = OptimiserPlot.NONE;
		else
			this.optimiserPlot = OptimiserPlot.values()[optimiserPlot];
	}
	
	public TraceMode getTraceMode()
	{
		if (traceMode == null)
			traceMode = TraceMode.LATEST_FORERUNNER;
		return traceMode;
	}
	
	public void setTraceMode(TraceMode traceMode)
	{
		this.traceMode = traceMode;
	}
	
	public void setTraceMode(int traceMode)
	{
		if (traceMode < 0 || traceMode >= TraceMode.values().length)
			this.traceMode = TraceMode.LATEST_FORERUNNER;
		else
			this.traceMode = TraceMode.values()[traceMode];
	}
	
	public ClusteringAlgorithm getClusteringAlgorithm()
	{
		if (clusteringAlgorithm == null)
			clusteringAlgorithm = ClusteringAlgorithm.PAIRWISE;
		return clusteringAlgorithm;
	}
	
	public void setClusteringAlgorithm(ClusteringAlgorithm clusteringAlgorithm)
	{
		this.clusteringAlgorithm = clusteringAlgorithm;
	}
	
	public void setClusteringAlgorithm(int clusteringAlgorithm)
	{
		if (clusteringAlgorithm < 0 || clusteringAlgorithm >= ClusteringAlgorithm.values().length)
			this.clusteringAlgorithm = ClusteringAlgorithm.PAIRWISE;
		else
			this.clusteringAlgorithm = ClusteringAlgorithm.values()[clusteringAlgorithm];
	}
}
