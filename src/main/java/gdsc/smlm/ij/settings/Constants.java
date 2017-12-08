package gdsc.smlm.ij.settings;

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
 * Define property constants. Enables plugins to share their properties using the ImageJ Prefs class.
 */
public class Constants
{
	public static final String settingsFilename = "gdsc.smlm.settingsFilename";
	public static final String settingsDirectory = "gdsc.smlm.settingsDirectory";

	// -=-=-=- 
	// Used in the GaussianFit class
	// -=-=-=-
	public static final String smooth = "gdsc.smlm.smooth";
	public static final String boxSize = "gdsc.smlm.boxSize";
	public static final String boxSize2 = "gdsc.smlm.boxSize2";
	public static final String background = "gdsc.smlm.background";
	public static final String peakHeight = "gdsc.smlm.peakHeight";
	public static final String fractionAboveBackground = "gdsc.smlm.fractionAboveBackground";
	public static final String peakWidth = "gdsc.smlm.peakWidth";
	public static final String topN = "gdsc.smlm.topN";
	public static final String blockFindAlgorithm = "gdsc.smlm.blockFindAlgorithm";
	public static final String neighbourCheck = "gdsc.smlm.neighbourCheck";
	public static final String border = "gdsc.smlm.border";
	public static final String fitFunction = "gdsc.smlm.fitFunction";
	public static final String fitBackground = "gdsc.smlm.fitBackground";
	public static final String fitCriteria = "gdsc.smlm.fitCriteria";
	public static final String logProgress = "gdsc.smlm.logProgress";
	public static final String maxIterations = "gdsc.smlm.maxIterations";
	public static final String relativeThreshold = "gdsc.smlm.relativeThreshold";
	public static final String absoluteThreshold = "gdsc.smlm.absoluteThreshold";
	public static final String singleFit = "gdsc.smlm.singleFit";
	public static final String singleRegionSize = "gdsc.smlm.singleRegionSize";
	public static final String initialPeakStdDev0 = "gdsc.smlm.initialPeakStdDev0";
	public static final String showDeviations = "gdsc.smlm.showDeviations";
	public static final String filterResults = "gdsc.smlm.filterResults";
	public static final String showFit = "gdsc.smlm.showFit";
	// -=-=-=-

	public static final String algorithm = "gdsc.smlm.algorithm";

	public static final String inputFilename = "gdsc.smlm.inputFilename";
	
	public static final String inputNmPerPixel = "gdsc.smlm.nmPerPixel";
	public static final String inputGain = "gdsc.smlm.gain";
	public static final String inputExposureTime = "gdsc.smlm.exposureTime";
	public static final String inputNoise = "gdsc.smlm.noise";
	
	public static final String tiffSeriesMode = "gdsc.smlm.tiffSeriesMode";
	public static final String tiffSeriesDirectory = "gdsc.smlm.tiffSeriesDirectory";
	public static final String tiffSeriesFile = "gdsc.smlm.tiffSeriesFile";
	public static final String tiffSeriesLogProgress = "gdsc.smlm.tiffSeriesLogProgress";
	public static final String tiffSeriesOutputMode = "gdsc.smlm.tiffSeriesOutputMode";
	public static final String tiffSeriesOutputNImages = "gdsc.smlm.tiffSeriesOutputNImages";
	public static final String tiffSeriesOutputDirectory = "gdsc.smlm.tiffSeriesOutputDirectory";
	
	public static final String sCMOSAnalysisDirectory = "gdsc.smlm.sCMOSAnalysisDirectory";
}
