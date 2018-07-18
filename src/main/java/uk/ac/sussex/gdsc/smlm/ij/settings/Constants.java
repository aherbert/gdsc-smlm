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
package uk.ac.sussex.gdsc.smlm.ij.settings;

/**
 * Define property constants. Enables plugins to share their properties using the ImageJ Prefs class.
 */
public class Constants
{
    /** uk.ac.sussex.gdsc.smlm.settingsFilename */
    public static final String settingsFilename = "uk.ac.sussex.gdsc.smlm.settingsFilename";
    /** uk.ac.sussex.gdsc.smlm.settingsDirectory */
    public static final String settingsDirectory = "uk.ac.sussex.gdsc.smlm.settingsDirectory";

    // -=-=-=-
    // Used in the GaussianFit class
    // -=-=-=-

    /** uk.ac.sussex.gdsc.smlm.smooth */
    public static final String smooth = "uk.ac.sussex.gdsc.smlm.smooth";
    /** uk.ac.sussex.gdsc.smlm.boxSize */
    public static final String boxSize = "uk.ac.sussex.gdsc.smlm.boxSize";
    /** uk.ac.sussex.gdsc.smlm.boxSize2 */
    public static final String boxSize2 = "uk.ac.sussex.gdsc.smlm.boxSize2";
    /** uk.ac.sussex.gdsc.smlm.background */
    public static final String background = "uk.ac.sussex.gdsc.smlm.background";
    /** uk.ac.sussex.gdsc.smlm.peakHeight */
    public static final String peakHeight = "uk.ac.sussex.gdsc.smlm.peakHeight";
    /** uk.ac.sussex.gdsc.smlm.fractionAboveBackground */
    public static final String fractionAboveBackground = "uk.ac.sussex.gdsc.smlm.fractionAboveBackground";
    /** uk.ac.sussex.gdsc.smlm.peakWidth */
    public static final String peakWidth = "uk.ac.sussex.gdsc.smlm.peakWidth";
    /** uk.ac.sussex.gdsc.smlm.topN */
    public static final String topN = "uk.ac.sussex.gdsc.smlm.topN";
    /** uk.ac.sussex.gdsc.smlm.blockFindAlgorithm */
    public static final String blockFindAlgorithm = "uk.ac.sussex.gdsc.smlm.blockFindAlgorithm";
    /** uk.ac.sussex.gdsc.smlm.neighbourCheck */
    public static final String neighbourCheck = "uk.ac.sussex.gdsc.smlm.neighbourCheck";
    /** uk.ac.sussex.gdsc.smlm.border */
    public static final String border = "uk.ac.sussex.gdsc.smlm.border";
    /** uk.ac.sussex.gdsc.smlm.fitFunction */
    public static final String fitFunction = "uk.ac.sussex.gdsc.smlm.fitFunction";
    /** uk.ac.sussex.gdsc.smlm.fitBackground */
    public static final String fitBackground = "uk.ac.sussex.gdsc.smlm.fitBackground";
    /** uk.ac.sussex.gdsc.smlm.fitCriteria */
    public static final String fitCriteria = "uk.ac.sussex.gdsc.smlm.fitCriteria";
    /** uk.ac.sussex.gdsc.smlm.logProgress */
    public static final String logProgress = "uk.ac.sussex.gdsc.smlm.logProgress";
    /** uk.ac.sussex.gdsc.smlm.maxIterations */
    public static final String maxIterations = "uk.ac.sussex.gdsc.smlm.maxIterations";
    /** uk.ac.sussex.gdsc.smlm.relativeThreshold */
    public static final String relativeThreshold = "uk.ac.sussex.gdsc.smlm.relativeThreshold";
    /** uk.ac.sussex.gdsc.smlm.absoluteThreshold */
    public static final String absoluteThreshold = "uk.ac.sussex.gdsc.smlm.absoluteThreshold";
    /** uk.ac.sussex.gdsc.smlm.singleFit */
    public static final String singleFit = "uk.ac.sussex.gdsc.smlm.singleFit";
    /** uk.ac.sussex.gdsc.smlm.singleRegionSize */
    public static final String singleRegionSize = "uk.ac.sussex.gdsc.smlm.singleRegionSize";
    /** uk.ac.sussex.gdsc.smlm.initialPeakStdDev0 */
    public static final String initialPeakStdDev0 = "uk.ac.sussex.gdsc.smlm.initialPeakStdDev0";
    /** uk.ac.sussex.gdsc.smlm.showDeviations */
    public static final String showDeviations = "uk.ac.sussex.gdsc.smlm.showDeviations";
    /** uk.ac.sussex.gdsc.smlm.filterResults */
    public static final String filterResults = "uk.ac.sussex.gdsc.smlm.filterResults";
    /** uk.ac.sussex.gdsc.smlm.showFit */
    public static final String showFit = "uk.ac.sussex.gdsc.smlm.showFit";
    // -=-=-=-

    /** uk.ac.sussex.gdsc.smlm.algorithm */
    public static final String algorithm = "uk.ac.sussex.gdsc.smlm.algorithm";

    /** uk.ac.sussex.gdsc.smlm.inputFilename */
    public static final String inputFilename = "uk.ac.sussex.gdsc.smlm.inputFilename";

    /** uk.ac.sussex.gdsc.smlm.nmPerPixel */
    public static final String inputNmPerPixel = "uk.ac.sussex.gdsc.smlm.nmPerPixel";
    /** uk.ac.sussex.gdsc.smlm.gain */
    public static final String inputGain = "uk.ac.sussex.gdsc.smlm.gain";
    /** uk.ac.sussex.gdsc.smlm.exposureTime */
    public static final String inputExposureTime = "uk.ac.sussex.gdsc.smlm.exposureTime";
    /** uk.ac.sussex.gdsc.smlm.noise */
    public static final String inputNoise = "uk.ac.sussex.gdsc.smlm.noise";

    /** uk.ac.sussex.gdsc.smlm.tiffSeriesMode */
    public static final String tiffSeriesMode = "uk.ac.sussex.gdsc.smlm.tiffSeriesMode";
    /** uk.ac.sussex.gdsc.smlm.tiffSeriesDirectory */
    public static final String tiffSeriesDirectory = "uk.ac.sussex.gdsc.smlm.tiffSeriesDirectory";
    /** uk.ac.sussex.gdsc.smlm.tiffSeriesFile */
    public static final String tiffSeriesFile = "uk.ac.sussex.gdsc.smlm.tiffSeriesFile";
    /** uk.ac.sussex.gdsc.smlm.tiffSeriesLogProgress */
    public static final String tiffSeriesLogProgress = "uk.ac.sussex.gdsc.smlm.tiffSeriesLogProgress";
    /** uk.ac.sussex.gdsc.smlm.tiffSeriesOutputMode */
    public static final String tiffSeriesOutputMode = "uk.ac.sussex.gdsc.smlm.tiffSeriesOutputMode";
    /** uk.ac.sussex.gdsc.smlm.tiffSeriesOutputNImages */
    public static final String tiffSeriesOutputNImages = "uk.ac.sussex.gdsc.smlm.tiffSeriesOutputNImages";
    /** uk.ac.sussex.gdsc.smlm.tiffSeriesOutputDirectory */
    public static final String tiffSeriesOutputDirectory = "uk.ac.sussex.gdsc.smlm.tiffSeriesOutputDirectory";

    /** uk.ac.sussex.gdsc.smlm.sCMOSAnalysisDirectory */
    public static final String sCMOSAnalysisDirectory = "uk.ac.sussex.gdsc.smlm.sCMOSAnalysisDirectory";
}