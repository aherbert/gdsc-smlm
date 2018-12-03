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
 * Define property constants. Enables plugins to share their properties using the ImageJ Prefs
 * class.
 */
public class Constants {
  /** gdsc.smlm.settingsFilename */
  public static final String settingsFilename = "gdsc.smlm.settingsFilename";
  /** gdsc.smlm.settingsDirectory */
  public static final String settingsDirectory = "gdsc.smlm.settingsDirectory";

  // -=-=-=-
  // Used in the GaussianFit class
  // -=-=-=-

  /** gdsc.smlm.smooth */
  public static final String smooth = "gdsc.smlm.smooth";
  /** gdsc.smlm.boxSize */
  public static final String boxSize = "gdsc.smlm.boxSize";
  /** gdsc.smlm.boxSize2 */
  public static final String boxSize2 = "gdsc.smlm.boxSize2";
  /** gdsc.smlm.background */
  public static final String background = "gdsc.smlm.background";
  /** gdsc.smlm.peakHeight */
  public static final String peakHeight = "gdsc.smlm.peakHeight";
  /** gdsc.smlm.fractionAboveBackground */
  public static final String fractionAboveBackground = "gdsc.smlm.fractionAboveBackground";
  /** gdsc.smlm.peakWidth */
  public static final String peakWidth = "gdsc.smlm.peakWidth";
  /** gdsc.smlm.topN */
  public static final String topN = "gdsc.smlm.topN";
  /** gdsc.smlm.blockFindAlgorithm */
  public static final String blockFindAlgorithm = "gdsc.smlm.blockFindAlgorithm";
  /** gdsc.smlm.neighbourCheck */
  public static final String neighbourCheck = "gdsc.smlm.neighbourCheck";
  /** gdsc.smlm.border */
  public static final String border = "gdsc.smlm.border";
  /** gdsc.smlm.fitFunction */
  public static final String fitFunction = "gdsc.smlm.fitFunction";
  /** gdsc.smlm.fitBackground */
  public static final String fitBackground = "gdsc.smlm.fitBackground";
  /** gdsc.smlm.fitCriteria */
  public static final String fitCriteria = "gdsc.smlm.fitCriteria";
  /** gdsc.smlm.logProgress */
  public static final String logProgress = "gdsc.smlm.logProgress";
  /** gdsc.smlm.maxIterations */
  public static final String maxIterations = "gdsc.smlm.maxIterations";
  /** gdsc.smlm.relativeThreshold */
  public static final String relativeThreshold = "gdsc.smlm.relativeThreshold";
  /** gdsc.smlm.absoluteThreshold */
  public static final String absoluteThreshold = "gdsc.smlm.absoluteThreshold";
  /** gdsc.smlm.singleFit */
  public static final String singleFit = "gdsc.smlm.singleFit";
  /** gdsc.smlm.singleRegionSize */
  public static final String singleRegionSize = "gdsc.smlm.singleRegionSize";
  /** gdsc.smlm.initialPeakStdDev0 */
  public static final String initialPeakStdDev0 = "gdsc.smlm.initialPeakStdDev0";
  /** gdsc.smlm.showDeviations */
  public static final String showDeviations = "gdsc.smlm.showDeviations";
  /** gdsc.smlm.filterResults */
  public static final String filterResults = "gdsc.smlm.filterResults";
  /** gdsc.smlm.showFit */
  public static final String showFit = "gdsc.smlm.showFit";
  // -=-=-=-

  /** gdsc.smlm.algorithm */
  public static final String algorithm = "gdsc.smlm.algorithm";

  /** gdsc.smlm.inputFilename */
  public static final String inputFilename = "gdsc.smlm.inputFilename";

  /** gdsc.smlm.nmPerPixel */
  public static final String inputNmPerPixel = "gdsc.smlm.nmPerPixel";
  /** gdsc.smlm.gain */
  public static final String inputGain = "gdsc.smlm.gain";
  /** gdsc.smlm.exposureTime */
  public static final String inputExposureTime = "gdsc.smlm.exposureTime";
  /** gdsc.smlm.noise */
  public static final String inputNoise = "gdsc.smlm.noise";

  /** gdsc.smlm.tiffSeriesMode */
  public static final String tiffSeriesMode = "gdsc.smlm.tiffSeriesMode";
  /** gdsc.smlm.tiffSeriesDirectory */
  public static final String tiffSeriesDirectory = "gdsc.smlm.tiffSeriesDirectory";
  /** gdsc.smlm.tiffSeriesFile */
  public static final String tiffSeriesFile = "gdsc.smlm.tiffSeriesFile";
  /** gdsc.smlm.tiffSeriesLogProgress */
  public static final String tiffSeriesLogProgress = "gdsc.smlm.tiffSeriesLogProgress";
  /** gdsc.smlm.tiffSeriesOutputMode */
  public static final String tiffSeriesOutputMode = "gdsc.smlm.tiffSeriesOutputMode";
  /** gdsc.smlm.tiffSeriesOutputNImages */
  public static final String tiffSeriesOutputNImages = "gdsc.smlm.tiffSeriesOutputNImages";
  /** gdsc.smlm.tiffSeriesOutputDirectory */
  public static final String tiffSeriesOutputDirectory = "gdsc.smlm.tiffSeriesOutputDirectory";

  /** gdsc.smlm.sCMOSAnalysisDirectory */
  public static final String sCMOSAnalysisDirectory = "gdsc.smlm.sCMOSAnalysisDirectory";
}
