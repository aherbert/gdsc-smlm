/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.ij.plugins;

/**
 * Define property constants for use with the ImageJ Prefs class.
 *
 * @see ij.Prefs Prefs
 */
class PrefsKey {
  // -=-=-=-
  // Used in the GaussianFit class
  // -=-=-=-

  /** gdsc.smlm.smooth */
  static final String smooth = "gdsc.smlm.smooth";
  /** gdsc.smlm.boxSize */
  static final String boxSize = "gdsc.smlm.boxSize";
  /** gdsc.smlm.boxSize2 */
  static final String boxSize2 = "gdsc.smlm.boxSize2";
  /** gdsc.smlm.background */
  static final String background = "gdsc.smlm.background";
  /** gdsc.smlm.peakHeight */
  static final String peakHeight = "gdsc.smlm.peakHeight";
  /** gdsc.smlm.fractionAboveBackground */
  static final String fractionAboveBackground = "gdsc.smlm.fractionAboveBackground";
  /** gdsc.smlm.peakWidth */
  static final String peakWidth = "gdsc.smlm.peakWidth";
  /** gdsc.smlm.topN */
  static final String topN = "gdsc.smlm.topN";
  /** gdsc.smlm.blockFindAlgorithm */
  static final String blockFindAlgorithm = "gdsc.smlm.blockFindAlgorithm";
  /** gdsc.smlm.neighbourCheck */
  static final String neighbourCheck = "gdsc.smlm.neighbourCheck";
  /** gdsc.smlm.border */
  static final String border = "gdsc.smlm.border";
  /** gdsc.smlm.fitFunction */
  static final String fitFunction = "gdsc.smlm.fitFunction";
  /** gdsc.smlm.fitBackground */
  static final String fitBackground = "gdsc.smlm.fitBackground";
  /** gdsc.smlm.fitCriteria */
  static final String fitCriteria = "gdsc.smlm.fitCriteria";
  /** gdsc.smlm.logProgress */
  static final String logProgress = "gdsc.smlm.logProgress";
  /** gdsc.smlm.maxIterations */
  static final String maxIterations = "gdsc.smlm.maxIterations";
  /** gdsc.smlm.relativeThreshold */
  static final String relativeThreshold = "gdsc.smlm.relativeThreshold";
  /** gdsc.smlm.absoluteThreshold */
  static final String absoluteThreshold = "gdsc.smlm.absoluteThreshold";
  /** gdsc.smlm.singleFit */
  static final String singleFit = "gdsc.smlm.singleFit";
  /** gdsc.smlm.singleRegionSize */
  static final String singleRegionSize = "gdsc.smlm.singleRegionSize";
  /** gdsc.smlm.initialPeakStdDev0 */
  static final String initialPeakStdDev0 = "gdsc.smlm.initialPeakStdDev0";
  /** gdsc.smlm.showDeviations */
  static final String showDeviations = "gdsc.smlm.showDeviations";
  /** gdsc.smlm.filterResults */
  static final String filterResults = "gdsc.smlm.filterResults";
  /** gdsc.smlm.showFit */
  static final String showFit = "gdsc.smlm.showFit";
  // -=-=-=-

  /** gdsc.smlm.algorithm */
  static final String algorithm = "gdsc.smlm.algorithm";

  /** gdsc.smlm.inputFilename */
  static final String inputFilename = "gdsc.smlm.inputFilename";

  /** gdsc.smlm.nmPerPixel */
  static final String inputNmPerPixel = "gdsc.smlm.nmPerPixel";
  /** gdsc.smlm.gain */
  static final String inputGain = "gdsc.smlm.gain";
  /** gdsc.smlm.exposureTime */
  static final String inputExposureTime = "gdsc.smlm.exposureTime";
  /** gdsc.smlm.noise */
  static final String inputNoise = "gdsc.smlm.noise";

  /** gdsc.smlm.tiffSeriesMode */
  static final String tiffSeriesMode = "gdsc.smlm.tiffSeriesMode";
  /** gdsc.smlm.tiffSeriesDirectory */
  static final String tiffSeriesDirectory = "gdsc.smlm.tiffSeriesDirectory";
  /** gdsc.smlm.tiffSeriesFile */
  static final String tiffSeriesFile = "gdsc.smlm.tiffSeriesFile";
  /** gdsc.smlm.tiffSeriesLogProgress */
  static final String tiffSeriesLogProgress = "gdsc.smlm.tiffSeriesLogProgress";
  /** gdsc.smlm.tiffSeriesOutputMode */
  static final String tiffSeriesOutputMode = "gdsc.smlm.tiffSeriesOutputMode";
  /** gdsc.smlm.tiffSeriesOutputNImages */
  static final String tiffSeriesOutputNImages = "gdsc.smlm.tiffSeriesOutputNImages";
  /** gdsc.smlm.tiffSeriesOutputDirectory */
  static final String tiffSeriesOutputDirectory = "gdsc.smlm.tiffSeriesOutputDirectory";

  /** gdsc.smlm.sCMOSAnalysisDirectory */
  static final String sCMOSAnalysisDirectory = "gdsc.smlm.sCMOSAnalysisDirectory";
}
