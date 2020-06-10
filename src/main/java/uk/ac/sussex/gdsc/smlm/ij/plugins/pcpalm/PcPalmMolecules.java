/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.ij.plugins.pcpalm;

import gnu.trove.list.array.TDoubleArrayList;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.plugin.frame.Recorder;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer.Optimum;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.InverseTransformDiscreteSampler;
import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import uk.ac.sussex.gdsc.core.annotation.Nullable;
import uk.ac.sussex.gdsc.core.clustering.Cluster;
import uk.ac.sussex.gdsc.core.clustering.ClusterPoint;
import uk.ac.sussex.gdsc.core.clustering.ClusteringAlgorithm;
import uk.ac.sussex.gdsc.core.clustering.ClusteringEngine;
import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.SimpleImageJTrackProgress;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.Plot2;
import uk.ac.sussex.gdsc.core.utils.DoubleData;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.StoredData;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.core.utils.rng.BinomialDiscreteInverseCumulativeProbabilityFunction;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;
import uk.ac.sussex.gdsc.core.utils.rng.UniformRandomProviders;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationHelper;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.data.config.PsfHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.function.SkewNormalFunction;
import uk.ac.sussex.gdsc.smlm.ij.plugins.HelpUrls;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ParameterUtils;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.plugins.SmlmUsageTracker;
import uk.ac.sussex.gdsc.smlm.model.MaskDistribution;
import uk.ac.sussex.gdsc.smlm.model.StandardFluorophoreSequenceModel;
import uk.ac.sussex.gdsc.smlm.model.UniformDistribution;
import uk.ac.sussex.gdsc.smlm.results.AttributePeakResult;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.NullSource;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.Trace;
import uk.ac.sussex.gdsc.smlm.results.TraceManager;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.PrecisionResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.StandardResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyrResultProcedure;

/**
 * Use the PC-PALM protocol to prepare a set of localisations into molecules. This can be used for
 * for clustering analysis.
 *
 * <p>See Sengupta, et al (2013). Quantifying spatial resolution in point-localisation
 * superresolution images using pair correlation analysis. Nature Protocols 8, pp345-354.
 *
 * <p>See also Veatch, et al (2012). Correlation Functions Quantify Super-Resolution Images and
 * Estimate Apparent Clustering Due to Over-Counting. PLoS One 7, Issue 2, e31457
 */
public class PcPalmMolecules implements PlugIn {
  /** The title. */
  static final String TITLE = "PC-PALM Molecules";

  private Rectangle roiBounds;
  private int roiImageWidth;
  private int roiImageHeight;
  private long start;

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    static final String[] RUN_MODE =
        {"PC-PALM", "Manual Tracing", "In-memory results", "Simulation"};
    static String[] singlesMode =
        new String[] {"Ignore", "Include in molecules histogram", "Include in final filtering"};
    static final String[] BLINKING_DISTRIBUTION =
        new String[] {"Poisson", "Geometric", "None", "Binomial"};
    static final String[] CLUSTER_SIMULATION =
        new String[] {"None", "Circles", "Non-overlapping circles", "Circles Mask"};

    String inputOption;
    boolean chooseRoi;
    double nmPerPixelLimit;

    int runMode;

    // Mode 0: PC-PALM protocol for estimating localisation precision and then tracing molecules
    String roiImage;
    int histogramBins;
    int singlesModeIndex;
    boolean simplexFitting;
    boolean showHistograms;
    boolean binaryImage;
    /** The blinking rate. */
    double blinkingRate;
    double pvalue;
    int blinkingDistribution;
    boolean clearResults;

    // Mode 1. Manual tracing of molecules
    double distanceThreshold;
    double timeThreshold;

    // Mode 2. Direct use of in-memory results
    // - No parameters needed

    // Mode 3. Random simulation of molecules
    int numberOfMolecules;
    double simulationSize;
    boolean distanceAnalysis;

    int clusterSimulation;
    double clusterNumber;
    double clusterNumberStdDev;
    double clusterRadius;
    boolean showClusterMask;

    // Low resolution image construction
    int lowResolutionImageSize;
    double roiSizeInUm;
    boolean showHighResolutionImage;

    /** The results. */
    MemoryPeakResults results;
    /** The minx. */
    double minx;
    /** The miny. */
    double miny;
    /** The maxx. */
    double maxx;
    /** The maxy. */
    double maxy;
    /** The nm per pixel. */
    double nmPerPixel;
    /** The molecules. */
    List<Molecule> molecules;
    /** The sigma S. */
    double sigmaS;
    /** The peak density. */
    double densityPeaks;
    /** The protein density. */
    double densityProtein;
    /** The seconds. */
    double seconds;
    /** The area. */
    double area;

    Settings() {
      // Set defaults
      inputOption = "";
      roiImage = "";
      histogramBins = 50;
      singlesModeIndex = 1;
      showHistograms = true;
      binaryImage = true;
      blinkingRate = 2;
      pvalue = 0.6;
      distanceThreshold = 150;
      timeThreshold = 1;
      numberOfMolecules = 2000;
      simulationSize = 16;
      clusterNumber = 3;
      clusterRadius = 50;
      lowResolutionImageSize = 1024;
      roiSizeInUm = 4;
      sigmaS = 20;
    }

    Settings(Settings source) {
      inputOption = source.inputOption;
      chooseRoi = source.chooseRoi;
      nmPerPixelLimit = source.nmPerPixelLimit;
      runMode = source.runMode;
      roiImage = source.roiImage;
      histogramBins = source.histogramBins;
      singlesModeIndex = source.singlesModeIndex;
      simplexFitting = source.simplexFitting;
      showHistograms = source.showHistograms;
      binaryImage = source.binaryImage;
      blinkingRate = source.blinkingRate;
      pvalue = source.pvalue;
      blinkingDistribution = source.blinkingDistribution;
      clearResults = source.clearResults;
      distanceThreshold = source.distanceThreshold;
      timeThreshold = source.timeThreshold;
      numberOfMolecules = source.numberOfMolecules;
      simulationSize = source.simulationSize;
      distanceAnalysis = source.distanceAnalysis;
      clusterSimulation = source.clusterSimulation;
      clusterNumber = source.clusterNumber;
      clusterNumberStdDev = source.clusterNumberStdDev;
      clusterRadius = source.clusterRadius;
      showClusterMask = source.showClusterMask;
      lowResolutionImageSize = source.lowResolutionImageSize;
      roiSizeInUm = source.roiSizeInUm;
      showHighResolutionImage = source.showHighResolutionImage;
      results = source.results;
      minx = source.minx;
      miny = source.miny;
      maxx = source.maxx;
      maxy = source.maxy;
      nmPerPixel = source.nmPerPixel;
      molecules = source.molecules;
      sigmaS = source.sigmaS;
      densityPeaks = source.densityPeaks;
      densityProtein = source.densityProtein;
      seconds = source.seconds;
      area = source.area;
    }

    Settings copy() {
      return new Settings(this);
    }

    /**
     * Load a copy of the settings.
     *
     * @return the settings
     */
    static Settings load() {
      return lastSettings.get().copy();
    }

    /**
     * Save the settings.
     */
    void save() {
      lastSettings.set(this);
    }
  }

  /**
   * Contains the results that are the used by other PC-PALM plugins.
   */
  static class MoleculesResults {
    /** The blinking rate. */
    final double blinkingRate;
    /** The results. */
    final MemoryPeakResults results;
    /** The minx. */
    final double minx;
    /** The miny. */
    final double miny;
    /** The maxx. */
    final double maxx;
    /** The maxy. */
    final double maxy;
    /** The nm per pixel. */
    final double nmPerPixel;
    /** The molecules. */
    final List<Molecule> molecules;
    /** The sigma S. */
    final double sigmaS;
    /** The peak density. */
    final double densityPeaks;
    /** The protein density. */
    final double densityProtein;
    /** The seconds. */
    final double seconds;
    /** The area. */
    final double area;

    /**
     * Instantiates a new molecules results.
     *
     * @param source the source
     */
    MoleculesResults(Settings source) {
      blinkingRate = source.blinkingRate;
      results = source.results;
      minx = source.minx;
      miny = source.miny;
      maxx = source.maxx;
      maxy = source.maxy;
      nmPerPixel = source.nmPerPixel;
      molecules = source.molecules;
      sigmaS = source.sigmaS;
      densityPeaks = source.densityPeaks;
      densityProtein = source.densityProtein;
      seconds = source.seconds;
      area = source.area;
    }
  }

  /**
   * Gets the last used analysis settings.
   *
   * @return the analysis settings
   */
  static MoleculesResults getMoleculesResults() {
    return new MoleculesResults(Settings.load());
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    // Require some fit results and selected regions
    final boolean resultsAvailable = MemoryPeakResults.countMemorySize() > 0;

    if (!getRunMode(resultsAvailable)) {
      return;
    }

    if (settings.runMode != 3) {
      settings.results = ResultsManager.loadInputResults(settings.inputOption, true, null, null);
      if (MemoryPeakResults.isEmpty(settings.results)) {
        IJ.error(TITLE, "No results could be loaded");
        return;
      }

      if (settings.results.getCalibration() == null) {
        IJ.error(TITLE, "Results are not calibrated");
        return;
      }

      // Get the lifetime before cropping as this is the true representation of the number of
      // frames.
      // This may truncate the lifetime if the first/last localisation are not near the end of the
      // acquisition lifetime
      getLifetime();

      settings.results = cropToRoi(settings.results);
      if (settings.results.size() == 0) {
        IJ.error(TITLE, "No results within the crop region");
        return;
      }
    }

    // Clear cached results
    settings.molecules = null;

    // Different run-modes for generating the set of molecules for analysis
    switch (settings.runMode) {
      case 0:
        runPcPalm();
        break;
      case 1:
        runManualTracing();
        break;
      case 2:
        runInMemoryResults();
        break;
      case 3:
      default:
        runSimulation(resultsAvailable);
        settings.area = settings.simulationSize * settings.simulationSize;
        settings.seconds = 100; // Use an arbitrary lifetime
        break;
    }

    if (settings.molecules == null) {
      return;
    }
    if (settings.molecules.size() < 2) {
      IJ.error(TITLE, "Not enough molecules to construct a binary image");
      return;
    }

    // Generate binary PALM image
    if (!createImage(settings.molecules)) {
      return;
    }

    // Density is required for the PC analysis
    settings.densityPeaks = calculatePeakDensity();
    if (settings.runMode == 0 || settings.runMode == 3) {
      // Blinking rate is mentioned in the PC-PALM protocol and so we include it here.
      // TODO - Add automated estimation of the blinking rate from the data using the method of
      // Annibale, et al (2011), Quantitative photo activated localization microscopy: unraveling
      // the effects of photoblinking. PLoS One, 6(7): e22678
      // (http://dx.doi.org/10.1371%2Fjournal.pone.0022678)
      settings.densityProtein = settings.densityPeaks / settings.blinkingRate;
      log("Peak Density = %s (um^-2). Protein Density = %s (um^-2)",
          MathUtils.rounded(settings.densityPeaks * 1e6),
          MathUtils.rounded(settings.densityProtein * 1e6));
    } else {
      // No blinking rate for non PC-PALM methods. This can be configured in later plugins if
      // required.
      settings.blinkingRate = 1;
      settings.densityProtein = settings.densityPeaks;
      log("Molecule Density = %s (um^-2)", MathUtils.rounded(settings.densityPeaks * 1e6));
    }

    log("Results lifetime = %s s", MathUtils.rounded(settings.seconds));

    // Use a second plugin filter that will work on a region drawn on the binary image
    // and compute the PALM analysis

    final double seconds = (System.currentTimeMillis() - start) / 1000.0;
    final String msg = TITLE + " complete : " + seconds + "s";
    IJ.showStatus(msg);
    log(msg);
  }

  private boolean getRunMode(boolean resultsAvailable) {
    ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

    // Build a list of all images with a region ROI
    final List<String> titles = new LinkedList<>();
    for (final int imageId : ImageJUtils.getIdList()) {
      final ImagePlus imp = WindowManager.getImage(imageId);
      if (imp != null && imp.getRoi() != null && imp.getRoi().isArea()) {
        titles.add(imp.getTitle());
      }
    }

    settings = Settings.load();

    if (!resultsAvailable) {
      settings.runMode = 3;

      gd.addMessage("Simulate molecules for cluster analysis.\n"
          + "Computes a binary image from localisation data");

      gd.addNumericField("Molecules", settings.numberOfMolecules, 0);
      gd.addNumericField("Simulation_size (um)", settings.simulationSize, 2);
      gd.addNumericField("Blinking_rate", settings.blinkingRate, 2);
      gd.addChoice("Blinking_distribution", Settings.BLINKING_DISTRIBUTION,
          settings.blinkingDistribution);
      gd.addNumericField("Average_precision (nm)", settings.sigmaS, 2);
      gd.addCheckbox("Show_histograms", settings.showHistograms);
      gd.addCheckbox("Distance_analysis", settings.distanceAnalysis);
      gd.addChoice("Cluster_simulation", Settings.CLUSTER_SIMULATION, settings.clusterSimulation);
      gd.addNumericField("Cluster_number", settings.clusterNumber, 2);
      gd.addNumericField("Cluster_variation (SD)", settings.clusterNumberStdDev, 2);
      gd.addNumericField("Cluster_radius", settings.clusterRadius, 2);
      gd.addCheckbox("Show_cluster_mask", settings.showClusterMask);

      Recorder.recordOption("Run_mode", Settings.RUN_MODE[settings.runMode]);
    } else {
      gd.addMessage("Prepare molecules for cluster analysis.\n"
          + "Computes a binary image from raw localisation data");
      ResultsManager.addInput(gd, settings.inputOption, InputSource.MEMORY);
      if (!titles.isEmpty()) {
        gd.addCheckbox((titles.size() == 1) ? "Use_ROI" : "Choose_ROI", settings.chooseRoi);
      }
      gd.addChoice("Run_mode", Settings.RUN_MODE, settings.runMode);
    }

    gd.addMessage("Select options for low resolution image:");
    gd.addSlider("Image_size (px)", 512, 2048, settings.lowResolutionImageSize);
    gd.addSlider("ROI_size (um)", 1.5, 4, settings.roiSizeInUm);
    gd.addMessage("Select options for high resolution image:");
    gd.addCheckbox("Show_high_res_image", settings.showHighResolutionImage);
    gd.addSlider("nm_per_pixel_limit", 0, 20, settings.nmPerPixelLimit);
    gd.addMessage("Optionally remove all analysis results from memory");
    gd.addCheckbox("Clear_results", settings.clearResults);

    gd.addHelp(HelpUrls.getUrl("pc-palm-molecules"));
    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    if (!resultsAvailable) {
      settings.numberOfMolecules = (int) Math.abs(gd.getNextNumber());
      settings.simulationSize = Math.abs(gd.getNextNumber());
      settings.blinkingRate = Math.abs(gd.getNextNumber());
      settings.blinkingDistribution = gd.getNextChoiceIndex();
      settings.sigmaS = Math.abs(gd.getNextNumber());
      settings.showHistograms = gd.getNextBoolean();
      settings.distanceAnalysis = gd.getNextBoolean();
      settings.clusterSimulation = gd.getNextChoiceIndex();
      settings.clusterNumber = Math.abs(gd.getNextNumber());
      settings.clusterNumberStdDev = Math.abs(gd.getNextNumber());
      settings.clusterRadius = Math.abs(gd.getNextNumber());
      settings.showClusterMask = gd.getNextBoolean();
    } else {
      settings.inputOption = ResultsManager.getInputSource(gd);
      if (!titles.isEmpty()) {
        settings.chooseRoi = gd.getNextBoolean();
      }
      settings.runMode = gd.getNextChoiceIndex();
    }
    settings.lowResolutionImageSize = (int) gd.getNextNumber();
    settings.roiSizeInUm = gd.getNextNumber();
    settings.showHighResolutionImage = gd.getNextBoolean();
    settings.nmPerPixelLimit = Math.abs(gd.getNextNumber());
    settings.clearResults = gd.getNextBoolean();

    settings.save();

    // Check arguments
    try {
      if (!resultsAvailable) {
        ParameterUtils.isAboveZero("Molecules", settings.numberOfMolecules);
        ParameterUtils.isAboveZero("Simulation size", settings.simulationSize);
        ParameterUtils.isEqualOrAbove("Blinking rate", settings.blinkingRate, 1);
        ParameterUtils.isEqualOrAbove("Cluster number", settings.clusterNumber, 1);
      }
      ParameterUtils.isAbove("Image scale", settings.lowResolutionImageSize, 1);
      ParameterUtils.isAboveZero("ROI size", settings.roiSizeInUm);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    if (!titles.isEmpty() && settings.chooseRoi && resultsAvailable) {
      if (titles.size() == 1) {
        settings.roiImage = titles.get(0);
        Recorder.recordOption("Image", settings.roiImage);
      } else {
        final String[] items = titles.toArray(new String[titles.size()]);
        gd = new ExtendedGenericDialog(TITLE);
        gd.addMessage("Select the source image for the ROI");
        gd.addChoice("Image", items, settings.roiImage);
        gd.showDialog();
        if (gd.wasCanceled()) {
          return false;
        }
        settings.roiImage = gd.getNextChoice();
      }
      final ImagePlus imp = WindowManager.getImage(settings.roiImage);

      roiBounds = imp.getRoi().getBounds();
      roiImageWidth = imp.getWidth();
      roiImageHeight = imp.getHeight();
    } else {
      roiBounds = null;
    }

    if (!resultsAvailable && !getPValue()) {
      return false;
    }

    if (settings.clearResults) {
      PcPalmAnalysis.clearResults();
      PcPalmFitting.clearResults();
    }

    return true;
  }

  private MemoryPeakResults cropToRoi(MemoryPeakResults results) {
    final Rectangle bounds = results.getBounds(true);
    settings.area =
        (bounds.width * bounds.height * results.getNmPerPixel() * results.getNmPerPixel()) / 1e6;
    if (roiBounds == null) {
      return results;
    }

    // Adjust bounds relative to input results image
    final double xscale = (double) roiImageWidth / bounds.width;
    final double yscale = (double) roiImageHeight / bounds.height;

    final float minX = (float) Math.floor(roiBounds.x / xscale);
    final float maxX = (float) Math.ceil((roiBounds.x + roiBounds.width) / xscale);
    final float minY = (float) Math.floor(roiBounds.y / yscale);
    final float maxY = (float) Math.ceil((roiBounds.y + roiBounds.height) / yscale);

    // Update the area with the cropped region
    settings.area *= (maxX - minX) / bounds.width;
    settings.area *= (maxY - minY) / bounds.height;

    // Create a new set of results within the bounds
    final MemoryPeakResults newResults = new MemoryPeakResults();
    newResults.begin();
    results.forEach(DistanceUnit.PIXEL, (XyrResultProcedure) (x, y, result) -> {
      if (x >= minX && x <= maxX && y >= minY && y <= maxY) {
        newResults.add(result);
      }
    });
    newResults.end();
    newResults.copySettings(results);
    newResults
        .setBounds(new Rectangle((int) minX, (int) minY, (int) (maxX - minX), (int) (maxY - minY)));
    return newResults;
  }

  private void runPcPalm() {
    if (!showPcPalmDialog()) {
      return;
    }

    startLog();

    // Follow the PC-PALM protocol
    log("Fitting localisation precision...");
    final List<Molecule> localisations = extractLocalisations(settings.results);
    final double sigmaRaw = calculateAveragePrecision(localisations, "Localisations");
    log("%d localisations with an average precision of %.2f", settings.results.size(), sigmaRaw);

    log("Fitting molecule precision...");
    final ArrayList<Molecule> singles = new ArrayList<>();
    final double[] precisions = localisations.stream().mapToDouble(m -> m.precision).toArray();
    settings.molecules = extractMolecules(settings.results, precisions, sigmaRaw, singles);
    if (settings.singlesModeIndex == 1) {
      settings.molecules.addAll(singles);
    }
    settings.sigmaS = calculateAveragePrecision(settings.molecules, "Molecules");
    log("%d molecules with an average precision of %.2f", settings.molecules.size(),
        settings.sigmaS);

    // Q. Should this filter the original localisations or just the grouped peaks?
    if (settings.singlesModeIndex == 2) {
      settings.molecules.addAll(singles);
    }
    settings.molecules = filterMolecules(settings.molecules, settings.sigmaS);
    log("%d molecules within precision %.2f", settings.molecules.size(), 3 * settings.sigmaS);
  }

  private void startLog() {
    logSpacer();
    log(TITLE);
    logSpacer();
    start = System.currentTimeMillis();
  }

  private boolean showPcPalmDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

    gd.addMessage("Estimate the average localisation precision by fitting histograms.\n"
        + "Use the precision to trace localisations into molecule pulses.");

    gd.addNumericField("Histogram_bins", settings.histogramBins, 0);
    gd.addChoice("Singles_mode", Settings.singlesMode, settings.singlesModeIndex);
    gd.addCheckbox("Simplex_fit", settings.simplexFitting);
    gd.addCheckbox("Show_histograms", settings.showHistograms);
    gd.addCheckbox("Binary_image", settings.binaryImage);
    gd.addNumericField("Blinking_rate", settings.blinkingRate, 2);

    gd.addHelp(HelpUrls.getUrl("pc-palm-molecules"));
    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.histogramBins = (int) gd.getNextNumber();
    settings.singlesModeIndex = gd.getNextChoiceIndex();
    settings.simplexFitting = gd.getNextBoolean();
    settings.showHistograms = gd.getNextBoolean();
    settings.binaryImage = gd.getNextBoolean();
    settings.blinkingRate = gd.getNextNumber();

    // Check arguments
    try {
      ParameterUtils.isEqualOrAbove("Histogram bins", settings.histogramBins, 0);
      ParameterUtils.isEqualOrAbove("Blinking rate", settings.blinkingRate, 1);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return true;
  }

  /**
   * Extract molecules for the PC-PALM analysis.
   *
   * <p>Estimate the localisation uncertainty (precision) of each molecule using the formula of
   * Mortensen, et al (2010), Nature Methods 7, 377-381. Store distance in nm and signal in photons
   * using the calibration
   *
   * @param results the results
   * @return the array list
   * @throws DataException If conversion to nm and photons with computed precision is not possible
   */
  public List<Molecule> extractLocalisations(MemoryPeakResults results) {
    final ArrayList<Molecule> list = new ArrayList<>(results.size());

    // Access calibrated data
    final StandardResultProcedure sp =
        new StandardResultProcedure(results, DistanceUnit.NM, IntensityUnit.PHOTON);
    sp.getIxy();
    final PrecisionResultProcedure pp = new PrecisionResultProcedure(results);
    pp.getPrecision();

    for (int i = 0, size = pp.size(); i < size; i++) {
      list.add(new Molecule(sp.x[i], sp.y[i], pp.precisions[i], sp.intensity[i]));
    }
    return list;
  }

  /**
   * Calculate the average precision by fitting a skewed Gaussian to the histogram of the precision
   * distribution.
   *
   * @param molecules the molecules
   * @param subTitle the sub title
   * @return The average precision
   */
  private double calculateAveragePrecision(List<Molecule> molecules, String subTitle) {
    final String title = (settings.showHistograms) ? TITLE + " Histogram " + subTitle : null;
    return calculateAveragePrecision(molecules, title, settings.histogramBins, true, true);
  }

  /**
   * Calculate the average precision by fitting a skewed Gaussian to the histogram of the precision
   * distribution.
   *
   * <p>A simple mean and SD of the histogram is computed. If the mean of the Skewed Gaussian does
   * not fit within 3 SDs of the simple mean then the simple mean is returned.
   *
   * @param molecules the molecules
   * @param title the plot title (null if no plot should be displayed)
   * @param histogramBins the histogram bins
   * @param logFitParameters Record the fit parameters to the ImageJ log
   * @param removeOutliers The distribution is created using all values within 1.5x the
   *        inter-quartile range (IQR) of the data
   * @return The average precision
   */
  public double calculateAveragePrecision(List<Molecule> molecules, String title, int histogramBins,
      boolean logFitParameters, boolean removeOutliers) {
    // Plot histogram of the precision
    final float[] data = new float[molecules.size()];
    final DescriptiveStatistics stats = new DescriptiveStatistics();
    double yMin = Double.NEGATIVE_INFINITY;
    double yMax = 0;
    for (int i = 0; i < data.length; i++) {
      data[i] = (float) molecules.get(i).precision;
      stats.addValue(data[i]);
    }

    // Set the min and max y-values using 1.5 x IQR
    if (removeOutliers) {
      final double lower = stats.getPercentile(25);
      final double upper = stats.getPercentile(75);
      if (Double.isNaN(lower) || Double.isNaN(upper)) {
        if (logFitParameters) {
          ImageJUtils.log("Error computing IQR: %f - %f", lower, upper);
        }
      } else {
        final double iqr = upper - lower;

        yMin = FastMath.max(lower - iqr, stats.getMin());
        yMax = FastMath.min(upper + iqr, stats.getMax());

        if (logFitParameters) {
          ImageJUtils.log("  Data range: %f - %f. Plotting 1.5x IQR: %f - %f", stats.getMin(),
              stats.getMax(), yMin, yMax);
        }
      }
    }

    if (yMin == Double.NEGATIVE_INFINITY) {
      yMin = stats.getMin();
      yMax = stats.getMax();

      if (logFitParameters) {
        ImageJUtils.log("  Data range: %f - %f", yMin, yMax);
      }
    }

    int bins;
    if (histogramBins <= 0) {
      bins = (int) Math.ceil((stats.getMax() - stats.getMin())
          / HistogramPlot.getBinWidthScottsRule(stats.getStandardDeviation(), (int) stats.getN()));
    } else {
      bins = histogramBins;
    }

    final float[][] hist = HistogramPlot.calcHistogram(data, yMin, yMax, bins);

    Plot2 plot = null;
    if (title != null) {
      plot = new Plot2(title, "Precision", "Frequency");
      final float[] xValues = hist[0];
      final float[] yValues = hist[1];
      if (xValues.length > 0) {
        final double xPadding = 0.05 * (xValues[xValues.length - 1] - xValues[0]);
        plot.setLimits(xValues[0] - xPadding, xValues[xValues.length - 1] + xPadding, 0,
            MathUtils.max(yValues) * 1.05);
      }
      plot.addPoints(xValues, yValues, Plot2.BAR);
      ImageJUtils.display(title, plot);
    }

    // Extract non-zero data
    float[] x = Arrays.copyOf(hist[0], hist[0].length);
    float[] y = hist[1];
    int count = 0;
    final float dx = (x[1] - x[0]) * 0.5f;
    for (int i = 0; i < y.length; i++) {
      if (y[i] > 0) {
        x[count] = x[i] + dx;
        y[count] = y[i];
        count++;
      }
    }
    x = Arrays.copyOf(x, count);
    y = Arrays.copyOf(y, count);

    // Sense check to fitted data. Get mean and SD of histogram
    final double[] stats2 = HistogramPlot.getHistogramStatistics(x, y);
    double mean = stats2[0];
    if (logFitParameters) {
      log("  Initial Statistics: %f +/- %f", stats2[0], stats2[1]);
    }

    // Standard Gaussian fit
    final double[] parameters = fitGaussian(x, y);
    if (parameters == null) {
      log("  Failed to fit initial Gaussian");
      return mean;
    }
    double newMean = parameters[1];
    double error = Math.abs(stats2[0] - newMean) / stats2[1];
    if (error > 3) {
      log("  Failed to fit Gaussian: %f standard deviations from histogram mean", error);
      return mean;
    }
    if (newMean < yMin || newMean > yMax) {
      log("  Failed to fit Gaussian: %f outside data range %f - %f", newMean, yMin, yMax);
      return mean;
    }

    mean = newMean;

    if (logFitParameters) {
      log("  Initial Gaussian: %f @ %f +/- %f", parameters[0], parameters[1], parameters[2]);
    }

    final double[] initialSolution = new double[] {parameters[0], parameters[1], parameters[2], -1};

    // Fit to a skewed Gaussian (or appropriate function)
    final double[] skewParameters = fitSkewGaussian(x, y, initialSolution);
    if (skewParameters == null) {
      log("  Failed to fit Skewed Gaussian");
      return mean;
    }

    final SkewNormalFunction sn = new SkewNormalFunction(skewParameters);
    if (logFitParameters) {
      log("  Skewed Gaussian: %f @ %f +/- %f (a = %f) => %f +/- %f", skewParameters[0],
          skewParameters[1], skewParameters[2], skewParameters[3], sn.getMean(),
          Math.sqrt(sn.getVariance()));
    }

    newMean = sn.getMean();
    error = Math.abs(stats2[0] - newMean) / stats2[1];
    if (error > 3) {
      log("  Failed to fit Skewed Gaussian: %f standard deviations from histogram mean", error);
      return mean;
    }
    if (newMean < yMin || newMean > yMax) {
      log("  Failed to fit Skewed Gaussian: %f outside data range %f - %f", newMean, yMin, yMax);
      return mean;
    }

    // Use original histogram x-axis to maintain all the bins
    if (plot != null) {
      x = hist[0];
      for (int i = 0; i < y.length; i++) {
        x[i] += dx;
      }
      plot.setColor(Color.red);
      addToPlot(plot, x, skewParameters, Plot.LINE);

      plot.setColor(Color.black);
      ImageJUtils.display(title, plot);
    }

    // Return the average precision from the fitted curve
    return newMean;
  }

  @Nullable
  private static double[] fitGaussian(float[] x, float[] y) {
    final WeightedObservedPoints obs = new WeightedObservedPoints();
    for (int i = 0; i < x.length; i++) {
      obs.add(x[i], y[i]);
    }

    final Collection<WeightedObservedPoint> observations = obs.toList();
    final GaussianCurveFitter fitter = GaussianCurveFitter.create().withMaxIterations(2000);
    final GaussianCurveFitter.ParameterGuesser guess =
        new GaussianCurveFitter.ParameterGuesser(observations);
    double[] initialGuess = null;
    try {
      initialGuess = guess.guess();
      return fitter.withStartPoint(initialGuess).fit(observations);
    } catch (final TooManyEvaluationsException ex) {
      // Use the initial estimate
      return initialGuess;
    } catch (final Exception ex) {
      // Just in case there is another exception type, or the initial estimate failed
      return null;
    }
  }

  @Nullable
  private double[] fitSkewGaussian(float[] x, float[] y, double[] initialSolution) {
    try {
      return (settings.simplexFitting) ? optimiseSimplex(x, y, initialSolution)
          : optimiseLeastSquares(x, y, initialSolution);
    } catch (final TooManyEvaluationsException ex) {
      return null;
    }
  }

  private double[] optimiseLeastSquares(float[] x, float[] y, double[] initialSolution) {
    // Least-squares optimisation using numerical gradients
    final SkewNormalDifferentiableFunction function =
        new SkewNormalDifferentiableFunction(initialSolution);
    function.addData(x, y);

    final LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();

    //@formatter:off
    final LeastSquaresProblem problem = new LeastSquaresBuilder()
        .maxEvaluations(Integer.MAX_VALUE)
        .maxIterations(3000)
        .start(initialSolution)
        .target(function.calculateTarget())
        .weight(new DiagonalMatrix(function.calculateWeights()))
        .model(function, function::jacobian)
        .build();
    //@formatter:on

    final Optimum optimum = optimizer.optimize(problem);

    return optimum.getPoint().toArray();
  }

  private double[] optimiseSimplex(float[] x, float[] y, double[] initialSolution) {
    // Simplex optimisation
    final SkewNormalMultivariateFunction sn2 = new SkewNormalMultivariateFunction(initialSolution);
    sn2.addData(x, y);
    final NelderMeadSimplex simplex = new NelderMeadSimplex(4);
    final SimplexOptimizer opt = new SimplexOptimizer(1e-6, 1e-10);
    final PointValuePair solution = opt.optimize(new MaxEval(1000),
        new InitialGuess(initialSolution), simplex, new ObjectiveFunction(sn2), GoalType.MINIMIZE);

    return solution.getPointRef();
  }

  /**
   * Add the skewed gaussian to the histogram plot.
   *
   * @param plot the plot
   * @param x the x
   * @param parameters Gaussian parameters
   * @param shape the shape
   */
  private static void addToPlot(Plot2 plot, float[] x, double[] parameters, int shape) {
    final SkewNormalFunction sn = new SkewNormalFunction(parameters);
    final float[] y = new float[x.length];
    for (int i = 0; i < x.length; i++) {
      y[i] = (float) sn.evaluate(x[i]);
    }
    plot.addPoints(x, y, shape);
  }

  /**
   * Group all localisations in successive frames within 2.5x of the initial precision estimate into
   * a single molecule
   *
   * @param results The results
   * @param precisions the precisions
   * @param sigmaRaw The initial precision estimate
   * @param singles a list of the singles (not grouped into molecules)
   * @return a list of molecules
   */
  private static ArrayList<Molecule> extractMolecules(MemoryPeakResults results,
      double[] precisions, double sigmaRaw, ArrayList<Molecule> singles) {
    return traceMolecules(results, precisions, sigmaRaw * 2.5, 1, singles);
  }

  /**
   * Trace localisations.
   *
   * @param results The results
   * @param precisions the precisions
   * @param distance The distance threshold (nm)
   * @param time The time threshold (frames)
   * @param singles a list of the singles (not grouped into molecules)
   * @return a list of molecules
   */
  private static ArrayList<Molecule> traceMolecules(MemoryPeakResults results, double[] precisions,
      double distance, int time, ArrayList<Molecule> singles) {
    // These plugins are not really supported so just leave them to throw an exception if
    // the data cannot be handled
    final TypeConverter<IntensityUnit> ic =
        results.getCalibrationReader().getIntensityConverter(IntensityUnit.PHOTON);
    final TypeConverter<DistanceUnit> dc =
        results.getCalibrationReader().getDistanceConverter(DistanceUnit.NM);

    // Create a new dataset with the precision
    final MemoryPeakResults results2 = new MemoryPeakResults(results.size());
    for (int i = 0, size = results.size(); i < size; i++) {
      final AttributePeakResult peak2 = new AttributePeakResult(results.get(i));
      peak2.setPrecision(precisions[i]);
      results2.add(peak2);
    }

    final TraceManager tm = new TraceManager(results2);
    final double distanceThreshold = dc.convertBack(distance);

    tm.traceMolecules(distanceThreshold, time);
    final Trace[] traces = tm.getTraces();
    final ArrayList<Molecule> molecules = new ArrayList<>(traces.length);

    for (final Trace t : traces) {
      final double p = t.getLocalisationPrecision(dc);
      final float[] centroid = t.getCentroid();
      final List<Molecule> list = t.size() == 1 ? singles : molecules;
      list.add(new Molecule(dc.convert(centroid[0]), dc.convert(centroid[1]), p,
          ic.convert(t.getSignal())));
    }
    log("  %d localisations traced to %d molecules (%d singles, %d traces) using d=%.2f nm,"
        + " t=%d frames (%s s)", results.size(), molecules.size() + singles.size(), singles.size(),
        molecules.size(), distance, time,
        MathUtils.rounded(time * results.getCalibrationReader().getExposureTime() / 1000.0));
    return molecules;
  }

  /**
   * Calculate the density of peaks in the original data.
   *
   * @return The peak density
   */
  private double calculatePeakDensity() {
    // Use the area from the source of the molecules
    return settings.molecules.size() / (settings.area * 1E6);
  }

  /**
   * Return a new list, removing all molecules with a precision over 3x of the precision estimate.
   *
   * @param molecules the molecules
   * @param sigmaS The precision estimate
   * @return the array list
   */
  private static ArrayList<Molecule> filterMolecules(List<Molecule> molecules, double sigmaS) {
    final ArrayList<Molecule> newMolecules = new ArrayList<>(molecules.size());
    final double limit = 3 * sigmaS;
    for (final Molecule m : molecules) {
      if (m.precision <= limit) {
        newMolecules.add(m);
      }
    }
    return newMolecules;
  }

  private void runManualTracing() {
    if (!showManualTracingDialog()) {
      return;
    }

    startLog();

    // Convert seconds to frames
    final int timeInFrames = FastMath.max(1, (int) Math.round(settings.timeThreshold * 1000.0
        / settings.results.getCalibrationReader().getExposureTime()));

    // Get precisions
    final PrecisionResultProcedure pp = new PrecisionResultProcedure(settings.results);
    pp.getPrecision();

    final ArrayList<Molecule> singles = new ArrayList<>();
    settings.molecules = traceMolecules(settings.results, pp.precisions, settings.distanceThreshold,
        timeInFrames, singles);
    settings.molecules.addAll(singles);
  }

  private boolean showManualTracingDialog() {
    final GenericDialog gd = new GenericDialog(TITLE);

    gd.addMessage("Use distance and time thresholds to trace localisations into molecules.");

    gd.addNumericField("Distance (nm)", settings.distanceThreshold, 0);
    gd.addNumericField("Time (seconds)", settings.timeThreshold, 2);

    gd.addHelp(HelpUrls.getUrl("pc-palm-molecules"));
    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.distanceThreshold = Math.abs(gd.getNextNumber());
    settings.timeThreshold = Math.abs(gd.getNextNumber());

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Distance threshold", settings.distanceThreshold);
      ParameterUtils.isAboveZero("Time threshold", settings.timeThreshold);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return true;
  }

  private void runInMemoryResults() {
    startLog();
    settings.molecules = extractLocalisations(settings.results);
  }

  private void runSimulation(boolean resultsAvailable) {
    if (resultsAvailable && !showSimulationDialog()) {
      return;
    }

    startLog();

    log("Simulation parameters");
    if (settings.blinkingDistribution == 3) {
      log("  - Clusters = %d", settings.numberOfMolecules);
      log("  - Simulation size = %s um", MathUtils.rounded(settings.simulationSize, 4));
      log("  - Molecules/cluster = %s", MathUtils.rounded(settings.blinkingRate, 4));
      log("  - Blinking distribution = %s",
          Settings.BLINKING_DISTRIBUTION[settings.blinkingDistribution]);
      log("  - p-Value = %s", MathUtils.rounded(settings.pvalue, 4));
    } else {
      log("  - Molecules = %d", settings.numberOfMolecules);
      log("  - Simulation size = %s um", MathUtils.rounded(settings.simulationSize, 4));
      log("  - Blinking rate = %s", MathUtils.rounded(settings.blinkingRate, 4));
      log("  - Blinking distribution = %s",
          Settings.BLINKING_DISTRIBUTION[settings.blinkingDistribution]);
    }
    log("  - Average precision = %s nm", MathUtils.rounded(settings.sigmaS, 4));
    log("  - Clusters simulation = " + Settings.CLUSTER_SIMULATION[settings.clusterSimulation]);
    if (settings.clusterSimulation > 0) {
      log("  - Cluster number = %s +/- %s", MathUtils.rounded(settings.clusterNumber, 4),
          MathUtils.rounded(settings.clusterNumberStdDev, 4));
      log("  - Cluster radius = %s nm", MathUtils.rounded(settings.clusterRadius, 4));
    }

    final double nmPerPixel = 100;
    double width = settings.simulationSize * 1000.0;
    // Allow a border of 3 x sigma for +/- precision
    width -= 3 * settings.sigmaS;
    final UniformRandomProvider rng = UniformRandomProviders.create();
    final UniformDistribution dist =
        new UniformDistribution(null, new double[] {width, width, 0}, rng.nextInt());
    final NormalizedGaussianSampler gauss = SamplerUtils.createNormalizedGaussianSampler(rng);

    settings.molecules = new ArrayList<>(settings.numberOfMolecules);
    // Create some dummy results since the calibration is required for later analysis
    settings.results = new MemoryPeakResults(PsfHelper.create(PSFType.CUSTOM));
    settings.results.setCalibration(CalibrationHelper.create(nmPerPixel, 1, 100));
    settings.results.setSource(new NullSource("Molecule Simulation"));
    settings.results.begin();
    int count = 0;

    // Generate a sequence of coordinates
    final ArrayList<double[]> xyz = new ArrayList<>((int) (settings.numberOfMolecules * 1.1));

    final Statistics statsRadius = new Statistics();
    final Statistics statsSize = new Statistics();
    final String maskTitle = TITLE + " Cluster Mask";
    ByteProcessor bp = null;
    double maskScale = 0;

    // TODO - Add a fluctuations model to this.

    if (settings.clusterSimulation > 0) {
      // Simulate clusters.

      // Note: In the Veatch et al. paper (Plos 1, e31457) correlation functions are built using
      // circles
      // with small radii of 4-8 Arbitrary Units (AU) or large radii of 10-30 AU. A fluctuations
      // model is
      // created at T = 1.075 Tc. It is not clear exactly how the particles are distributed.
      // It may be that a mask is created first using the model. The particles are placed on the
      // mask using
      // a specified density. This simulation produces a figure to show either a damped cosine
      // function
      // (circles) or an exponential (fluctuations). The number of particles in each circle may be
      // randomly
      // determined just by density. The figure does not discuss the derivation of the cluster size
      // statistic.
      //
      // If this plugin simulation is run with a uniform distribution and blinking rate of 1 then
      // the damped
      // cosine function is reproduced. The curve crosses g(r)=1 at a value equivalent to the
      // average
      // distance to the centre-of-mass of each drawn cluster, not the input cluster radius
      // parameter (which
      // is a hard upper limit on the distance to centre).

      final int maskSize = settings.lowResolutionImageSize;
      int[] mask = null;
      maskScale = width / maskSize; // scale is in nm/pixel

      final ArrayList<double[]> clusterCentres = new ArrayList<>();
      int totalSteps = 1 + (int) Math.ceil(settings.numberOfMolecules / settings.clusterNumber);
      if (settings.clusterSimulation == 2 || settings.clusterSimulation == 3) {
        // Clusters are non-overlapping circles

        // Ensure the circles do not overlap by using an exclusion mask that accumulates
        // out-of-bounds pixels by drawing the last cluster (plus some border) on an image. When no
        // more pixels are available then stop generating molecules.
        // This is done by cumulatively filling a mask and using the MaskDistribution to select
        // a new point. This may be slow but it works.

        // TODO - Allow clusters of different sizes...

        mask = new int[maskSize * maskSize];
        Arrays.fill(mask, 255);
        MaskDistribution maskDistribution =
            new MaskDistribution(mask, maskSize, maskSize, 0, maskScale, maskScale, rng);
        double[] centre;
        IJ.showStatus("Computing clusters mask");
        final int roiRadius = (int) Math.round((settings.clusterRadius * 2) / maskScale);

        if (settings.clusterSimulation == 3) {
          // Generate a mask of circles then sample from that.
          // If we want to fill the mask completely then adjust the total steps to be the number of
          // circles that can fit inside the mask.
          totalSteps = (int) (maskSize * maskSize
              / (Math.PI * MathUtils.pow2(settings.clusterRadius / maskScale)));
        }

        while ((centre = maskDistribution.next()) != null && clusterCentres.size() < totalSteps) {
          IJ.showProgress(clusterCentres.size(), totalSteps);
          // The mask returns the coordinates with the centre of the image at 0,0
          centre[0] += width / 2;
          centre[1] += width / 2;
          clusterCentres.add(centre);

          // Fill in the mask around the centre to exclude any more circles that could overlap
          final double cx = centre[0] / maskScale;
          final double cy = centre[1] / maskScale;
          fillMask(mask, maskSize, (int) cx, (int) cy, roiRadius, 0);
          try {
            maskDistribution =
                new MaskDistribution(mask, maskSize, maskSize, 0, maskScale, maskScale, rng);
          } catch (final IllegalArgumentException ex) {
            // This can happen when there are no more non-zero pixels
            log("WARNING: No more room for clusters on the mask area (created %d of estimated %d)",
                clusterCentres.size(), totalSteps);
            break;
          }
        }
        ImageJUtils.finished();
      } else {
        // Pick centres randomly from the distribution
        while (clusterCentres.size() < totalSteps) {
          clusterCentres.add(dist.next());
        }
      }

      if (settings.showClusterMask || settings.clusterSimulation == 3) {
        // Show the mask for the clusters
        if (mask == null) {
          mask = new int[maskSize * maskSize];
        } else {
          Arrays.fill(mask, 0);
        }
        final int roiRadius = (int) Math.round((settings.clusterRadius) / maskScale);
        for (final double[] c : clusterCentres) {
          final double cx = c[0] / maskScale;
          final double cy = c[1] / maskScale;
          fillMask(mask, maskSize, (int) cx, (int) cy, roiRadius, 1);
        }

        if (settings.clusterSimulation == 3) {
          // We have the mask. Now pick points at random from the mask.
          final MaskDistribution maskDistribution =
              new MaskDistribution(mask, maskSize, maskSize, 0, maskScale, maskScale, rng);

          // Allocate each molecule position to a parent circle so defining clusters.
          final int[][] clusters = new int[clusterCentres.size()][];
          final int[] clusterSize = new int[clusters.length];

          for (int i = 0; i < settings.numberOfMolecules; i++) {
            final double[] centre = maskDistribution.next();
            // The mask returns the coordinates with the centre of the image at 0,0
            centre[0] += width / 2;
            centre[1] += width / 2;
            xyz.add(centre);

            // Output statistics on cluster size and number.
            // TODO - Finding the closest cluster could be done better than an all-vs-all comparison
            double max = distance2(centre, clusterCentres.get(0));
            int cluster = 0;
            for (int j = 1; j < clusterCentres.size(); j++) {
              final double d2 = distance2(centre, clusterCentres.get(j));
              if (d2 < max) {
                max = d2;
                cluster = j;
              }
            }

            // Assign point i to cluster
            centre[2] = cluster;

            if (clusterSize[cluster] == 0) {
              clusters[cluster] = new int[10];
            }
            if (clusters[cluster].length <= clusterSize[cluster]) {
              clusters[cluster] =
                  Arrays.copyOf(clusters[cluster], (int) (clusters[cluster].length * 1.5));
            }
            clusters[cluster][clusterSize[cluster]++] = i;
          }

          // Generate real cluster size statistics
          for (int j = 0; j < clusterSize.length; j++) {
            final int size = clusterSize[j];
            if (size == 0) {
              continue;
            }

            statsSize.add(size);

            if (size == 1) {
              statsRadius.add(0);
              continue;
            }

            // Find centre of cluster and add the distance to each point
            final double[] com = new double[2];
            for (int n = 0; n < size; n++) {
              final double[] xy = xyz.get(clusters[j][n]);
              for (int k = 0; k < 2; k++) {
                com[k] += xy[k];
              }
            }
            for (int k = 0; k < 2; k++) {
              com[k] /= size;
            }
            for (int n = 0; n < size; n++) {
              final double dx = xyz.get(clusters[j][n])[0] - com[0];
              final double dy = xyz.get(clusters[j][n])[1] - com[1];
              statsRadius.add(Math.sqrt(dx * dx + dy * dy));
            }
          }
        }

        if (settings.showClusterMask) {
          bp = new ByteProcessor(maskSize, maskSize);
          for (int i = 0; i < mask.length; i++) {
            if (mask[i] != 0) {
              bp.set(i, 128);
            }
          }
          ImageJUtils.display(maskTitle, bp);
        }
      }

      // Use the simulated cluster centres to create clusters of the desired size
      if (settings.clusterSimulation == 1 || settings.clusterSimulation == 2) {
        for (final double[] clusterCentre : clusterCentres) {
          final int clusterN = (int) Math.round((settings.clusterNumberStdDev > 0)
              ? settings.clusterNumber + gauss.sample() * settings.clusterNumberStdDev
              : settings.clusterNumber);
          if (clusterN < 1) {
            continue;
          }
          if (clusterN == 1) {
            // No need for a cluster around a point
            xyz.add(clusterCentre);
            statsRadius.add(0);
            statsSize.add(1);
          } else {
            // Generate N random points within a circle of the chosen cluster radius.
            // Locate the centre-of-mass and the average distance to the centre.
            final double[] com = new double[3];
            int size = 0;
            while (size < clusterN) {
              // Generate a random point within a circle uniformly
              // http://stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly
              final double t = 2.0 * Math.PI * rng.nextDouble();
              final double u = rng.nextDouble() + rng.nextDouble();
              final double r = settings.clusterRadius * ((u > 1) ? 2 - u : u);
              final double x = r * Math.cos(t);
              final double y = r * Math.sin(t);
              final double[] xy = new double[] {clusterCentre[0] + x, clusterCentre[1] + y};
              xyz.add(xy);
              for (int k = 0; k < 2; k++) {
                com[k] += xy[k];
              }
              size++;
            }
            // Add the distance of the points from the centre of the cluster.
            // Note this does not account for the movement due to precision.
            statsSize.add(size);
            if (size == 1) {
              statsRadius.add(0);
            } else {
              for (int k = 0; k < 2; k++) {
                com[k] /= size;
              }
              while (size > 0) {
                final double dx = xyz.get(xyz.size() - size)[0] - com[0];
                final double dy = xyz.get(xyz.size() - size)[1] - com[1];
                statsRadius.add(Math.sqrt(dx * dx + dy * dy));
                size--;
              }
            }
          }
        }
      }
    } else {
      // Random distribution
      for (int i = 0; i < settings.numberOfMolecules; i++) {
        xyz.add(dist.next());
      }
    }

    // The Gaussian sigma should be applied so the overall distance from the centre
    // ( sqrt(x^2+y^2) ) has a standard deviation of sigmaS?
    final double sigma1D = settings.sigmaS / Math.sqrt(2);

    // Show optional histograms
    StoredDataStatistics intraDistances = null;
    StoredData blinks = null;
    if (settings.showHistograms) {
      final int capacity = (int) (xyz.size() * settings.blinkingRate);
      intraDistances = new StoredDataStatistics(capacity);
      blinks = new StoredData(capacity);
    }

    final Statistics statsSigma = new Statistics();
    for (int i = 0; i < xyz.size(); i++) {
      int occurrences = getBlinks(rng, settings.blinkingRate);
      if (blinks != null) {
        blinks.add(occurrences);
      }

      final int size = settings.molecules.size();

      // Get coordinates in nm
      final double[] moleculeXyz = xyz.get(i);

      if (bp != null && occurrences > 0) {
        bp.putPixel((int) Math.round(moleculeXyz[0] / maskScale),
            (int) Math.round(moleculeXyz[1] / maskScale), 255);
      }

      while (occurrences-- > 0) {
        final double[] localisationXy = Arrays.copyOf(moleculeXyz, 2);
        // Add random precision
        if (sigma1D > 0) {
          final double dx = gauss.sample() * sigma1D;
          final double dy = gauss.sample() * sigma1D;
          localisationXy[0] += dx;
          localisationXy[1] += dy;
          if (!dist.isWithinXy(localisationXy)) {
            continue;
          }
          // Calculate mean-squared displacement
          statsSigma.add(dx * dx + dy * dy);
        }
        final double x = localisationXy[0];
        final double y = localisationXy[1];
        settings.molecules.add(new Molecule(x, y, i, 1));

        // Store in pixels
        final float xx = (float) (x / nmPerPixel);
        final float yy = (float) (y / nmPerPixel);
        final float[] params = PeakResult.createParams(0, 0, xx, yy, 0);
        settings.results.add(i + 1, (int) xx, (int) yy, 0, 0, 0, 0, params, null);
      }

      if (settings.molecules.size() > size) {
        count++;
        if (intraDistances != null) {
          final int newCount = settings.molecules.size() - size;
          if (newCount == 1) {
            // No intra-molecule distances
            continue;
          }

          // Get the distance matrix between these molecules
          final double[][] matrix = new double[newCount][newCount];
          for (int ii = size, x = 0; ii < settings.molecules.size(); ii++, x++) {
            for (int jj = size + 1, y = 1; jj < settings.molecules.size(); jj++, y++) {
              final double d2 = settings.molecules.get(ii).distance2(settings.molecules.get(jj));
              matrix[x][y] = matrix[y][x] = d2;
            }
          }

          // Get the maximum distance for particle linkage clustering of this molecule
          double max = 0;
          for (int x = 0; x < newCount; x++) {
            // Compare to all-other molecules and get the minimum distance
            // needed to join at least one
            double linkDistance = Double.POSITIVE_INFINITY;
            for (int y = 0; y < newCount; y++) {
              if (x == y) {
                continue;
              }
              if (matrix[x][y] < linkDistance) {
                linkDistance = matrix[x][y];
              }
            }
            // Check if this is larger
            if (max < linkDistance) {
              max = linkDistance;
            }
          }
          intraDistances.add(Math.sqrt(max));
        }
      }
    }
    settings.results.end();

    if (bp != null) {
      ImageJUtils.display(maskTitle, bp);
    }

    log("Simulation results");
    log("  * Molecules = %d (%d activated)", xyz.size(), count);
    log("  * Blinking rate = %s",
        MathUtils.rounded((double) settings.molecules.size() / xyz.size(), 4));
    log("  * Precision (Mean-displacement) = %s nm",
        (statsSigma.getN() > 0) ? MathUtils.rounded(Math.sqrt(statsSigma.getMean()), 4) : "0");
    if (intraDistances != null) {
      if (intraDistances.getN() == 0) {
        log("  * Mean Intra-Molecule particle linkage distance = 0 nm");
        log("  * Fraction of inter-molecule particle linkage @ 0 nm = 0 %%");
      } else {
        plot(blinks, "Blinks/Molecule", true);
        final double[][] intraHist =
            plot(intraDistances, "Intra-molecule particle linkage distance", false);

        // Determine 95th and 99th percentile
        // Will not be null as we requested a non-integer histogram.
        int p99 = intraHist[0].length - 1;
        final double limit1 = 0.99 * intraHist[1][p99];
        final double limit2 = 0.95 * intraHist[1][p99];
        while (intraHist[1][p99] > limit1 && p99 > 0) {
          p99--;
        }
        int p95 = p99;
        while (intraHist[1][p95] > limit2 && p95 > 0) {
          p95--;
        }

        log("  * Mean Intra-Molecule particle linkage distance = %s nm"
            + " (95%% = %s, 99%% = %s, 100%% = %s)", MathUtils.rounded(intraDistances.getMean(), 4),
            MathUtils.rounded(intraHist[0][p95], 4), MathUtils.rounded(intraHist[0][p99], 4),
            MathUtils.rounded(intraHist[0][intraHist[0].length - 1], 4));

        if (settings.distanceAnalysis) {
          performDistanceAnalysis(intraHist, p99);
        }
      }
    }
    if (settings.clusterSimulation > 0) {
      log("  * Cluster number = %s +/- %s", MathUtils.rounded(statsSize.getMean(), 4),
          MathUtils.rounded(statsSize.getStandardDeviation(), 4));
      log("  * Cluster radius = %s +/- %s nm (mean distance to centre-of-mass)",
          MathUtils.rounded(statsRadius.getMean(), 4),
          MathUtils.rounded(statsRadius.getStandardDeviation(), 4));
    }
  }

  @Nullable
  private static double[][] plot(DoubleData stats, String label, boolean integerBins) {
    final String title = TITLE + " " + label;

    if (integerBins) {
      // The histogram is not need for the return statement
      new HistogramPlotBuilder(title, stats, label).setMinBinWidth(1).show();
      return null;
    }

    // Show a cumulative histogram so that the bin size is not relevant
    final double[][] hist = MathUtils.cumulativeHistogram(stats.values(), false);

    // Create the axes
    final double[] xValues = hist[0];
    final double[] yValues = hist[1];

    // Plot
    final Plot2 plot = new Plot2(title, label, "Frequency", xValues, yValues);
    ImageJUtils.display(title, plot);

    return hist;
  }

  private void performDistanceAnalysis(double[][] intraHist, int p99) {
    // We want to know the fraction of distances between molecules at the 99th percentile
    // that are intra- rather than inter-molecule.
    // Do single linkage clustering of closest pair at this distance and count the number of
    // links that are inter and intra.

    // Convert molecules for clustering
    final ArrayList<ClusterPoint> points = new ArrayList<>(settings.molecules.size());
    for (final Molecule m : settings.molecules) {
      // Precision was used to store the molecule ID
      points.add(ClusterPoint.newClusterPoint((int) m.precision, m.x, m.y, m.photons));
    }
    final ClusteringEngine engine = new ClusteringEngine(Prefs.getThreads(),
        ClusteringAlgorithm.PARTICLE_SINGLE_LINKAGE, SimpleImageJTrackProgress.getInstance());
    IJ.showStatus("Clustering to check inter-molecule distances");
    engine.setTrackJoins(true);
    final List<Cluster> clusters = engine.findClusters(points, intraHist[0][p99]);
    IJ.showStatus("");
    if (clusters != null) {
      final double[] intraIdDistances = engine.getIntraIdDistances();
      final double[] interIdDistances = engine.getInterIdDistances();

      final int all = interIdDistances.length + intraIdDistances.length;

      log("  * Fraction of inter-molecule particle linkage @ %s nm = %s %%",
          MathUtils.rounded(intraHist[0][p99], 4),
          (all > 0) ? MathUtils.rounded(100.0 * interIdDistances.length / all, 4) : "0");

      // Show a double cumulative histogram plot
      final double[][] intraIdHist = MathUtils.cumulativeHistogram(intraIdDistances, false);
      final double[][] interIdHist = MathUtils.cumulativeHistogram(interIdDistances, false);

      // Plot
      final String title = TITLE + " molecule linkage distance";
      final Plot2 plot = new Plot2(title, "Distance", "Frequency", intraIdHist[0], intraIdHist[1]);
      double max = (intraIdHist[1].length > 0) ? intraIdHist[1][intraIdHist[1].length - 1] : 0;
      if (interIdHist[1].length > 0) {
        max = FastMath.max(max, interIdHist[1][interIdHist[1].length - 1]);
      }
      plot.setLimits(0, intraIdHist[0][intraIdHist[0].length - 1], 0, max);
      plot.setColor(Color.blue);
      plot.addPoints(interIdHist[0], interIdHist[1], Plot.LINE);
      plot.setColor(Color.black);
      ImageJUtils.display(title, plot);
    } else {
      log("Aborted clustering to check inter-molecule distances");
    }
  }

  private static double distance2(double[] centre1, double[] centre2) {
    final double dx = centre1[0] - centre2[0];
    final double dy = centre1[1] - centre2[1];
    return dx * dx + dy * dy;
  }

  /**
   * Fill the given mask with the fill value using a circle with the specified centre and radius.
   *
   * @param mask the mask
   * @param maskSize the mask size
   * @param cx the cx
   * @param cy the cy
   * @param roiRadius the roi radius
   * @param fill the fill
   */
  private static void fillMask(int[] mask, final int maskSize, final int cx, final int cy,
      final int roiRadius, final int fill) {
    int minx = cx - roiRadius;
    int maxx = cx + roiRadius;
    final int miny = cy - roiRadius;
    final int maxy = cy + roiRadius;
    final int r2 = roiRadius * roiRadius;

    // Pre-calculate x range
    if (minx < 0) {
      minx = 0;
    }
    if (maxx >= maskSize) {
      maxx = maskSize - 1;
    }
    if (minx > maxx) {
      return;
    }

    int count = 0;
    final int[] dx2 = new int[roiRadius * 2 + 1];
    for (int x = minx; x <= maxx; x++) {
      dx2[count++] = (cx - x) * (cx - x);
    }

    if (miny < 0 || maxy >= maskSize) {
      for (int y = miny, dy = -roiRadius; y <= maxy; y++, dy++) {
        if (y < 0) {
          continue;
        }
        if (y >= maskSize) {
          break;
        }
        final int limit = r2 - (dy * dy);
        for (int i = y * maskSize + minx, nn = 0; nn < count; i++, nn++) {
          if (dx2[nn] <= limit) {
            mask[i] = fill;
          }
        }
      }
    } else {
      for (int y = miny, dy = -roiRadius; y <= maxy; y++, dy++) {
        final int limit = r2 - (dy * dy);
        for (int i = y * maskSize + minx, nn = 0; nn < count; i++, nn++) {
          if (dx2[nn] <= limit) {
            mask[i] = fill;
          }
        }
      }
    }
  }

  private int getBlinks(UniformRandomProvider rng, double averageBlinks) {
    switch (settings.blinkingDistribution) {
      case 3:
        // Binomial distribution
        final int trials = (int) Math.round(averageBlinks);
        return new InverseTransformDiscreteSampler(rng,
            new BinomialDiscreteInverseCumulativeProbabilityFunction(trials, settings.pvalue))
                .sample();

      case 2:
        return (int) Math.round(averageBlinks);
      case 1:
        return StandardFluorophoreSequenceModel.getBlinks(true, rng, averageBlinks);
      default:
        return StandardFluorophoreSequenceModel.getBlinks(false, rng, averageBlinks);
    }
  }

  private boolean showSimulationDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

    gd.addMessage("Simulate a random distribution of molecules.");

    gd.addNumericField("Molecules", settings.numberOfMolecules, 0);
    gd.addNumericField("Simulation_size (um)", settings.simulationSize, 2);
    gd.addNumericField("Blinking_rate", settings.blinkingRate, 2);
    gd.addChoice("Blinking_distribution", Settings.BLINKING_DISTRIBUTION,
        settings.blinkingDistribution);
    gd.addNumericField("Average_precision (nm)", settings.sigmaS, 2);
    gd.addCheckbox("Show_histograms", settings.showHistograms);
    gd.addCheckbox("Distance_analysis", settings.distanceAnalysis);

    gd.addChoice("Cluster_simulation", Settings.CLUSTER_SIMULATION, settings.clusterSimulation);
    gd.addNumericField("Cluster_number", settings.clusterNumber, 2);
    gd.addNumericField("Cluster_variation (SD)", settings.clusterNumberStdDev, 2);
    gd.addNumericField("Cluster_radius", settings.clusterRadius, 2);
    gd.addCheckbox("Show_cluster_mask", settings.showClusterMask);

    gd.addHelp(HelpUrls.getUrl("pc-palm-molecules"));
    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.numberOfMolecules = (int) Math.abs(gd.getNextNumber());
    settings.simulationSize = Math.abs(gd.getNextNumber());
    settings.blinkingRate = Math.abs(gd.getNextNumber());
    settings.blinkingDistribution = gd.getNextChoiceIndex();
    settings.sigmaS = Math.abs(gd.getNextNumber());
    settings.showHistograms = gd.getNextBoolean();
    settings.distanceAnalysis = gd.getNextBoolean();
    settings.clusterSimulation = gd.getNextChoiceIndex();
    settings.clusterNumber = Math.abs(gd.getNextNumber());
    settings.clusterNumberStdDev = Math.abs(gd.getNextNumber());
    settings.clusterRadius = Math.abs(gd.getNextNumber());
    settings.showClusterMask = gd.getNextBoolean();

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Molecules", settings.numberOfMolecules);
      ParameterUtils.isAboveZero("Simulation size", settings.simulationSize);
      ParameterUtils.isEqualOrAbove("Blinking rate", settings.blinkingRate, 1);
      ParameterUtils.isEqualOrAbove("Cluster number", settings.clusterNumber, 1);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return getPValue();
  }

  private boolean getPValue() {
    if (settings.blinkingDistribution == 3) {
      final GenericDialog gd = new GenericDialog(TITLE);
      gd.addMessage("Binomial distribution requires a p-value");
      gd.addSlider("p-Value (%)", 0, 100, 100 * settings.pvalue);
      gd.showDialog();
      if (gd.wasCanceled()) {
        return false;
      }
      settings.pvalue = FastMath.max(FastMath.min(gd.getNextNumber(), 100), 0) / 100;
    }
    return true;
  }

  private static class FrameProcedure implements PeakResultProcedure {
    int start;
    int end;

    public FrameProcedure(int start, int end) {
      this.start = start;
      this.end = end;
    }

    @Override
    public void execute(PeakResult result) {
      if (start > result.getFrame()) {
        start = result.getFrame();
      }
      if (end < result.getEndFrame()) {
        end = result.getEndFrame();
      }
    }
  }

  /**
   * Get the lifetime of the results using the earliest and latest frames and the calibrated
   * exposure time.
   */
  private void getLifetime() {
    if (settings.results.isEmpty()) {
      settings.seconds = 0;
      return;
    }
    int start = settings.results.getFirstFrame();
    int end = start;
    final FrameProcedure p = new FrameProcedure(start, end);
    settings.results.forEach(p);
    start = p.start;
    end = p.end;
    settings.seconds =
        (end - start + 1) * settings.results.getCalibrationReader().getExposureTime() / 1000;
  }

  private boolean createImage(List<Molecule> molecules) {
    if (molecules.isEmpty()) {
      return false;
    }

    // Find the limits of the image
    settings.minx = settings.maxx = molecules.get(0).x;
    settings.miny = settings.maxy = molecules.get(0).y;

    // Compute limits
    for (int i = molecules.size(); i-- > 0;) {
      final Molecule m1 = molecules.get(i);
      if (settings.minx > m1.x) {
        settings.minx = m1.x;
      } else if (settings.maxx < m1.x) {
        settings.maxx = m1.x;
      }
      if (settings.miny > m1.y) {
        settings.miny = m1.y;
      } else if (settings.maxy < m1.y) {
        settings.maxy = m1.y;
      }
    }

    // Assign to a grid
    final int gridSize = 500;
    final double xBinSize = (settings.maxx - settings.minx) / gridSize;
    final double yBinSize = (settings.maxy - settings.miny) / gridSize;
    final int nXBins = 1 + (int) ((settings.maxx - settings.minx) / xBinSize);
    final int nybins = 1 + (int) ((settings.maxy - settings.miny) / yBinSize);
    final Molecule[][] grid = new Molecule[nXBins][nybins];
    for (final Molecule m : molecules) {
      final int xbin = (int) ((m.x - settings.minx) / xBinSize);
      final int ybin = (int) ((m.y - settings.miny) / yBinSize);
      // Build a single linked list
      m.next = grid[xbin][ybin];
      grid[xbin][ybin] = m;
    }

    // Find the minimum distance between molecules.
    double distanceMin = Double.POSITIVE_INFINITY;

    IJ.showStatus("Computing minimum distance ...");
    IJ.showProgress(0);
    final Molecule[] neighbours = new Molecule[5];
    for (int ybin = 0, currentIndex = 0, finalIndex = nXBins * nybins; ybin < nybins; ybin++) {
      for (int xbin = 0; xbin < nXBins; xbin++) {
        IJ.showProgress(currentIndex, finalIndex);

        for (Molecule m1 = grid[xbin][ybin]; m1 != null; m1 = m1.next) {
          // Build a list of which cells to compare up to a maximum of 4
          // @formatter:off
          //      | 0,0 | 1,0
          // ------------+-----
          // -1,1 | 0,1 | 1,1
          // @formatter:on

          int count = 0;
          neighbours[count++] = m1.next;

          if (ybin < nybins - 1) {
            neighbours[count++] = grid[xbin][ybin + 1];
            if (xbin > 0) {
              neighbours[count++] = grid[xbin - 1][ybin + 1];
            }
          }
          if (xbin < nXBins - 1) {
            neighbours[count++] = grid[xbin + 1][ybin];
            if (ybin < nybins - 1) {
              neighbours[count++] = grid[xbin + 1][ybin + 1];
            }
          }

          // Compare to neighbours
          while (count-- > 0) {
            for (Molecule m2 = neighbours[count]; m2 != null; m2 = m2.next) {
              if (distanceMin > m1.distance2(m2)) {
                distanceMin = m1.distance2(m2);
              }
            }
          }
        }
      }
    }
    ImageJUtils.finished();

    settings.nmPerPixel = Math.sqrt(distanceMin);
    log("Minimum distance between molecules = %g nm", settings.nmPerPixel);
    if (settings.nmPerPixel == 0 && settings.nmPerPixelLimit == 0) {
      IJ.error(TITLE, "Zero minimum distance between molecules - please enter a nm/pixel limit "
          + "for image reconstruction");
      return false;
    }

    if (settings.nmPerPixel < settings.nmPerPixelLimit) {
      log("Minimum distance adjusted to user defined limit %.2f nm", settings.nmPerPixelLimit);
      settings.nmPerPixel = settings.nmPerPixelLimit;
    }

    // Compute the minimum size we can use and stay within memory.
    // Assume a 4um x 4um analysis section with 800nm feature radius.
    final double limit = (4000.0 + 800.0 + 1) / 4096;
    if (settings.nmPerPixel < limit) {
      log("Minimum distance adjusted to %.2f nm to fit in memory", limit);
      settings.nmPerPixel = limit;
    }
    log("X-range %.2f - %.2f : Y-range %.2f - %.2f (nm)", settings.minx, settings.maxx,
        settings.miny, settings.maxy);

    // Construct a binary representation

    final String namePrefix =
        settings.results.getName() + " " + ((settings.binaryImage) ? "Binary" : "Count") + " Image";

    double lowResNmPerPixel;
    final double xrange = settings.maxx - settings.minx;
    final double yrange = settings.maxy - settings.miny;
    if (xrange > 0 || yrange > 0) {
      lowResNmPerPixel = FastMath.max(xrange, yrange) / settings.lowResolutionImageSize;
    } else {
      // The resolution does not matter
      lowResNmPerPixel = 100;
    }
    final ImagePlus imp = displayImage(namePrefix + " (low res)", molecules, settings.minx,
        settings.miny, settings.maxx, settings.maxy, lowResNmPerPixel, false, settings.binaryImage);

    // Add an ROI to allow the user to select regions. PC-PALM recommends 2x2 to 4x4 um^2
    final int size = (int) (settings.roiSizeInUm * 1000.0 / lowResNmPerPixel);
    imp.setRoi(new Rectangle(0, 0, size, size));

    if (settings.showHighResolutionImage) {
      displayImage(namePrefix + " (high res)", molecules, settings.minx, settings.miny,
          settings.maxx, settings.maxy, settings.nmPerPixel, false, settings.binaryImage);
    }

    // Store the molecules, the data range and the dMin.
    // This will be used by a filter plugin that crops sections from the image for PC analysis
    this.settings.molecules = molecules;

    return true;
  }

  /**
   * Draw an image of the molecules.
   *
   * @param molecules the molecules
   * @param minx the minx
   * @param miny the miny
   * @param maxx the maxx
   * @param maxy the maxy
   * @param nmPerPixel the nm per pixel
   * @param checkBounds Set to true to check the molecules is within the bounds
   * @param binary the binary
   * @return the image
   */
  static ImageProcessor drawImage(List<Molecule> molecules, double minx, double miny, double maxx,
      double maxy, double nmPerPixel, boolean checkBounds, boolean binary) {
    final double scalex = maxx - minx;
    final double scaley = maxy - miny;
    final int width = (int) Math.round(scalex / nmPerPixel) + 1;
    final int height = (int) Math.round(scaley / nmPerPixel) + 1;

    // ***
    // The PC-PALM + PLoS One papers describe using a binary image.
    // However both papers provide MatLab code where the number of particles is
    // calculated using sum(sum(I)). This indicates a non binary image could be input
    // to the routine to calculate the correlation function g(r).
    // ***
    if (binary) {
      final byte[] data = new byte[width * height];
      for (final Molecule m : molecules) {
        if (checkBounds) {
          if (m.x < minx || m.x >= maxx || m.y < miny || m.y >= maxy) {
            continue;
          }
        }

        // Shift to the origin. This makes the image more memory efficient.
        final int x = (int) Math.round((m.x - minx) / nmPerPixel);
        final int y = (int) Math.round((m.y - miny) / nmPerPixel);
        final int index = y * width + x;

        // Construct a binary image
        data[index] = (byte) 1;
      }

      final ByteProcessor ip = new ByteProcessor(width, height, data, null);
      ip.setMinAndMax(0, 1);
      return ip;
    }

    final short[] data = new short[width * height];
    for (final Molecule m : molecules) {
      if (checkBounds && (m.x < minx || m.x >= maxx || m.y < miny || m.y >= maxy)) {
        continue;
      }

      // Shift to the origin. This makes the image more memory efficient.
      final int x = (int) Math.round((m.x - minx) / nmPerPixel);
      final int y = (int) Math.round((m.y - miny) / nmPerPixel);
      final int index = y * width + x;

      // Construct a count image
      data[index]++;
    }

    final ShortProcessor ip = new ShortProcessor(width, height, data, null);
    ip.setMinAndMax(0, MathUtils.max(data));
    return ip;
  }

  /**
   * Display the image.
   *
   * @param title the title
   * @param ip the image processor
   * @param nmPerPixel the nm per pixel
   * @return the image
   */
  static ImagePlus displayImage(String title, ImageProcessor ip, double nmPerPixel) {
    final ImagePlus imp = ImageJUtils.display(title, ip);
    final Calibration cal = new Calibration();
    cal.setUnit("um");
    cal.pixelWidth = cal.pixelHeight = nmPerPixel / 1000;
    imp.setCalibration(cal);
    return imp;
  }

  /**
   * Display an image of the molecules.
   *
   * @param title the title
   * @param molecules the molecules
   * @param minx the minx
   * @param miny the miny
   * @param maxx the maxx
   * @param maxy the maxy
   * @param nmPerPixel the nm per pixel
   * @param checkBounds Set to true to check the molecules is within the bounds
   * @param binary the binary
   * @return the image
   */
  static ImagePlus displayImage(String title, List<Molecule> molecules, double minx, double miny,
      double maxx, double maxy, double nmPerPixel, boolean checkBounds, boolean binary) {
    final ImageProcessor ip =
        drawImage(molecules, minx, miny, maxx, maxy, nmPerPixel, checkBounds, binary);
    return displayImage(title, ip, nmPerPixel);
  }

  /**
   * Allow optimisation using Apache Commons Math 3 Optimiser.
   */
  private abstract class SkewNormalOptimiserFunction extends SkewNormalFunction {
    public SkewNormalOptimiserFunction(double[] parameters) {
      super(parameters);
    }

    protected TDoubleArrayList x;
    protected TDoubleArrayList y;

    public void addData(float[] x, float[] y) {
      this.x = new TDoubleArrayList();
      this.y = new TDoubleArrayList();
      for (int i = 0; i < x.length; i++) {
        this.x.add(x[i]);
        this.y.add(y[i]);
      }
    }

    public double[] calculateTarget() {
      return y.toArray();
    }

    public double[] calculateWeights() {
      final double[] w = new double[y.size()];
      for (int i = 0; i < y.size(); i++) {
        w[i] = 1;
      }
      return w;
    }
  }

  /**
   * Allow optimisation using Apache Commons Math 3 Gradient Optimiser.
   */
  private class SkewNormalDifferentiableFunction extends SkewNormalOptimiserFunction
      implements MultivariateVectorFunction {
    // Adapted from http://commons.apache.org/proper/commons-math/userguide/optimization.html
    // Use the deprecated API since the new one is not yet documented.

    public SkewNormalDifferentiableFunction(double[] parameters) {
      super(parameters);
    }

    private double[][] jacobian(double[] variables) {
      // Compute the gradients using numerical differentiation

      final double[][] jacobian = new double[x.size()][4];
      final double delta = 0.001;
      final double[][] d = new double[variables.length][variables.length];
      for (int i = 0; i < variables.length; i++) {
        // Should the delta be changed for each parameter?
        d[i][i] = delta * Math.abs(variables[i]);
      }
      for (int i = 0; i < jacobian.length; ++i) {
        final double x = this.x.getQuick(i);
        final double value = evaluate(x, variables);
        for (int j = 0; j < variables.length; j++) {
          final double value2 = evaluate(x, variables[0] + d[0][j], variables[1] + d[1][j],
              variables[2] + d[2][j], variables[3] + d[3][j]);
          jacobian[i][j] = (value2 - value) / d[j][j];
        }
      }
      return jacobian;
    }

    @Override
    public double[] value(double[] variables) {
      final double[] values = new double[x.size()];
      for (int i = 0; i < values.length; i++) {
        values[i] = evaluate(x.getQuick(i), variables);
      }
      return values;
    }
  }

  /**
   * Allow optimisation using Apache Commons Math 3 Simplex.
   */
  private class SkewNormalMultivariateFunction extends SkewNormalOptimiserFunction
      implements MultivariateFunction {
    public SkewNormalMultivariateFunction(double[] parameters) {
      super(parameters);
    }

    @Override
    public double value(double[] point) {
      // Objective function is to minimise sum-of-squares
      double ss = 0;
      for (int i = x.size(); i-- > 0;) {
        final double dx = y.get(i) - evaluate(x.getQuick(i), point);
        ss += dx * dx;
      }
      return ss;
    }
  }

  /**
   * Log a message to the IJ log window.
   *
   * @param format the format
   * @param args the args
   */
  private static void log(String format, Object... args) {
    ImageJUtils.log(format, args);
  }

  /**
   * Output a spacer to the ImageJ log.
   */
  static void logSpacer() {
    log("-=-=-=-=-=-=-");
  }
}
