/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.RandomUtils;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.PSFEstimatorSettings;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.engine.FitConfiguration;
import uk.ac.sussex.gdsc.smlm.engine.FitEngine;
import uk.ac.sussex.gdsc.smlm.engine.FitEngineConfiguration;
import uk.ac.sussex.gdsc.smlm.engine.FitJob;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.ij.IJImageSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.ij.utils.ImageJImageConverter;
import uk.ac.sussex.gdsc.smlm.results.AggregatedImageSource;
import uk.ac.sussex.gdsc.smlm.results.ImageSource;
import uk.ac.sussex.gdsc.smlm.results.InterlacedImageSource;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResultStore;
import uk.ac.sussex.gdsc.smlm.results.PeakResults;
import uk.ac.sussex.gdsc.smlm.results.ThreadSafePeakResults;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Roi;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.text.TextWindow;

import org.apache.commons.lang3.concurrent.ConcurrentRuntimeException;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.StatisticalSummary;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.stat.inference.TestUtils;
import org.apache.commons.rng.simple.RandomSource;

import java.awt.Rectangle;
import java.util.Collection;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Iteratively fits local maxima using a 2D Gaussian until the PSF converges.
 */
public class PsfEstimator implements PlugInFilter, ThreadSafePeakResults {
  private static final String TITLE = "PSF Estimator";
  private static final int FLAGS = DOES_16 | DOES_8G | DOES_32 | NO_CHANGES;

  private static final int TRY_AGAIN = 0;
  private static final int COMPLETE = 1;
  private static final int INSUFFICIENT_PEAKS = 2;
  private static final int ABORTED = 3;
  private static final int EXCEPTION = 4;
  private static final int BAD_ESTIMATE = 5;

  private static TextWindow resultsWindow;

  private double initialPeakStdDev0 = 1;
  private double initialPeakStdDev1 = 1;
  private double initialPeakAngle;

  private boolean extraOptions;
  private static int optionIntegrateFrames = 1;
  private int integrateFrames = 1;
  private static boolean optionInterlacedData;
  private static int optionDataStart = 1;
  private static int optionDataBlock = 1;
  private static int optionDataSkip;
  private boolean interlacedData;
  private int dataStart = 1;
  private int dataBlock = 1;
  private int dataSkip;

  private FitEngineConfiguration config;
  private PSFEstimatorSettings.Builder settings;

  private ImagePlus imp;

  // Required for the significance tests
  private static final int ANGLE = 0;
  private static final int X = 1;
  private static final int Y = 2;
  private static final int XY = 3;
  private static final String[] NAMES = {"Angle", "X SD", "Y SD"};
  private final DescriptiveStatistics[] sampleNew = new DescriptiveStatistics[3];
  private final DescriptiveStatistics[] sampleOld = new DescriptiveStatistics[3];
  private final boolean[] ignore = new boolean[3];

  @Override
  public int setup(String arg, ImagePlus imp) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    extraOptions = ImageJUtils.isExtraOptions();
    if (imp == null) {
      IJ.noImage();
      return DONE;
    }

    settings = SettingsManager.readPsfEstimatorSettings(0).toBuilder();
    // Reset
    if (IJ.controlKeyDown()) {
      config = new FitEngineConfiguration();
      final Calibration calibration = SettingsManager.readCalibration(0);
      config.getFitConfiguration().setCalibration(calibration);
    } else {
      config = SettingsManager.readFitEngineConfiguration(0);
    }

    final Roi roi = imp.getRoi();
    if (roi != null && roi.getType() != Roi.RECTANGLE) {
      IJ.error("Rectangular ROI required");
      return DONE;
    }

    return showDialog(imp);
  }

  /**
   * Show dialog.
   *
   * @param imp the imp
   * @return the int
   */
  private int showDialog(ImagePlus imp) {
    // Keep class variables for the parameters we are fitting
    final FitConfiguration fitConfig = config.getFitConfiguration();
    initialPeakStdDev0 = 1;
    initialPeakStdDev1 = 1;
    initialPeakAngle = 0;
    try {
      initialPeakStdDev0 = fitConfig.getInitialXSd();
      initialPeakStdDev1 = fitConfig.getInitialYSd();
      initialPeakAngle = fitConfig.getInitialAngle();
    } catch (final IllegalStateException ex) {
      // Ignore this as the current PSF is not a 2 axis and theta Gaussian PSF
    }

    if (!extraOptions) {
      interlacedData = false;
      integrateFrames = 1;
    }

    this.imp = imp;

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);
    gd.addMessage("Estimate 2D Gaussian to fit maxima");

    gd.addNumericField("Initial_StdDev0", initialPeakStdDev0, 3);
    gd.addNumericField("Initial_StdDev1", initialPeakStdDev1, 3);
    gd.addNumericField("Initial_Angle", initialPeakAngle, 3);
    gd.addNumericField("Number_of_peaks", settings.getNumberOfPeaks(), 0);

    // pValue sets the smallest significance level probability level at which they are said to be
    // different.
    // i.e. p <= pValue they are different

    // lower pValue means harder to be found different.
    // lower pValue means easier to be found the same.

    gd.addNumericField("p-Value", settings.getPValue(), 4);
    gd.addCheckbox("Update_preferences", settings.getUpdatePreferences());
    gd.addCheckbox("Log_progress", settings.getDebugPsfEstimator());
    gd.addCheckbox("Iterate", settings.getIterate());
    gd.addCheckbox("Show_histograms", settings.getShowHistograms());
    gd.addNumericField("Histogram_bins", settings.getHistogramBins(), 0);

    PeakFit.addCameraOptions(gd, fitConfig);
    PeakFit.addPsfOptions(gd, fitConfig);
    final PeakFit.SimpleFitEngineConfigurationProvider provider =
        new PeakFit.SimpleFitEngineConfigurationProvider(config);
    PeakFit.addDataFilterOptions(gd, provider);
    PeakFit.addSearchOptions(gd, provider);
    PeakFit.addBorderOptions(gd, provider);
    PeakFit.addFittingOptions(gd, provider);

    if (extraOptions) {
      gd.addCheckbox("Interlaced_data", optionInterlacedData);
      gd.addSlider("Integrate_frames", 1, 5, optionIntegrateFrames);
    }

    gd.addMessage("--- Gaussian fitting ---");

    gd.addChoice("Fit_solver", SettingsManager.getFitSolverNames(),
        fitConfig.getFitSolver().ordinal());

    // Parameters specific to each Fit solver are collected in a second dialog

    gd.addNumericField("Fail_limit", config.getFailuresLimit(), 0);
    gd.addNumericField("Pass_rate", config.getPassRate(), 2);
    gd.addCheckbox("Include_neighbours", config.isIncludeNeighbours());
    gd.addSlider("Neighbour_height", 0.01, 1, config.getNeighbourHeightThreshold());
    gd.addSlider("Residuals_threshold", 0.01, 1, config.getResidualsThreshold());

    gd.addMessage(
        "--- Peak filtering ---\nDiscard fits that shift; are too low; or expand/contract");
    gd.addCheckbox("Smart_filter", fitConfig.isSmartFilter());
    gd.addCheckbox("Disable_simple_filter", fitConfig.isDisableSimpleFilter());
    gd.addSlider("Shift_factor", 0.01, 2, fitConfig.getCoordinateShiftFactor());
    gd.addNumericField("Signal_strength", fitConfig.getSignalStrength(), 2);
    gd.addNumericField("Min_photons", fitConfig.getMinPhotons(), 0);
    gd.addSlider("Min_width_factor", 0, 0.99, fitConfig.getMinWidthFactor());
    gd.addSlider("Width_factor", 1, 4.5, fitConfig.getMaxWidthFactor());
    PeakFit.addPrecisionOptions(gd, new PeakFit.SimpleFitConfigurationProvider(fitConfig));

    gd.showDialog();

    if (gd.wasCanceled() || !readDialog(gd)) {
      return DONE;
    }

    return FLAGS;
  }

  private boolean readDialog(ExtendedGenericDialog gd) {
    initialPeakStdDev0 = gd.getNextNumber();
    initialPeakStdDev1 = gd.getNextNumber();
    initialPeakAngle = gd.getNextNumber();

    settings.setNumberOfPeaks((int) gd.getNextNumber());
    settings.setPValue(gd.getNextNumber());
    settings.setUpdatePreferences(gd.getNextBoolean());
    settings.setDebugPsfEstimator(gd.getNextBoolean());
    settings.setIterate(gd.getNextBoolean());
    settings.setShowHistograms(gd.getNextBoolean());
    settings.setHistogramBins((int) gd.getNextNumber());

    final FitConfiguration fitConfig = config.getFitConfiguration();
    fitConfig.setPsfType(PeakFit.getPsfTypeValues()[gd.getNextChoiceIndex()]);
    config.setDataFilterType(gd.getNextChoiceIndex());
    config.setDataFilter(gd.getNextChoiceIndex(), Math.abs(gd.getNextNumber()), false, 0);
    config.setSearch(gd.getNextNumber());
    config.setBorder(gd.getNextNumber());
    config.setFitting(gd.getNextNumber());

    if (extraOptions) {
      interlacedData = optionInterlacedData = gd.getNextBoolean();
      integrateFrames = optionIntegrateFrames = (int) gd.getNextNumber();
    }

    fitConfig.setFitSolver(gd.getNextChoiceIndex());
    config.setFailuresLimit((int) gd.getNextNumber());
    config.setPassRate(gd.getNextNumber());
    config.setIncludeNeighbours(gd.getNextBoolean());
    config.setNeighbourHeightThreshold(gd.getNextNumber());
    config.setResidualsThreshold(gd.getNextNumber());

    fitConfig.setSmartFilter(gd.getNextBoolean());
    fitConfig.setDisableSimpleFilter(gd.getNextBoolean());
    fitConfig.setCoordinateShiftFactor(gd.getNextNumber());
    fitConfig.setSignalStrength(gd.getNextNumber());
    fitConfig.setMinPhotons(gd.getNextNumber());
    fitConfig.setMinWidthFactor(gd.getNextNumber());
    fitConfig.setMaxWidthFactor(gd.getNextNumber());
    fitConfig.setPrecisionThreshold(gd.getNextNumber());

    gd.collectOptions();

    if (gd.invalidNumber()) {
      return false;
    }

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Initial SD0", initialPeakStdDev0);
      ParameterUtils.isAboveZero("Initial SD1", initialPeakStdDev1);
      ParameterUtils.isPositive("Initial angle", initialPeakAngle);
      ParameterUtils.isPositive("Number of peaks", settings.getNumberOfPeaks());
      ParameterUtils.isAboveZero("P-value", settings.getPValue());
      ParameterUtils.isEqualOrBelow("P-value", settings.getPValue(), 0.5);
      if (settings.getShowHistograms()) {
        ParameterUtils.isAboveZero("Histogram bins", settings.getHistogramBins());
      }
      ParameterUtils.isAboveZero("Search width", config.getSearch());
      ParameterUtils.isAboveZero("Fitting width", config.getFitting());
      // Can be negative to disable
      // Parameters.isPositive("Failures limit", config.getFailuresLimit());
      ParameterUtils.isPositive("Neighbour height threshold", config.getNeighbourHeightThreshold());
      ParameterUtils.isPositive("Residuals threshold", config.getResidualsThreshold());
      ParameterUtils.isPositive("Coordinate Shift factor", fitConfig.getCoordinateShiftFactor());
      ParameterUtils.isPositive("Signal strength", fitConfig.getSignalStrength());
      ParameterUtils.isPositive("Min photons", fitConfig.getMinPhotons());
      ParameterUtils.isPositive("Min width factor", fitConfig.getMinWidthFactor());
      ParameterUtils.isPositive("Width factor", fitConfig.getMaxWidthFactor());
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    if (fitConfig.getPsfType() == PSFType.ONE_AXIS_GAUSSIAN_2D && fitConfig.isFixedPsf()) {
      final String msg =
          "ERROR: A width-fitting function must be selected (i.e. not fixed-width fitting)";
      IJ.error(TITLE, msg);
      log(msg);
      return false;
    }

    SettingsManager.writeSettings(config, 0);

    if (!PeakFit.configureSmartFilter(config, 0)) {
      return false;
    }
    if (!PeakFit.configureDataFilter(config, 0)) {
      return false;
    }
    if (!PeakFit.configureFitSolver(config, IJImageSource.getBounds(imp), null, 0)) {
      return false;
    }

    // Extra parameters are needed for interlaced data
    if (interlacedData) {
      gd = new ExtendedGenericDialog(TITLE);
      gd.addMessage("Interlaced data requires a repeating pattern of frames to process.\n"
          + "Describe the regular repeat of the data:\n \n"
          + "Start = The first frame that contains data\n"
          + "Block = The number of continuous frames containing data\n"
          + "Skip = The number of continuous frames to ignore before the next data\n \n"
          + "E.G. 2:9:1 = Data was imaged from frame 2 for 9 frames, 1 frame to ignore,"
          + " then repeat.");
      gd.addNumericField("Start", optionDataStart, 0);
      gd.addNumericField("Block", optionDataBlock, 0);
      gd.addNumericField("Skip", optionDataSkip, 0);

      gd.showDialog();
      if (gd.wasCanceled()) {
        return false;
      }

      if (!gd.wasCanceled()) {
        dataStart = (int) gd.getNextNumber();
        dataBlock = (int) gd.getNextNumber();
        dataSkip = (int) gd.getNextNumber();

        if (dataStart > 0 && dataBlock > 0 && dataSkip > 0) {
          // Store options for next time
          optionInterlacedData = true;
          optionDataStart = dataStart;
          optionDataBlock = dataBlock;
          optionDataSkip = dataSkip;
        }
      } else {
        interlacedData = false;
      }
    }

    return true;
  }

  @Override
  public void run(ImageProcessor ip) {
    int result;
    while (true) {
      result = estimatePsf();
      if (settings.getIterate() && result == TRY_AGAIN) {
        continue;
      }
      break;
    }
    if (result < INSUFFICIENT_PEAKS) {
      log("Finished. Check the table for final parameters");

      // Only save if successful
      if (settings.getUpdatePreferences()) {
        SettingsManager.writeSettings(settings.build());
        SettingsManager.writeSettings(config, 0);
      }
    }
  }

  private int estimatePsf() {
    log("Estimating PSF ... Press escape to abort");

    final PeakFit fitter = createFitter();

    // Use the fit configuration to generate a Gaussian function to test what is being evaluated
    final Gaussian2DFunction gf = config.getFitConfiguration().createGaussianFunction(1, 1, 1);
    createResultsWindow();
    ignore[ANGLE] = !gf.evaluatesAngle();
    ignore[X] = !gf.evaluatesSD0();
    ignore[Y] = !gf.evaluatesSD1();

    final double[] params = new double[] {gf.evaluatesAngle() ? initialPeakAngle : 0,
        gf.evaluatesSD0() ? initialPeakStdDev0 : 0, gf.evaluatesSD1() ? initialPeakStdDev1 : 0, 0,
        0};
    final double[] paramsDev = new double[3];
    final boolean[] identical = new boolean[4];
    final double[] p = new double[] {Double.NaN, Double.NaN, Double.NaN, Double.NaN};

    addToResultTable(0, 0, params, paramsDev, p);

    if (!calculateStatistics(fitter, params, paramsDev)) {
      return (ImageJUtils.isInterrupted()) ? ABORTED : INSUFFICIENT_PEAKS;
    }

    if (!addToResultTable(1, size(), params, paramsDev, p)) {
      return BAD_ESTIMATE;
    }

    boolean tryAgain = false;

    int iteration = 2;
    do {
      if (!calculateStatistics(fitter, params, paramsDev)) {
        return (ImageJUtils.isInterrupted()) ? ABORTED : INSUFFICIENT_PEAKS;
      }

      try {
        for (int i = 0; i < 3; i++) {
          getP(i, p, identical);
        }

        if (!ignore[Y]) {
          getPairedP(sampleNew[X], sampleNew[Y], XY, p, identical);
        }

        if (!addToResultTable(iteration++, size(), params, paramsDev, p)) {
          return BAD_ESTIMATE;
        }

        if ((ignore[ANGLE] || identical[ANGLE] || identical[XY]) && (ignore[X] || identical[X])
            && (ignore[Y] || identical[Y])) {
          tryAgain = checkAngleSignificance() || checkXySignificance(identical);

          // Update recommended values. Only use if significant
          params[X] = sampleNew[X].getMean();
          params[Y] = (!ignore[Y] && !identical[XY]) ? sampleNew[Y].getMean() : params[X];
          params[ANGLE] = (!ignore[ANGLE]) ? sampleNew[ANGLE].getMean() : 0;

          // update starting configuration
          initialPeakAngle = (float) params[ANGLE];
          initialPeakStdDev0 = (float) params[X];
          initialPeakStdDev1 = (float) params[Y];

          if (settings.getUpdatePreferences()) {
            config.getFitConfiguration().setInitialPeakStdDev0((float) params[X]);
            config.getFitConfiguration().setInitialPeakStdDev1((float) params[Y]);
            config.getFitConfiguration().setInitialAngle((float) params[ANGLE]);
          }

          break;
        }

        if (IJ.escapePressed()) {
          IJ.beep();
          IJ.showStatus("Aborted");
          return ABORTED;
        }
      } catch (final Exception ex) {
        Logger.getLogger(getClass().getName()).log(Level.WARNING, "Error while estimating the PSF",
            ex);
        return EXCEPTION;
      }
    }
    while (true);

    return (tryAgain) ? TRY_AGAIN : COMPLETE;
  }

  private boolean checkAngleSignificance() {
    boolean tryAgain = false;
    if (ignore[ANGLE]) {
      return tryAgain;
    }

    // The angle is relative to the major axis (X).
    // It could be close to 0, 90 or 180 to allow it to be ignored in favour of a free circular
    // function.

    final double[] angles = sampleNew[ANGLE].getValues();

    for (final double testAngle : new double[] {90, 0, 180}) {
      // The angle will be in the 0-180 domain.
      // We need to compute the Statistical summary around the testAngle.
      StatisticalSummary sampleStats;
      if (testAngle == 0 || testAngle == 180) {
        final SummaryStatistics stats = new SummaryStatistics();
        final boolean zeroAngle = (testAngle == 0);
        for (double a : angles) {
          if (zeroAngle) {
            // Convert to -90-90 domain
            if (a > 90) {
              a -= 180;
            }
            // Convert to 90-270 domain
          } else if (a < 90) {
            a += 180;
          }
          stats.addValue(a);
        }
        sampleStats = stats;
      } else {
        // Already in the 0-180 domain around the angle 90
        sampleStats = sampleNew[ANGLE];
      }

      final double p = TestUtils.tTest(testAngle, sampleStats);
      if (p > settings.getPValue()) {
        log("NOTE: Angle is not significant: %g ~ %g (p=%g) => Re-run with fixed zero angle",
            sampleStats.getMean(), testAngle, p);
        ignore[ANGLE] = true;
        config.getFitConfiguration().setPsfType(PSFType.TWO_AXIS_GAUSSIAN_2D);
        tryAgain = true;
        break;
      }
      debug("  NOTE: Angle is significant: %g !~ %g (p=%g)", sampleNew[ANGLE].getMean(), testAngle,
          p);
    }
    return tryAgain;
  }

  private boolean checkXySignificance(boolean[] identical) {
    boolean tryAgain = false;
    if (identical[XY]) {
      log("NOTE: X-width and Y-width are not significantly different: %g ~ %g => "
          + "Re-run with circular function", sampleNew[X].getMean(), sampleNew[Y].getMean());
      config.getFitConfiguration().setPsfType(PSFType.ONE_AXIS_GAUSSIAN_2D);
      tryAgain = true;
    }
    return tryAgain;
  }

  private void getP(int index, double[] pvalue, boolean[] identical) {
    getP(sampleNew[index], sampleOld[index], index, pvalue, identical);
  }

  private void getP(StatisticalSummary sample1, StatisticalSummary sample2, int index,
      double[] pvalue, boolean[] identical) {
    if (sample1.getN() < 2) {
      return;
    }

    // The number returned is the smallest significance level at which one can reject the null
    // hypothesis that the mean of the paired differences is 0 in favor of the two-sided alternative
    // that the mean paired difference is not equal to 0. For a one-sided test, divide the returned
    // value by 2
    pvalue[index] = TestUtils.tTest(sample1, sample2);
    identical[index] = (pvalue[index] > settings.getPValue());
  }

  private void getPairedP(DescriptiveStatistics sample1, DescriptiveStatistics sample2, int index,
      double[] pvalue, boolean[] identical) {
    if (sample1.getN() < 2) {
      return;
    }

    // The number returned is the smallest significance level at which one can reject the null
    // hypothesis that the mean of the paired differences is 0 in favor of the two-sided alternative
    // that the mean paired difference is not equal to 0. For a one-sided test, divide the returned
    // value by 2
    pvalue[index] = TestUtils.pairedTTest(sample1.getValues(), sample2.getValues());
    identical[index] = (pvalue[index] > settings.getPValue());
  }

  private boolean calculateStatistics(PeakFit fitter, double[] params, double[] paramsDev) {
    debug("  Fitting PSF");

    swapStatistics();

    // Create the fit engine using the PeakFit plugin
    final FitConfiguration fitConfig = config.getFitConfiguration();
    fitConfig.setInitialAngle((float) params[0]);
    fitConfig.setInitialPeakStdDev0((float) params[1]);
    fitConfig.setInitialPeakStdDev1((float) params[2]);

    final ImageStack stack = imp.getImageStack();
    final Rectangle roi = stack.getProcessor(1).getRoi();

    ImageSource source = new IJImageSource(imp);
    // Allow interlaced data by wrapping the image source
    if (interlacedData) {
      source = new InterlacedImageSource(source, dataStart, dataBlock, dataSkip);
    }
    // Allow frame aggregation by wrapping the image source
    if (integrateFrames > 1) {
      source = new AggregatedImageSource(source, integrateFrames);
    }
    fitter.initialiseImage(source, roi, true);

    fitter.addPeakResults(this);
    fitter.initialiseFitting();

    final FitEngine engine = fitter.createFitEngine();

    // Use random slices
    final int[] slices = new int[stack.getSize()];
    for (int i = 0; i < slices.length; i++) {
      slices[i] = i + 1;
    }
    RandomUtils.shuffle(slices, RandomSource.create(RandomSource.SPLIT_MIX_64));

    IJ.showStatus("Fitting ...");

    // Use multi-threaded code for speed
    int sliceIndex;
    for (sliceIndex = 0; sliceIndex < slices.length; sliceIndex++) {
      final int slice = slices[sliceIndex];
      IJ.showProgress(size(), settings.getNumberOfPeaks());

      final ImageProcessor ip = stack.getProcessor(slice);
      ip.setRoi(roi); // stack processor does not set the bounds required by ImageConverter
      final FitJob job = new FitJob(slice, ImageJImageConverter.getData(ip), roi);
      engine.run(job);

      if (sampleSizeReached() || ImageJUtils.isInterrupted()) {
        break;
      }
    }

    if (ImageJUtils.isInterrupted()) {
      IJ.showProgress(1);
      engine.end(true);
      return false;
    }

    // Wait until we have enough results
    while (!sampleSizeReached() && !engine.isQueueEmpty()) {
      IJ.showProgress(size(), settings.getNumberOfPeaks());
      try {
        Thread.sleep(50);
      } catch (final InterruptedException ex) {
        Thread.currentThread().interrupt();
        throw new ConcurrentRuntimeException("Unexpected interruption", ex);
      }
    }
    // End now if we have enough samples
    engine.end(sampleSizeReached());

    IJ.showStatus("");
    IJ.showProgress(1);

    // This count will be an over-estimate given that the provider is ahead of the consumer
    // in this multi-threaded system
    debug("  Processed %d/%d slices (%d peaks)", sliceIndex, slices.length, size());

    setParams(ANGLE, params, paramsDev, sampleNew[ANGLE]);
    setParams(X, params, paramsDev, sampleNew[X]);
    setParams(Y, params, paramsDev, sampleNew[Y]);

    if (settings.getShowHistograms()) {
      final HistogramPlotBuilder builder =
          new HistogramPlotBuilder(TITLE).setNumberOfBins(settings.getHistogramBins());
      final WindowOrganiser wo = new WindowOrganiser();
      for (int ii = 0; ii < 3; ii++) {
        if (sampleNew[ii].getN() == 0) {
          continue;
        }
        final StoredDataStatistics stats = StoredDataStatistics.create(sampleNew[ii].getValues());
        builder.setData(stats).setName(NAMES[ii])
            .setPlotLabel("Mean = " + MathUtils.rounded(stats.getMean()) + ". Median = "
                + MathUtils.rounded(sampleNew[ii].getPercentile(50)))
            .show(wo);
      }
      wo.tile();
    }

    if (size() < 2) {
      log("ERROR: Insufficient number of fitted peaks, terminating ...");
      return false;
    }
    return true;
  }

  private static void setParams(int index, double[] params, double[] paramsDev,
      DescriptiveStatistics sample) {
    if (sample.getN() > 0) {
      params[index] = sample.getMean();
      paramsDev[index] = sample.getStandardDeviation();
    }
  }

  private void swapStatistics() {
    sampleOld[ANGLE] = sampleNew[ANGLE];
    sampleOld[X] = sampleNew[X];
    sampleOld[Y] = sampleNew[Y];
  }

  private PeakFit createFitter() {
    return new PeakFit(config, null);
  }

  /**
   * Create the result window (if it is not available).
   */
  private static void createResultsWindow() {
    if (!ImageJUtils.isShowing(resultsWindow)) {
      resultsWindow = new TextWindow(TITLE + " Results", createResultsHeader(), "", 900, 300);
    }
  }

  private static String createResultsHeader() {
    final StringBuilder sb = new StringBuilder();
    sb.append("Iteration\t");
    sb.append("N-peaks\t");
    sb.append("Angle\t");
    sb.append("+/-\t");
    sb.append("p(Angle same)\t");
    sb.append("X SD\t");
    sb.append("+/-\t");
    sb.append("p(X same)\t");
    sb.append("Y SD\t");
    sb.append("+/-\t");
    sb.append("p(Y same)\t");
    sb.append("p(XY same)\t");
    return sb.toString();
  }

  private boolean addToResultTable(int iteration, int n, double[] params, double[] paramsDev,
      double[] pvalue) {
    final StringBuilder sb = new StringBuilder();
    sb.append(iteration).append('\t').append(n).append('\t');
    for (int i = 0; i < 3; i++) {
      sb.append(params[i]).append('\t');
      sb.append(paramsDev[i]).append('\t');
      sb.append(pvalue[i]).append('\t');
    }
    sb.append(pvalue[XY]).append('\t');
    resultsWindow.append(sb.toString());

    if (params[X] > imp.getWidth() || params[Y] > imp.getWidth()) {
      log("ERROR: Bad width estimation (try altering the peak validation parameters),"
          + " terminating ...");
      return false;
    }
    return true;
  }

  private void debug(String format, Object... args) {
    if (settings.getDebugPsfEstimator()) {
      log(format, args);
    }
  }

  private static void log(String format, Object... args) {
    IJ.log(String.format(format, args));
  }

  @Override
  public void begin() {
    sampleNew[ANGLE] = new DescriptiveStatistics();
    sampleNew[X] = new DescriptiveStatistics();
    sampleNew[Y] = new DescriptiveStatistics();
  }

  @Override
  public void add(PeakResult result) {
    add(result.getFrame(), result.getOrigX(), result.getOrigY(), result.getOrigValue(),
        result.getError(), result.getNoise(), result.getMeanIntensity(), result.getParameters(),
        result.getParameterDeviations());
  }

  @Override
  public synchronized void add(int peak, int origX, int origY, float origValue, double chiSquared,
      float noise, float meanIntensity, float[] params, float[] paramsStdDev) {
    if (!sampleSizeReached()) {
      if (!ignore[ANGLE]) {
        // Assume fitting a 2D Gaussian
        sampleNew[ANGLE].addValue(params[Gaussian2DFunction.ANGLE]);
      }
      sampleNew[X].addValue(params[PeakResult.X]);
      if (!ignore[Y]) {
        sampleNew[Y].addValue(params[PeakResult.Y]);
      }
    }
  }

  private boolean sampleSizeReached() {
    return size() >= settings.getNumberOfPeaks();
  }

  @Override
  public synchronized void addAll(PeakResult[] results) {
    for (final PeakResult result : results) {
      add(result.getFrame(), result.getOrigX(), result.getOrigY(), result.getOrigValue(),
          result.getError(), result.getNoise(), result.getMeanIntensity(), result.getParameters(),
          result.getParameterDeviations());
    }
  }

  @Override
  public synchronized void addAll(Collection<PeakResult> results) {
    for (final PeakResult result : results) {
      add(result.getFrame(), result.getOrigX(), result.getOrigY(), result.getOrigValue(),
          result.getError(), result.getNoise(), result.getMeanIntensity(), result.getParameters(),
          result.getParameterDeviations());
    }
  }

  @Override
  public void addAll(PeakResultStore results) {
    addAll(results.toArray());
  }

  @Override
  public int size() {
    return (int) sampleNew[X].getN();
  }

  @Override
  public void end() {
    // Do nothing
  }

  @Override
  public boolean isActive() {
    return true;
  }

  @Override
  public ImageSource getSource() {
    // Ignored
    return null;
  }

  @Override
  public void setBounds(Rectangle bounds) {
    // Ignored
  }

  @Override
  public Rectangle getBounds() {
    // Ignored
    return null;
  }

  @Override
  public void setCalibration(Calibration calibration) {
    // Ignored
  }

  @Override
  public Calibration getCalibration() {
    // Ignored
    return null;
  }

  @Override
  public void setConfiguration(String configuration) {
    // Ignored
  }

  @Override
  public String getConfiguration() {
    // Ignored
    return null;
  }

  @Override
  public void copySettings(PeakResults peakResults) {
    // Ignored
  }

  @Override
  public void setSource(ImageSource source) {
    // Ignored
  }

  @Override
  public String getName() {
    // Ignored
    return null;
  }

  @Override
  public void setName(String name) {
    // Ignored
  }

  @Override
  public PSF getPsf() {
    // Ignored
    return null;
  }

  @Override
  public void setPsf(PSF psf) {
    // Ignored
  }
}
