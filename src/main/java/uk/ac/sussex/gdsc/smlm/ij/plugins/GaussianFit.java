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

package uk.ac.sussex.gdsc.smlm.ij.plugins;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.measure.Measurements;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.Blitter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.process.ShortProcessor;
import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Rectangle;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicReference;
import uk.ac.sussex.gdsc.core.annotation.Nullable;
import uk.ac.sussex.gdsc.core.filters.FilteredNonMaximumSuppression;
import uk.ac.sussex.gdsc.core.ij.ImageJPluginLoggerHelper;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.OffsetPointRoi;
import uk.ac.sussex.gdsc.core.utils.ImageExtractor;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.SortUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitSolver;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.data.config.PsfProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.AngleUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.engine.FitConfiguration;
import uk.ac.sussex.gdsc.smlm.filters.BlockMeanFilter;
import uk.ac.sussex.gdsc.smlm.fitting.FitResult;
import uk.ac.sussex.gdsc.smlm.fitting.FitStatus;
import uk.ac.sussex.gdsc.smlm.fitting.Gaussian2DFitter;
import uk.ac.sussex.gdsc.smlm.function.gaussian.EllipticalGaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.ij.results.ImageJTablePeakResults;
import uk.ac.sussex.gdsc.smlm.ij.settings.Constants;
import uk.ac.sussex.gdsc.smlm.ij.utils.ImageJImageConverter;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;

/**
 * Fits the selected rectangular ROI using a 2D Gaussian.
 */
public class GaussianFit implements ExtendedPlugInFilter {
  private static final String TITLE = "Gaussian Fit";
  private static final int FLAGS = DOES_16 | DOES_8G | DOES_32 | FINAL_PROCESSING | SNAPSHOT;

  private ImagePlus imp;

  private int[] maxIndices;
  private FitResult fitResult;
  private double chiSquared;

  private ImageJTablePeakResults results;
  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    double smooth;
    int boxSize;
    float background;
    float peakHeight;
    float fractionAboveBackground;
    float peakWidth;
    int topN;
    boolean blockFindAlgorithm;
    boolean neighbourCheck;
    int border;
    int fitFunction;
    boolean fitBackground;
    boolean logProgress;
    int maxIterations;
    double relativeThreshold;
    double absoluteThreshold;
    boolean singleFit;
    int singleRegionSize;
    double initialPeakStdDev;
    boolean showDeviations;
    boolean filterResults;
    boolean showFit;

    Settings() {
      // Set defaults
      smooth = Prefs.get(Constants.smooth, 0);
      boxSize = (int) Prefs.get(Constants.boxSize, 1);
      background = (float) Prefs.get(Constants.background, 0);
      peakHeight = (float) Prefs.get(Constants.peakHeight, 0);
      fractionAboveBackground = (float) Prefs.get(Constants.fractionAboveBackground, 0);
      peakWidth = (float) Prefs.get(Constants.peakWidth, 0);
      topN = (int) Prefs.get(Constants.topN, 0);
      blockFindAlgorithm = Prefs.get(Constants.blockFindAlgorithm, true);
      neighbourCheck = Prefs.get(Constants.neighbourCheck, false);
      border = (int) Prefs.get(Constants.border, 0);
      fitFunction = (int) Prefs.get(Constants.fitFunction, 0);
      fitBackground = Prefs.get(Constants.fitBackground, true);
      logProgress = Prefs.get(Constants.logProgress, false);
      maxIterations = (int) Prefs.get(Constants.maxIterations, 20);
      relativeThreshold = Prefs.get(Constants.relativeThreshold, 1e-5);
      absoluteThreshold = Prefs.get(Constants.absoluteThreshold, 1e-10);
      singleFit = Prefs.get(Constants.singleFit, false);
      singleRegionSize = (int) Prefs.get(Constants.singleRegionSize, 10);
      initialPeakStdDev = Prefs.get(Constants.initialPeakStdDev0, 0);
      showDeviations = Prefs.get(Constants.showDeviations, false);
      filterResults = Prefs.get(Constants.filterResults, false);
      showFit = Prefs.get(Constants.showFit, false);
    }

    Settings(Settings source) {
      smooth = source.smooth;
      boxSize = source.boxSize;
      background = source.background;
      peakHeight = source.peakHeight;
      fractionAboveBackground = source.fractionAboveBackground;
      peakWidth = source.peakWidth;
      topN = source.topN;
      blockFindAlgorithm = source.blockFindAlgorithm;
      neighbourCheck = source.neighbourCheck;
      border = source.border;
      fitFunction = source.fitFunction;
      fitBackground = source.fitBackground;
      logProgress = source.logProgress;
      maxIterations = source.maxIterations;
      relativeThreshold = source.relativeThreshold;
      absoluteThreshold = source.absoluteThreshold;
      singleFit = source.singleFit;
      singleRegionSize = source.singleRegionSize;
      initialPeakStdDev = source.initialPeakStdDev;
      showDeviations = source.showDeviations;
      filterResults = source.filterResults;
      showFit = source.showFit;
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
      Prefs.set(Constants.smooth, smooth);
      Prefs.set(Constants.boxSize, boxSize);
      Prefs.set(Constants.background, background);
      Prefs.set(Constants.peakHeight, peakHeight);
      Prefs.set(Constants.fractionAboveBackground, fractionAboveBackground);
      Prefs.set(Constants.peakWidth, peakWidth);
      Prefs.set(Constants.topN, topN);
      Prefs.set(Constants.blockFindAlgorithm, blockFindAlgorithm);
      Prefs.set(Constants.neighbourCheck, neighbourCheck);
      Prefs.set(Constants.border, border);
      Prefs.set(Constants.fitFunction, fitFunction);
      Prefs.set(Constants.fitBackground, fitBackground);
      Prefs.set(Constants.logProgress, logProgress);
      Prefs.set(Constants.showDeviations, showDeviations);
      Prefs.set(Constants.filterResults, filterResults);
      Prefs.set(Constants.showFit, showFit);
      Prefs.set(Constants.maxIterations, maxIterations);
      Prefs.set(Constants.relativeThreshold, relativeThreshold);
      Prefs.set(Constants.absoluteThreshold, absoluteThreshold);
      Prefs.set(Constants.singleFit, singleFit);
      Prefs.set(Constants.singleRegionSize, singleRegionSize);
      Prefs.set(Constants.initialPeakStdDev0, initialPeakStdDev);
    }
  }

  private static class PsfTypeLoader {
    private static final PSFType[] psfTypeValues;
    private static final String[] psfTypeNames;

    static {
      //@formatter:off
      final EnumSet<PSFType> set = EnumSet.of(
          PSFType.ONE_AXIS_GAUSSIAN_2D,
          PSFType.TWO_AXIS_GAUSSIAN_2D,
          PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D);
      //@formatter:on
      psfTypeValues = set.toArray(new PSFType[set.size()]);
      psfTypeNames = new String[psfTypeValues.length];
      for (int i = 0; i < psfTypeValues.length; i++) {
        psfTypeNames[i] = PsfProtosHelper.getName(psfTypeValues[i]);
      }
    }
  }

  @Override
  public int setup(String arg, ImagePlus imp) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (imp == null) {
      IJ.noImage();
      return DONE;
    }

    final Roi roi = imp.getRoi();
    if (roi != null && roi.getType() != Roi.RECTANGLE) {
      IJ.error("Rectangular ROI required");
      return DONE;
    }

    this.imp = imp;
    if ("final".equals(arg)) {
      runFinal(imp.getProcessor());
      return DONE;
    }

    return FLAGS;
  }

  /**
   * Gets the PSF type values.
   *
   * @return the PSF type values
   */
  public static PSFType[] getPsfTypeValues() {
    return PsfTypeLoader.psfTypeValues;
  }


  /**
   * Gets the PSF type names.
   *
   * @return the PSF type names
   */
  public static String[] getPsfTypeNames() {
    return PsfTypeLoader.psfTypeNames;
  }

  @Override
  public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
    final double[] limits = getLimits(imp.getProcessor());
    final double minValue = limits[0];
    final double maxValue = limits[1];

    settings = Settings.load();

    if (settings.background > maxValue) {
      settings.background = (int) maxValue;
    }
    if (settings.background < minValue) {
      settings.background = (int) minValue;
    }

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("gaussian-fit"));

    gd.addMessage("Fit 2D Gaussian to identified maxima");

    gd.addMessage("--- Image smoothing ---\n" + "- Within a 2n+1 box\n");
    gd.addSlider("Smoothing", 0, 4.5, settings.smooth);

    gd.addMessage("--- Maxima identification ---\n" + "- Within a 2n+1 box\n");
    gd.addSlider("Box_size", 1, 15, settings.boxSize);
    gd.addSlider("Background", minValue, maxValue, settings.background);
    gd.addSlider("Min_height", 0, maxValue, settings.peakHeight);
    gd.addSlider("Fraction_above_background", 0, 1.01, settings.fractionAboveBackground);
    gd.addSlider("Min_width", 0, 20, settings.peakWidth);
    gd.addSlider("Top_N", 0, 20, settings.topN);
    gd.addCheckbox("Block_find_algorithm", settings.blockFindAlgorithm);
    gd.addCheckbox("Neighbour_check", settings.neighbourCheck);
    gd.addSlider("Border", 0, 15, settings.border);

    gd.addMessage("--- Gaussian fitting ---");
    gd.addChoice("PSF", getPsfTypeNames(), PsfProtosHelper.getName(getPsfType()));
    gd.addCheckbox("Fit_background", settings.fitBackground);
    gd.addNumericField("Max_iterations", settings.maxIterations, 0);
    gd.addNumericField("Relative_threshold", settings.relativeThreshold, -3);
    gd.addNumericField("Absolute_threshold", settings.absoluteThreshold, -3);
    gd.addCheckbox("Single_fit", settings.singleFit);
    gd.addNumericField("Single_region_size", settings.singleRegionSize, 0);
    gd.addNumericField("Initial_StdDev", settings.initialPeakStdDev, 3);
    gd.addCheckbox("Log_progress", settings.logProgress);
    gd.addCheckbox("Show_deviations", settings.showDeviations);
    gd.addCheckbox("Filter_results", settings.filterResults);
    gd.addCheckbox("Show_fit", settings.showFit);

    gd.addPreviewCheckbox(pfr);
    gd.addDialogListener(this::dialogItemChanged);

    gd.showDialog();

    if (gd.wasCanceled() || !dialogItemChanged(gd, null)) {
      imp.setOverlay(null);
      return DONE;
    }

    return FLAGS;
  }

  /**
   * Calculate the min/max limits for the image. The max is set at the 99th percentile of the data.
   *
   * @param ip the ip
   * @return The limits
   */
  private static double[] getLimits(ImageProcessor ip) {
    final ImageStatistics stats = ImageStatistics.getStatistics(ip, Measurements.MIN_MAX, null);
    final double[] limits = new double[] {stats.min, stats.max};

    // Use histogram to cover x% of the data
    final int[] data = ip.getHistogram();

    if (data == null) {
      return limits;
    }

    // 8/16 bit greyscale image. Set upper limit to the height of the 99% percentile.
    final int limit = (int) (99.0 * ip.getPixelCount() / 100.0);
    int count = 0;
    int index = 0;
    while (index < data.length) {
      count += data[index];
      if (count > limit) {
        break;
      }
      index++;
    }

    limits[1] = index;

    return limits;
  }

  private boolean dialogItemChanged(GenericDialog gd, @SuppressWarnings("unused") AWTEvent event) {
    settings.smooth = gd.getNextNumber();
    settings.boxSize = (int) gd.getNextNumber();
    settings.background = (float) gd.getNextNumber();
    settings.peakHeight = (float) gd.getNextNumber();
    settings.fractionAboveBackground = (float) gd.getNextNumber();
    settings.peakWidth = (float) gd.getNextNumber();
    settings.topN = (int) gd.getNextNumber();
    settings.blockFindAlgorithm = gd.getNextBoolean();
    settings.neighbourCheck = gd.getNextBoolean();
    settings.border = (int) gd.getNextNumber();

    settings.fitFunction = gd.getNextChoiceIndex();
    settings.fitBackground = gd.getNextBoolean();
    settings.maxIterations = (int) gd.getNextNumber();
    settings.relativeThreshold = gd.getNextNumber();
    settings.absoluteThreshold = gd.getNextNumber();
    settings.singleFit = gd.getNextBoolean();
    settings.singleRegionSize = (int) gd.getNextNumber();
    settings.initialPeakStdDev = gd.getNextNumber();
    settings.logProgress = gd.getNextBoolean();
    settings.showDeviations = gd.getNextBoolean();
    settings.filterResults = gd.getNextBoolean();
    settings.showFit = gd.getNextBoolean();

    if (gd.invalidNumber()) {
      return false;
    }

    // Check arguments
    try {
      ParameterUtils.isPositive("Smoothing", settings.smooth);
      ParameterUtils.isAboveZero("Box size", settings.boxSize);
      ParameterUtils.isPositive("Peak height", settings.peakHeight);
      ParameterUtils.isPositive("Fraction above background", settings.fractionAboveBackground);
      ParameterUtils.isPositive("Peak width", settings.peakWidth);
      ParameterUtils.isPositive("Top N", settings.topN);
      ParameterUtils.isPositive("Border", settings.border);
      ParameterUtils.isAboveZero("Relative threshold", settings.relativeThreshold);
      ParameterUtils.isAboveZero("Absolute threshold", settings.absoluteThreshold);
      ParameterUtils.isAboveZero("Max iterations", settings.maxIterations);
      ParameterUtils.isAboveZero("Single region size", settings.singleRegionSize);
      ParameterUtils.isPositive("Initial peak StdDev", settings.initialPeakStdDev);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    settings.save();

    if (!gd.getPreviewCheckbox().getState()) {
      imp.setOverlay(null);
    }

    return true;
  }

  @Override
  public void run(ImageProcessor ip) {
    final Rectangle bounds = ip.getRoi();

    // Crop to the ROI
    final float[] data = ImageJImageConverter.getData(ip);

    final int width = bounds.width;
    final int height = bounds.height;

    if (getSmooth() > 0) {
      // No need for a copy since we are using a snapshot buffer
      final BlockMeanFilter filter = new BlockMeanFilter();
      filter.stripedBlockFilter(data, width, height, (float) getSmooth());
    }

    maxIndices = getMaxima(data, width, height);

    if (settings.topN > 0 && maxIndices.length > settings.topN) {
      SortUtils.sortIndices(maxIndices, data, true);
      maxIndices = Arrays.copyOf(maxIndices, settings.topN);
    }

    // Show an overlay of the indices
    if (maxIndices.length > 0) {
      final int nMaxima = maxIndices.length;
      final float[] xpoints = new float[nMaxima];
      final float[] ypoints = new float[nMaxima];
      int count = 0;
      for (final int index : maxIndices) {
        final int x = index % width;
        final int y = index / width;
        xpoints[count] = 0.5f + bounds.x + x;
        ypoints[count] = 0.5f + bounds.y + y;
        count++;
      }

      setOverlay(nMaxima, xpoints, ypoints);
    } else {
      imp.setOverlay(null);
    }

    for (int index = data.length; index-- > 0;) {
      ip.setf(bounds.x + index % width, bounds.y + index / width, data[index]);
    }
  }

  /**
   * Show the points as an overlay.
   *
   * @param npoints the number of points
   * @param xpoints the x points
   * @param ypoints the y points
   */
  private void setOverlay(int npoints, float[] xpoints, float[] ypoints) {
    final PointRoi roi = new OffsetPointRoi(xpoints, ypoints, npoints);

    final Color strokeColor = Color.yellow;
    final Color fillColor = Color.green;

    roi.setStrokeColor(strokeColor);
    roi.setFillColor(fillColor);
    roi.setShowLabels(false);

    imp.setOverlay(roi, strokeColor, 2, fillColor);
  }

  /**
   * Perform fitting using the chosen maxima. Update the overlay if successful.
   *
   * @param ip The input image
   */
  private void runFinal(ImageProcessor ip) {
    ip.reset();
    final Rectangle bounds = ip.getRoi();

    // Crop to the ROI
    final float[] data = ImageJImageConverter.getData(ip);

    final int width = bounds.width;
    final int height = bounds.height;

    // Sort the maxima
    float[] smoothData = data;
    if (getSmooth() > 0) {
      // Smoothing destructively modifies the data so create a copy
      smoothData = Arrays.copyOf(data, width * height);
      final BlockMeanFilter filter = new BlockMeanFilter();
      if (settings.smooth <= settings.border) {
        filter.stripedBlockFilterInternal(smoothData, width, height, (float) settings.smooth);
      } else {
        filter.stripedBlockFilter(smoothData, width, height, (float) settings.smooth);
      }
    }
    SortUtils.sortIndices(maxIndices, smoothData, true);

    // Show the candidate peaks
    if (maxIndices.length > 0) {
      final String message = String.format("Identified %d peaks", maxIndices.length);
      if (isLogProgress()) {
        IJ.log(message);
        for (final int index : maxIndices) {
          IJ.log(String.format("  %.2f @ [%d,%d]", data[index], bounds.x + index % width,
              bounds.y + index / width));
        }
      }

      // Check whether to run if the number of peaks is large
      if (maxIndices.length > 10) {
        final GenericDialog gd = new GenericDialog("Warning");
        gd.addMessage(message + "\nDo you want to fit?");
        gd.showDialog();
        if (gd.wasCanceled()) {
          return;
        }
      }
    } else {
      IJ.log("No maxima identified");
      return;
    }

    results = new ImageJTablePeakResults(settings.showDeviations,
        imp.getTitle() + " [" + imp.getCurrentSlice() + "]");
    final CalibrationWriter cw = new CalibrationWriter();
    cw.setIntensityUnit(IntensityUnit.COUNT);
    cw.setDistanceUnit(DistanceUnit.PIXEL);
    cw.setAngleUnit(AngleUnit.RADIAN);
    results.setCalibration(cw.getCalibration());
    results.setPsf(PsfProtosHelper.getDefaultPsf(getPsfType()));
    results.setShowFittingData(true);
    results.setAngleUnit(AngleUnit.DEGREE);
    results.begin();

    // Perform the Gaussian fit
    long ellapsed = 0;

    final FloatProcessor renderedImage =
        settings.showFit ? new FloatProcessor(ip.getWidth(), ip.getHeight()) : null;

    if (!settings.singleFit) {
      if (isLogProgress()) {
        IJ.log("Combined fit");
      }

      // Estimate height from smoothed data
      final double[] estimatedHeights = new double[maxIndices.length];
      for (int i = 0; i < estimatedHeights.length; i++) {
        estimatedHeights[i] = smoothData[maxIndices[i]];
      }

      final FitConfiguration config = new FitConfiguration();
      setupPeakFiltering(config);

      final long time = System.nanoTime();
      final double[] params = fitMultiple(data, width, height, maxIndices, estimatedHeights);
      ellapsed = System.nanoTime() - time;

      if (params != null) {
        // Copy all the valid parameters into a new array
        final double[] validParams = new double[params.length];
        int count = 0;
        int validPeaks = 0;
        validParams[count++] = params[0];

        final double[] initialParams = convertParameters(fitResult.getInitialParameters());
        final double[] paramsDev = convertParameters(fitResult.getParameterDeviations());
        final Rectangle regionBounds = new Rectangle();

        final float[] xpoints = new float[maxIndices.length];
        final float[] ypoints = new float[maxIndices.length];
        int npoints = 0;

        for (int i = 1, n = 0; i < params.length;
            i += Gaussian2DFunction.PARAMETERS_PER_PEAK, n++) {
          final int y = maxIndices[n] / width;
          final int x = maxIndices[n] % width;

          // Check the peak is a good fit
          if (settings.filterResults
              && config.validatePeak(n, initialParams, params, paramsDev) != FitStatus.OK) {
            continue;
          }

          if (settings.showFit) {
            // Copy the valid parameters before there are adjusted to global bounds
            validPeaks++;
            for (int ii = i, j = 0; j < Gaussian2DFunction.PARAMETERS_PER_PEAK; ii++, j++) {
              validParams[count++] = params[ii];
            }
          }

          final double[] peakParams = extractParams(params, i);
          final double[] peakParamsDev = extractParams(paramsDev, i);

          addResult(bounds, regionBounds, peakParams, peakParamsDev, npoints, x, y,
              data[maxIndices[n]]);

          // Add fit result to the overlay - Coords are updated with the region offsets in addResult
          final double xf = peakParams[Gaussian2DFunction.X_POSITION];
          final double yf = peakParams[Gaussian2DFunction.Y_POSITION];
          xpoints[npoints] = (float) xf;
          ypoints[npoints] = (float) yf;
          npoints++;
        }

        setOverlay(npoints, xpoints, ypoints);

        // Draw the fit
        if (validPeaks != 0) {
          addToImage(bounds.x, bounds.y, renderedImage, validParams, validPeaks, width, height);
        }
      } else {
        if (isLogProgress()) {
          IJ.log("Failed to fit " + TextUtils.pleural(maxIndices.length, "peak") + ": "
              + getReason(fitResult));
        }
        imp.setOverlay(null);
      }
    } else {
      if (isLogProgress()) {
        IJ.log("Individual fit");
      }

      int npoints = 0;
      final float[] xpoints = new float[maxIndices.length];
      final float[] ypoints = new float[maxIndices.length];

      // Extract each peak and fit individually
      final ImageExtractor ie = ImageExtractor.wrap(data, width, height);
      float[] region = null;
      final Gaussian2DFitter gf = createGaussianFitter(settings.filterResults);
      double[] validParams = null;

      final ShortProcessor renderedImageCount =
          settings.showFit ? new ShortProcessor(ip.getWidth(), ip.getHeight()) : null;

      for (int n = 0; n < maxIndices.length; n++) {
        final int y = maxIndices[n] / width;
        final int x = maxIndices[n] % width;

        final long time = System.nanoTime();
        final Rectangle regionBounds = ie.getBoxRegionBounds(x, y, settings.singleRegionSize);
        region = ie.crop(regionBounds, region);

        final int newIndex = (y - regionBounds.y) * regionBounds.width + x - regionBounds.x;

        if (isLogProgress()) {
          IJ.log("Fitting peak " + (n + 1));
        }

        final double[] peakParams = fitSingle(gf, region, regionBounds.width, regionBounds.height,
            newIndex, smoothData[maxIndices[n]]);
        ellapsed += System.nanoTime() - time;

        // Output fit result
        if (peakParams != null) {
          if (settings.showFit) {
            // Copy the valid parameters before there are adjusted to global bounds
            validParams = peakParams.clone();
          }

          double[] peakParamsDev = null;
          if (settings.showDeviations) {
            peakParamsDev = convertParameters(fitResult.getParameterDeviations());
          }

          addResult(bounds, regionBounds, peakParams, peakParamsDev, n, x, y, data[maxIndices[n]]);

          // Add fit result to the overlay - Coords are updated with the region offsets in addResult
          final double xf = peakParams[Gaussian2DFunction.X_POSITION];
          final double yf = peakParams[Gaussian2DFunction.Y_POSITION];
          xpoints[npoints] = (float) xf;
          ypoints[npoints] = (float) yf;
          npoints++;

          // Draw the fit
          if (settings.showDeviations) {
            final int ox = bounds.x + regionBounds.x;
            final int oy = bounds.y + regionBounds.y;
            addToImage(ox, oy, renderedImage, validParams, 1, regionBounds.width,
                regionBounds.height);
            addCount(ox, oy, renderedImageCount, regionBounds.width, regionBounds.height);
          }
        } else if (isLogProgress()) {
          IJ.log("Failed to fit peak " + (n + 1) + ": " + getReason(fitResult));
        }
      }

      // Update the overlay
      if (npoints > 0) {
        setOverlay(npoints, xpoints, ypoints);
      } else {
        imp.setOverlay(null);
      }
      // Create the mean
      if (settings.showFit) {
        for (int i = renderedImageCount.getPixelCount(); i-- > 0;) {
          final int count = renderedImageCount.get(i);
          if (count > 1) {
            renderedImage.setf(i, renderedImage.getf(i) / count);
          }
        }
      }
    }

    results.end();

    if (renderedImage != null) {
      ImageJUtils.display(TITLE, renderedImage);
    }

    if (isLogProgress()) {
      IJ.log("Time = " + (ellapsed / 1000000.0) + "ms");
    }
  }

  /**
   * Adds the function to the image at the specified origin.
   *
   * @param ox the x-origin
   * @param oy the y-origin
   * @param renderedImage the rendered image
   * @param params the function parameters
   * @param npeaks the number of peaks
   * @param width the function width
   * @param height the function height
   */
  private void addToImage(int ox, int oy, final FloatProcessor renderedImage, double[] params,
      int npeaks, int width, int height) {
    final EllipticalGaussian2DFunction f = new EllipticalGaussian2DFunction(npeaks, width, height);
    invertParameters(params);
    final FloatProcessor fp = new FloatProcessor(width, height, f.computeValues(params));
    // Insert into a full size image
    renderedImage.copyBits(fp, ox, oy, Blitter.ADD);
  }

  /**
   * Adds a count of 1 to the image at the specified origin.
   *
   * @param ox the x-origin
   * @param oy the y-origin
   * @param image the image
   * @param width the function width
   * @param height the function height
   */
  private static void addCount(int ox, int oy, final ShortProcessor image, int width, int height) {
    image.setRoi(ox, oy, width, height);
    image.add(1);
  }

  /**
   * Gets the reason.
   *
   * @param fitResult the fit result
   * @return the reason
   */
  private static String getReason(FitResult fitResult) {
    if (fitResult == null || fitResult.getStatus() == null) {
      return "unknown reason";
    }
    final FitStatus status = fitResult.getStatus();
    return status.toString().toLowerCase(Locale.US).replace("_", " ");
  }

  private static double[] extractParams(double[] params, int index) {
    if (params == null) {
      return null;
    }

    // 0 is the background. Then the peaks are packed.
    final double[] p = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    p[0] = params[0];
    System.arraycopy(params, index, p, 1, Gaussian2DFunction.PARAMETERS_PER_PEAK);
    return p;
  }

  private void addResult(Rectangle bounds, Rectangle regionBounds, double[] params,
      double[] paramsDev, int n, int x, int y, float value) {
    x += bounds.x;
    y += bounds.y;
    params[Gaussian2DFunction.X_POSITION] += bounds.x + regionBounds.x;
    params[Gaussian2DFunction.Y_POSITION] += bounds.y + regionBounds.y;
    results.add(n + 1, x, y, value, chiSquared, 0f, 0f, SimpleArrayUtils.toFloat(params),
        SimpleArrayUtils.toFloat(paramsDev));
  }

  /**
   * Find the indices of the maxima using the currently configured parameters.
   *
   * <p>Data must be arranged in yx block order, i.e. height rows of width.
   *
   * @param data the data
   * @param width the width
   * @param height the height
   * @return Indices of the maxima
   */
  public int[] getMaxima(float[] data, int width, int height) {
    // Find maxima
    final FilteredNonMaximumSuppression nms = new FilteredNonMaximumSuppression();
    nms.setBackground(getBackground());
    nms.setFractionAboveBackground(getFractionAboveBackground());
    nms.setMinimumHeight(getPeakHeight());
    nms.setMinimumWidth(getPeakWidth());
    nms.setNeighbourCheck(isNeighbourCheck());

    int[] indices;
    if (isBlockFindAlgorithm()) {
      if (getBorder() > 0) {
        indices = nms.blockFindInternal(data, width, height, getBoxSize(), getBorder());
      } else {
        indices = nms.blockFind(data, width, height, getBoxSize());
      }
    } else if (getBorder() > 0) {
      indices = nms.maxFindInternal(data, width, height, getBoxSize(), getBorder());
    } else {
      indices = nms.maxFind(data, width, height, getBoxSize());
    }
    return indices;
  }

  /**
   * Fits a 2D Gaussian to the given data. Fits all the specified peaks.
   *
   * <p>Data must be arranged in yx block order, i.e. height rows of width.
   *
   * <p>Note: The fit coordinates should be offset by 0.5 if the input data represents pixels
   *
   * @param data the data
   * @param width the width
   * @param height the height
   * @param maxIndices Indices of the data to fit
   * @return Array containing the fitted curve data: The first value is the Background. The
   *         remaining values are Amplitude, PosX, PosY, StdDevX, StdDevY for each fitted peak. If
   *         elliptical fitting is performed the values are Amplitude, Angle, PosX, PosY, StdDevX,
   *         StdDevY for each fitted peak Null if no fit is possible.
   */
  public @Nullable double[] fitMultiple(float[] data, int width, int height, int[] maxIndices) {
    return fitMultiple(data, width, height, maxIndices, null);
  }

  /**
   * Fits a 2D Gaussian to the given data. Fits all the specified peaks.
   *
   * <p>Data must be arranged in yx block order, i.e. height rows of width.
   *
   * <p>Note: The fit coordinates should be offset by 0.5 if the input data represents pixels
   *
   * @param data the data
   * @param width the width
   * @param height the height
   * @param maxIndices Indices of the data to fit
   * @param estimatedHeights Estimated heights for the peaks (input from smoothed data)
   * @return Array containing the fitted curve data: The first value is the Background. The
   *         remaining values are Amplitude, PosX, PosY, StdDevX, StdDevY for each fitted peak. If
   *         elliptical fitting is performed the values are Amplitude, Angle, PosX, PosY, StdDevX,
   *         StdDevY for each fitted peak Null if no fit is possible.
   */
  @Nullable
  private double[] fitMultiple(float[] data, int width, int height, int[] maxIndices,
      double[] estimatedHeights) {
    if (data == null || data.length != width * height) {
      return null;
    }

    if (maxIndices == null || maxIndices.length == 0) {
      return null;
    }

    final Gaussian2DFitter gf = createGaussianFitter(false);
    this.fitResult =
        gf.fit(SimpleArrayUtils.toDouble(data), width, height, maxIndices, estimatedHeights);
    if (fitResult.getStatus() == FitStatus.OK) {
      chiSquared = fitResult.getError();
      final double[] params = fitResult.getParameters();
      convertParameters(params);
      return params;
    }

    return null;
  }

  @Nullable
  private double[] convertParameters(double[] params) {
    if (params == null) {
      return null;
    }

    // Convert coordinates with +0.5 pixel offset
    // Convert radians to degrees (if elliptical fitting)
    final int n = params.length / Gaussian2DFunction.PARAMETERS_PER_PEAK;
    for (int i = 0, j = 0; i < n; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK) {
      params[j + Gaussian2DFunction.X_POSITION] += 0.5;
      params[j + Gaussian2DFunction.Y_POSITION] += 0.5;
      if (isEllipticalFitting()) {
        params[j + Gaussian2DFunction.ANGLE] *= 180.0 / Math.PI;
      }
    }
    return params;
  }

  @Nullable
  private double[] invertParameters(double[] params) {
    if (params == null) {
      return null;
    }

    // Convert coordinates with -0.5 pixel offset
    // Convert degrees to radians (if elliptical fitting)
    final int n = params.length / Gaussian2DFunction.PARAMETERS_PER_PEAK;
    for (int i = 0, j = 0; i < n; i++, j += Gaussian2DFunction.PARAMETERS_PER_PEAK) {
      params[j + Gaussian2DFunction.X_POSITION] -= 0.5;
      params[j + Gaussian2DFunction.Y_POSITION] -= 0.5;
      if (isEllipticalFitting()) {
        params[j + Gaussian2DFunction.ANGLE] *= Math.PI / 180.0;
      }
    }
    return params;
  }

  /**
   * Fits a 2D Gaussian to the given data. Fits all the specified peaks.
   *
   * <p>Data must be arranged in yx block order, i.e. height rows of width.
   *
   * @param gf the gf
   * @param data the data
   * @param width the width
   * @param height the height
   * @param index Index of the data to fit
   * @param estimatedHeight Estimated height for the peak (input from smoothed data)
   * @return Array containing the fitted curve data: The first value is the Background. The
   *         remaining values are Amplitude, PosX, PosY, StdDevX, StdDevY for each fitted peak. Null
   *         if no fit is possible.
   */
  @Nullable
  private double[] fitSingle(Gaussian2DFitter gf, float[] data, int width, int height, int index,
      double estimatedHeight) {
    this.fitResult = gf.fit(SimpleArrayUtils.toDouble(data), width, height, new int[] {index},
        new double[] {estimatedHeight});
    if (fitResult.getStatus() == FitStatus.OK) {
      chiSquared = fitResult.getError();
      final double[] params = fitResult.getParameters();
      convertParameters(params);
      // Check the fit is within the data
      if (params[Gaussian2DFunction.X_POSITION] < 0 || params[Gaussian2DFunction.X_POSITION] > width
          || params[Gaussian2DFunction.Y_POSITION] < 0
          || params[Gaussian2DFunction.Y_POSITION] > height) {
        fitResult = new FitResult(FitStatus.OUTSIDE_FIT_REGION, fitResult.getDegreesOfFreedom(),
            fitResult.getError(), fitResult.getInitialParameters(), fitResult.getParameters(),
            fitResult.getParameterDeviations(), fitResult.getNumberOfPeaks(),
            fitResult.getNumberOfFittedParameters(), fitResult.getStatusData(),
            fitResult.getIterations(), fitResult.getEvaluations());
        return null;
      }
      return params;
    }

    return null;
  }

  private Gaussian2DFitter createGaussianFitter(boolean simpleFiltering) {
    final FitConfiguration config = new FitConfiguration();
    config.setFitSolver(FitSolver.LVM_LSE);
    config.setPsf(PsfProtosHelper.getDefaultPsf(getPsfType()));
    config.setMaxIterations(getMaxIterations());
    config.setRelativeThreshold(settings.relativeThreshold);
    config.setAbsoluteThreshold(settings.absoluteThreshold);
    config.setInitialPeakStdDev(getInitialPeakStdDev());
    config.setComputeDeviations(settings.showDeviations);

    // Set-up peak filtering only for single fitting
    config.setDisableSimpleFilter(!simpleFiltering);
    setupPeakFiltering(config);

    if (isLogProgress()) {
      config.setLog(ImageJPluginLoggerHelper.getLogger(getClass()));
    }

    config.setBackgroundFitting(settings.fitBackground);

    return new Gaussian2DFitter(config);
  }

  /**
   * Sets up peak filtering.
   *
   * @param config the configuration
   */
  protected void setupPeakFiltering(FitConfiguration config) {
    final double mk = getSmooth() * 2 + 1;
    final double halfMk = 0.5 * mk;
    config.setCoordinateShift(halfMk);
    config.setSignalStrength(0);
    config.setMinWidthFactor(0.5);
    config.setMaxWidthFactor(3);
    if (settings.logProgress) {
      config.setLog(ImageJPluginLoggerHelper.getLogger(getClass()));
    }
  }

  /**
   * Fits a single 2D Gaussian to the image within the image processor. The fit is initialised at
   * the highest pixel value and then optimised.
   *
   * <p>The angle parameter is only set if using elliptical Gaussian fitting.
   *
   * <p>Note: The fitted coordinates are offset by 0.5, i.e. using the middle of the pixel. This
   * equates to input data 0,0 representing 0.5,0.5.
   *
   * @param ip the ip
   * @return Array containing the fitted curve data: Background, Amplitude, PosX, PosY, StdDevX,
   *         StdDevY, Angle. Null if no fit is possible.
   */
  public double[] fit(ImageProcessor ip) {
    // Note: This is a library function used in e.g. GDSC ImageJ plugins: FindFoci

    final float[] data = (float[]) ip.toFloat(0, null).getPixels();

    final double[] result = fit(data, ip.getWidth(), ip.getHeight());
    if (result != null) {
      result[Gaussian2DFunction.X_POSITION] += 0.5;
      result[Gaussian2DFunction.Y_POSITION] += 0.5;
    }
    return result;
  }

  /**
   * Fits a single 2D Gaussian to the data. The fit is initialised at the highest value and then
   * optimised.
   *
   * <p>Data must be arranged in yx block order, i.e. height rows of width.
   *
   * <p>The angle parameter is only set if using elliptical Gaussian fitting.
   *
   * <p>Note: The returned fit coordinates should be offset by 0.5 if the input data represents
   * pixels
   *
   * @param data the data
   * @param width the width
   * @param height the height
   * @return Array containing the fitted curve data: Background, Amplitude, PosX, PosY, StdDevX,
   *         StdDevY, Angle. Null if no fit is possible.
   */
  @Nullable
  public double[] fit(float[] data, int width, int height) {
    // Note: This is a library function used in e.g. GDSC ImageJ plugins: FindFoci

    if (data == null || data.length != width * height) {
      return null;
    }

    // Get the limits
    float max = Float.MIN_VALUE;
    int maxIndex = -1;
    for (int i = data.length; i-- > 0;) {
      final float f = data[i];
      if (max < f) {
        max = f;
        maxIndex = i;
      }
    }

    if (maxIndex < 0) {
      return null;
    }

    final Gaussian2DFitter gf = createGaussianFitter(false);
    final FitResult fitResult =
        gf.fit(SimpleArrayUtils.toDouble(data), width, height, new int[] {maxIndex});
    if (fitResult.getStatus() == FitStatus.OK) {
      chiSquared = fitResult.getError();
      final double[] params = fitResult.getParameters();

      // Check bounds
      final double x = params[Gaussian2DFunction.X_POSITION];
      final double y = params[Gaussian2DFunction.Y_POSITION];
      if (x < 0 || x >= width || y < 0 || y >= height) {
        return null;
      }

      // Re-arrange order for backwards compatibility with old code.
      final double background = params[Gaussian2DFunction.BACKGROUND];
      final double intensity = params[Gaussian2DFunction.SIGNAL];
      final double sx = params[Gaussian2DFunction.X_SD];
      final double sy = params[Gaussian2DFunction.Y_SD];
      final double angle = params[Gaussian2DFunction.ANGLE];
      final double amplitude = Gaussian2DPeakResultHelper.getAmplitude(intensity, sx, sy);

      return new double[] {background, amplitude, x, y, sx, sy, angle};
    }

    return null;
  }

  @Override
  public void setNPasses(int passes) {
    // Nothing to do
  }

  /**
   * Sets the smooth.
   *
   * @param smooth the new smooth
   */
  public void setSmooth(double smooth) {
    this.settings.smooth = smooth;
  }

  /**
   * Gets the smooth.
   *
   * @return the smooth
   */
  public double getSmooth() {
    return settings.smooth;
  }

  /**
   * Sets the box size.
   *
   * @param boxSize the new box size
   */
  public void setBoxSize(int boxSize) {
    this.settings.boxSize = boxSize;
  }

  /**
   * Gets the box size.
   *
   * @return the box size
   */
  public int getBoxSize() {
    return settings.boxSize;
  }

  /**
   * Sets the background.
   *
   * @param background the new background
   */
  public void setBackground(float background) {
    this.settings.background = background;
  }

  /**
   * Gets the background.
   *
   * @return the background
   */
  public float getBackground() {
    return settings.background;
  }

  /**
   * Sets the peak height.
   *
   * @param peakHeight the new peak height
   */
  public void setPeakHeight(float peakHeight) {
    this.settings.peakHeight = peakHeight;
  }

  /**
   * Gets the peak height.
   *
   * @return the peak height
   */
  public float getPeakHeight() {
    return settings.peakHeight;
  }

  /**
   * Sets the fraction above background.
   *
   * @param fractionAboveBackground the new fraction above background
   */
  public void setFractionAboveBackground(float fractionAboveBackground) {
    this.settings.fractionAboveBackground = fractionAboveBackground;
  }

  /**
   * Gets the fraction above background.
   *
   * @return the fraction above background
   */
  public float getFractionAboveBackground() {
    return settings.fractionAboveBackground;
  }

  /**
   * Sets the peak width.
   *
   * @param peakWidth the new peak width
   */
  public void setPeakWidth(float peakWidth) {
    this.settings.peakWidth = peakWidth;
  }

  /**
   * Gets the peak width.
   *
   * @return the peak width
   */
  public float getPeakWidth() {
    return settings.peakWidth;
  }

  /**
   * Gets the top N.
   *
   * @return the top N
   */
  public int getTopN() {
    return settings.topN;
  }

  /**
   * Sets the top N.
   *
   * @param topN the new top N
   */
  public void setTopN(int topN) {
    this.settings.topN = topN;
  }

  /**
   * Sets the block find algorithm.
   *
   * @param blockFindAlgorithm the new block find algorithm
   */
  public void setBlockFindAlgorithm(boolean blockFindAlgorithm) {
    this.settings.blockFindAlgorithm = blockFindAlgorithm;
  }

  /**
   * Checks if is block find algorithm.
   *
   * @return true, if is block find algorithm
   */
  public boolean isBlockFindAlgorithm() {
    return settings.blockFindAlgorithm;
  }

  /**
   * Sets the neighbour check.
   *
   * @param neighbourCheck the new neighbour check
   */
  public void setNeighbourCheck(boolean neighbourCheck) {
    this.settings.neighbourCheck = neighbourCheck;
  }

  /**
   * Checks if is neighbour check.
   *
   * @return true, if is neighbour check
   */
  public boolean isNeighbourCheck() {
    return settings.neighbourCheck;
  }

  /**
   * Sets the border.
   *
   * @param border the new border
   */
  public void setBorder(int border) {
    this.settings.border = border;
  }

  /**
   * Gets the border.
   *
   * @return the border
   */
  public int getBorder() {
    return settings.border;
  }

  /**
   * Sets the fit function.
   *
   * @param fitFunction the new fit function
   */
  public void setFitFunction(int fitFunction) {
    this.settings.fitFunction = fitFunction;
  }

  /**
   * Gets the fit function.
   *
   * @return the fit function
   */
  public int getFitFunction() {
    return settings.fitFunction;
  }

  /**
   * Gets the PSF type.
   *
   * @return the PSF type
   */
  public PSFType getPsfType() {
    if (settings.fitFunction >= 0 && settings.fitFunction < getPsfTypeValues().length) {
      return getPsfTypeValues()[settings.fitFunction];
    }
    return PSFType.ONE_AXIS_GAUSSIAN_2D;
  }

  /**
   * Sets the fit background.
   *
   * @param fitBackground the new fit background
   */
  public void setFitBackground(boolean fitBackground) {
    this.settings.fitBackground = fitBackground;
  }

  /**
   * Checks if is fit background.
   *
   * @return true, if is fit background
   */
  public boolean isFitBackground() {
    return settings.fitBackground;
  }

  /**
   * Sets the log progress.
   *
   * @param logProgress the new log progress
   */
  public void setLogProgress(boolean logProgress) {
    this.settings.logProgress = logProgress;
  }

  /**
   * Checks if is log progress.
   *
   * @return true, if is log progress
   */
  public boolean isLogProgress() {
    return settings.logProgress;
  }

  /**
   * Sets the max iterations.
   *
   * @param maxIterations the new max iterations
   */
  public void setMaxIterations(int maxIterations) {
    this.settings.maxIterations = maxIterations;
  }

  /**
   * Gets the max iterations.
   *
   * @return the max iterations
   */
  public int getMaxIterations() {
    return settings.maxIterations;
  }

  /**
   * Checks if is elliptical fitting.
   *
   * @return true, if is elliptical fitting
   */
  public boolean isEllipticalFitting() {
    return settings.fitFunction == 3;
  }

  /**
   * Sets the initial peak standard deviation. This is estimated from the data if zero
   *
   * @param initialPeakStdDev the new initial peak standard deviation
   */
  public void setInitialPeakStdDev(double initialPeakStdDev) {
    this.settings.initialPeakStdDev = initialPeakStdDev;
  }

  /**
   * Gets the initial peak standard deviation.
   *
   * @return the initial peak standard deviation
   */
  public double getInitialPeakStdDev() {
    return settings.initialPeakStdDev;
  }
}
