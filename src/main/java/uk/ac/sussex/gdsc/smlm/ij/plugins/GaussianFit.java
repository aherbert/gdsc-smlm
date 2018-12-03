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
package uk.ac.sussex.gdsc.smlm.ij.plugins;

import uk.ac.sussex.gdsc.core.filters.FilteredNonMaximumSuppression;
import uk.ac.sussex.gdsc.core.ij.ImageJPluginLoggerHelper;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.utils.ImageExtractor;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.SortUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.FitSolver;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtosHelper;
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
import uk.ac.sussex.gdsc.smlm.ij.results.IJTablePeakResults;
import uk.ac.sussex.gdsc.smlm.ij.settings.Constants;
import uk.ac.sussex.gdsc.smlm.ij.utils.IJImageConverter;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.measure.Measurements;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;

import java.awt.AWTEvent;
import java.awt.Color;
import java.awt.Rectangle;
import java.util.Arrays;
import java.util.EnumSet;

/**
 * Fits the selected rectangular ROI using a 2D Gaussian.
 */
public class GaussianFit implements ExtendedPlugInFilter, DialogListener {
  private static final String TITLE = "Gaussian Fit";
  private double smooth = Prefs.get(Constants.smooth, 0);
  private int boxSize = (int) Prefs.get(Constants.boxSize, 1);
  private float background = (float) Prefs.get(Constants.background, 0);
  private float peakHeight = (float) Prefs.get(Constants.peakHeight, 0);
  private float fractionAboveBackground = (float) Prefs.get(Constants.fractionAboveBackground, 0);
  private float peakWidth = (float) Prefs.get(Constants.peakWidth, 0);
  private int topN = (int) Prefs.get(Constants.topN, 0);
  private boolean blockFindAlgorithm = Prefs.get(Constants.blockFindAlgorithm, true);
  private boolean neighbourCheck = Prefs.get(Constants.neighbourCheck, false);
  private int border = (int) Prefs.get(Constants.border, 0);
  private int fitFunction = (int) Prefs.get(Constants.fitFunction, 0);
  private boolean fitBackground = Prefs.get(Constants.fitBackground, true);
  private boolean logProgress = Prefs.get(Constants.logProgress, false);
  private int maxIterations = (int) Prefs.get(Constants.maxIterations, 20);
  private double relativeThreshold = Prefs.get(Constants.relativeThreshold, 1e-5);
  private double absoluteThreshold = Prefs.get(Constants.absoluteThreshold, 1e-10);
  private boolean singleFit = Prefs.get(Constants.singleFit, false);
  private int singleRegionSize = (int) Prefs.get(Constants.singleRegionSize, 10);
  private double initialPeakStdDev = Prefs.get(Constants.initialPeakStdDev0, 0);
  private boolean showDeviations = Prefs.get(Constants.showDeviations, false);
  private boolean filterResults = Prefs.get(Constants.filterResults, false);
  private boolean showFit = Prefs.get(Constants.showFit, false);

  private final int flags = DOES_16 | DOES_8G | DOES_32 | FINAL_PROCESSING | SNAPSHOT;
  private ImagePlus imp;

  private int[] maxIndices;
  private FitResult fitResult;
  private double chiSquared;

  private IJTablePeakResults results;

  /** {@inheritDoc} */
  @Override
  public int setup(String arg, ImagePlus imp) {
    SMLMUsageTracker.recordPlugin(this.getClass(), arg);

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
    if (arg.equals("final")) {
      runFinal(imp.getProcessor());
      return DONE;
    }

    return flags;
  }

  private static PSFType[] _PSFTypeValues;

  /**
   * Gets the PSF type values.
   *
   * @return the PSF type values
   */
  public static PSFType[] getPSFTypeValues() {
    if (_PSFTypeValues == null) {
      initPSFType();
    }
    return _PSFTypeValues;
  }

  private static String[] _PSFTypeNames;

  /**
   * Gets the PSF type names.
   *
   * @return the PSF type names
   */
  public static String[] getPSFTypeNames() {
    if (_PSFTypeNames == null) {
      initPSFType();
    }
    return _PSFTypeNames;
  }

  private static void initPSFType() {
    //@formatter:off
    final EnumSet<PSFType> d = EnumSet.of(
        PSFType.ONE_AXIS_GAUSSIAN_2D,
        PSFType.TWO_AXIS_GAUSSIAN_2D,
        PSFType.TWO_AXIS_AND_THETA_GAUSSIAN_2D);
    //@formatter:on
    _PSFTypeValues = d.toArray(new PSFType[d.size()]);
    _PSFTypeNames = new String[_PSFTypeValues.length];
    for (int i = 0; i < _PSFTypeValues.length; i++) {
      _PSFTypeNames[i] = PSFProtosHelper.getName(_PSFTypeValues[i]);
    }
  }

  /** {@inheritDoc} */
  @Override
  public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
    final double[] limits = getLimits(imp.getProcessor());
    final double minValue = limits[0];
    final double maxValue = limits[1];

    if (background > maxValue) {
      background = (int) maxValue;
    }
    if (background < minValue) {
      background = (int) minValue;
    }

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    gd.addMessage("Fit 2D Gaussian to identified maxima");

    gd.addMessage("--- Image smoothing ---\n" + "- Within a 2n+1 box\n");
    gd.addSlider("Smoothing", 0, 4.5, smooth);

    gd.addMessage("--- Maxima identification ---\n" + "- Within a 2n+1 box\n");
    gd.addSlider("Box_size", 1, 15, boxSize);
    gd.addSlider("Background", minValue, maxValue, background);
    gd.addSlider("Min_height", 0, maxValue, peakHeight);
    gd.addSlider("Fraction_above_background", 0, 1.01, fractionAboveBackground);
    gd.addSlider("Min_width", 0, 20, peakWidth);
    gd.addSlider("Top_N", 0, 20, topN);
    gd.addCheckbox("Block_find_algorithm", blockFindAlgorithm);
    gd.addCheckbox("Neighbour_check", neighbourCheck);
    gd.addSlider("Border", 0, 15, border);

    gd.addMessage("--- Gaussian fitting ---");
    gd.addChoice("PSF", getPSFTypeNames(), PSFProtosHelper.getName(getPSFType()));
    gd.addCheckbox("Fit_background", fitBackground);
    gd.addNumericField("Max_iterations", maxIterations, 0);
    gd.addNumericField("Relative_threshold", relativeThreshold, -3);
    gd.addNumericField("Absolute_threshold", absoluteThreshold, -3);
    gd.addCheckbox("Single_fit", singleFit);
    gd.addNumericField("Single_region_size", singleRegionSize, 0);
    gd.addNumericField("Initial_StdDev", initialPeakStdDev, 3);
    gd.addCheckbox("Log_progress", logProgress);
    gd.addCheckbox("Show_deviations", showDeviations);
    gd.addCheckbox("Filter_results", filterResults);
    gd.addCheckbox("Show_fit", showFit);

    gd.addPreviewCheckbox(pfr);
    gd.addDialogListener(this);

    // // Initialise preview
    // gd.getPreviewCheckbox().setState(true);
    // setProperties();
    // this.run(imp.getProcessor());

    gd.showDialog();

    if (gd.wasCanceled() || !dialogItemChanged(gd, null)) {
      // imp.getProcessor().reset();
      imp.setOverlay(null);
      return DONE;
    }

    return flags;
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
    int i = 0;
    while (i < data.length) {
      count += data[i];
      if (count > limit) {
        break;
      }
      i++;
    }

    limits[1] = i;

    return limits;
  }

  /** {@inheritDoc} */
  @Override
  public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
    smooth = gd.getNextNumber();
    boxSize = (int) gd.getNextNumber();
    background = (float) gd.getNextNumber();
    peakHeight = (float) gd.getNextNumber();
    fractionAboveBackground = (float) gd.getNextNumber();
    peakWidth = (float) gd.getNextNumber();
    topN = (int) gd.getNextNumber();
    blockFindAlgorithm = gd.getNextBoolean();
    neighbourCheck = gd.getNextBoolean();
    border = (int) gd.getNextNumber();

    fitFunction = gd.getNextChoiceIndex();
    fitBackground = gd.getNextBoolean();
    maxIterations = (int) gd.getNextNumber();
    relativeThreshold = gd.getNextNumber();
    absoluteThreshold = gd.getNextNumber();
    singleFit = gd.getNextBoolean();
    singleRegionSize = (int) gd.getNextNumber();
    initialPeakStdDev = gd.getNextNumber();
    logProgress = gd.getNextBoolean();
    showDeviations = gd.getNextBoolean();
    filterResults = gd.getNextBoolean();
    showFit = gd.getNextBoolean();

    if (gd.invalidNumber()) {
      return false;
    }

    // Check arguments
    try {
      Parameters.isPositive("Smoothing", smooth);
      Parameters.isAboveZero("Box size", boxSize);
      Parameters.isPositive("Peak height", peakHeight);
      Parameters.isPositive("Fraction above background", fractionAboveBackground);
      Parameters.isPositive("Peak width", peakWidth);
      Parameters.isPositive("Top N", topN);
      Parameters.isPositive("Border", border);
      Parameters.isAboveZero("Relative threshold", relativeThreshold);
      Parameters.isAboveZero("Absolute threshold", absoluteThreshold);
      Parameters.isAboveZero("Max iterations", maxIterations);
      Parameters.isAboveZero("Single region size", singleRegionSize);
      Parameters.isPositive("Initial peak StdDev", initialPeakStdDev);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    setProperties();

    if (!gd.getPreviewCheckbox().getState()) {
      imp.setOverlay(null);
    }

    return true;
  }

  private void setProperties() {
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

  /** {@inheritDoc} */
  @Override
  public void run(ImageProcessor ip) {
    final Rectangle bounds = ip.getRoi();

    // Crop to the ROI
    final float[] data = IJImageConverter.getData(ip);

    final int width = bounds.width;
    final int height = bounds.height;

    if (getSmooth() > 0) {
      // No need for a copy since we are using a snapshot buffer
      final BlockMeanFilter filter = new BlockMeanFilter();
      filter.stripedBlockFilter(data, width, height, (float) getSmooth());
    }

    maxIndices = getMaxima(data, width, height);

    if (topN > 0 && maxIndices.length > topN) {
      SortUtils.sortIndices(maxIndices, data, true);
      maxIndices = Arrays.copyOf(maxIndices, topN);
    }

    // Show an overlay of the indices
    if (maxIndices.length > 0) {
      final int nMaxima = maxIndices.length;
      final float[] xpoints = new float[nMaxima];
      final float[] ypoints = new float[nMaxima];
      int n = 0;
      for (final int index : maxIndices) {
        xpoints[n] = 0.5f + bounds.x + index % width;
        ypoints[n] = 0.5f + bounds.y + index / width;
        n++;
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
   * @param nMaxima the number of maxima
   * @param xpoints the xpoints
   * @param ypoints the ypoints
   */
  private void setOverlay(int nMaxima, float[] xpoints, float[] ypoints) {
    final PointRoi roi = new PointRoi(xpoints, ypoints, nMaxima);

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
    final float[] data = IJImageConverter.getData(ip);

    final int width = bounds.width;
    final int height = bounds.height;

    // Sort the maxima
    float[] smoothData = data;
    if (getSmooth() > 0) {
      // Smoothing destructively modifies the data so create a copy
      smoothData = Arrays.copyOf(data, width * height);
      final BlockMeanFilter filter = new BlockMeanFilter();
      // filter.blockAverage(smoothData, width, height, smooth);
      if (smooth <= border) {
        filter.stripedBlockFilterInternal(smoothData, width, height, (float) smooth);
      } else {
        filter.stripedBlockFilter(smoothData, width, height, (float) smooth);
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

    results =
        new IJTablePeakResults(showDeviations, imp.getTitle() + " [" + imp.getCurrentSlice() + "]");
    final CalibrationWriter cw = new CalibrationWriter();
    cw.setIntensityUnit(IntensityUnit.COUNT);
    cw.setDistanceUnit(DistanceUnit.PIXEL);
    cw.setAngleUnit(AngleUnit.RADIAN);
    results.setCalibration(cw.getCalibration());
    results.setPSF(PSFProtosHelper.getDefaultPSF(getPSFType()));
    results.setShowFittingData(true);
    results.setAngleUnit(AngleUnit.DEGREE);
    results.begin();

    // Perform the Gaussian fit
    long ellapsed = 0;

    if (!singleFit) {
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
        int c = 0;
        int validPeaks = 0;
        validParams[c++] = params[0];

        final double[] initialParams = convertParameters(fitResult.getInitialParameters());
        final double[] paramsDev = convertParameters(fitResult.getParameterDeviations());
        final Rectangle regionBounds = new Rectangle();

        final float[] xpoints = new float[maxIndices.length];
        final float[] ypoints = new float[maxIndices.length];
        int nMaxima = 0;

        for (int i = 1, n = 0; i < params.length; i +=
            Gaussian2DFunction.PARAMETERS_PER_PEAK, n++) {
          final int y = maxIndices[n] / width;
          final int x = maxIndices[n] % width;

          // Check the peak is a good fit
          if (filterResults
              && config.validatePeak(n, initialParams, params, paramsDev) != FitStatus.OK) {
            continue;
          }

          if (showFit) {
            // Copy the valid parameters
            validPeaks++;
            for (int ii = i, j = 0; j < Gaussian2DFunction.PARAMETERS_PER_PEAK; ii++, j++) {
              validParams[c++] = params[ii];
            }
          }

          final double[] peakParams = extractParams(params, i);
          final double[] peakParamsDev = extractParams(paramsDev, i);

          addResult(bounds, regionBounds, peakParams, peakParamsDev, nMaxima, x, y,
              data[maxIndices[n]]);

          // Add fit result to the overlay - Coords are updated with the region offsets in addResult
          final double xf = peakParams[Gaussian2DFunction.X_POSITION];
          final double yf = peakParams[Gaussian2DFunction.Y_POSITION];
          xpoints[nMaxima] = (float) xf;
          ypoints[nMaxima] = (float) yf;
          nMaxima++;
        }

        setOverlay(nMaxima, xpoints, ypoints);

        // Draw the fit
        if (showFit && validPeaks != 0) {
          final double[] pixels = new double[data.length];
          final EllipticalGaussian2DFunction f =
              new EllipticalGaussian2DFunction(validPeaks, width, height);
          invertParameters(validParams);
          f.initialise(validParams);
          for (int x = 0; x < pixels.length; x++) {
            pixels[x] = f.eval(x);
          }
          final FloatProcessor fp = new FloatProcessor(width, height, pixels);
          // Insert into a full size image
          final FloatProcessor fp2 = new FloatProcessor(ip.getWidth(), ip.getHeight());
          fp2.insert(fp, bounds.x, bounds.y);
          ImageJUtils.display(TITLE, fp2);
        }
      } else {
        if (isLogProgress()) {
          IJ.log("Failed to fit " + TextUtils.pleural(maxIndices.length, "peak")
              + getReason(fitResult));
        }
        imp.setOverlay(null);
      }
    } else {
      if (isLogProgress()) {
        IJ.log("Individual fit");
      }

      int nMaxima = 0;
      final float[] xpoints = new float[maxIndices.length];
      final float[] ypoints = new float[maxIndices.length];

      // Extract each peak and fit individually
      final ImageExtractor ie = ImageExtractor.wrap(data, width, height);
      float[] region = null;
      final Gaussian2DFitter gf = createGaussianFitter(filterResults);

      for (int n = 0; n < maxIndices.length; n++) {
        final int y = maxIndices[n] / width;
        final int x = maxIndices[n] % width;

        final long time = System.nanoTime();
        final Rectangle regionBounds = ie.getBoxRegionBounds(x, y, singleRegionSize);
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
          double[] peakParamsDev = null;
          if (showDeviations) {
            peakParamsDev = convertParameters(fitResult.getParameterDeviations());
          }

          addResult(bounds, regionBounds, peakParams, peakParamsDev, n, x, y, data[maxIndices[n]]);

          // Add fit result to the overlay - Coords are updated with the region offsets in addResult
          final double xf = peakParams[Gaussian2DFunction.X_POSITION];
          final double yf = peakParams[Gaussian2DFunction.Y_POSITION];
          xpoints[nMaxima] = (float) xf;
          ypoints[nMaxima] = (float) yf;
          nMaxima++;
        } else if (isLogProgress()) {
          IJ.log("Failed to fit peak " + (n + 1) + getReason(fitResult));
        }
      }

      // Update the overlay
      if (nMaxima > 0) {
        setOverlay(nMaxima, xpoints, ypoints);
      } else {
        imp.setOverlay(null);
      }
    }

    results.end();

    if (isLogProgress()) {
      IJ.log("Time = " + (ellapsed / 1000000.0) + "ms");
    }
  }

  /**
   * Gets the reason.
   *
   * @param fitResult the fit result
   * @return the reason
   */
  private static String getReason(FitResult fitResult) {
    if (fitResult == null || fitResult.getStatus() == null) {
      return "";
    }
    final FitStatus status = fitResult.getStatus();
    return status.toString().toLowerCase().replace("_", " ");
  }

  private static double[] extractParams(double[] params, int i) {
    if (params == null) {
      return null;
    }

    // 0 is the background. Then the peaks are packed.
    final double[] p = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
    p[0] = params[0];
    System.arraycopy(params, i, p, 1, Gaussian2DFunction.PARAMETERS_PER_PEAK);
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
   * Find the indices of the maxima using the currently configured parameters <p> Data must be
   * arranged in yx block order, i.e. height rows of width.
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

    int[] maxIndices;
    if (isBlockFindAlgorithm()) {
      if (getBorder() > 0) {
        maxIndices = nms.blockFindInternal(data, width, height, getBoxSize(), getBorder());
      } else {
        maxIndices = nms.blockFind(data, width, height, getBoxSize());
      }
    } else if (getBorder() > 0) {
      maxIndices = nms.maxFindInternal(data, width, height, getBoxSize(), getBorder());
    } else {
      maxIndices = nms.maxFind(data, width, height, getBoxSize());
    }
    return maxIndices;
  }

  /**
   * Fits a 2D Gaussian to the given data. Fits all the specified peaks. <p> Data must be arranged
   * in yx block order, i.e. height rows of width. <p> Note: The fit coordinates should be offset by
   * 0.5 if the input data represents pixels
   *
   * @param data the data
   * @param width the width
   * @param height the height
   * @param maxIndices Indices of the data to fit
   * @return Array containing the fitted curve data: The first value is the Background. The
   *         remaining values are Amplitude, PosX, PosY, StdDevX, StdDevY for each fitted peak. If
   *         elliptical fitting is performed the values are Amplitude, Angle, PosX, PosY, StdDevX,
   *         StdDevY for each fitted peak <p> Null if no fit is possible.
   */
  public double[] fitMultiple(float[] data, int width, int height, int[] maxIndices) {
    return fitMultiple(data, width, height, maxIndices, null);
  }

  /**
   * Fits a 2D Gaussian to the given data. Fits all the specified peaks. <p> Data must be arranged
   * in yx block order, i.e. height rows of width. <p> Note: The fit coordinates should be offset by
   * 0.5 if the input data represents pixels
   *
   * @param data the data
   * @param width the width
   * @param height the height
   * @param maxIndices Indices of the data to fit
   * @param estimatedHeights Estimated heights for the peaks (input from smoothed data)
   * @return Array containing the fitted curve data: The first value is the Background. The
   *         remaining values are Amplitude, PosX, PosY, StdDevX, StdDevY for each fitted peak. If
   *         elliptical fitting is performed the values are Amplitude, Angle, PosX, PosY, StdDevX,
   *         StdDevY for each fitted peak <p> Null if no fit is possible.
   */
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

  private double[] convertParameters(double[] params) {
    if (params == null) {
      return null;
    }

    // Convert coordinates with 0.5 pixel offset
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

  private double[] invertParameters(double[] params) {
    if (params == null) {
      return null;
    }

    // Convert coordinates with 0.5 pixel offset
    // Convert radians to degrees (if elliptical fitting)
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
   * Fits a 2D Gaussian to the given data. Fits all the specified peaks. <p> Data must be arranged
   * in yx block order, i.e. height rows of width.
   *
   * @param gf the gf
   * @param data the data
   * @param width the width
   * @param height the height
   * @param index Index of the data to fit
   * @param estimatedHeight Estimated height for the peak (input from smoothed data)
   * @return Array containing the fitted curve data: The first value is the Background. The
   *         remaining values are Amplitude, PosX, PosY, StdDevX, StdDevY for each fitted peak. <p>
   *         Null if no fit is possible.
   */
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
    config.setPSF(PSFProtosHelper.getDefaultPSF(getPSFType()));
    config.setMaxIterations(getMaxIterations());
    config.setRelativeThreshold(relativeThreshold);
    config.setAbsoluteThreshold(absoluteThreshold);
    config.setInitialPeakStdDev(getInitialPeakStdDev());
    config.setComputeDeviations(showDeviations);

    // Set-up peak filtering only for single fitting
    config.setDisableSimpleFilter(!simpleFiltering);
    setupPeakFiltering(config);

    if (isLogProgress()) {
      config.setLog(ImageJPluginLoggerHelper.getLogger(getClass()));
    }

    config.setBackgroundFitting(fitBackground);

    return new Gaussian2DFitter(config);
  }

  /**
   * Sets up peak filtering.
   *
   * @param config the configuration
   */
  protected void setupPeakFiltering(FitConfiguration config) {
    final double Mk = getSmooth() * 2 + 1;
    final double halfMk = 0.5f * Mk;
    config.setCoordinateShift(halfMk);
    config.setSignalStrength(0);
    config.setMinWidthFactor(0.5);
    config.setWidthFactor(3);
    if (logProgress) {
      config.setLog(ImageJPluginLoggerHelper.getLogger(getClass()));
    }
  }

  /**
   * Fits a single 2D Gaussian to the image within the image processor. The fit is initialised at
   * the highest pixel value and then optimised. <p> The angle parameter is only set if using
   * elliptical Gaussian fitting. <p> Note: The fitted coordinates are offset by 0.5, i.e. using the
   * middle of the pixel. This equates to input data 0,0 representing 0.5,0.5.
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
   * optimised. <p> Data must be arranged in yx block order, i.e. height rows of width. <p> The
   * angle parameter is only set if using elliptical Gaussian fitting. <p> Note: The returned fit
   * coordinates should be offset by 0.5 if the input data represents pixels
   *
   * @param data the data
   * @param width the width
   * @param height the height
   * @return Array containing the fitted curve data: Background, Amplitude, PosX, PosY, StdDevX,
   *         StdDevY, Angle. Null if no fit is possible.
   */
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
  public void setNPasses(int nPasses) {
    // Nothing to do
  }

  /**
   * @param smooth the smooth to set
   */
  public void setSmooth(double smooth) {
    this.smooth = smooth;
  }

  /**
   * @return the smooth.
   */
  public double getSmooth() {
    return smooth;
  }

  /**
   * @param boxSize the boxSize to set
   */
  public void setBoxSize(int boxSize) {
    this.boxSize = boxSize;
  }

  /**
   * @return the boxSize.
   */
  public int getBoxSize() {
    return boxSize;
  }

  /**
   * @param background the background to set
   */
  public void setBackground(float background) {
    this.background = background;
  }

  /**
   * @return the background.
   */
  public float getBackground() {
    return background;
  }

  /**
   * @param peakHeight the peakHeight to set
   */
  public void setPeakHeight(float peakHeight) {
    this.peakHeight = peakHeight;
  }

  /**
   * @return the peakHeight.
   */
  public float getPeakHeight() {
    return peakHeight;
  }

  /**
   * @param fractionAboveBackground the fractionAboveBackground to set
   */
  public void setFractionAboveBackground(float fractionAboveBackground) {
    this.fractionAboveBackground = fractionAboveBackground;
  }

  /**
   * @return the fractionAboveBackground.
   */
  public float getFractionAboveBackground() {
    return fractionAboveBackground;
  }

  /**
   * @param peakWidth the peakWidth to set
   */
  public void setPeakWidth(float peakWidth) {
    this.peakWidth = peakWidth;
  }

  /**
   * @return the peakWidth.
   */
  public float getPeakWidth() {
    return peakWidth;
  }

  /**
   * @return the topN.
   */
  public int getTopN() {
    return topN;
  }

  /**
   * @param topN the topN to set
   */
  public void setTopN(int topN) {
    this.topN = topN;
  }

  /**
   * @param blockFindAlgorithm the blockFindAlgorithm to set
   */
  public void setBlockFindAlgorithm(boolean blockFindAlgorithm) {
    this.blockFindAlgorithm = blockFindAlgorithm;
  }

  /**
   * @return the blockFindAlgorithm.
   */
  public boolean isBlockFindAlgorithm() {
    return blockFindAlgorithm;
  }

  /**
   * @param neighbourCheck the neighbourCheck to set
   */
  public void setNeighbourCheck(boolean neighbourCheck) {
    this.neighbourCheck = neighbourCheck;
  }

  /**
   * @return the neighbourCheck.
   */
  public boolean isNeighbourCheck() {
    return neighbourCheck;
  }

  /**
   * @param border the border to set
   */
  public void setBorder(int border) {
    this.border = border;
  }

  /**
   * @return the border.
   */
  public int getBorder() {
    return border;
  }

  /**
   * @param fitFunction the fitFunction to set
   */
  public void setFitFunction(int fitFunction) {
    this.fitFunction = fitFunction;
  }

  /**
   * @return the fitFunction.
   */
  public int getFitFunction() {
    return fitFunction;
  }

  /**
   * Gets the PSF type.
   *
   * @return the PSF type
   */
  public PSFType getPSFType() {
    if (fitFunction >= 0 && fitFunction < getPSFTypeValues().length) {
      return getPSFTypeValues()[fitFunction];
    }
    return PSFType.ONE_AXIS_GAUSSIAN_2D;
  }

  /**
   * @param fitBackground the fitBackground to set
   */
  public void setFitBackground(boolean fitBackground) {
    this.fitBackground = fitBackground;
  }

  /**
   * @return the fitBackground.
   */
  public boolean isFitBackground() {
    return fitBackground;
  }

  /**
   * @param logProgress the logProgress to set
   */
  public void setLogProgress(boolean logProgress) {
    this.logProgress = logProgress;
  }

  /**
   * @return the logProgress.
   */
  public boolean isLogProgress() {
    return logProgress;
  }

  /**
   * @param maxIterations the maxIterations to set
   */
  public void setMaxIterations(int maxIterations) {
    this.maxIterations = maxIterations;
  }

  /**
   * @return the maxIterations.
   */
  public int getMaxIterations() {
    return maxIterations;
  }

  /**
   * @return True if fitting an elliptical Gaussian.
   */
  public boolean isEllipticalFitting() {
    return fitFunction == 3;
  }

  /**
   * @param initialPeakStdDev the initial peak standard deviation. This is estimated from the data
   *        if zero
   */
  public void setInitialPeakStdDev(double initialPeakStdDev) {
    this.initialPeakStdDev = initialPeakStdDev;
  }

  /**
   * @return the the initial peak standard deviation.
   */
  public double getInitialPeakStdDev() {
    return initialPeakStdDev;
  }
}
