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

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.Roi;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import ij.process.Blitter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Point;
import java.awt.Rectangle;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicReference;
import java.util.regex.Pattern;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.LoessInterpolator;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import uk.ac.sussex.gdsc.core.annotation.Nullable;
import uk.ac.sussex.gdsc.core.ij.AlignImagesFft;
import uk.ac.sussex.gdsc.core.ij.AlignImagesFft.SubPixelMethod;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.SimpleImageJTrackProgress;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.ij.process.Fht;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.logging.TrackProgress;
import uk.ac.sussex.gdsc.core.utils.FileUtils;
import uk.ac.sussex.gdsc.core.utils.ImageWindow.WindowMethod;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.concurrent.ConcurrencyUtils;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageMode;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsImageType;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.results.ImageJImagePeakResults;
import uk.ac.sussex.gdsc.smlm.ij.results.ImagePeakResultsFactory;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.count.Counter;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.StandardResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyrResultProcedure;

/**
 * Calculates drift in localisation results. Can use the feducial markers within ROI added to the
 * ROI manager or by aligning N consecutive frames with the overall image.
 */
public class DriftCalculator implements PlugIn {
  private static final String TITLE = "Drift Calculator";

  private static AtomicReference<PlotWindow> plotx = new AtomicReference<>();
  private static AtomicReference<PlotWindow> ploty = new AtomicReference<>();

  private int interpolationStart;
  private int interpolationEnd;
  private double[] calculatedTimepoints;
  private double[] lastdx;
  private double[] lastdy;

  private final TrackProgress tracker = SimpleImageJTrackProgress.getInstance();

  // Used to multi-thread the image alignment
  private ExecutorService executor;

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    static final String SUB_IMAGE_ALIGNMENT = "Localisation Sub-Images";
    static final String DRIFT_FILE = "Drift File";
    static final String STACK_ALIGNMENT = "Reference Stack Alignment";
    static final String MARKED_ROIS = "Marked ROIs";
    static final String[] UPDATE_METHODS =
        new String[] {"None", "Update", "New dataset", "New truncated dataset"};
    static final String[] SIZES = new String[] {"128", "256", "512", "1024", "2048"};

    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    String driftFilename;
    String method;
    int updateMethod;

    String inputOption;
    int maxIterations;
    double relativeError;
    double smoothing;
    boolean limitSmoothing;
    int minSmoothingPoints;
    int maxSmoothingPoints;
    int iterations;
    boolean plotDrift;
    boolean saveDrift;

    // Parameters to control the image alignment algorithm
    int frames;
    int minimimLocalisations;
    String reconstructionSize;
    String stackTitle;
    int startFrame;
    int frameSpacing;
    int interpolationMethod;
    SubPixelMethod subPixelMethod;

    Settings() {
      // Set defaults
      driftFilename = "";
      method = "";
      inputOption = "";
      maxIterations = 50;
      relativeError = 0.01;
      smoothing = 0.25;
      limitSmoothing = true;
      minSmoothingPoints = 10;
      maxSmoothingPoints = 50;
      iterations = 1;
      plotDrift = true;
      frames = 2000;
      minimimLocalisations = 50;
      reconstructionSize = SIZES[1];
      stackTitle = "";
      startFrame = 1;
      frameSpacing = 1;
      interpolationMethod = ImageProcessor.BILINEAR;
      subPixelMethod = AlignImagesFft.SubPixelMethod.CUBIC;
    }

    Settings(Settings source) {
      driftFilename = source.driftFilename;
      method = source.method;
      updateMethod = source.updateMethod;
      inputOption = source.inputOption;
      maxIterations = source.maxIterations;
      relativeError = source.relativeError;
      smoothing = source.smoothing;
      limitSmoothing = source.limitSmoothing;
      minSmoothingPoints = source.minSmoothingPoints;
      maxSmoothingPoints = source.maxSmoothingPoints;
      iterations = source.iterations;
      plotDrift = source.plotDrift;
      saveDrift = source.saveDrift;
      frames = source.frames;
      minimimLocalisations = source.minimimLocalisations;
      reconstructionSize = source.reconstructionSize;
      stackTitle = source.stackTitle;
      startFrame = source.startFrame;
      frameSpacing = source.frameSpacing;
      interpolationMethod = source.interpolationMethod;
      subPixelMethod = source.subPixelMethod;
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
   * Align images to the reference initialised in the given aligner.
   */
  private static class ImageAligner implements Runnable {
    final AlignImagesFft aligner;
    final ImageProcessor[] ip;
    final int[] time;
    final Rectangle alignBounds;
    final List<double[]> alignments;
    final int from;
    final int to;
    final Ticker ticker;
    final SubPixelMethod subPixelMethod;

    ImageAligner(AlignImagesFft aligner, ImageProcessor[] ip, int[] time, Rectangle alignBounds,
        List<double[]> alignments, int from, int to, Ticker ticker, SubPixelMethod subPixelMethod) {
      this.aligner = aligner;
      this.ip = ip;
      this.time = time;
      this.alignBounds = alignBounds;
      this.alignments = alignments;
      this.from = from;
      this.to = to;
      this.ticker = ticker;
      this.subPixelMethod = subPixelMethod;
    }

    @Override
    public void run() {
      for (int i = from; i < to && i < ip.length; i++) {
        // Window method is ignored since the image processor is already an FHT image
        double[] result = aligner.align(ip[i], WindowMethod.TUKEY, alignBounds, subPixelMethod);
        // Create a result for failures
        if (result == null) {
          result = new double[] {Double.NaN, Double.NaN, time[i]};
        }
        // Store the time point with the result
        result[2] = time[i];
        alignments.add(result);
        ticker.tick();
      }
    }
  }

  /**
   * Duplicate and translate images.
   */
  private static class ImageTranslator implements Runnable {
    final ImageProcessor[] images;
    final ImageProcessor[] ip;
    final double[] dx;
    final double[] dy;
    final int from;
    final int to;
    final Ticker ticker;
    final int interpolationMethod;

    ImageTranslator(ImageProcessor[] images, ImageProcessor[] ip, double[] dx, double[] dy,
        int from, int to, Ticker ticker, int interpolationMethod) {
      this.images = images;
      this.ip = ip;
      this.dx = dx;
      this.dy = dy;
      this.from = from;
      this.to = to;
      this.ticker = ticker;
      this.interpolationMethod = interpolationMethod;
    }

    @Override
    public void run() {
      for (int i = from; i < to && i < ip.length; i++) {
        ip[i] = images[i].duplicate();
        if (dx[i] != 0 || dy[i] != 0) {
          ip[i].setInterpolationMethod(interpolationMethod);
          ip[i].translate(dx[i], dy[i]);
        }
        ticker.tick();
      }
    }
  }

  /**
   * Creates an image reconstruction from the provided localisations.
   */
  private static class ImageBuilder implements Runnable {
    final List<Localisation> localisations;
    final ImageProcessor[] images;
    final int image;
    final Rectangle bounds;
    final float scale;
    final double[] dx;
    final double[] dy;
    final Ticker ticker;

    ImageBuilder(List<Localisation> localisations, ImageProcessor[] images, int image,
        Rectangle bounds, float scale, double[] dx, double[] dy, Ticker ticker) {
      this.localisations = localisations;
      this.images = images;
      this.image = image;
      this.bounds = bounds;
      this.scale = scale;
      this.dx = dx;
      this.dy = dy;
      this.ticker = ticker;
    }

    @Override
    public void run() {
      final ImageJImagePeakResults blockImage = newImage(bounds, scale);
      for (final Localisation r : localisations) {
        blockImage.add(r.time, (float) (r.x + dx[r.time]), (float) (r.y + dy[r.time]), r.signal);
      }
      images[image] = getImage(blockImage);
      ticker.tick();
    }
  }

  /**
   * Prepare the slices in a stack for image correlation.
   */
  private static class ImageFhtInitialiser implements Runnable {
    final ImageStack stack;
    final ImageProcessor[] images;
    final AlignImagesFft aligner;
    final Fht[] fhtImages;
    final int from;
    final int to;
    final Ticker ticker;

    public ImageFhtInitialiser(ImageStack stack, ImageProcessor[] images, AlignImagesFft aligner,
        Fht[] fhtImages, int from, int to, final Ticker ticker) {
      this.stack = stack;
      this.images = images;
      this.aligner = aligner;
      this.fhtImages = fhtImages;
      this.from = from;
      this.to = to;
      this.ticker = ticker;
    }

    @Override
    public void run() {
      for (int i = from; i < to && i < images.length; i++) {
        images[i] = stack.getProcessor(i + 1);
        AlignImagesFft.applyWindowSeparable(images[i], WindowMethod.TUKEY);
        fhtImages[i] = aligner.transformTarget(images[i], WindowMethod.NONE);
        ticker.tick();
      }
    }
  }

  /**
   * Used to precalculate the localisation signal and store it with T,X,Y values with double
   * precision.
   */
  private static class Spot {
    final int time;
    final double x;
    final double y;
    final double signal;

    Spot(int time, double x, double y, double signal) {
      this.time = time;
      this.x = x;
      this.y = y;
      this.signal = signal;
    }

    static int compare(Spot r1, Spot r2) {
      // Sort in time order
      if (r1.time == r2.time) {
        // ... then signal
        return Double.compare(r2.signal, r1.signal);
      }
      return (r1.time < r2.time) ? -1 : 1;
    }
  }

  /**
   * Used to precalculate the localisation signal and store it with T,X,Y values.
   */
  private static class Localisation {
    final int time;
    final float x;
    final float y;
    final float signal;

    Localisation(int time, float x, float y, float signal) {
      this.time = time;
      this.x = x;
      this.y = y;
      this.signal = signal;
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    // Require some fit results and selected regions
    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "There are no fitting results in memory");
      return;
    }
    final Roi[] rois = getRois();
    final String[] stackTitles = createStackImageList();

    if (!showDialog(rois, stackTitles)) {
      return;
    }

    final MemoryPeakResults results =
        ResultsManager.loadInputResults(settings.inputOption, false, DistanceUnit.PIXEL, null);
    if (results == null || results.size() < 2) {
      IJ.error(TITLE, "There are not enough fitting results for drift correction");
      return;
    }
    double[][] drift = null;
    final int[] limits = findTimeLimits(results);
    try {
      if (settings.method.equals(Settings.MARKED_ROIS)) {
        drift = calculateUsingMarkers(results, limits, rois);
      } else if (settings.method.equals(Settings.STACK_ALIGNMENT)) {
        final ImageStack stack = showStackDialog(stackTitles);
        if (stack == null) {
          return;
        }
        drift = calculateUsingImageStack(stack, limits);
      } else if (settings.method.equals(Settings.DRIFT_FILE)) {
        drift = calculateUsingDriftFile(limits);
      } else {
        if (!showSubImageDialog()) {
          return;
        }
        drift =
            calculateUsingFrames(results, limits, Integer.parseInt(settings.reconstructionSize));
      }
    } finally {
      if (executor != null) {
        executor.shutdown();
      }
    }

    if (drift == null) {
      return;
    }

    ImageJUtils.log("Drift correction interpolated for frames [%d - %d] of [%d - %d] (%s%%)",
        interpolationStart, interpolationEnd, limits[0], limits[1], MathUtils.rounded(
            (100.0 * (interpolationEnd - interpolationStart + 1)) / (limits[1] - limits[0] + 1)));

    applyDriftCorrection(results, drift);
  }

  @Nullable
  private static Roi[] getRois() {
    final RoiManager rmanager = RoiManager.getInstance();
    if (rmanager == null || rmanager.getCount() == 0) {
      IJ.log("To use feducial markers for drift correction, add ROIs to the RoiManager"
          + " (select a region then press [t]).");
      return null;
    }
    return rmanager.getRoisAsArray();
  }

  private boolean showDialog(Roi[] rois, String[] stackTitles) {
    settings = Settings.load();

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("drift-calculator"));

    gd.addMessage("Correct the drift in localisation results");
    ResultsManager.addInput(gd, settings.inputOption, InputSource.MEMORY);
    final ArrayList<String> methods = new ArrayList<>(4);
    methods.add(Settings.SUB_IMAGE_ALIGNMENT);
    methods.add(Settings.DRIFT_FILE);
    if (rois != null) {
      methods.add(Settings.MARKED_ROIS);
    }
    if (stackTitles != null) {
      methods.add(Settings.STACK_ALIGNMENT);
    }
    final String[] items = methods.toArray(new String[0]);
    gd.addChoice("Method", items, settings.method);
    gd.addMessage("Stopping criteria");
    gd.addSlider("Max_iterations", 0, 100, settings.maxIterations);
    gd.addNumericField("Relative_error", settings.relativeError, 3);
    gd.addMessage("LOESS smoothing parameters");
    gd.addSlider("Smoothing", 0.001, 1, settings.smoothing);
    gd.addCheckbox("Limit_smoothing", settings.limitSmoothing);
    gd.addSlider("Min_smoothing_points", 5, 50, settings.minSmoothingPoints);
    gd.addSlider("Max_smoothing_points", 5, 50, settings.maxSmoothingPoints);
    gd.addSlider("Smoothing_iterations", 1, 10, settings.iterations);
    gd.addCheckbox("Plot_drift", settings.plotDrift);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.inputOption = ResultsManager.getInputSource(gd);
    settings.method = gd.getNextChoice();
    settings.maxIterations = (int) gd.getNextNumber();
    settings.relativeError = gd.getNextNumber();
    settings.smoothing = gd.getNextNumber();
    settings.limitSmoothing = gd.getNextBoolean();
    settings.minSmoothingPoints = (int) gd.getNextNumber();
    settings.maxSmoothingPoints = (int) gd.getNextNumber();
    settings.iterations = (int) gd.getNextNumber();
    settings.plotDrift = gd.getNextBoolean();
    settings.save();

    // Check arguments
    try {
      ParameterUtils.isPositive("Max iterations", settings.maxIterations);
      ParameterUtils.isAboveZero("Relative error", settings.relativeError);
      ParameterUtils.isPositive("Smoothing", settings.smoothing);
      if (settings.limitSmoothing) {
        ParameterUtils.isEqualOrAbove("Min smoothing points", settings.minSmoothingPoints, 3);
        ParameterUtils.isEqualOrAbove("Max smoothing points", settings.maxSmoothingPoints, 3);
        ParameterUtils.isEqualOrAbove("Max smoothing points", settings.maxSmoothingPoints,
            settings.minSmoothingPoints);
      }
      ParameterUtils.isEqualOrBelow("Smoothing", settings.smoothing, 1);
      ParameterUtils.isPositive("Smoothing iterations", settings.iterations);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return true;
  }

  private boolean showSubImageDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("drift-calculator"));

    gd.addMessage("Compute the drift using localisation sub-image alignment");
    gd.addNumericField("Frames", settings.frames, 0);
    gd.addSlider("Minimum_localisations", 10, 50, settings.minimimLocalisations);
    gd.addChoice("FFT size", Settings.SIZES, settings.reconstructionSize);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.frames = (int) gd.getNextNumber();
    settings.minimimLocalisations = (int) gd.getNextNumber();
    settings.reconstructionSize = gd.getNextChoice();

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Frames", settings.frames);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return true;
  }

  private ImageStack showStackDialog(String[] stackTitles) {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("drift-calculator"));

    gd.addMessage("Compute the drift using a reference stack alignment");

    gd.addChoice("Stack_image", stackTitles, settings.stackTitle);
    gd.addMessage("Frame = previous + spacing");
    gd.addNumericField("Start_frame", settings.startFrame, 0);
    gd.addSlider("Frame_spacing", 1, 20, settings.frameSpacing);
    final String[] methods = ImageProcessor.getInterpolationMethods();
    gd.addChoice("Interpolation_method", methods, methods[settings.interpolationMethod]);
    gd.showDialog();

    if (gd.wasCanceled()) {
      return null;
    }

    settings.stackTitle = gd.getNextChoice();
    settings.startFrame = (int) gd.getNextNumber();
    settings.frameSpacing = (int) gd.getNextNumber();
    settings.interpolationMethod = gd.getNextChoiceIndex();

    try {
      ParameterUtils.isAboveZero("Start frame", settings.startFrame);
      ParameterUtils.isAboveZero("Frame spacing", settings.frameSpacing);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return null;
    }

    final ImagePlus imp = WindowManager.getImage(settings.stackTitle);
    if (imp != null && imp.getStackSize() > 1) {
      return imp.getImageStack();
    }
    return null;
  }

  /**
   * Build a list of suitable stack images.
   *
   * @return the list of suitable stack images.
   */
  @Nullable
  private static String[] createStackImageList() {
    final int[] idList = WindowManager.getIDList();
    if (idList != null) {
      final String[] list = new String[idList.length];
      int count = 0;
      for (final int id : idList) {
        final ImagePlus imp = WindowManager.getImage(id);
        if (imp != null && imp.getStackSize() > 1) {
          list[count++] = imp.getTitle();
        }
      }
      return Arrays.copyOf(list, count);
    }
    return null;
  }

  private void applyDriftCorrection(MemoryPeakResults results, double[][] drift) {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage("Apply drift correction to in-memory results?");
    gd.addChoice("Update_method", Settings.UPDATE_METHODS, settings.updateMethod);
    // Option to save the drift unless it was loaded from file
    if (!Settings.DRIFT_FILE.equals(settings.method)) {
      gd.addCheckbox("Save_drift", settings.saveDrift);
    }
    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    settings.updateMethod = gd.getNextChoiceIndex();
    if (!Settings.DRIFT_FILE.equals(settings.method)) {
      settings.saveDrift = gd.getNextBoolean();
      saveDrift(calculatedTimepoints, lastdx, lastdy);
    }
    if (settings.updateMethod == 0) {
      return;
    }

    final double[] dx = drift[0];
    final double[] dy = drift[1];

    // Note: We can use the raw procedure on the results because we requested
    // the results were in pixels

    if (settings.updateMethod == 1) {
      // Update the results in memory
      ImageJUtils.log("Applying drift correction to the results set: " + results.getName());
      results.forEach((PeakResultProcedure) result -> {
        result.setXPosition((float) (result.getXPosition() + dx[result.getFrame()]));
        result.setYPosition((float) (result.getYPosition() + dy[result.getFrame()]));
      });
    } else {
      // Create a new set of results
      final MemoryPeakResults newResults = new MemoryPeakResults(results.size());
      newResults.copySettings(results);
      newResults.setName(results.getName() + " (Corrected)");
      MemoryPeakResults.addResults(newResults);
      final boolean truncate = settings.updateMethod == 3;
      ImageJUtils.log("Creating %sdrift corrected results set: " + newResults.getName(),
          (truncate) ? "truncated " : "");
      results.forEach((PeakResultProcedure) result -> {
        if (truncate
            && (result.getFrame() < interpolationStart || result.getFrame() > interpolationEnd)) {
          return;
        }
        result.setXPosition((float) (result.getXPosition() + dx[result.getFrame()]));
        result.setYPosition((float) (result.getYPosition() + dy[result.getFrame()]));
        newResults.add(result);
      });
    }
  }

  /**
   * Calculates drift using the feducial markers within ROI.
   *
   * <p>Adapted from the drift calculation method in QuickPALM.
   *
   * @param results the results
   * @param limits the limits
   * @param rois the rois
   * @return the drift { dx[], dy[] }
   */
  @Nullable
  private double[][] calculateUsingMarkers(MemoryPeakResults results, int[] limits, Roi[] rois) {
    final Spot[][] roiSpots = findSpots(results, rois, limits);

    // Check we have enough data
    if (roiSpots.length == 0) {
      IJ.error("No peak fit results in the selected ROIs");
      return null;
    }

    final double[] dx = new double[limits[1] + 1];
    final double[] dy = new double[dx.length];

    final double[] sum = new double[roiSpots.length];
    final double[] weights = calculateWeights(roiSpots, dx.length, sum);

    final double smoothing = updateSmoothingParameter(weights);

    lastdx = null;
    double change =
        calculateDriftUsingMarkers(roiSpots, weights, sum, dx, dy, smoothing, settings.iterations);
    if (Double.isNaN(change) || tracker.isEnded()) {
      return null;
    }
    ImageJUtils.log("Drift Calculator : Initial drift " + MathUtils.rounded(change));

    for (int i = 1; i <= settings.maxIterations; i++) {
      change = calculateDriftUsingMarkers(roiSpots, weights, sum, dx, dy, smoothing,
          settings.iterations);
      if (Double.isNaN(change)) {
        return null;
      }

      if (converged(i, change, getTotalDrift(dx, dy, weights))) {
        break;
      }
    }

    if (tracker.isEnded()) {
      return null;
    }

    interpolate(dx, dy, weights);

    plotDrift(limits, dx, dy);
    saveDrift(weights, dx, dy);

    return new double[][] {dx, dy};
  }

  /**
   * Update the smoothing parameter using the upper and lower limits for the number of points to use
   * for smoothing.
   *
   * @param data The data to be smoothed
   * @return The updated smoothing parameter
   */
  private double updateSmoothingParameter(double[] data) {
    if (!settings.limitSmoothing) {
      return settings.smoothing;
    }

    final int n = countNonZeroValues(data);

    int bandwidthInpoints = (int) (settings.smoothing * n);

    // Check the bounds for the smoothing
    final int original = bandwidthInpoints;
    if (settings.minSmoothingPoints > 0) {
      bandwidthInpoints = Math.max(bandwidthInpoints, settings.minSmoothingPoints);
    }
    if (settings.maxSmoothingPoints > 0) {
      bandwidthInpoints = Math.min(bandwidthInpoints, settings.maxSmoothingPoints);
    }

    final double newSmoothing = (double) bandwidthInpoints / n;
    if (original != bandwidthInpoints) {
      ImageJUtils.log("Updated smoothing parameter for %d data points to %s (%d smoothing points)",
          n, MathUtils.rounded(newSmoothing), bandwidthInpoints);
    }

    return newSmoothing;
  }

  /**
   * Count the number of points where the data array is not zero.
   *
   * @param data the data
   * @return the number of points
   */
  private static int countNonZeroValues(double[] data) {
    int count = 0;
    for (final double d : data) {
      if (d != 0) {
        count++;
      }
    }
    return count;
  }

  private static double getTotalDrift(double[] dx, double[] dy, double[] originalDriftTimePoints) {
    double totalDrift = 0;
    for (int t = 0; t < dx.length; t++) {
      if (originalDriftTimePoints[t] != 0) {
        totalDrift += Math.sqrt(dx[t] * dx[t] + dy[t] * dy[t]);
      }
    }
    return totalDrift;
  }

  private boolean converged(int iteration, double change, double totalDrift) {
    final double error = change / totalDrift;
    ImageJUtils.log("Iteration %d : Drift %s : Total change %s : Relative change %s", iteration,
        MathUtils.rounded(totalDrift), MathUtils.rounded(change), MathUtils.rounded(error));
    if (error < settings.relativeError || change < 1e-16) {
      return true;
    }
    if (tracker.isEnded()) {
      ImageJUtils.log("WARNING : Drift calculation was interrupted");
      return true;
    }
    return false;
  }

  private static boolean smooth(double[] newDx, double[] newDy, double[] originalDriftTimePoints,
      double smoothing, int iterations) {
    final double[][] values =
        extractValues(originalDriftTimePoints, 0, newDx.length - 1, newDx, newDy);

    // Smooth
    final LoessInterpolator loess = new LoessInterpolator(smoothing, iterations);
    values[1] = loess.smooth(values[0], values[1]);
    values[2] = loess.smooth(values[0], values[2]);

    // Add back
    int count = 0;
    for (int t = 0; t < newDx.length; t++) {
      if (originalDriftTimePoints[t] != 0) {
        newDx[t] = values[1][count];
        newDy[t] = values[2][count];
        count++;

        if (Double.isNaN(newDx[t])) {
          ImageJUtils.log("ERROR : Loess smoothing created bad X-estimate at point %d/%d", t,
              newDx.length);
          return false;
        }
        if (Double.isNaN(newDy[t])) {
          ImageJUtils.log("ERROR : Loess smoothing created bad Y-estimate at point %d/%d", t,
              newDx.length);
          return false;
        }
      }
    }

    return true;
  }

  private void interpolate(double[] dx, double[] dy, double[] originalDriftTimePoints) {
    // Interpolator can only create missing values within the range provided by the input values.
    // The two ends have to be extrapolated.
    // TODO: Perform extrapolation. Currently the end values are used.

    // Find end points
    int startT = 0;
    while (originalDriftTimePoints[startT] == 0) {
      startT++;
    }
    int endT = originalDriftTimePoints.length - 1;
    while (originalDriftTimePoints[endT] == 0) {
      endT--;
    }

    // Extrapolate using a constant value
    for (int t = startT; t-- > 0;) {
      dx[t] = dx[startT];
      dy[t] = dy[startT];
    }
    for (int t = endT; ++t < dx.length;) {
      dx[t] = dx[endT];
      dy[t] = dy[endT];
    }

    final double[][] values = extractValues(originalDriftTimePoints, startT, endT, dx, dy);

    PolynomialSplineFunction fx;
    PolynomialSplineFunction fy;
    if (values[0].length < 3) {
      fx = new LinearInterpolator().interpolate(values[0], values[1]);
      fy = new LinearInterpolator().interpolate(values[0], values[2]);
    } else {
      fx = new SplineInterpolator().interpolate(values[0], values[1]);
      fy = new SplineInterpolator().interpolate(values[0], values[2]);
    }

    for (int t = startT; t <= endT; t++) {
      if (originalDriftTimePoints[t] == 0) {
        dx[t] = fx.value(t);
        dy[t] = fy.value(t);
      }
    }

    this.interpolationStart = startT;
    this.interpolationEnd = endT;
  }

  private static int[] findTimeLimits(MemoryPeakResults results) {
    final StandardResultProcedure sp = new StandardResultProcedure(results);
    sp.getT();
    return MathUtils.limits(sp.frame);
  }

  /**
   * Build a list of the points that are within each roi.
   *
   * @param results the results
   * @param rois the rois
   * @param limits the limits
   * @return the spots
   */
  private static Spot[][] findSpots(MemoryPeakResults results, Roi[] rois, int[] limits) {
    final ArrayList<Spot[]> roiSpots = new ArrayList<>(rois.length);
    for (int i = 0; i < rois.length; i++) {
      final Spot[] spots = findSpots(results, rois[i].getBounds(), limits);
      if (spots.length > 0) {
        roiSpots.add(spots);
      }
    }
    return roiSpots.toArray(new Spot[0][]);
  }

  private static Spot[] findSpots(MemoryPeakResults results, Rectangle bounds, int[] limits) {
    final LocalList<Spot> list = new LocalList<>(limits[1] - limits[0] + 1);
    final float minx = bounds.x;
    final float miny = bounds.y;
    final float maxx = (float) bounds.x + bounds.width;
    final float maxy = (float) bounds.y + bounds.height;

    // Find spots within the ROI
    results.forEach(DistanceUnit.PIXEL, (XyrResultProcedure) (x, y, result) -> {
      if (x > minx && x < maxx && y > miny && y < maxy) {
        list.add(new Spot(result.getFrame(), x, y, result.getIntensity()));
      }
    });

    // For each frame pick the strongest spot
    Collections.sort(list, Spot::compare);

    final LocalList<Spot> newList = new LocalList<>(list.size());

    int currentT = -1;
    for (final Spot spot : list) {
      if (currentT != spot.time) {
        newList.add(spot);
        currentT = spot.time;
      }
    }

    return newList.toArray(new Spot[0]);
  }

  /**
   * For each ROI calculate the sum of the spot intensity. Also compute the sum of the intensity for
   * each time point.
   *
   * @param roiSpots the roi spots
   * @param timepoints the timepoints
   * @param sum The sum of the intensity for each ROI
   * @return The sum of the intensity for each time point.
   */
  private static double[] calculateWeights(Spot[][] roiSpots, int timepoints, double[] sum) {
    final double[] weights = new double[timepoints];
    for (int i = 0; i < roiSpots.length; i++) {
      for (final Spot s : roiSpots[i]) {
        weights[s.time] += s.signal;
        sum[i] += s.signal;
      }
    }
    return weights;
  }

  /**
   * Calculate the drift as displacement of each spot from the centre-of-mass. Update the current
   * drift parameters.
   *
   * @param roiSpots the roi spots
   * @param weights the weights
   * @param sum the sum
   * @param dx the dx
   * @param dy the dy
   * @param smoothing loess smoothing fraction
   * @param iterations Iterations for loess smoothing
   * @return The total update to the drift parameters (Euclidian distance)
   */
  private double calculateDriftUsingMarkers(Spot[][] roiSpots, double[] weights, double[] sum,
      double[] dx, double[] dy, double smoothing, int iterations) {
    final double[] newDx = new double[dx.length];
    final double[] newDy = new double[dy.length];

    // For each ROI
    for (int i = 0; i < roiSpots.length; i++) {
      // Calculate centre-of-mass using the current position (coord + drift)
      double cx = 0;
      double cy = 0;
      for (final Spot s : roiSpots[i]) {
        cx += s.signal * (s.x + dx[s.time]);
        cy += s.signal * (s.y + dy[s.time]);
      }
      cx /= sum[i];
      cy /= sum[i];

      // Calculate update to the drift as centre-of-mass minus the current position (coord + drift)
      for (final Spot s : roiSpots[i]) {
        newDx[s.time] += s.signal * (cx - (s.x + dx[s.time]));
        newDy[s.time] += s.signal * (cy - (s.y + dy[s.time]));
      }
    }

    // Normalise
    for (int t = 0; t < dx.length; t++) {
      if (weights[t] != 0) {
        newDx[t] /= weights[t];
        newDy[t] /= weights[t];
      }

      // New drift = previous drift + update
      newDx[t] += dx[t];
      newDy[t] += dy[t];
    }

    // Store the pure drift values for plotting
    calculatedTimepoints = Arrays.copyOf(weights, weights.length);
    lastdx = Arrays.copyOf(newDx, newDx.length);
    lastdy = Arrays.copyOf(newDy, newDy.length);

    // Perform smoothing
    if (smoothing > 0 && !smooth(newDx, newDy, weights, smoothing, iterations)) {
      return Double.NaN;
    }

    // Average drift correction for the calculated points should be zero to allow change comparison
    normalise(newDx, weights);
    normalise(newDy, weights);

    // Calculate change and update the input drift parameters
    double change = 0;
    for (int t = 0; t < dx.length; t++) {
      if (weights[t] != 0) {
        final double d1 = dx[t] - newDx[t];
        final double d2 = dy[t] - newDy[t];
        change += Math.sqrt(d1 * d1 + d2 * d2);
        dx[t] = newDx[t];
        dy[t] = newDy[t];
      }
    }
    return change;
  }

  /**
   * Normalise the data so that the points identified by non-zeros in the toProcess array have a
   * centre of mass of zero. The shift is calculated on a subset of the points but applied to all
   * points.
   *
   * @param data the data
   * @param toProcess the to process
   */
  private static void normalise(double[] data, double[] toProcess) {
    double av1 = 0;
    int count = 0;
    for (int i = 0; i < data.length; i++) {
      if (toProcess[i] != 0) {
        av1 += data[i];
        count++;
      }
    }
    av1 /= count;

    for (int i = 0; i < data.length; i++) {
      data[i] -= av1;
    }
  }

  /**
   * For all indices between min and max, if the data array is not zero then add the index and the
   * values from array 1 and 2 to the output.
   *
   * @param data the data
   * @param minT the min T
   * @param maxT the max T
   * @param array1 the array 1
   * @param array2 the array 2
   * @return Array of [index][array1][array2]
   */
  private static double[][] extractValues(double[] data, int minT, int maxT, double[] array1,
      double[] array2) {
    // Extract data points for smoothing
    final int timepoints = maxT - minT + 1;
    final double[][] values = new double[3][timepoints];
    int count = 0;
    for (int t = minT; t <= maxT; t++) {
      if (data[t] != 0) {
        values[0][count] = t;
        values[1][count] = array1[t];
        values[2][count] = array2[t];
        count++;
      }
    }
    values[0] = Arrays.copyOf(values[0], count);
    values[1] = Arrays.copyOf(values[1], count);
    values[2] = Arrays.copyOf(values[2], count);
    return values;
  }

  private void plotDrift(int[] limits, double[] dx, double[] dy) {
    if (!settings.plotDrift) {
      return;
    }

    // Build an array of timepoints from the min to the max
    final double[] completeT = new double[limits[1] + 1];
    for (int i = limits[0]; i < completeT.length; i++) {
      completeT[i] = i;
    }

    // Drift should be centred around zero for the calculated points to produce a fair plot
    normalise(lastdx, calculatedTimepoints);
    normalise(lastdy, calculatedTimepoints);

    // Extract the interpolated points and the original drift
    final double[][] interpolated = extractValues(dx, limits[0], limits[1], dx, dy);
    final double[][] original =
        extractValues(calculatedTimepoints, limits[0], limits[1], lastdx, lastdy);

    final PlotWindow window = plotDrift(null, interpolated, original, "Drift X", 1);
    ploty.set(plotDrift(window, interpolated, original, "Drift Y", 2));
    plotx.set(window);
  }

  private static PlotWindow plotDrift(PlotWindow parent, double[][] interpolated,
      double[][] original, String name, int index) {
    // Create plot
    final double[] xlimits = MathUtils.limits(interpolated[0]);
    double[] ylimits = MathUtils.limits(original[index]);
    ylimits = MathUtils.limits(ylimits, interpolated[index]);

    final Plot plot = new Plot(name, "Frame", "Drift (px)");
    plot.setLimits(xlimits[0], xlimits[1], ylimits[0], ylimits[1]);
    plot.setColor(new Color(0, 0, 155)); // De-saturated blue
    plot.addPoints(original[0], original[index], Plot.CROSS);
    plot.setColor(java.awt.Color.RED);
    plot.addPoints(interpolated[0], interpolated[index], Plot.LINE);
    final WindowOrganiser wo = new WindowOrganiser();
    final PlotWindow window = ImageJUtils.display(name, plot, wo);

    if (wo.isNotEmpty() && parent != null) {
      final Point location = parent.getLocation();
      location.y += parent.getHeight();
      window.setLocation(location);
    }

    return window;
  }

  /**
   * Saves the T,X,Y values to file for all t in the originalDriftTimePoints array which are not
   * zero.
   *
   * @param originalDriftTimePoints the original drift time points
   * @param dx the dx
   * @param dy the dy
   */
  private void saveDrift(double[] originalDriftTimePoints, double[] dx, double[] dy) {
    if (!settings.saveDrift) {
      return;
    }
    if (!getDriftFilename()) {
      return;
    }
    try (BufferedWriter out = Files.newBufferedWriter(Paths.get(settings.driftFilename))) {
      out.write("Time\tX\tY");
      out.newLine();
      for (int t = 0; t < dx.length; t++) {
        if (originalDriftTimePoints[t] != 0) {
          out.write(String.format("%d\t%f\t%f%n", t, dx[t], dy[t]));
        }
      }
      ImageJUtils.log("Saved calculated drift to file: " + settings.driftFilename);
    } catch (final IOException ex) {
      // Do nothing
    }
  }

  /**
   * Calculates drift using T,X,Y records read from a file.
   *
   * @param limits the limits
   * @return the drift { dx[], dy[] }
   */
  @Nullable
  private double[][] calculateUsingDriftFile(int[] limits) {
    // Read drift TXY from file
    calculatedTimepoints = new double[limits[1] + 1];
    lastdx = new double[calculatedTimepoints.length];
    lastdy = new double[calculatedTimepoints.length];

    if (!getDriftFilename()) {
      return null;
    }

    if (readDriftFile(limits) < 2) {
      ImageJUtils.log("ERROR : Not enough drift points within the time limits %d - %d", limits[0],
          limits[1]);
      return null;
    }

    final double[] dx = Arrays.copyOf(lastdx, lastdx.length);
    final double[] dy = Arrays.copyOf(lastdy, lastdy.length);

    final double smoothing = updateSmoothingParameter(calculatedTimepoints);

    // Perform smoothing
    if (smoothing > 0 && !smooth(dx, dy, calculatedTimepoints, smoothing, settings.iterations)) {
      return null;
    }

    // Average drift correction for the calculated points should be zero
    normalise(dx, calculatedTimepoints);
    normalise(dy, calculatedTimepoints);

    interpolate(dx, dy, calculatedTimepoints);

    plotDrift(limits, dx, dy);

    return new double[][] {dx, dy};
  }

  private boolean getDriftFilename() {
    final String[] path = ImageJUtils.decodePath(settings.driftFilename);
    final OpenDialog chooser = new OpenDialog("Drift_file", path[0], path[1]);
    if (chooser.getFileName() == null) {
      return false;
    }
    settings.driftFilename = chooser.getDirectory() + chooser.getFileName();
    FileUtils.replaceExtension(settings.driftFilename, "tsv");
    return true;
  }

  /**
   * Read the drift file storing the T,X,Y into the class level calculatedTimepoints, lastdx and
   * lastdy arrays. Ignore any records where T is outside the limits.
   *
   * @param limits the limits
   * @return The number of records read
   */
  private int readDriftFile(int[] limits) {
    int ok = 0;
    try (BufferedReader input = Files.newBufferedReader(Paths.get(settings.driftFilename))) {
      String line;
      final Pattern pattern = Pattern.compile("[\t, ]+");
      while ((line = input.readLine()) != null) {
        if (line.length() == 0) {
          continue;
        }
        if (Character.isDigit(line.charAt(0))) {
          try (Scanner scanner = new Scanner(line)) {
            scanner.useDelimiter(pattern);
            scanner.useLocale(Locale.US);
            final int t = scanner.nextInt();
            if (t < limits[0] || t > limits[1]) {
              continue;
            }
            final double x = scanner.nextDouble();
            final double y = scanner.nextDouble();
            calculatedTimepoints[t] = ++ok;
            lastdx[t] = x;
            lastdy[t] = y;
          } catch (final NoSuchElementException ex) {
            // Do nothing
          }
        }
      }
    } catch (final IOException ex) {
      // ignore
    }
    return ok;
  }

  private static class BlockPeakResultProcedure implements PeakResultProcedure {
    final ArrayList<ArrayList<Localisation>> blocks = new ArrayList<>();
    ArrayList<Localisation> nextBlock;
    final Counter counter = new Counter();
    final int frames;
    final int minimimLocalisations;

    BlockPeakResultProcedure(Settings settings) {
      this.frames = settings.frames;
      this.minimimLocalisations = settings.minimimLocalisations;
    }

    @Override
    public void execute(PeakResult result) {
      if (result.getFrame() > counter.getCount()) {
        while (result.getFrame() > counter.getCount()) {
          counter.increment(frames);
        }
        // To avoid blocks without many results only create a new block if the min size has been met
        if (nextBlock == null || nextBlock.size() >= minimimLocalisations) {
          nextBlock = new ArrayList<>();
        }
        blocks.add(nextBlock);
      }
      nextBlock.add(new Localisation(result.getFrame(), result.getXPosition(),
          result.getYPosition(), result.getIntensity()));
    }
  }

  /**
   * Calculates drift using images from N consecutive frames aligned to the overall image.
   *
   * @param results the results
   * @param limits the limits
   * @param reconstructionSize the reconstruction size
   * @return the drift { dx[], dy[] }
   */
  @Nullable
  private double[][] calculateUsingFrames(MemoryPeakResults results, int[] limits,
      int reconstructionSize) {
    // Extract the localisations into blocks of N consecutive frames
    final BlockPeakResultProcedure p = new BlockPeakResultProcedure(settings);
    results.sort();
    results.forEach(p);

    final ArrayList<ArrayList<Localisation>> blocks = p.blocks;
    final ArrayList<Localisation> nextBlock = p.nextBlock;

    if (blocks.size() < 2) {
      tracker.log("ERROR : Require at least 2 images for drift calculation");
      return null;
    }

    // Check the final block has enough localisations
    if (nextBlock.size() < settings.minimimLocalisations) {
      blocks.remove(blocks.size() - 1);
      if (blocks.size() < 2) {
        tracker.log("ERROR : Require at least 2 images for drift calculation");
        return null;
      }

      final ArrayList<Localisation> combinedBlock = blocks.get(blocks.size() - 1);
      combinedBlock.addAll(nextBlock);
    }

    // Find the average time point for each block
    final int[] blockT = new int[blocks.size()];
    int time = 0;
    for (final ArrayList<Localisation> block : blocks) {
      long sum = 0;
      for (final Localisation r : block) {
        sum += r.time;
      }
      blockT[time++] = (int) (sum / block.size());
    }

    // Calculate a scale to use when constructing the images for alignment
    final Rectangle bounds = results.getBounds(true);
    final float scale = (reconstructionSize - 1f) / Math.max(bounds.width, bounds.height);

    executor = Executors.newFixedThreadPool(Prefs.getThreads());

    final double[] dx = new double[limits[1] + 1];
    final double[] dy = new double[dx.length];

    final double[] originalDriftTimePoints = getOriginalDriftTimePoints(dx, blockT);
    lastdx = null;

    final double smoothing = updateSmoothingParameter(originalDriftTimePoints);

    double change = calculateDriftUsingFrames(blocks, blockT, bounds, scale, dx, dy,
        originalDriftTimePoints, smoothing, settings.iterations);
    if (Double.isNaN(change) || tracker.isEnded()) {
      return null;
    }

    plotDrift(limits, dx, dy);
    ImageJUtils.log("Drift Calculator : Initial drift " + MathUtils.rounded(change));

    for (int i = 1; i <= settings.maxIterations; i++) {
      change = calculateDriftUsingFrames(blocks, blockT, bounds, scale, dx, dy,
          originalDriftTimePoints, smoothing, settings.iterations);
      if (Double.isNaN(change)) {
        return null;
      }

      plotDrift(limits, dx, dy);

      if (converged(i, change, getTotalDrift(dx, dy, originalDriftTimePoints))) {
        break;
      }
    }

    if (tracker.isEnded()) {
      return null;
    }

    plotDrift(limits, dx, dy);

    return new double[][] {dx, dy};
  }

  /**
   * Create an array to show the time-point of the original calculated drift alignment.
   *
   * @param dx The drift array
   * @param timepoints Array of timepoints for which there is a drift calculation
   * @return array matching dx length with non-zero values for each identified timepoint
   */
  private static double[] getOriginalDriftTimePoints(double[] dx, int[] timepoints) {
    final double[] originalDriftTimePoints = new double[dx.length];
    for (int i = 0; i < timepoints.length; i++) {
      originalDriftTimePoints[timepoints[i]] = 1;
    }
    return originalDriftTimePoints;
  }

  /**
   * Calculate the drift by aligning N consecutive frames with the overall image. Update the current
   * drift parameters.
   *
   * @param blocks the blocks
   * @param blockT the block T
   * @param bounds the bounds
   * @param scale the scale
   * @param dx the dx
   * @param dy the dy
   * @param originalDriftTimePoints the original drift time points
   * @param smoothing the smoothing
   * @param iterations the iterations
   * @return the double
   */
  private double calculateDriftUsingFrames(ArrayList<ArrayList<Localisation>> blocks, int[] blockT,
      Rectangle bounds, float scale, double[] dx, double[] dy, double[] originalDriftTimePoints,
      double smoothing, int iterations) {
    // Construct images using the current drift
    tracker.status("Constructing images");

    // Built an image for each block of results.
    final ImageProcessor[] images = new ImageProcessor[blocks.size()];

    final List<Future<?>> futures = new LinkedList<>();
    final Ticker ticker = Ticker.createStarted(tracker, images.length * 2L, true);

    for (int i = 0; i < images.length; i++) {
      futures.add(executor
          .submit(new ImageBuilder(blocks.get(i), images, i, bounds, scale, dx, dy, ticker)));
    }
    ConcurrencyUtils.waitForCompletionUnchecked(futures);

    for (int i = 0; i < blocks.size(); i++) {
      tracker.progress(i, blocks.size());
      final ImageJImagePeakResults blockImage = newImage(bounds, scale);
      for (final Localisation r : blocks.get(i)) {
        blockImage.add(r.time, (float) (r.x + dx[r.time]), (float) (r.y + dy[r.time]), r.signal);
      }
      images[i] = getImage(blockImage);
    }

    // Build an image with all results.
    final FloatProcessor allIp = new FloatProcessor(images[0].getWidth(), images[0].getHeight());
    for (final ImageProcessor ip : images) {
      allIp.copyBits(ip, 0, 0, Blitter.ADD);
    }

    return calculateDrift(blockT, scale, dx, dy, originalDriftTimePoints, smoothing, iterations,
        images, allIp, true, ticker);
  }

  /**
   * Calculate the drift of images to the reference image. Update the current drift parameters.
   *
   * @param imageT The frame number for each image
   * @param scale The image scale (used to adjust the drift to the correct size)
   * @param dx The X drift
   * @param dy The Y drift
   * @param originalDriftTimePoints Non-zero when the frame number refers to an aligned image frame
   * @param smoothing LOESS smoothing parameter
   * @param iterations LOESS iterations parameter
   * @param images The images to align
   * @param reference The reference image
   * @param includeCurrentDrift Set to true if the input images already have the current drift
   *        applied. The new drift will be added to the current drift.
   * @param ticker the ticker
   * @return the change in the drift
   */
  private double calculateDrift(int[] imageT, float scale, double[] dx, double[] dy,
      double[] originalDriftTimePoints, double smoothing, int iterations,
      final ImageProcessor[] images, FloatProcessor reference, boolean includeCurrentDrift,
      Ticker ticker) {
    // Align
    tracker.status("Aligning images");
    final AlignImagesFft aligner = new AlignImagesFft();
    aligner.initialiseReference(reference, WindowMethod.NONE, false);
    final Rectangle alignBounds = AlignImagesFft.createHalfMaxBounds(reference.getWidth(),
        reference.getHeight(), reference.getWidth(), reference.getHeight());

    final List<double[]> alignments =
        Collections.synchronizedList(new ArrayList<double[]>(images.length));
    final List<Future<?>> futures = new LinkedList<>();

    final int imagesPerThread = getImagesPerThread(images);
    for (int i = 0; i < images.length; i += imagesPerThread) {
      futures.add(executor.submit(new ImageAligner(aligner, images, imageT, alignBounds, alignments,
          i, i + imagesPerThread, ticker, settings.subPixelMethod)));
    }
    ConcurrencyUtils.waitForCompletionUnchecked(futures);
    tracker.progress(1);

    // Used to flag when an alignment has failed
    originalDriftTimePoints =
        Arrays.copyOf(originalDriftTimePoints, originalDriftTimePoints.length);
    final double[] newDx = new double[dx.length];
    final double[] newDy = new double[dy.length];
    int ok = 0;
    for (final double[] result : alignments) {
      final int t = (int) result[2];
      if (Double.isNaN(result[0])) {
        // Q: How to ignore bad alignments?
        // Only do smoothing where there was an alignment?
        originalDriftTimePoints[t] = 0;
        tracker.log("WARNING : Unable to align image for time %d to the overall projection", t);
      } else {
        ok++;
        newDx[t] = result[0] / scale;
        newDy[t] = result[1] / scale;
        if (includeCurrentDrift) {
          // New drift = update + previous drift
          newDx[t] += dx[t];
          newDy[t] += dy[t];
        }
      }
    }

    if (ok < 2) {
      tracker.log("ERROR : Unable to align more than 1 image to the overall projection");
      return Double.NaN;
    }

    // Store the pure drift values for plotting
    calculatedTimepoints = Arrays.copyOf(originalDriftTimePoints, originalDriftTimePoints.length);
    lastdx = Arrays.copyOf(newDx, newDx.length);
    lastdy = Arrays.copyOf(newDy, newDy.length);

    // Perform smoothing
    if (smoothing > 0) {
      tracker.status("Smoothing drift");
      if (!smooth(newDx, newDy, originalDriftTimePoints, smoothing, iterations)) {
        return Double.NaN;
      }
    }

    // Interpolate values for all time limits
    tracker.status("Interpolating drift");
    interpolate(newDx, newDy, originalDriftTimePoints);

    // Average drift correction for the calculated points should be zero to allow change comparison
    normalise(newDx, originalDriftTimePoints);
    normalise(newDy, originalDriftTimePoints);

    // Calculate change and update the input drift parameters
    double change = 0;
    for (int t = 0; t < dx.length; t++) {
      if (originalDriftTimePoints[t] != 0) {
        final double d1 = dx[t] - newDx[t];
        final double d2 = dy[t] - newDy[t];
        change += Math.sqrt(d1 * d1 + d2 * d2);
      }
      // Update all points since interpolation has already been done
      dx[t] = newDx[t];
      dy[t] = newDy[t];
    }

    tracker.status("");
    return change;
  }

  /**
   * Get the number of images that should be processed on each thread.
   *
   * @param images The list of images
   * @return The images per thread
   */
  private static int getImagesPerThread(final ImageProcessor[] images) {
    return Math.max(1, (int) Math.round((double) images.length / Prefs.getThreads()));
  }

  private static ImageJImagePeakResults newImage(Rectangle bounds, float imageScale) {
    final ImageJImagePeakResults image =
        ImagePeakResultsFactory.createPeakResultsImage(ResultsImageType.DRAW_INTENSITY, true, false,
            "", bounds, 100, 1, imageScale, 0, ResultsImageMode.IMAGE_ADD);
    image.setDisplayImage(false);
    image.begin();
    return image;
  }

  private static ImageProcessor getImage(ImageJImagePeakResults imageResults) {
    imageResults.end();
    return imageResults.getImagePlus().getProcessor();
  }

  /**
   * Calculates drift using images from a reference stack aligned to the overall z-projection.
   *
   * @param stack the stack
   * @param limits the limits
   * @return the drift { dx[], dy[] }
   */
  @Nullable
  private double[][] calculateUsingImageStack(ImageStack stack, int[] limits) {
    // Update the limits using the stack size
    final int upperT = settings.startFrame + settings.frameSpacing * (stack.getSize() - 1);
    limits[1] = Math.max(limits[1], upperT);

    // TODO - Truncate the stack if there are far too many frames for the localisation limits

    tracker.status("Constructing images");

    executor = Executors.newFixedThreadPool(Prefs.getThreads());

    // Built an image and FHT image for each slice
    final ImageProcessor[] images = new ImageProcessor[stack.getSize()];
    final Fht[] fhtImages = new Fht[stack.getSize()];

    final List<Future<?>> futures = new LinkedList<>();
    final Ticker ticker = Ticker.createStarted(tracker, images.length, true);

    final int imagesPerThread = getImagesPerThread(images);
    final AlignImagesFft aligner = new AlignImagesFft();
    final FloatProcessor referenceIp = stack.getProcessor(1).toFloat(0, null);
    // We do not care about the window method because this processor will not
    // actually be used for alignment, it is a reference for the FHT size
    aligner.initialiseReference(referenceIp, WindowMethod.NONE, false);
    for (int i = 0; i < images.length; i += imagesPerThread) {
      futures.add(executor.submit(new ImageFhtInitialiser(stack, images, aligner, fhtImages, i,
          i + imagesPerThread, ticker)));
    }
    ConcurrencyUtils.waitForCompletionUnchecked(futures);
    tracker.progress(1);

    if (tracker.isEnded()) {
      return null;
    }

    final double[] dx = new double[limits[1] + 1];
    final double[] dy = new double[dx.length];

    final double[] originalDriftTimePoints = new double[dx.length];
    final int[] blockT = new int[stack.getSize()];
    for (int i = 0, t = settings.startFrame; i < stack.getSize(); i++, t += settings.frameSpacing) {
      originalDriftTimePoints[t] = 1;
      blockT[i] = t;
    }

    final double smoothing = updateSmoothingParameter(originalDriftTimePoints);

    lastdx = null;
    // For the first iteration calculate drift to the first image in the stack
    // (since the average projection may have a large drift blurring the image)
    double change = calculateDriftUsingImageStack(referenceIp, images, fhtImages, blockT, dx, dy,
        originalDriftTimePoints, smoothing, settings.iterations);
    if (Double.isNaN(change) || tracker.isEnded()) {
      return null;
    }

    plotDrift(limits, dx, dy);
    ImageJUtils.log("Drift Calculator : Initial drift " + MathUtils.rounded(change));

    for (int i = 1; i <= settings.maxIterations; i++) {
      change = calculateDriftUsingImageStack(null, images, fhtImages, blockT, dx, dy,
          originalDriftTimePoints, smoothing, settings.iterations);
      if (Double.isNaN(change)) {
        return null;
      }

      plotDrift(limits, dx, dy);

      if (converged(i, change, getTotalDrift(dx, dy, originalDriftTimePoints))) {
        break;
      }
    }

    if (tracker.isEnded()) {
      return null;
    }

    plotDrift(limits, dx, dy);

    return new double[][] {dx, dy};
  }

  /**
   * Calculate the drift of images to the reference image. If no reference is provided then produce
   * a combined z-projection. Update the current drift parameters.
   *
   * @param reference the reference
   * @param images The images to align
   * @param fhtImages The images to align (pre-transformed to a FHT)
   * @param blockT The frame number for each image
   * @param dx The X drift
   * @param dy The Y drift
   * @param originalDriftTimePoints Non-zero when the frame number refers to an aligned image frame
   * @param smoothing the smoothing
   * @param iterations the iterations
   * @return The change in the drift (NaN is an error occurred)
   */
  private double calculateDriftUsingImageStack(FloatProcessor reference, ImageProcessor[] images,
      Fht[] fhtImages, int[] blockT, double[] dx, double[] dy, double[] originalDriftTimePoints,
      double smoothing, int iterations) {
    Ticker ticker = Ticker.createStarted(tracker, images.length, true);

    if (reference == null) {
      // Construct images using the current drift
      tracker.status("Constructing reference image");

      // Build an image using the current drift
      final List<Future<?>> futures = new LinkedList<>();
      ticker = Ticker.createStarted(tracker, images.length * 2L, true);

      final ImageProcessor[] blockIp = new ImageProcessor[images.length];
      final double[] threadDx = new double[images.length];
      final double[] threadDy = new double[images.length];
      for (int i = 0; i < images.length; i++) {
        threadDx[i] = dx[blockT[i]];
        threadDy[i] = dy[blockT[i]];
      }

      final int imagesPerThread = getImagesPerThread(images);
      for (int i = 0; i < images.length; i += imagesPerThread) {
        futures.add(executor.submit(new ImageTranslator(images, blockIp, threadDx, threadDy, i,
            i + imagesPerThread, ticker, settings.interpolationMethod)));
      }
      ConcurrencyUtils.waitForCompletionUnchecked(futures);

      // Build an image with all results.
      reference = new FloatProcessor(blockIp[0].getWidth(), blockIp[0].getHeight());
      for (final ImageProcessor ip : blockIp) {
        reference.copyBits(ip, 0, 0, Blitter.ADD);
      }
    }

    // Ensure the reference is windowed
    AlignImagesFft.applyWindowSeparable(reference, WindowMethod.TUKEY);

    return calculateDrift(blockT, 1f, dx, dy, originalDriftTimePoints, smoothing, iterations,
        fhtImages, reference, false, ticker);
  }
}
