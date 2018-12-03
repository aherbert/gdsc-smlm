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

import java.awt.Rectangle;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.math3.random.HaltonSequenceGenerator;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.plugin.PlugIn;
import ij.text.TextWindow;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog.OptionListener;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.TurboList;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CameraType;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSFType;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtosHelper;
import uk.ac.sussex.gdsc.smlm.engine.FitConfiguration;
import uk.ac.sussex.gdsc.smlm.engine.FitEngineConfiguration;
import uk.ac.sussex.gdsc.smlm.fitting.FitStatus;
import uk.ac.sussex.gdsc.smlm.fitting.FunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.Gaussian2DFitter;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.ij.IJImageSource;
import uk.ac.sussex.gdsc.smlm.ij.plugins.CreateData.BenchmarkParameters;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.ij.utils.IJImageConverter;
import uk.ac.sussex.gdsc.smlm.model.camera.CameraModel;

/**
 * Fits the benchmark image created by CreateData plugin.
 */
public class BenchmarkFit implements PlugIn {
  private static final String TITLE = "Benchmark Fit";

  private static int regionSize = 4;
  private static double lastId = 0;

  private static final String[] ORIGIN_XY = {"Origin", "Centre-of-Mass", "Offset"};
  private static final String[] ORIGIN_Z = {"Origin", "0", "Offset"};
  private static int originXY = 0;
  private static int originZ = 0;
  private static double offsetX = 0;
  private static double offsetY = 0;
  private static double offsetZ = 0;

  private static boolean zeroOffset = true;
  private static double offsetPoints = 0;
  private static double offsetRangeX = 0.5;
  private static double offsetRangeY = 0.5;
  private static double offsetRangeZ = 0.5;

  private static boolean backgroundFitting = true;
  private static boolean estimateBackground = true;
  private static boolean signalFitting = true;
  private static boolean estimateSignal = true;
  private static boolean showHistograms = false;
  private static boolean saveRawData = false;
  private static String rawDataDirectory = "";
  private static int histogramBins = 100;

  private static TextWindow summaryTable = null, analysisTable = null;

  //@formatter:off
  // These are assuming a Gaussian 2D PSF
  private static final String[] NAMES = new String[] {
      "dB (photons)",
      "dSignal (photons)",
      "dX (nm)",
      "dY (nm)",
      "dZ (nm)",
      "dSx (nm)",
      "dSy (nm)",
      "dAngle (deg)",
      "Time (ms)",
      "dActualSignal (photons)",
      "dSax (nm)",
      "dSay (nm)" };
  //@formatter:on
  private static final int TIME = 8;
  private static final int ACTUAL_SIGNAL = 9;
  private static final int ADJUSTED_X_SD = 10;
  private static final int ADJUSTED_Y_SD = 11;
  private static boolean[] displayHistograms = new boolean[NAMES.length];
  static {
    for (int i = 0; i < displayHistograms.length; i++) {
      displayHistograms[i] = true;
    }
  }

  private FitConfiguration fitConfig;
  private ImagePlus imp;
  private CreateData.BenchmarkParameters benchmarkParameters;
  private final double[] answer = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
  private Rectangle region = null;
  private final AtomicInteger comValid = new AtomicInteger();

  // Used to store all the results for cross-method comparison
  private double[][] results;
  private long[] resultsTime;

  /**
   * Store the benchmark result.
   */
  public class BenchmarkResult {
    /**
     * The parameters used to create the data.
     */
    final BenchmarkParameters benchmarkParameters;
    /**
     * The actual parameters (with XY position adjusted for the region size).
     */
    final double[] answer;
    /**
     * The string description of the parameters used to create and then fit the data.
     */
    final String parameters;
    /**
     * Results conversion factors.
     */
    final double[] convert;
    /**
     * The results of fitting the data. Results are only stored if fitting was successful.
     */
    final double[][] results;
    /**
     * The time for fitting the data.
     */
    final long[] resultsTime;

    /**
     * Instantiates a new benchmark result.
     *
     * @param benchmarkParameters the benchmark parameters
     * @param answer the answer
     * @param parameters the parameters
     * @param convert the convert
     * @param results the results
     * @param resultsTime the results time
     */
    public BenchmarkResult(BenchmarkParameters benchmarkParameters, double[] answer,
        String parameters, double[] convert, double[][] results, long[] resultsTime) {
      this.benchmarkParameters = benchmarkParameters;
      this.answer = answer;
      this.parameters = parameters;
      this.convert = convert;
      this.results = results;
      this.resultsTime = resultsTime;
    }
  }

  /**
   * Store all the results from fitting on the same benchmark dataset.
   */
  public static LinkedList<BenchmarkResult> benchmarkResults = new LinkedList<>();

  /**
   * Used to allow multi-threading of the fitting method.
   */
  private class Worker implements Runnable {
    volatile boolean finished = false;
    final BlockingQueue<Integer> jobs;
    final Statistics[] stats = new Statistics[NAMES.length];
    final ImageStack stack;
    final Rectangle region;
    final double[][] offsets;
    final FitConfiguration fitConfig;
    final CameraModel cameraModel;
    final double sa;
    final int size;
    final int totalFrames;
    final double[] origin;

    float[] data = null;
    private double[] lb, ub = null;
    private double[] lc, uc = null;

    public Worker(BlockingQueue<Integer> jobs, ImageStack stack, Rectangle region,
        FitConfiguration fitConfig, CameraModel cameraModel) {
      this.jobs = jobs;
      this.stack = stack;
      this.region = region;
      this.fitConfig = fitConfig.clone();
      this.cameraModel = cameraModel;
      this.offsets = startPoints;

      for (int i = 0; i < stats.length; i++) {
        stats[i] = (showHistograms || saveRawData) ? new StoredDataStatistics() : new Statistics();
      }
      sa = getSa();

      createBounds();

      size = region.height;
      totalFrames = benchmarkParameters.frames;
      origin = new double[] {answer[Gaussian2DFunction.X_POSITION],
          answer[Gaussian2DFunction.Y_POSITION], answer[Gaussian2DFunction.Z_POSITION]};
    }

    /** {@inheritDoc} */
    @Override
    public void run() {
      try {
        while (true) {
          final Integer job = jobs.take();
          if (job == null || job.intValue() < 0) {
            break;
          }
          if (!finished) {
            // Only run if not finished to allow queue to be emptied
            run(job.intValue());
          }
        }
      } catch (final InterruptedException ex) {
        System.out.println(ex.toString());
        throw new RuntimeException(ex);
      } finally {
        finished = true;
      }
    }

    private void run(int frame) {
      if (ImageJUtils.isInterrupted()) {
        finished = true;
        return;
      }

      showProgress();

      // Extract the data
      data = IJImageConverter.getData(stack.getPixels(frame + 1), stack.getWidth(),
          stack.getHeight(), region, data);

      // Use a camera model to pre-process the data.
      // The camera model is cropped to the correct region size.
      if (fitConfig.isFitCameraCounts()) {
        cameraModel.removeBias(data);
      } else {
        cameraModel.removeBiasAndGain(data);
      }

      final double[] data = SimpleArrayUtils.toDouble(this.data);

      // Get the background and signal estimate for fitting in the correct units
      final double b = (backgroundFitting && estimateBackground) ? getBackground(data, size, size)
          // Convert the answer to the correct units
          : answer[Gaussian2DFunction.BACKGROUND]
              * ((fitConfig.isFitCameraCounts()) ? benchmarkParameters.gain : 1);
      final double signal = (signalFitting && estimateSignal) ? getSignal(data, b)
          // Convert the answer to the correct units
          : answer[Gaussian2DFunction.SIGNAL]
              * ((fitConfig.isFitCameraCounts()) ? benchmarkParameters.gain : 1);

      // Find centre-of-mass estimate
      final double[] com = new double[2];
      getCentreOfMass(data, size, size, com);

      // Update the origin from the answer
      if (originXY == 1) {
        // Use Centre of mass as the origin
        origin[0] = com[0];
        origin[1] = com[1];
      }
      if (originZ == 1) {
        // Use zero as the origin
        origin[2] = 0;
      }

      final double dx = com[0] - answer[Gaussian2DFunction.X_POSITION];
      final double dy = com[1] - answer[Gaussian2DFunction.Y_POSITION];
      if (Math.abs(dx) < offsetRangeX && Math.abs(dy) < offsetRangeY) {
        comValid.getAndIncrement();
      }

      final double[] initialParams = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
      initialParams[Gaussian2DFunction.BACKGROUND] = b;
      initialParams[Gaussian2DFunction.SIGNAL] = signal;
      initialParams[Gaussian2DFunction.X_SD] = fitConfig.getInitialXSD();
      initialParams[Gaussian2DFunction.Y_SD] = fitConfig.getInitialYSD();

      double[][] bounds = null;
      final double[][] result = new double[offsets.length][];
      final long[] time = new long[offsets.length];
      int c = 0;
      int resultPosition = frame;
      for (final double[] offset : offsets) {
        final long start = System.nanoTime();

        // Do fitting
        final double[] params = initialParams.clone();
        params[Gaussian2DFunction.X_POSITION] = origin[0] + offset[0];
        params[Gaussian2DFunction.Y_POSITION] = origin[1] + offset[1];
        params[Gaussian2DFunction.Z_POSITION] = origin[2] + offset[2];
        fitConfig.initialise(1, size, size, params);
        final FunctionSolver solver = fitConfig.getFunctionSolver();
        if (solver.isBounded()) {
          bounds = setBounds(solver, initialParams, bounds);
        } else if (solver.isConstrained()) {
          setConstraints(solver);
        }

        final FitStatus status = solver.fit(data, null, params, null);
        if (isValid(status, params, size)) {
          // TODO - Check this is OK for the MLE camera model
          // That estimates counts. We require an estimate of photons.
          if (fitConfig.isFitCameraCounts()) {
            // Update all the parameters to be in photons
            params[Gaussian2DFunction.BACKGROUND] /= benchmarkParameters.gain;
            params[Gaussian2DFunction.SIGNAL] /= benchmarkParameters.gain;
          }
          result[c] = params;
          time[c] = System.nanoTime() - start;
          // Store all the results for later analysis
          results[resultPosition] = params;
          resultsTime[resultPosition] = time[c];
          c++;
        } else {
          // System.out.println(status);
        }
        resultPosition += totalFrames;
      }

      addResults(stats, answer, benchmarkParameters.p[frame], sa, time, result, c);
    }

    /**
     * Set background using the average value of the edge in the data.
     *
     * @param data the data
     * @param maxx the maxx
     * @param maxy the maxy
     * @return The background
     */
    private double getBackground(double[] data, int maxx, int maxy) {
      return Gaussian2DFitter.getBackground(data, maxx, maxy, 1);
    }

    private double[][] setBounds(FunctionSolver solver, double[] params, double[][] bounds) {
      if (bounds == null) {
        double[] lower = null;
        double[] upper = null;
        // Check the bounds
        if (params[Gaussian2DFunction.BACKGROUND] < lb[Gaussian2DFunction.BACKGROUND]) {
          lower = lb.clone();
          lower[Gaussian2DFunction.BACKGROUND] = params[Gaussian2DFunction.BACKGROUND]
              - Math.abs(lb[Gaussian2DFunction.BACKGROUND] - params[Gaussian2DFunction.BACKGROUND]);
          if (fitConfig.requireStrictlyPositiveFunction()
              && lower[Gaussian2DFunction.BACKGROUND] < 0) {
            lower[Gaussian2DFunction.BACKGROUND] = 0;
          }
        }
        if (params[Gaussian2DFunction.BACKGROUND] > ub[Gaussian2DFunction.BACKGROUND]) {
          upper = ub.clone();
          upper[Gaussian2DFunction.BACKGROUND] = params[Gaussian2DFunction.BACKGROUND]
              + Math.abs(ub[Gaussian2DFunction.BACKGROUND] - params[Gaussian2DFunction.BACKGROUND]);
        }
        if (params[Gaussian2DFunction.SIGNAL] < lb[Gaussian2DFunction.SIGNAL]) {
          if (lower == null) {
            lower = lb.clone();
          }
          lower[Gaussian2DFunction.SIGNAL] = params[Gaussian2DFunction.SIGNAL]
              - Math.abs(lb[Gaussian2DFunction.SIGNAL] - params[Gaussian2DFunction.SIGNAL]);
          if (fitConfig.requireStrictlyPositiveFunction() && lower[Gaussian2DFunction.SIGNAL] < 0) {
            lower[Gaussian2DFunction.SIGNAL] = 0;
          }
        }
        if (params[Gaussian2DFunction.SIGNAL] > ub[Gaussian2DFunction.SIGNAL]) {
          if (upper == null) {
            upper = ub.clone();
          }
          upper[Gaussian2DFunction.SIGNAL] = params[Gaussian2DFunction.SIGNAL]
              + Math.abs(ub[Gaussian2DFunction.SIGNAL] - params[Gaussian2DFunction.SIGNAL]);
        }
        if (lower == null) {
          lower = lb;
        }
        if (upper == null) {
          upper = ub;
        }
        bounds = new double[][] {lower, upper};
      }
      solver.setBounds(bounds[0], bounds[1]);
      return bounds;
    }

    private void createBounds() {
      if (ub == null) {
        ub = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
        lb = new double[ub.length];

        // Background could be zero so always have an upper limit
        ub[Gaussian2DFunction.BACKGROUND] = Math.max(0, 2 * benchmarkParameters.getBackground());
        final double signal = benchmarkParameters.getSignal();
        lb[Gaussian2DFunction.SIGNAL] = signal * 0.5;
        ub[Gaussian2DFunction.SIGNAL] = signal * 2;
        ub[Gaussian2DFunction.X_POSITION] = 2 * regionSize + 1;
        ub[Gaussian2DFunction.Y_POSITION] = 2 * regionSize + 1;
        lb[Gaussian2DFunction.ANGLE] = -Math.PI;
        ub[Gaussian2DFunction.ANGLE] = Math.PI;
        lb[Gaussian2DFunction.Z_POSITION] = Double.NEGATIVE_INFINITY;
        ub[Gaussian2DFunction.Z_POSITION] = Double.POSITIVE_INFINITY;
        final double wf = 1.5;
        final double s = benchmarkParameters.s / benchmarkParameters.a;
        lb[Gaussian2DFunction.X_SD] = s / wf;
        ub[Gaussian2DFunction.X_SD] = s * wf;
        lb[Gaussian2DFunction.Y_SD] = s / wf;
        ub[Gaussian2DFunction.Y_SD] = s * wf;
        if (fitConfig.isFitCameraCounts()) {
          final double gain = benchmarkParameters.gain;
          // Update all the parameters affected by gain
          lb[Gaussian2DFunction.BACKGROUND] *= gain;
          lb[Gaussian2DFunction.SIGNAL] *= gain;
          ub[Gaussian2DFunction.BACKGROUND] *= gain;
          ub[Gaussian2DFunction.SIGNAL] *= gain;
        }
      }
    }

    private void setConstraints(FunctionSolver solver) {
      if (uc == null) {
        lc = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
        uc = new double[lc.length];
        Arrays.fill(lc, Float.NEGATIVE_INFINITY);
        Arrays.fill(uc, Float.POSITIVE_INFINITY);
        lc[Gaussian2DFunction.BACKGROUND] = 0;
        lc[Gaussian2DFunction.SIGNAL] = 0;
      }
      solver.setConstraints(lc, uc);
    }

    private boolean isValid(FitStatus status, double[] params, int size) {
      if (status != FitStatus.OK) {
        return false;
      }

      // Reject fits that are outside the bounds of the data
      if (params[Gaussian2DFunction.SIGNAL] < 0 || params[Gaussian2DFunction.X_POSITION] < 0
          || params[Gaussian2DFunction.Y_POSITION] < 0
          || params[Gaussian2DFunction.X_POSITION] > size
          || params[Gaussian2DFunction.Y_POSITION] > size) {
        return false;
      }

      // Q. Should we do width bounds checking?
      if (fitConfig.isXSDFitting()) {
        if (params[Gaussian2DFunction.X_SD] < lb[Gaussian2DFunction.X_SD]
            || params[Gaussian2DFunction.X_SD] > ub[Gaussian2DFunction.X_SD]) {
          return false;
        }
      }
      if (fitConfig.isYSDFitting()) {
        if (params[Gaussian2DFunction.Y_SD] < lb[Gaussian2DFunction.Y_SD]
            || params[Gaussian2DFunction.Y_SD] > ub[Gaussian2DFunction.Y_SD]) {
          return false;
        }
      }

      return true;
    }
  }

  /**
   * Add the results to the statistics.
   *
   * @param stats the stats
   * @param answer the answer
   * @param photons the photons
   * @param sa the sa
   * @param time the time
   * @param result the result
   * @param c Count of the number of results
   */
  private static void addResults(Statistics[] stats, double[] answer, double photons, double sa,
      long[] time, double[][] result, int c) {
    // Store the results from each run
    for (int i = 0; i < c; i++) {
      addResult(stats, answer, photons, sa, result[i], time[i]);
    }
  }

  /**
   * Add the given results to the statistics.
   *
   * @param stats the stats
   * @param answer the answer
   * @param photons the photons
   * @param sa the sa
   * @param result the result
   * @param time the time
   */
  private static void addResult(Statistics[] stats, double[] answer, double photons, double sa,
      double[] result, long time) {
    for (int j = 0; j < result.length; j++) {
      stats[j].add(result[j] - answer[j]);
    }
    stats[TIME].add(time);
    stats[ACTUAL_SIGNAL].add(result[Gaussian2DFunction.SIGNAL] - photons);
    stats[ADJUSTED_X_SD].add(result[Gaussian2DFunction.X_SD] - sa);
    stats[ADJUSTED_Y_SD].add(result[Gaussian2DFunction.Y_SD] - sa);
  }

  @Override
  public void run(String arg) {
    SMLMUsageTracker.recordPlugin(this.getClass(), arg);

    if ("analysis".equals(arg)) {
      if (benchmarkResults.isEmpty()) {
        IJ.error(TITLE, "No benchmark results in memory.\n \n" + TextUtils.wrap(
            "Run the Fit Benchmark Data plugin and results will be stored for comparison analysis.",
            60));
        return;
      }
      runAnalysis();
    } else {
      if (CreateData.benchmarkParameters == null) {
        IJ.error(TITLE,
            "No benchmark parameters in memory.\n \n" + TextUtils.wrap(
                "Run the " + CreateData.TITLE
                    + " plugin in benchmark mode with a fixed number of photons per localisation.",
                60));
        return;
      }
      benchmarkParameters = CreateData.benchmarkParameters;
      imp = CreateData.getImage();
      if (imp == null || imp.getStackSize() != benchmarkParameters.frames) {
        IJ.error(TITLE, "No benchmark image to match the parameters in memory");
        return;
      }

      if (!showDialog()) {
        return;
      }

      run();
    }
  }

  private boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    final double sa = getSa();
    gd.addMessage(String.format(
        "Fits the benchmark image created by CreateData plugin.\nPSF width = %s, adjusted = %s",
        MathUtils.rounded(benchmarkParameters.s / benchmarkParameters.a), MathUtils.rounded(sa)));

    final FitEngineConfiguration config = SettingsManager.readFitEngineConfiguration(0);
    fitConfig = config.getFitConfiguration();
    fitConfig.setNmPerPixel(benchmarkParameters.a);

    // For each new benchmark width, reset the PSF width to the square pixel adjustment
    if (lastId != benchmarkParameters.id) {
      lastId = benchmarkParameters.id;
      fitConfig.setInitialPeakStdDev(benchmarkParameters.s / benchmarkParameters.a);
      // The adjusted width is only relevant when using a single point approximation
      // for a Gaussian over the pixel. Using the ERF function computes the actual
      // integral over the pixel.
      // fitConfig.setInitialPeakStdDev(sa);

      // Set the PSF. This requires the CreateData plugin to store the most appropriate
      // PSF used for the simulation.
      fitConfig.setPSF(benchmarkParameters.psf);
    }

    gd.addSlider("Region_size", 2, 20, regionSize);
    PeakFit.addPSFOptions(gd, fitConfig);
    gd.addChoice("Fit_solver", SettingsManager.getFitSolverNames(),
        fitConfig.getFitSolver().ordinal());
    gd.addChoice("Origin_XY", ORIGIN_XY, originXY, new OptionListener<Integer>() {
      @Override
      public boolean collectOptions(Integer value) {
        originXY = value;
        return collectOptions(false);
      }

      @Override
      public boolean collectOptions() {
        return collectOptions(true);
      }

      private boolean collectOptions(boolean silent) {
        if (originXY != 2) {
          return false;
        }
        final ExtendedGenericDialog egd = new ExtendedGenericDialog("Origin XY");
        egd.addNumericField("Offset_X", offsetX, 2, 6, "nm");
        egd.addNumericField("Offset_Y", offsetY, 2, 6, "nm");
        egd.setSilent(silent);
        egd.showDialog(true, gd);
        if (egd.wasCanceled()) {
          return false;
        }
        offsetX = gd.getNextNumber();
        offsetY = gd.getNextNumber();
        return true;
      }
    });
    gd.addChoice("Origin_Z", ORIGIN_Z, originZ, new OptionListener<Integer>() {
      @Override
      public boolean collectOptions(Integer value) {
        originZ = value;
        return collectOptions(false);
      }

      @Override
      public boolean collectOptions() {
        return collectOptions(true);
      }

      private boolean collectOptions(boolean silent) {
        if (originZ != 2) {
          return false;
        }
        final ExtendedGenericDialog egd = new ExtendedGenericDialog("Origin Z");
        egd.addNumericField("Offset_Z", offsetZ, 2, 6, "nm");
        egd.setSilent(silent);
        egd.showDialog(true, gd);
        if (egd.wasCanceled()) {
          return false;
        }
        offsetZ = gd.getNextNumber();
        return true;
      }
    });

    gd.addCheckbox("Zero_offset", zeroOffset);
    gd.addNumericField("Offset_points", offsetPoints, 0, new OptionListener<Double>() {
      @Override
      public boolean collectOptions(Double value) {
        offsetPoints = Math.max(0, value);
        return collectOptions(false);
      }

      @Override
      public boolean collectOptions() {
        return collectOptions(true);
      }

      private boolean collectOptions(boolean silent) {
        if (offsetPoints == 0) {
          return false;
        }
        final ExtendedGenericDialog egd = new ExtendedGenericDialog("Offset range");
        egd.addSlider("Offset_range_x", 0, 2.5, offsetRangeX);
        egd.addSlider("Offset_range_y", 0, 2.5, offsetRangeY);
        egd.addSlider("Offset_range_z", 0, 2.5, offsetRangeZ);
        egd.setSilent(silent);
        egd.showDialog(true, gd);
        if (egd.wasCanceled()) {
          return false;
        }
        offsetRangeX = Math.max(0, egd.getNextNumber());
        offsetRangeY = Math.max(0, egd.getNextNumber());
        offsetRangeZ = Math.max(0, egd.getNextNumber());
        return true;
      }
    });
    gd.addCheckbox("Background_fitting", backgroundFitting, new OptionListener<Boolean>() {
      @Override
      public boolean collectOptions(Boolean value) {
        backgroundFitting = value;
        return collectOptions(false);
      }

      @Override
      public boolean collectOptions() {
        return collectOptions(true);
      }

      private boolean collectOptions(boolean silent) {
        if (!backgroundFitting) {
          return false;
        }
        final ExtendedGenericDialog egd = new ExtendedGenericDialog("Background fitting");
        egd.addCheckbox("Estimate_background", estimateBackground);
        egd.setSilent(silent);
        egd.showDialog(true, gd);
        if (egd.wasCanceled()) {
          return false;
        }
        estimateBackground = egd.getNextBoolean();
        return true;
      }
    });
    gd.addMessage("Signal fitting can be disabled for "
        + PSFProtosHelper.getName(PSFType.ONE_AXIS_GAUSSIAN_2D) + " function");
    gd.addCheckbox("Signal_fitting", signalFitting, new OptionListener<Boolean>() {
      @Override
      public boolean collectOptions(Boolean value) {
        signalFitting = value;
        return collectOptions(false);
      }

      @Override
      public boolean collectOptions() {
        return collectOptions(true);
      }

      private boolean collectOptions(boolean silent) {
        if (!signalFitting) {
          return false;
        }
        final ExtendedGenericDialog egd = new ExtendedGenericDialog("Signal fitting");
        egd.addCheckbox("Estimate_signal", estimateSignal);
        egd.setSilent(silent);
        egd.showDialog(true, gd);
        if (egd.wasCanceled()) {
          return false;
        }
        estimateSignal = egd.getNextBoolean();
        return true;
      }
    });
    gd.addCheckbox("Show_histograms", showHistograms);
    gd.addCheckbox("Save_raw_data", saveRawData);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    regionSize = (int) Math.abs(gd.getNextNumber());
    fitConfig.setPSFType(PeakFit.getPSFTypeValues()[gd.getNextChoiceIndex()]);
    fitConfig.setFitSolver(gd.getNextChoiceIndex());
    originXY = gd.getNextChoiceIndex();
    originZ = gd.getNextChoiceIndex();
    zeroOffset = gd.getNextBoolean();
    offsetPoints = Math.max(0, gd.getNextNumber());
    backgroundFitting = gd.getNextBoolean();
    signalFitting = gd.getNextBoolean();
    showHistograms = gd.getNextBoolean();
    saveRawData = gd.getNextBoolean();

    gd.collectOptions();

    // Do this before the call to is3D()
    if (!PeakFit.configurePSFModel(config)) {
      return false;
    }

    getStartPoints(fitConfig.is3D());

    if (startPoints.length == 0) {
      IJ.error(TITLE, "No initial fitting positions");
      return false;
    }

    if (regionSize < 1) {
      regionSize = 1;
    }

    if (gd.invalidNumber()) {
      return false;
    }

    // Initialise the correct calibration
    final CalibrationWriter calibration = new CalibrationWriter(fitConfig.getCalibration());
    calibration.setNmPerPixel(benchmarkParameters.a);
    calibration.setCountPerPhoton(benchmarkParameters.gain);
    calibration.setQuantumEfficiency(benchmarkParameters.qe);
    calibration.setBias(benchmarkParameters.bias);
    calibration.setCameraType(benchmarkParameters.cameraType);
    calibration.setReadNoise(benchmarkParameters.readNoise);
    calibration.setExposureTime(1000);
    fitConfig.setCalibration(calibration.getCalibration());

    fitConfig.setCameraModelName(benchmarkParameters.cameraModelName);

    if (!PeakFit.configureFitSolver(config, IJImageSource.getBounds(imp), null, 0)) {
      return false;
    }

    if (showHistograms) {
      final ExtendedGenericDialog gd2 = new ExtendedGenericDialog(TITLE);
      gd2.addMessage("Select the histograms to display");
      gd2.addNumericField("Histogram_bins", histogramBins, 0);

      final double[] convert = getConversionFactors();

      for (int i = 0; i < displayHistograms.length; i++) {
        if (convert[i] != 0) {
          gd2.addCheckbox(NAMES[i].replace(' ', '_'), displayHistograms[i]);
        }
      }
      gd2.showDialog();
      if (gd2.wasCanceled()) {
        return false;
      }
      histogramBins = (int) Math.abs(gd2.getNextNumber());
      for (int i = 0; i < displayHistograms.length; i++) {
        if (convert[i] != 0) {
          displayHistograms[i] = gd2.getNextBoolean();
        }
      }
    }

    return true;
  }

  private double getSa() {
    final double sa =
        PSFCalculator.squarePixelAdjustment(benchmarkParameters.s, benchmarkParameters.a)
            / benchmarkParameters.a;
    return sa;
  }

  private int progress, stepProgress, totalProgress;

  /**
   * Show progress.
   */
  private synchronized void showProgress() {
    if (progress % stepProgress == 0) {
      if (ImageJUtils.showStatus("Frame: " + progress + " / " + totalProgress)) {
        IJ.showProgress(progress, totalProgress);
      }
    }
    progress++;
  }

  private void run() {
    // Initialise the answer.
    answer[Gaussian2DFunction.BACKGROUND] = benchmarkParameters.getBackground();
    answer[Gaussian2DFunction.SIGNAL] = benchmarkParameters.getSignal();
    answer[Gaussian2DFunction.X_POSITION] = benchmarkParameters.x;
    answer[Gaussian2DFunction.Y_POSITION] = benchmarkParameters.y;
    answer[Gaussian2DFunction.Z_POSITION] = benchmarkParameters.z;
    answer[Gaussian2DFunction.X_SD] = benchmarkParameters.s / benchmarkParameters.a;
    answer[Gaussian2DFunction.Y_SD] = benchmarkParameters.s / benchmarkParameters.a;

    // Set up the fit region. Always round down since 0.5 is the centre of the pixel.
    int x = (int) benchmarkParameters.x;
    int y = (int) benchmarkParameters.y;
    region = new Rectangle(x - regionSize, y - regionSize, 2 * regionSize + 1, 2 * regionSize + 1);
    if (!new Rectangle(0, 0, imp.getWidth(), imp.getHeight()).contains(region)) {
      // Check if it is incorrect by only 1 pixel
      if (region.width <= imp.getWidth() + 1 && region.height <= imp.getHeight() + 1) {
        ImageJUtils.log("Adjusting region %s to fit within image bounds (%dx%d)", region.toString(),
            imp.getWidth(), imp.getHeight());
        region = new Rectangle(0, 0, imp.getWidth(), imp.getHeight());
      } else {
        IJ.error(TITLE, "Fit region does not fit within the image");
        return;
      }
    }

    // Adjust the centre & account for 0.5 pixel offset during fitting
    x -= region.x;
    y -= region.y;
    answer[Gaussian2DFunction.X_POSITION] -= (region.x + 0.5);
    answer[Gaussian2DFunction.Y_POSITION] -= (region.y + 0.5);

    // Configure for fitting
    fitConfig.setBackgroundFitting(backgroundFitting);
    fitConfig.setNotSignalFitting(!signalFitting);
    fitConfig.setComputeDeviations(false);

    // Create the camera model
    CameraModel cameraModel = fitConfig.getCameraModel();
    // Crop for speed. Reset origin first so the region is within the model
    cameraModel.setOrigin(0, 0);
    cameraModel = cameraModel.crop(region, false);

    final ImageStack stack = imp.getImageStack();

    // Create a pool of workers
    final int nThreads = Prefs.getThreads();
    final BlockingQueue<Integer> jobs = new ArrayBlockingQueue<>(nThreads * 2);
    final List<Worker> workers = new LinkedList<>();
    final List<Thread> threads = new LinkedList<>();
    for (int i = 0; i < nThreads; i++) {
      final Worker worker = new Worker(jobs, stack, region, fitConfig, cameraModel);
      final Thread t = new Thread(worker);
      workers.add(worker);
      threads.add(t);
      t.start();
    }

    final int totalFrames = benchmarkParameters.frames;

    // Store all the fitting results
    results = new double[totalFrames * startPoints.length][];
    resultsTime = new long[results.length];

    // Fit the frames
    totalProgress = totalFrames;
    stepProgress = ImageJUtils.getProgressInterval(totalProgress);
    progress = 0;
    for (int i = 0; i < totalFrames; i++) {
      // Only fit if there were simulated photons
      if (benchmarkParameters.p[i] > 0) {
        put(jobs, i);
      }
    }
    // Finish all the worker threads by passing in a null job
    for (int i = 0; i < threads.size(); i++) {
      put(jobs, -1);
    }

    // Wait for all to finish
    for (int i = 0; i < threads.size(); i++) {
      try {
        threads.get(i).join();
      } catch (final InterruptedException ex) {
        ex.printStackTrace();
      }
    }
    threads.clear();

    if (hasOffsetXY()) {
      ImageJUtils.log(TITLE + ": CoM within start offset = %d / %d (%s%%)", comValid.intValue(),
          totalFrames, MathUtils.rounded((100.0 * comValid.intValue()) / totalFrames));
    }

    IJ.showProgress(1);
    IJ.showStatus("Collecting results ...");

    // Collect the results
    Statistics[] stats = null;
    for (int i = 0; i < workers.size(); i++) {
      final Statistics[] next = workers.get(i).stats;
      if (stats == null) {
        stats = next;
        continue;
      }
      for (int j = 0; j < next.length; j++) {
        stats[j].add(next[j]);
      }
    }
    workers.clear();

    // Show a table of the results
    summariseResults(stats, cameraModel);

    // Optionally show histograms
    if (showHistograms && stats != null) {
      IJ.showStatus("Calculating histograms ...");

      final WindowOrganiser windowOrganiser = new WindowOrganiser();
      final double[] convert = getConversionFactors();

      final HistogramPlotBuilder builder =
          new HistogramPlotBuilder(TITLE).setNumberOfBins(histogramBins);
      for (int i = 0; i < NAMES.length; i++) {
        if (displayHistograms[i] && convert[i] != 0) {
          // We will have to convert the values...
          final double[] tmp = ((StoredDataStatistics) stats[i]).getValues();
          for (int j = 0; j < tmp.length; j++) {
            tmp[j] *= convert[i];
          }
          final StoredDataStatistics tmpStats = StoredDataStatistics.create(tmp);
          builder.setData(tmpStats).setName(NAMES[i])
              .setPlotLabel(String.format("%s +/- %s", MathUtils.rounded(tmpStats.getMean()),
                  MathUtils.rounded(tmpStats.getStandardDeviation())))
              .show(windowOrganiser);
        }
      }

      windowOrganiser.tile();
    }

    if (saveRawData) {
      final String dir = ImageJUtils.getDirectory("Data_directory", rawDataDirectory);
      if (dir != null) {
        saveData(stats, dir);
      }
    }

    IJ.showStatus("");
  }

  private static void saveData(Statistics[] stats, String dir) {
    rawDataDirectory = dir;
    for (int i = 0; i < NAMES.length; i++) {
      saveStatistics((StoredDataStatistics) stats[i], NAMES[i]);
    }
  }

  private static void saveStatistics(StoredDataStatistics stats, String title) {
    final String filename = rawDataDirectory + title.replace(" ", "_") + ".txt";

    try (BufferedWriter out =
        new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename), "UTF-8"))) {
      // out.write(title);
      // out.newLine();
      final double[] data = stats.getValues();
      Arrays.sort(data);
      for (final double d : data) {
        // out.write(MathUtils.rounded(d, 4)); // rounded
        out.write(Double.toString(d));
        out.newLine();
      }
    } catch (final Exception ex) {
      ex.printStackTrace();
    }
  }

  private static void put(BlockingQueue<Integer> jobs, int i) {
    try {
      jobs.put(i);
    } catch (final InterruptedException ex) {
      throw new RuntimeException("Unexpected interruption", ex);
    }
  }

  private double[][] startPoints = null;

  /**
   * Gets the start points.
   *
   * @param is3D Set to true if 3D
   * @return The starting points for the fitting
   */
  private double[][] getStartPoints(boolean is3D) {
    if (startPoints != null) {
      return startPoints;
    }

    final TurboList<double[]> list = new TurboList<>();

    // Set up origin with an offset
    final double[] origin = new double[3];
    switch (originXY) {
      case 2:
        // Offset from the origin in pixels
        origin[0] = offsetX / benchmarkParameters.a;
        origin[1] = offsetY / benchmarkParameters.a;
        break;
      case 0:
      case 1:
        // No offset
        break;
      default:
        throw new IllegalStateException();
    }
    switch (originZ) {
      case 2:
        // Offset from the origin in pixels
        origin[2] = offsetZ / benchmarkParameters.a;
        break;
      case 0:
      case 1:
        // No offset
        break;
      default:
        throw new IllegalStateException();
    }

    if (zeroOffset) {
      list.add(origin.clone());
    }

    if (offsetPoints > 0 && ((offsetRangeX > 0 || offsetRangeY > 0) || is3D && offsetRangeZ > 0)) {
      final double[] min = new double[] {-Math.max(0, offsetRangeX), -Math.max(0, offsetRangeY),
          -Math.max(0, offsetRangeZ)};
      final double[] range = new double[] {2 * min[0], 2 * min[1], 2 * min[2]};
      final HaltonSequenceGenerator halton = new HaltonSequenceGenerator((is3D) ? 3 : 2);
      for (int i = 0; i < offsetPoints; i++) {
        final double[] offset = origin.clone();
        final double[] v = halton.nextVector();
        for (int j = 0; j < v.length; j++) {
          offset[j] += v[j] * range[j] + min[j];
        }
        list.add(offset);
      }
    }

    startPoints = list.toArray(new double[list.size()][]);
    return startPoints;
  }

  private static boolean hasOffsetXY() {
    return offsetRangeX > 0 && offsetRangeY > 0;
  }

  /**
   * Sum the intensity above background to estimate the signal.
   *
   * @param data the data
   * @param b background
   * @return The signal
   */
  public static double getSignal(double[] data, double b) {
    double s = 0;
    for (final double d : data) {
      s += d;
    }
    // Subtract the background per pixel and ensure at least 1 photon in the signal
    return Math.max(1, s - b * data.length);
  }

  /**
   * Get the centre of mass of the data.
   *
   * @param data the data
   * @param maxx the maxx
   * @param maxy the maxy
   * @param com The centre-of-mass
   */
  public static void getCentreOfMass(double[] data, int maxx, int maxy, double[] com) {
    com[0] = com[1] = 0;
    double sum = 0;
    for (int y = 0, index = 0; y < maxy; y++) {
      double sum1 = 0;
      for (int x = 0; x < maxx; x++, index++) {
        final double value = data[index];
        sum1 += value;
        com[0] += x * value;
      }
      com[1] += y * sum1;
      sum += sum1;
    }

    for (int i = 2; i-- > 0;) {
      com[i] /= sum;
    }
  }

  private void summariseResults(Statistics[] stats, CameraModel cameraModel) {
    createTable();

    final StringBuilder sb = new StringBuilder();

    // Create the benchmark settings and the fitting settings
    sb.append(benchmarkParameters.getMolecules()).append('\t');
    sb.append(MathUtils.rounded(benchmarkParameters.getSignal())).append('\t');
    sb.append(MathUtils.rounded(benchmarkParameters.s)).append('\t');
    sb.append(MathUtils.rounded(benchmarkParameters.a)).append('\t');
    sb.append(MathUtils.rounded(getSa() * benchmarkParameters.a)).append('\t');
    // Report XY in nm from the pixel centre
    sb.append(MathUtils.rounded(distanceFromCentre(benchmarkParameters.x))).append('\t');
    sb.append(MathUtils.rounded(distanceFromCentre(benchmarkParameters.y))).append('\t');
    sb.append(MathUtils.rounded(benchmarkParameters.a * benchmarkParameters.z)).append('\t');

    final CameraType cameraType = benchmarkParameters.cameraType;
    if (cameraType == CameraType.SCMOS) {
      sb.append("sCMOS (").append(benchmarkParameters.cameraModelName).append(") ");
      final Rectangle bounds = benchmarkParameters.cameraBounds;
      final Rectangle cropBounds = cameraModel.getBounds();
      sb.append(" ").append(bounds.x + cropBounds.x).append(",").append(bounds.y + cropBounds.y);
      sb.append(" ").append(region.width).append("x").append(region.width);
    } else {
      sb.append(CalibrationProtosHelper.getName(cameraType));
      sb.append(" Gain=").append(benchmarkParameters.gain);
      sb.append(" B=").append(benchmarkParameters.bias);
    }
    sb.append('\t');

    sb.append(MathUtils.rounded(benchmarkParameters.getBackground())).append('\t');
    sb.append(MathUtils.rounded(benchmarkParameters.noise)).append('\t');

    sb.append(MathUtils.rounded(benchmarkParameters.getSignal() / benchmarkParameters.noise))
        .append('\t');
    sb.append(MathUtils.rounded(benchmarkParameters.precisionN)).append('\t');
    sb.append(MathUtils.rounded(benchmarkParameters.precisionX)).append('\t');
    sb.append(MathUtils.rounded(benchmarkParameters.precisionXML)).append('\t');
    sb.append(region.width).append("x");
    sb.append(region.height).append('\t');
    sb.append(MathUtils.rounded(fitConfig.getInitialPeakStdDev() * benchmarkParameters.a))
        .append('\t');
    sb.append(PSFProtosHelper.getName(fitConfig.getPSF().getPsfType()));
    if (fitConfig.isFixedPSF()) {
      // Only fixed fitting can ignore the signal
      if (!signalFitting) {
        sb.append("NS");
      }
    }
    if (!backgroundFitting) {
      sb.append("NB");
    }
    sb.append(":").append(PeakFit.getSolverName(fitConfig));
    if (fitConfig.isModelCameraMLE()) {
      sb.append(":Camera\t");

      // Add details of the noise model for the MLE
      final CalibrationReader r = new CalibrationReader(fitConfig.getCalibration());
      sb.append("EM=").append(r.isEMCCD());
      sb.append(":G=").append(r.getCountPerPhoton());
      sb.append(":N=").append(r.getReadNoise());
    } else {
      sb.append('\t');
    }

    // Convert to units of the image (ADUs and pixels)
    final double[] convert = getConversionFactors();

    // Store the results for fitting on this benchmark dataset
    final BenchmarkResult benchmarkResult = new BenchmarkResult(benchmarkParameters, answer,
        sb.toString(), convert, this.results, this.resultsTime);
    if (!benchmarkResults.isEmpty()) {
      // Clear the results if the benchmark has changed
      if (benchmarkResults.getFirst().benchmarkParameters.id != benchmarkParameters.id) {
        benchmarkResults.clear();
      }
    }
    benchmarkResults.add(benchmarkResult);

    // Now output the actual results ...
    sb.append('\t');
    final double recall =
        (stats[0].getN() / (double) startPoints.length) / benchmarkParameters.getMolecules();
    sb.append(MathUtils.rounded(recall));

    for (int i = 0; i < stats.length; i++) {
      if (convert[i] != 0) {
        sb.append('\t').append(MathUtils.rounded(stats[i].getMean() * convert[i], 6)).append('\t')
            .append(MathUtils.rounded(stats[i].getStandardDeviation() * convert[i]));
      } else {
        sb.append("\t0\t0");
      }
    }
    summaryTable.append(sb.toString());
  }

  /**
   * Get the factors to convert the fitted units into calibrated photons and nm units. Set the
   * conversion to zero if the function does not fit the specified statistic.
   *
   * @return The conversion factors
   */
  private double[] getConversionFactors() {
    final double[] convert = new double[NAMES.length];
    convert[Gaussian2DFunction.BACKGROUND] = (fitConfig.isBackgroundFitting()) ? 1 : 0;
    convert[Gaussian2DFunction.SIGNAL] =
        (fitConfig.isNotSignalFitting() && fitConfig.isFixedPSF()) ? 0 : 1;
    convert[Gaussian2DFunction.X_POSITION] = benchmarkParameters.a;
    convert[Gaussian2DFunction.Y_POSITION] = benchmarkParameters.a;
    convert[Gaussian2DFunction.Z_POSITION] = (fitConfig.isZFitting()) ? benchmarkParameters.a : 0;
    convert[Gaussian2DFunction.X_SD] = (fitConfig.isXSDFitting()) ? benchmarkParameters.a : 0;
    convert[Gaussian2DFunction.Y_SD] = (fitConfig.isYSDFitting()) ? benchmarkParameters.a : 0;
    convert[Gaussian2DFunction.ANGLE] = (fitConfig.isAngleFitting()) ? 180.0 / Math.PI : 0;
    convert[TIME] = 1e-6;
    convert[ACTUAL_SIGNAL] = convert[Gaussian2DFunction.SIGNAL];
    convert[ADJUSTED_X_SD] = convert[Gaussian2DFunction.X_SD];
    convert[ADJUSTED_Y_SD] = convert[Gaussian2DFunction.Y_SD];
    return convert;
  }

  private double distanceFromCentre(double x) {
    // This assumes a
    x -= 0.5;
    final int i = (int) Math.round(x);
    x = x - i;
    return x * benchmarkParameters.a;
  }

  private static void createTable() {
    if (summaryTable == null || !summaryTable.isVisible()) {
      summaryTable = new TextWindow(TITLE, createHeader(false), "", 1000, 300);
      summaryTable.setVisible(true);
    }
  }

  private static void createAnalysisTable() {
    if (analysisTable == null || !analysisTable.isVisible()) {
      analysisTable =
          new TextWindow(TITLE + " Combined Analysis", createHeader(true), "", 1000, 300);
      analysisTable.setVisible(true);
    }
  }

  private static String createHeader(boolean extraRecall) {
    final StringBuilder sb = new StringBuilder(createParameterHeader() + "\tRecall");
    if (extraRecall) {
      sb.append("\tOrigRecall");
    }
    for (int i = 0; i < NAMES.length; i++) {
      sb.append('\t').append(NAMES[i]).append("\t+/-");
    }
    return sb.toString();
  }

  private static String createParameterHeader() {
    return "Molecules\tN\ts (nm)\ta (nm)\tsa (nm)\tX (nm)\tY (nm)\tZ (nm)\tCamera model\tB (photons)\tNoise (photons)\tSNR\tLimit N\tLimit X\tLimit X ML\tRegion\tWidth\tMethod\tOptions";
  }

  private void runAnalysis() {
    benchmarkParameters = benchmarkResults.getFirst().benchmarkParameters;
    final double sa = getSa();

    // The fitting could have used centre-of-mass or not making the number of points different.
    // Find the shortest array (this will be the one where the centre-of-mass was not used)
    int length = Integer.MAX_VALUE;
    for (final BenchmarkResult benchmarkResult : benchmarkResults) {
      if (length > benchmarkResult.results.length) {
        length = benchmarkResult.results.length;
      }
    }

    // Build a list of all the frames which have results
    final int[] valid = new int[length];
    int j = 0;
    final int[] count = new int[benchmarkResults.size()];
    for (final BenchmarkResult benchmarkResult : benchmarkResults) {
      int c = 0;
      for (int i = 0; i < valid.length; i++) {
        if (benchmarkResult.results[i] != null) {
          c++;
          valid[i]++;
        }
      }
      count[j++] = c;
    }

    final int target = benchmarkResults.size();

    // Check that we have data
    if (!validData(valid, target)) {
      IJ.error(TITLE, "No frames have fitting results from all methods");
      return;
    }

    // Get the number of start points valid for all the results
    final int totalFrames = benchmarkParameters.frames;
    final double numberOfStartPoints = length / totalFrames;

    createAnalysisTable();

    // Create the results using only frames where all the fitting methods were successful
    j = 0;
    for (final BenchmarkResult benchmarkResult : benchmarkResults) {
      final double[] answer = benchmarkResult.answer;

      final Statistics[] stats = new Statistics[NAMES.length];
      for (int i = 0; i < stats.length; i++) {
        stats[i] = new Statistics();
      }

      for (int i = 0; i < valid.length; i++) {
        if (valid[i] < target) {
          continue;
        }

        addResult(stats, answer, benchmarkParameters.p[i % totalFrames], sa,
            benchmarkResult.results[i], benchmarkResult.resultsTime[i]);
      }

      final StringBuilder sb = new StringBuilder(benchmarkResult.parameters);

      // Now output the actual results ...
      sb.append('\t');
      final double recall =
          (stats[0].getN() / numberOfStartPoints) / benchmarkParameters.getMolecules();
      sb.append(MathUtils.rounded(recall));
      // Add the original recall
      sb.append('\t');
      final double recall2 =
          (count[j++] / numberOfStartPoints) / benchmarkParameters.getMolecules();
      sb.append(MathUtils.rounded(recall2));

      // Convert to units of the image (ADUs and pixels)
      final double[] convert = benchmarkResult.convert;

      for (int i = 0; i < stats.length; i++) {
        if (convert[i] != 0) {
          sb.append('\t').append(MathUtils.rounded(stats[i].getMean() * convert[i], 6)).append('\t')
              .append(MathUtils.rounded(stats[i].getStandardDeviation() * convert[i]));
        } else {
          sb.append("\t0\t0");
        }
      }
      analysisTable.append(sb.toString());
    }
  }

  private static boolean validData(int[] valid, int target) {
    for (int i = 0; i < valid.length; i++) {
      if (valid[i] == target) {
        return true;
      }
    }
    return false;
  }
}
