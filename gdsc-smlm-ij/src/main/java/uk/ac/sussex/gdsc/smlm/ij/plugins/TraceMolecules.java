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
import ij.Prefs;
import ij.WindowManager;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.LutLoader;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.text.TextWindow;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import java.awt.Checkbox;
import java.awt.TextField;
import java.io.BufferedWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Consumer;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import uk.ac.sussex.gdsc.core.clustering.Cluster;
import uk.ac.sussex.gdsc.core.clustering.ClusterPoint;
import uk.ac.sussex.gdsc.core.clustering.ClusteringAlgorithm;
import uk.ac.sussex.gdsc.core.clustering.ClusteringEngine;
import uk.ac.sussex.gdsc.core.data.utils.Converter;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.SimpleImageJTrackProgress;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.utils.FileUtils;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationHelper;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.CalibrationOrBuilder;
import uk.ac.sussex.gdsc.smlm.data.config.UnitConverterUtils;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.TimeUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.ClusteringSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.DynamicMultipleTargetTracing;
import uk.ac.sussex.gdsc.smlm.results.DynamicMultipleTargetTracing.DmttConfiguration;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.TextFilePeakResults;
import uk.ac.sussex.gdsc.smlm.results.Trace;
import uk.ac.sussex.gdsc.smlm.results.TraceManager;
import uk.ac.sussex.gdsc.smlm.results.TraceManager.TraceMode;
import uk.ac.sussex.gdsc.smlm.results.count.Counter;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.PrecisionResultProcedure;

/**
 * Run a tracing algorithm on the peak results to trace molecules across the frames.
 */

public class TraceMolecules implements PlugIn {
  private static final double MIN_BLINKING_RATE = 1; // Should never be <= 0

  private static final AtomicReference<TextWindow> SUMMARY_TABLE = new AtomicReference<>();
  private static final AtomicBoolean LOG_HEADER = new AtomicBoolean(true);

  private static final String[] NAMES = new String[] {"Total Signal", "Signal/Frame", "Blinks",
      "t-On (s)", "t-Off (s)", "Total t-On (s)", "Total t-Off (s)", "Dwell time (s)"};
  private static final String[] FILENAMES = new String[] {"total_signal", "signal_per_frame",
      "blinks", "t_on", "t_off", "total_t_on", "total_t_off", "dwell_time"};

  private static final int TOTAL_SIGNAL = 0;
  private static final int SIGNAL_PER_FRAME = 1;
  private static final int BLINKS = 2;
  private static final int T_ON = 3;
  private static final int T_OFF = 4;
  private static final int TOTAL_T_ON = 5;
  private static final int TOTAL_T_OFF = 6;
  private static final int DWELL_TIME = 7;

  private static final boolean[] INTEGER_DISPLAY;

  static {
    INTEGER_DISPLAY = new boolean[NAMES.length];
    INTEGER_DISPLAY[BLINKS] = true;
  }

  private static final boolean[] ALWAYS_REMOVE_OUTLIERS;

  static {
    ALWAYS_REMOVE_OUTLIERS = new boolean[NAMES.length];
    ALWAYS_REMOVE_OUTLIERS[TOTAL_SIGNAL] = false;
  }

  private String pluginTitle = "Trace or Cluster Molecules";
  /** The output name set using either "Cluster" or "Trace" depending on the run mode. */
  private String outputName;

  /** The plugin settings. */
  private Settings pluginSettings;
  /** The clustering settings. */
  private ClusteringSettings.Builder settings;
  private MemoryPeakResults results;
  /** Store exposure time in seconds. */
  private double exposureTime;

  // Used for the plotting
  private double[] ddistanceThresholds;
  private int[] timeThresholds;
  private ArrayList<double[]> zeroCrossingPoints;
  private FloatProcessor fp;
  private Calibration cal;
  // Store the pixel value for the first plotted result
  private int origX;
  private int origY;
  private boolean debugMode;
  private boolean altKeyDown;
  private boolean optimiseBlinkingRate;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> INSTANCE = new AtomicReference<>(new Settings());

    private static final String KEY_INPUT_OPTION = "gdsc.smlm.tracemolecules.inputOption";
    private static final String KEY_INPUT_DEBUG_MODE = "gdsc.smlm.tracemolecules.inputDebugMode";
    private static final String KEY_INPUT_OPTIMISE =
        "gdsc.smlm.tracemolecules.inputOptimiseBlinkingRate";
    private static final String KEY_DISPLAY_HISTOGRAMS =
        "gdsc.smlm.tracemolecules.displayHistograms";
    private static final String KEY_FILENAME = "gdsc.smlm.tracemolecules.filename";
    private static final char ON = '1';
    private static final char OFF = '0';

    String inputOption;
    boolean inputDebugMode;
    boolean inputOptimiseBlinkingRate;
    boolean[] displayHistograms;
    String filename;

    Settings() {
      inputOption = Prefs.get(KEY_INPUT_OPTION, "");
      inputDebugMode = Prefs.getBoolean(KEY_INPUT_DEBUG_MODE, true);
      inputOptimiseBlinkingRate = Prefs.getBoolean(KEY_INPUT_OPTIMISE, false);
      displayHistograms = new boolean[NAMES.length];
      Arrays.fill(displayHistograms, true);
      final String display = Prefs.get(KEY_DISPLAY_HISTOGRAMS, "");
      for (int i = Math.min(display.length(), displayHistograms.length); i-- != 0;) {
        displayHistograms[i] = display.charAt(i) == ON;
      }
      filename = Prefs.get(KEY_FILENAME, "");
    }

    Settings(Settings source) {
      inputOption = source.inputOption;
      inputDebugMode = source.inputDebugMode;
      inputOptimiseBlinkingRate = source.inputOptimiseBlinkingRate;
      displayHistograms = source.displayHistograms.clone();
      filename = source.filename;
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
      return INSTANCE.get().copy();
    }

    /**
     * Save the settings.
     */
    void save() {
      INSTANCE.set(this);
      Prefs.set(KEY_INPUT_OPTION, inputOption);
      Prefs.set(KEY_INPUT_DEBUG_MODE, inputDebugMode);
      Prefs.set(KEY_INPUT_OPTIMISE, inputOptimiseBlinkingRate);
      final StringBuilder sb = new StringBuilder(displayHistograms.length);
      for (int i = 0; i < displayHistograms.length; i++) {
        sb.append(displayHistograms[i] ? ON : OFF);
      }
      Prefs.set(KEY_DISPLAY_HISTOGRAMS, sb.toString());
      Prefs.set(KEY_FILENAME, filename);
    }
  }

  private enum OptimiserPlot {
    //@formatter:off
    NONE{ @Override
    public String getName() { return "None"; }},
    NEAREST_NEIGHBOUR{ @Override
    public String getName() { return "Nearest neighbour"; }},
    BILINEAR{ @Override
    public String getName() { return "Bi-linear"; }};
    //@formatter:on

    @Override
    public String toString() {
      return getName();
    }

    /**
     * Gets the name.
     *
     * @return the name
     */
    public abstract String getName();

    public static OptimiserPlot get(int ordinal) {
      if (ordinal < 0 || ordinal >= values().length) {
        ordinal = 0;
      }
      return values()[ordinal];
    }
  }

  private static TraceMode getTraceMode(int traceMode) {
    if (traceMode < 0 || traceMode >= TraceMode.values().length) {
      return TraceMode.LATEST_FORERUNNER;
    }
    return TraceMode.values()[traceMode];
  }

  private static ClusteringAlgorithm getClusteringAlgorithm(int clusteringAlgorithm) {
    if (clusteringAlgorithm < 0 || clusteringAlgorithm >= ClusteringAlgorithm.values().length) {
      return ClusteringAlgorithm.PAIRWISE;
    }
    return ClusteringAlgorithm.values()[clusteringAlgorithm];
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(pluginTitle, "No localisations in memory");
      return;
    }
    altKeyDown = ImageJUtils.isExtraOptions();

    Trace[] traces = null;
    int totalFiltered = 0;
    if ("dynamic".equals(arg)) {
      // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      // Dynamic Mutliple Target Tracing
      // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      outputName = "Dynamic Trace";

      if (!showDynamicTraceDialog()) {
        return;
      }

      final DmttConfiguration config = createDmttConfiguration();
      traces =
          new DynamicMultipleTargetTracing(results).traceMolecules(config).toArray(new Trace[0]);
    } else if ("cluster".equals(arg)) {
      // -=-=-=-=-=
      // Clustering
      // -=-=-=-=-=
      outputName = "Cluster";

      if (!showClusterDialog()) {
        return;
      }

      final ClusteringEngine engine = new ClusteringEngine(Prefs.getThreads(),
          getClusteringAlgorithm(settings.getClusteringAlgorithm()),
          SimpleImageJTrackProgress.getInstance());

      if (settings.getSplitPulses()) {
        engine.setPulseInterval(settings.getPulseInterval());
        limitTimeThreshold(settings.getPulseInterval());
      }

      final List<Cluster> clusters = engine.findClusters(convertToClusterPoints(),
          getDistance(settings.getDistanceThreshold(), results.getCalibration()),
          timeThresholdInFrames());

      if (clusters == null) {
        ImageJUtils.log("Aborted");
        return;
      }

      traces = convertToTraces(clusters);
    } else {
      // -=-=-=-
      // Tracing
      // -=-=-=-
      outputName = "Trace";

      if (!showDialog()) {
        return;
      }

      final TraceManager manager = new TraceManager(results);
      manager.setTraceMode(getTraceMode(settings.getTraceMode()));
      manager.setActivationFrameInterval(settings.getPulseInterval());
      manager.setActivationFrameWindow(settings.getPulseWindow());
      manager.setDistanceExclusion(
          getDistance(settings.getDistanceExclusion(), results.getCalibration()));

      if (settings.getOptimise()) {
        // Optimise before configuring for a pulse interval
        runOptimiser(manager);
      }

      if (settings.getSplitPulses()) {
        manager.setPulseInterval(settings.getPulseInterval());
        limitTimeThreshold(settings.getPulseInterval());
      }

      manager.setTracker(SimpleImageJTrackProgress.getInstance());
      manager.traceMolecules(getDistance(settings.getDistanceThreshold(), results.getCalibration()),
          timeThresholdInFrames());
      traces = manager.getTraces();
      totalFiltered = manager.getTotalFiltered();
    }

    // --=-=-=-=-=-
    // Results processing
    // --=-=-=-=-=-

    outputName += (outputName.endsWith("e") ? "" : "e") + "d";
    saveResults(results, traces, outputName);

    // Save singles + single localisations in a trace
    saveCentroidResults(results, getSingles(traces), outputName + " Singles");
    final Trace[] multiTraces = getTraces(traces);
    saveResults(results, multiTraces, outputName + " Multi");

    // Save centroids
    outputName += " Centroids";
    final MemoryPeakResults tracedResults = saveCentroidResults(results, traces, outputName);

    // Save traces separately
    saveCentroidResults(results, multiTraces, outputName + " Multi");

    // Sort traces by time to assist the results source in extracting frames sequentially.
    // Do this before saving to assist in debugging using the saved traces file.
    sortByTime(traces);

    if (settings.getSaveTraces()) {
      saveTraces(traces);
    }

    summarise(createSummaryTable(), traces, totalFiltered, settings.getDistanceThreshold(),
        timeThresholdInSeconds());

    IJ.showStatus(String.format("%d localisations => %d traces (%d filtered)", results.size(),
        tracedResults.size(), totalFiltered));
  }

  private static double getDistance(double distanceThreshold, CalibrationOrBuilder calibration) {
    // Convert from NM to native units
    final Converter c = CalibrationHelper.getDistanceConverter(calibration, DistanceUnit.NM);
    return c.convertBack(distanceThreshold);
  }

  /**
   * Limit the time threshold to the pulse interval duration.
   *
   * @param pulseInterval the pulse interval
   */
  private void limitTimeThreshold(int pulseInterval) {
    // The pulse interval is in frames.
    // Convert the interval to the correct units.

    double limit;
    if (settings.getTimeUnit() == TimeUnit.FRAME) {
      limit = settings.getPulseInterval();
    } else {
      final TypeConverter<TimeUnit> convert = UnitConverterUtils.createConverter(TimeUnit.FRAME,
          settings.getTimeUnit(), getExposureTimeInMilliSeconds());
      limit = convert.convert(pulseInterval);
    }

    if (settings.getTimeThreshold() > limit) {
      settings.setTimeThreshold(limit);
    }
  }

  private List<ClusterPoint> convertToClusterPoints() {
    return convertToClusterPoints(results);
  }

  /**
   * Convert a list of peak results into points for the clustering engine.
   *
   * @param results the results
   * @return the list of clusters
   */
  public static List<ClusterPoint> convertToClusterPoints(MemoryPeakResults results) {
    final ArrayList<ClusterPoint> points = new ArrayList<>(results.size());
    final Counter counter = new Counter();
    results.forEach((PeakResultProcedure) result -> points.add(ClusterPoint.newTimeClusterPoint(
        counter.getAndIncrement(), result.getXPosition(), result.getYPosition(),
        result.getIntensity(), result.getFrame(), result.getEndFrame())));
    return points;
  }

  private Trace[] convertToTraces(List<Cluster> clusters) {
    return convertToTraces(results, clusters);
  }

  /**
   * Convert the clusters from the clustering engine into traces composed of the original list of
   * peak results.
   *
   * @param results the results
   * @param clusters the clusters
   * @return the traces
   */
  public static Trace[] convertToTraces(MemoryPeakResults results, List<Cluster> clusters) {
    final Trace[] traces = new Trace[clusters.size()];
    int index = 0;
    for (final Cluster cluster : clusters) {
      final Trace trace = new Trace();
      trace.setId(index + 1);
      for (ClusterPoint point = cluster.getHeadClusterPoint(); point != null;
          point = point.getNext()) {
        // The point Id was the position in the original results array
        trace.add(results.get(point.getId()));
      }
      trace.sort();
      traces[index++] = trace;
    }
    return traces;
  }

  /**
   * Sort traces by time.
   *
   * @param traces the traces
   */
  static void sortByTime(Trace[] traces) {
    for (final Trace t : traces) {
      t.sort();
    }
    Arrays.sort(traces,
        (o1, o2) -> Integer.compare(o1.getHead().getFrame(), o2.getHead().getFrame()));
  }

  /**
   * Convert the traces to results.
   *
   * @param sourceResults the source results
   * @param traces the traces
   * @param name the name
   * @return the memory peak results
   */
  static MemoryPeakResults saveResults(MemoryPeakResults sourceResults, Trace[] traces,
      String name) {
    final MemoryPeakResults tracedResults =
        TraceManager.convertToPeakResults(sourceResults, traces);
    tracedResults.setName(sourceResults.getName() + " " + name);
    MemoryPeakResults.addResults(tracedResults);
    return tracedResults;
  }

  private static MemoryPeakResults saveCentroidResults(MemoryPeakResults sourceResults,
      Trace[] traces, String name) {
    final MemoryPeakResults tracedResults =
        TraceManager.convertToCentroidPeakResults(sourceResults, traces);
    tracedResults.setName(sourceResults.getName() + " " + name);
    MemoryPeakResults.addResults(tracedResults);
    return tracedResults;
  }

  private static Trace[] getSingles(Trace[] traces) {
    final ArrayList<Trace> result = new ArrayList<>();
    for (final Trace t : traces) {
      if (t.size() == 1) {
        result.add(t);
      }
    }
    return result.toArray(new Trace[0]);
  }

  private static Trace[] getTraces(Trace[] traces) {
    final ArrayList<Trace> result = new ArrayList<>();
    for (final Trace t : traces) {
      if (t.size() != 1) {
        result.add(t);
      }
    }
    return result.toArray(new Trace[0]);
  }

  private void saveTraces(Trace[] traces) {
    pluginSettings.filename =
        saveTraces(results, traces, createSettingsComment(), pluginSettings.filename, 0);
  }

  /**
   * Save the traces to the file. A File open dialog is presented and the selected filename
   * returned.
   *
   * <p>If the id is above zero then the file open dialog title will have the id appended and the
   * filename is searched for .[0-9]+. and it is replaced with .id.
   *
   * @param sourceResults the source results
   * @param traces the traces
   * @param comment the comment
   * @param filename The initial filename
   * @param id The traces id (used if above 0)
   * @return The select filename (or null)
   */
  static String saveTraces(MemoryPeakResults sourceResults, Trace[] traces, String comment,
      String filename, int id) {
    IJ.showStatus("Saving traces");
    String title = "Traces_File";
    if (id > 0) {
      title += id;
      final String regex = "\\.[0-9]+\\.";
      if (filename.matches(regex)) {
        filename = filename.replaceAll(regex, "." + (id) + ".");
      } else {
        filename = FileUtils.replaceExtension(filename, id + ".xls");
      }
    }
    filename = ImageJUtils.getFilename(title, filename);
    if (filename != null) {
      filename = FileUtils.replaceExtension(filename, "xls");

      final boolean showDeviations = sourceResults.hasDeviations();
      // Assume that are results are from a single frame but store the trace ID
      final TextFilePeakResults traceResults =
          new TextFilePeakResults(filename, showDeviations, false, true);
      traceResults.copySettings(sourceResults);
      traceResults.begin();
      if (!traceResults.isActive()) {
        IJ.error("Failed to write to file: " + filename);
      } else {
        traceResults.addComment(comment);
        for (final Trace trace : traces) {
          traceResults.addTrace(trace);
        }
        traceResults.end();
      }
    }
    IJ.showStatus("");
    return filename;
  }

  private String createSettingsComment() {
    return String.format(
        "Molecule tracing : distance-threshold = %f nm : time-threshold = %f s (%d frames)",
        settings.getDistanceThreshold(), timeThresholdInSeconds(), timeThresholdInFrames());
  }

  private void summarise(Consumer<String> output, Trace[] traces, int filtered,
      double distanceThreshold, double timeThreshold) {
    IJ.showStatus("Calculating summary ...");

    final Statistics[] stats = new Statistics[NAMES.length];
    for (int i = 0; i < stats.length; i++) {
      stats[i] =
          (settings.getShowHistograms() || settings.getSaveTraceData()) ? new StoredDataStatistics()
              : new Statistics();
    }
    int singles = 0;
    for (final Trace trace : traces) {
      final int nBlinks = trace.getBlinks() - 1;
      stats[BLINKS].add(nBlinks);
      final int[] onTimes = trace.getOnTimes();
      final int[] offTimes = trace.getOffTimes();
      double timeOn = 0;
      for (final int t : onTimes) {
        stats[T_ON].add(t * exposureTime);
        timeOn += t * exposureTime;
      }
      stats[TOTAL_T_ON].add(timeOn);
      if (offTimes != null) {
        double timeOff = 0;
        for (final int t : offTimes) {
          stats[T_OFF].add(t * exposureTime);
          timeOff += t * exposureTime;
        }
        stats[TOTAL_T_OFF].add(timeOff);
      }
      final double signal = trace.getSignal() / results.getGain();
      stats[TOTAL_SIGNAL].add(signal);
      stats[SIGNAL_PER_FRAME].add(signal / trace.size());
      stats[DWELL_TIME]
          .add((trace.getTail().getEndFrame() - trace.getHead().getFrame() + 1) * exposureTime);
      if (trace.size() == 1) {
        singles++;
      }
    }

    // Add to the summary table
    final StringBuilder sb = new StringBuilder();
    sb.append(results.getName()).append('\t');
    sb.append(
        outputName.equals("Cluster") ? getClusteringAlgorithm(settings.getClusteringAlgorithm())
            : getTraceMode(settings.getTraceMode()))
        .append('\t');
    sb.append(MathUtils.rounded(getExposureTimeInMilliSeconds(), 3)).append('\t');
    sb.append(MathUtils.rounded(distanceThreshold, 3)).append('\t');
    sb.append(MathUtils.rounded(timeThreshold, 3));
    if (settings.getSplitPulses()) {
      sb.append(" *");
    }
    sb.append('\t');
    sb.append(convertSecondsToFrames(timeThreshold)).append('\t');
    sb.append(traces.length).append('\t');
    sb.append(filtered).append('\t');
    sb.append(singles).append('\t');
    sb.append(traces.length - singles).append('\t');
    for (int i = 0; i < stats.length; i++) {
      sb.append(MathUtils.rounded(stats[i].getMean(), 3)).append('\t');
    }
    if (java.awt.GraphicsEnvironment.isHeadless()) {
      IJ.log(sb.toString());
      return;
    }
    output.accept(sb.toString());

    if (settings.getShowHistograms()) {
      IJ.showStatus("Calculating histograms ...");

      final WindowOrganiser windowOrganiser = new WindowOrganiser();
      final HistogramPlotBuilder builder =
          new HistogramPlotBuilder(pluginTitle).setNumberOfBins(settings.getHistogramBins());
      for (int i = 0; i < NAMES.length; i++) {
        if (pluginSettings.displayHistograms[i]) {
          builder.setData((StoredDataStatistics) stats[i]).setName(NAMES[i])
              .setIntegerBins(INTEGER_DISPLAY[i])
              .setRemoveOutliersOption(
                  (settings.getRemoveOutliers() || ALWAYS_REMOVE_OUTLIERS[i]) ? 2 : 0)
              .show(windowOrganiser);
        }
      }

      windowOrganiser.tile();
    }

    if (settings.getSaveTraceData()) {
      saveTraceData(stats);
    }

    IJ.showStatus("");
  }

  private Consumer<String> createSummaryTable() {
    if (java.awt.GraphicsEnvironment.isHeadless()) {
      if (LOG_HEADER.compareAndSet(true, false)) {
        IJ.log(createHeader());
      }
      return IJ::log;
    }
    TextWindow window = SUMMARY_TABLE.get();
    if (window == null || !window.isVisible()) {
      window = new TextWindow(pluginTitle + " Data Summary", createHeader(), "", 800, 300);
      window.setVisible(true);
      SUMMARY_TABLE.set(window);
    }
    return window::append;
  }

  private static String createHeader() {
    final StringBuilder sb = new StringBuilder(
        "Dataset\tAlgorithm\tExposure time (ms)\tD-threshold (nm)\tT-threshold (s)\t"
            + "(Frames)\tMolecules\tFiltered\tSingles\tClusters");
    for (int i = 0; i < NAMES.length; i++) {
      sb.append('\t').append(NAMES[i]);
    }
    return sb.toString();
  }

  private void saveTraceData(Statistics[] stats) {
    // Get the directory
    IJ.showStatus("Saving trace data");
    final String directory =
        ImageJUtils.getDirectory("Trace_data_directory", settings.getTraceDataDirectory());
    if (directory != null) {
      settings.setTraceDataDirectory(directory);
      writeSettings();
      for (int i = 0; i < NAMES.length; i++) {
        saveTraceData((StoredDataStatistics) stats[i], NAMES[i], FILENAMES[i]);
      }
    }
    IJ.showStatus("");
  }

  private void saveTraceData(StoredDataStatistics stats, String name, String fileSuffix) {
    try (BufferedWriter file = Files.newBufferedWriter(
        Paths.get(settings.getTraceDataDirectory(), pluginTitle + "." + fileSuffix + ".txt"))) {
      file.append(name);
      file.newLine();

      for (final double d : stats.getValues()) {
        file.append(MathUtils.rounded(d, 4));
        file.newLine();
      }
    } catch (final Exception ex) {
      // Q. Add better handling of errors?
      Logger.getLogger(getClass().getName()).log(Level.WARNING, ex,
          () -> "Failed to save trace data to results directory: "
              + settings.getTraceDataDirectory());
    }
  }

  private boolean showDialog() {
    pluginTitle = outputName + " Molecules";
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(pluginTitle);
    gd.addHelp(HelpUrls.getUrl("trace-molecules"));

    readSettings();

    ResultsManager.addInput(gd, pluginSettings.inputOption, InputSource.MEMORY);

    gd.addNumericField("Distance_Threshold", settings.getDistanceThreshold(), 2, 6, "nm");
    gd.addNumericField("Distance_Exclusion", settings.getDistanceExclusion(), 2, 6, "nm");
    gd.addNumericField("Time_Threshold", settings.getTimeThreshold(), 2);
    gd.addChoice("Time_unit", SettingsManager.getTimeUnitNames(), settings.getTimeUnit().ordinal());
    final String[] traceModes =
        SettingsManager.getNames((Object[]) TraceManager.TraceMode.values());
    gd.addChoice("Trace_mode", traceModes,
        traceModes[getTraceMode(settings.getTraceMode()).ordinal()]);
    gd.addNumericField("Pulse_interval", settings.getPulseInterval(), 0, 6, "Frames");
    gd.addNumericField("Pulse_window", settings.getPulseWindow(), 0, 6, "Frames");
    gd.addCheckbox("Split_pulses", settings.getSplitPulses());
    gd.addCheckbox("Optimise", settings.getOptimise());
    gd.addCheckbox("Save_traces", settings.getSaveTraces());
    gd.addCheckbox("Show_histograms", settings.getShowHistograms());
    gd.addCheckbox("Save_trace_data", settings.getSaveTraceData());
    if (altKeyDown) {
      gd.addCheckbox("Debug", pluginSettings.inputDebugMode);
    }

    gd.showDialog();

    if (gd.wasCanceled() || !readDialog(gd)) {
      return false;
    }

    // Load the results
    results = ResultsManager.loadInputResults(pluginSettings.inputOption, true, null, null);
    if (MemoryPeakResults.isEmpty(results)) {
      IJ.error(pluginTitle, "No results could be loaded");
      IJ.showStatus("");
      return false;
    }

    // Store exposure time in seconds
    exposureTime = results.getCalibrationReader().getExposureTime() / 1000;

    return true;
  }

  private boolean readDialog(ExtendedGenericDialog gd) {
    pluginSettings.inputOption = ResultsManager.getInputSource(gd);
    settings.setDistanceThreshold(gd.getNextNumber());
    settings.setDistanceExclusion(gd.getNextNumber());
    settings.setTimeThreshold(gd.getNextNumber());
    settings.setTimeUnit(SettingsManager.getTimeUnitValues()[gd.getNextChoiceIndex()]);
    settings.setTraceMode(gd.getNextChoiceIndex());
    settings.setPulseInterval((int) gd.getNextNumber());
    settings.setPulseWindow((int) gd.getNextNumber());
    settings.setSplitPulses(gd.getNextBoolean());
    settings.setOptimise(gd.getNextBoolean());
    settings.setSaveTraces(gd.getNextBoolean());
    settings.setShowHistograms(gd.getNextBoolean());
    settings.setSaveTraceData(gd.getNextBoolean());
    if (altKeyDown) {
      debugMode = pluginSettings.inputDebugMode = gd.getNextBoolean();
    }

    writeSettings();

    if (gd.invalidNumber() || !showHistogramsDialog()) {
      return false;
    }

    // Check arguments
    try {
      // Allow distance of 0 for exact clustering
      ParameterUtils.isPositive("Distance threshold", settings.getDistanceThreshold());
      ParameterUtils.isAboveZero("Time threshold", settings.getTimeThreshold());
      ParameterUtils.isPositive("Pulse interval", settings.getPulseInterval());
      ParameterUtils.isPositive("Pulse window", settings.getPulseWindow());
    } catch (final IllegalArgumentException ex) {
      IJ.error(pluginTitle, ex.getMessage());
      return false;
    }

    return true;
  }

  private boolean showHistogramsDialog() {
    if (settings.getShowHistograms()) {
      final ExtendedGenericDialog gd = new ExtendedGenericDialog(pluginTitle);
      gd.addMessage("Select the histograms to display");
      gd.addCheckbox("Remove_outliers", settings.getRemoveOutliers());
      gd.addNumericField("Histogram_bins", settings.getHistogramBins(), 0);
      for (int i = 0; i < pluginSettings.displayHistograms.length; i++) {
        gd.addCheckbox(NAMES[i].replace(' ', '_'), pluginSettings.displayHistograms[i]);
      }
      gd.showDialog();
      if (gd.wasCanceled()) {
        return false;
      }
      settings.setRemoveOutliers(gd.getNextBoolean());
      settings.setHistogramBins((int) Math.abs(gd.getNextNumber()));
      for (int i = 0; i < pluginSettings.displayHistograms.length; i++) {
        pluginSettings.displayHistograms[i] = gd.getNextBoolean();
      }
      writeSettings();
    }
    return true;
  }

  private boolean showClusterDialog() {
    pluginTitle = outputName + " Molecules";
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(pluginTitle);
    gd.addHelp(HelpUrls.getUrl("cluster-molecules"));

    readSettings();

    ResultsManager.addInput(gd, pluginSettings.inputOption, InputSource.MEMORY);

    gd.addNumericField("Distance_Threshold", settings.getDistanceThreshold(), 2, 6, "nm");
    gd.addNumericField("Time_Threshold", settings.getTimeThreshold(), 2);
    gd.addChoice("Time_unit", SettingsManager.getTimeUnitNames(), settings.getTimeUnit().ordinal());
    final String[] algorithm = SettingsManager.getNames((Object[]) ClusteringAlgorithm.values());
    gd.addChoice("Clustering_algorithm", algorithm,
        algorithm[getClusteringAlgorithm(settings.getClusteringAlgorithm()).ordinal()]);
    gd.addNumericField("Pulse_interval", settings.getPulseInterval(), 0, 6, "Frames");
    gd.addCheckbox("Split_pulses", settings.getSplitPulses());
    gd.addCheckbox("Save_clusters", settings.getSaveTraces());
    gd.addCheckbox("Show_histograms", settings.getShowHistograms());
    gd.addCheckbox("Save_cluster_data", settings.getSaveTraceData());
    if (altKeyDown) {
      gd.addCheckbox("Debug", pluginSettings.inputDebugMode);
    }

    gd.showDialog();

    if (gd.wasCanceled() || !readClusterDialog(gd)) {
      return false;
    }

    // Load the results
    results = ResultsManager.loadInputResults(pluginSettings.inputOption, true, null, null);
    if (MemoryPeakResults.isEmpty(results)) {
      IJ.error(pluginTitle, "No results could be loaded");
      IJ.showStatus("");
      return false;
    }

    // Store exposure time in seconds
    exposureTime = results.getCalibrationReader().getExposureTime() / 1000;

    return true;
  }

  private void readSettings() {
    settings = SettingsManager.readClusteringSettings(0).toBuilder();
    pluginSettings = Settings.load();
  }

  private void writeSettings() {
    SettingsManager.writeSettings(settings.build());
    pluginSettings.save();
  }

  private boolean readClusterDialog(ExtendedGenericDialog gd) {
    pluginSettings.inputOption = ResultsManager.getInputSource(gd);
    settings.setDistanceThreshold(gd.getNextNumber());
    settings.setTimeThreshold(gd.getNextNumber());
    settings.setTimeUnit(SettingsManager.getTimeUnitValues()[gd.getNextChoiceIndex()]);
    settings.setClusteringAlgorithm(gd.getNextChoiceIndex());
    settings.setPulseInterval((int) gd.getNextNumber());
    settings.setSplitPulses(gd.getNextBoolean());
    settings.setSaveTraces(gd.getNextBoolean());
    settings.setShowHistograms(gd.getNextBoolean());
    settings.setSaveTraceData(gd.getNextBoolean());
    if (altKeyDown) {
      debugMode = pluginSettings.inputDebugMode = gd.getNextBoolean();
    }

    writeSettings();

    if (gd.invalidNumber() || !showHistogramsDialog()) {
      return false;
    }

    // Check arguments
    try {
      // Allow distance of 0 for exact clustering
      ParameterUtils.isPositive("Distance threshold", settings.getDistanceThreshold());
      final ClusteringAlgorithm clusteringAlgorithm =
          getClusteringAlgorithm(settings.getClusteringAlgorithm());
      if (clusteringAlgorithm == ClusteringAlgorithm.CENTROID_LINKAGE_DISTANCE_PRIORITY
          || clusteringAlgorithm == ClusteringAlgorithm.CENTROID_LINKAGE_TIME_PRIORITY) {
        ParameterUtils.isAboveZero("Time threshold", settings.getTimeThreshold());
        ParameterUtils.isPositive("Pulse interval", settings.getPulseInterval());
      }
    } catch (final IllegalArgumentException ex) {
      IJ.error(pluginTitle, ex.getMessage());
      return false;
    }

    return true;
  }

  private boolean showDynamicTraceDialog() {
    pluginTitle = outputName + " Molecules";
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(pluginTitle);
    gd.addHelp(HelpUrls.getUrl("dynamic-trace-molecules"));

    readSettings();

    ResultsManager.addInput(gd, pluginSettings.inputOption, InputSource.MEMORY);

    // Dynamic Multiple Target Tracing
    final TextField tfD = gd.addAndGetNumericField("Diffusion_coefficient",
        settings.getDiffusionCoefficentMaximum(), 3, 6, "um^2/s");
    final TextField tfW =
        gd.addAndGetNumericField("Temporal_window", settings.getTemporalWindow(), 0, 6, "frames");
    final TextField tfLdw =
        gd.addAndGetNumericField("Local_diffusion_weight", settings.getLocalDiffusionWeight(), 2);
    final TextField tfOiw =
        gd.addAndGetNumericField("On_intensity_weight", settings.getOnIntensityWeight(), 2);
    final TextField tfDdf = gd.addAndGetNumericField("Disappearance_decay_factor",
        settings.getDisappearanceDecayFactor(), 0, 6, "frames");
    final TextField tfDt = gd.addAndGetNumericField("Disappearance_threshold",
        settings.getDisappearanceThreshold(), 0, 6, "frames");
    final Checkbox cbDld = gd.addAndGetCheckbox("Disable_local_diffusion_model",
        settings.getDisableLocalDiffusionModel());
    final Checkbox cbDim =
        gd.addAndGetCheckbox("Disable_intensity_model", settings.getDisableIntensityModel());
    // Allow reset to defaults
    gd.addAndGetButton("Defaults", e -> {
      final DmttConfiguration config = DmttConfiguration.newBuilder(1).build();
      tfD.setText(String.valueOf(settings.getDiffusionCoefficentMaximum()));
      tfW.setText(String.valueOf(config.getTemporalWindow()));
      tfLdw.setText(String.valueOf(config.getLocalDiffusionWeight()));
      tfOiw.setText(String.valueOf(config.getOnIntensityWeight()));
      tfDdf.setText(String.valueOf(config.getDisappearanceDecayFactor()));
      tfDt.setText(String.valueOf(config.getDisappearanceThreshold()));
      cbDld.setState(config.isDisableLocalDiffusionModel());
      cbDim.setState(config.isDisableIntensityModel());
    });
    gd.addCheckbox("Save_traces", settings.getSaveTraces());
    gd.addCheckbox("Show_histograms", settings.getShowHistograms());
    gd.addCheckbox("Save_trace_data", settings.getSaveTraceData());

    gd.showDialog();

    if (gd.wasCanceled() || !readDynamicTraceDialog(gd)) {
      return false;
    }

    // Load the results
    results = ResultsManager.loadInputResults(pluginSettings.inputOption, true, null, null);
    if (MemoryPeakResults.isEmpty(results)) {
      IJ.error(pluginTitle, "No results could be loaded");
      IJ.showStatus("");
      return false;
    }

    // Store exposure time in seconds
    exposureTime = results.getCalibrationReader().getExposureTime() / 1000;

    return true;
  }

  private boolean readDynamicTraceDialog(ExtendedGenericDialog gd) {
    pluginSettings.inputOption = ResultsManager.getInputSource(gd);
    settings.setDiffusionCoefficentMaximum(gd.getNextNumber());
    settings.setTemporalWindow((int) gd.getNextNumber());
    settings.setLocalDiffusionWeight(gd.getNextNumber());
    settings.setOnIntensityWeight(gd.getNextNumber());
    settings.setDisappearanceDecayFactor(gd.getNextNumber());
    settings.setDisappearanceThreshold((int) gd.getNextNumber());
    settings.setDisableLocalDiffusionModel(gd.getNextBoolean());
    settings.setDisableIntensityModel(gd.getNextBoolean());
    settings.setSaveTraces(gd.getNextBoolean());
    settings.setShowHistograms(gd.getNextBoolean());
    settings.setSaveTraceData(gd.getNextBoolean());

    writeSettings();

    if (gd.invalidNumber() || !showHistogramsDialog()) {
      return false;
    }

    // Check arguments
    try {
      // Validate using the Builder
      createDmttConfiguration();
    } catch (final IllegalArgumentException ex) {
      IJ.error(pluginTitle, ex.getMessage());
      return false;
    }

    return true;
  }

  /**
   * Creates the DMTT configuration from the clustering settings.
   *
   * @return the DMTTconfiguration
   */
  private DmttConfiguration createDmttConfiguration() {
    return DmttConfiguration.newBuilder(settings.getDiffusionCoefficentMaximum())
        .setTemporalWindow(settings.getTemporalWindow())
        .setLocalDiffusionWeight(settings.getLocalDiffusionWeight())
        .setOnIntensityWeight(settings.getOnIntensityWeight())
        .setDisappearanceDecayFactor(settings.getDisappearanceDecayFactor())
        .setDisappearanceThreshold(settings.getDisappearanceThreshold())
        .setDisableLocalDiffusionModel(settings.getDisableLocalDiffusionModel())
        .setDisableIntensityModel(settings.getDisableIntensityModel()).build();
  }

  private void runOptimiser(TraceManager manager) {
    // Get an estimate of the number of molecules without blinking
    final Statistics stats = new Statistics();
    final double nmPerPixel = this.results.getNmPerPixel();
    final PrecisionResultProcedure pp = new PrecisionResultProcedure(results);
    pp.getPrecision();
    stats.add(pp.precisions);
    // Use twice the precision to get the initial distance threshold

    // Use 2.5x sigma as per the PC-PALM protocol in Sengupta, et al (2013) Nature Protocols 8, 345
    final double dEstimate = stats.getMean() * 2.5 / nmPerPixel;
    final int traceCount = manager.traceMolecules(dEstimate, 1);

    if (!getParameters(traceCount, dEstimate)) {
      return;
    }

    // TODO - Convert the distance threshold to use nm instead of pixels?
    final List<double[]> results = runTracing(manager, settings.getMinDistanceThreshold(),
        settings.getMaxDistanceThreshold(), settings.getMinTimeThreshold(),
        settings.getMaxTimeThreshold(), settings.getOptimiserSteps());

    // Compute fractional difference from the true value:
    // Use blinking rate directly or the estimated number of molecules
    double reference;
    int statistic;
    if (optimiseBlinkingRate) {
      reference = settings.getBlinkingRate();
      statistic = 3;
      IJ.log(String.format("Estimating blinking rate: %.2f", reference));
    } else {
      reference = traceCount / settings.getBlinkingRate();
      statistic = 2;
      IJ.log(String.format("Estimating number of molecules: %d / %.2f = %.2f", traceCount,
          settings.getBlinkingRate(), reference));
    }

    for (final double[] result : results) {
      if (optimiseBlinkingRate) {
        result[2] = (reference - result[statistic]) / reference;
      } else {
        result[2] = (result[statistic] - reference) / reference;
      }
    }

    // Locate the optimal parameters with a fit of the zero contour
    final boolean found = findOptimalParameters(results);

    createPlotResults(results);

    if (!found) {
      return;
    }

    // Make fractional difference absolute so that lowest is best
    for (final double[] result : results) {
      result[2] = Math.abs(result[2]);
    }

    // Set the optimal thresholds using the lowest value
    double[] best = new double[] {0, 0, Double.MAX_VALUE};
    for (final double[] result : results) {
      if (best[2] > result[2]) {
        best = result;
      }
    }

    settings.setDistanceThreshold(best[0]);

    // The optimiser works using frames so convert back to the correct units
    final TypeConverter<TimeUnit> convert = UnitConverterUtils.createConverter(TimeUnit.FRAME,
        settings.getTimeUnit(), getExposureTimeInMilliSeconds());
    settings.setTimeThreshold(convert.convert(best[1]));

    IJ.log(String.format(
        "Optimal fractional difference @ D-threshold=%g nm, T-threshold=%f s (%d frames)",
        settings.getDistanceThreshold(), timeThresholdInSeconds(), timeThresholdInFrames()));
    writeSettings();
  }

  private boolean getParameters(int traceCount, double distance) {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(pluginTitle + " Optimiser");

    final String msg = String.format("Estimate %d molecules at d=%f, t=1", traceCount, distance);
    IJ.log(msg);
    gd.addMessage(msg);
    gd.addNumericField("Min_Distance_Threshold (px)", settings.getMinDistanceThreshold(), 2);
    gd.addNumericField("Max_Distance_Threshold (px)", settings.getMaxDistanceThreshold(), 2);
    gd.addNumericField("Min_Time_Threshold (frames)", settings.getMinTimeThreshold(), 0);
    gd.addNumericField("Max_Time_Threshold (frames)", settings.getMaxTimeThreshold(), 0);
    gd.addSlider("Steps", 1, 20, settings.getOptimiserSteps());
    gd.addNumericField("Blinking_rate", settings.getBlinkingRate(), 2);
    final String[] plotNames = SettingsManager.getNames((Object[]) OptimiserPlot.values());
    gd.addChoice("Plot", plotNames,
        plotNames[OptimiserPlot.get(settings.getOptimiserPlot()).ordinal()]);
    if (altKeyDown) {
      gd.addCheckbox("Optimise_blinking", pluginSettings.inputOptimiseBlinkingRate);
    }

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.setMinDistanceThreshold(gd.getNextNumber());
    settings.setMaxDistanceThreshold(gd.getNextNumber());
    settings.setMinTimeThreshold((int) gd.getNextNumber());
    settings.setMaxTimeThreshold((int) gd.getNextNumber());
    settings.setOptimiserSteps((int) gd.getNextNumber());
    settings.setBlinkingRate(gd.getNextNumber());
    settings.setOptimiserPlot(gd.getNextChoiceIndex());
    if (altKeyDown) {
      optimiseBlinkingRate = pluginSettings.inputOptimiseBlinkingRate = gd.getNextBoolean();
    }

    writeSettings();

    if (gd.invalidNumber()) {
      return false;
    }

    if (settings.getMinDistanceThreshold() < 0) {
      settings.setMinDistanceThreshold(0);
    }
    if (settings.getMaxDistanceThreshold() < settings.getMinDistanceThreshold()) {
      settings.setMaxDistanceThreshold(settings.getMinDistanceThreshold());
    }
    if (settings.getMinTimeThreshold() < 0) {
      settings.setMinTimeThreshold(0);
    }
    if (settings.getMaxTimeThreshold() < settings.getMinTimeThreshold()) {
      settings.setMaxTimeThreshold(settings.getMinTimeThreshold());
    }
    if (settings.getOptimiserSteps() < 0) {
      settings.setOptimiserSteps(1);
    }
    if (settings.getBlinkingRate() < MIN_BLINKING_RATE) {
      IJ.error(gd.getTitle(), "Blinking rate must be above " + MIN_BLINKING_RATE);
      return false;
    }

    if (settings.getMinDistanceThreshold() == settings.getMaxDistanceThreshold()
        && settings.getMinTimeThreshold() == settings.getMaxTimeThreshold()) {
      IJ.error(gd.getTitle(), "Nothing to optimise");
      return false;
    }

    writeSettings();

    return true;
  }

  private int convertSecondsToFrames(double timeInSeconds) {
    return (int) Math.round(timeInSeconds / exposureTime);
  }

  private int timeThresholdInFrames() {
    return (int) Math.round(timeThresholdIn(TimeUnit.FRAME));
  }

  private double timeThresholdInSeconds() {
    return timeThresholdIn(TimeUnit.SECOND);
  }

  private double timeThresholdIn(TimeUnit timeUnit) {
    return UnitConverterUtils
        .createConverter(settings.getTimeUnit(), timeUnit, getExposureTimeInMilliSeconds())
        .convert(settings.getTimeThreshold());
  }

  private double getExposureTimeInMilliSeconds() {
    return exposureTime * 1000;
  }

  /**
   * Runs the tracing algorithm using distances and time thresholds between min and max with the
   * configured number of steps. Steps are spaced using a logarithmic scale.
   *
   * <p>Returns a list of [distance,time,N traces]
   *
   * @param peakResults the peak results
   * @param minDistanceThreshold the min distance threshold
   * @param maxDistanceThreshold the max distance threshold
   * @param minTimeThreshold the min time threshold
   * @param maxTimeThreshold the max time threshold
   * @param optimiserSteps the optimiser steps
   * @return a list of [distance,time,N traces,blinking rate]
   */
  public List<double[]> runTracing(MemoryPeakResults peakResults, double minDistanceThreshold,
      double maxDistanceThreshold, int minTimeThreshold, int maxTimeThreshold, int optimiserSteps) {
    return runTracing(new TraceManager(peakResults), minDistanceThreshold, maxDistanceThreshold,
        minTimeThreshold, maxTimeThreshold, optimiserSteps);
  }

  /**
   * Runs the tracing algorithm using distances and time thresholds between min and max with the
   * configured number of steps. Steps are spaced using a logarithmic scale.
   *
   * <p>Returns a list of [distance,time,N traces]
   *
   * @param manager the manager
   * @param minDistanceThreshold the min distance threshold
   * @param maxDistanceThreshold the max distance threshold
   * @param minTimeThreshold the min time threshold
   * @param maxTimeThreshold the max time threshold
   * @param optimiserSteps the optimiser steps
   * @return a list of [distance,time,N traces,blinking rate]
   */
  public List<double[]> runTracing(TraceManager manager, double minDistanceThreshold,
      double maxDistanceThreshold, int minTimeThreshold, int maxTimeThreshold, int optimiserSteps) {
    ddistanceThresholds = getIntervals(minDistanceThreshold, maxDistanceThreshold, optimiserSteps);
    timeThresholds = convert(getIntervals(minTimeThreshold, maxTimeThreshold, optimiserSteps));

    final int total = ddistanceThresholds.length * timeThresholds.length;
    final ArrayList<double[]> resultsList = new ArrayList<>(total);

    IJ.showStatus("Optimising tracing (" + total + " steps) ...");

    if (debugMode) {
      IJ.log("Optimising tracing ...");
    }

    int step = 0;
    for (final double d : ddistanceThresholds) {
      for (final int t : timeThresholds) {
        IJ.showProgress(step++, total);
        final int n = manager.traceMolecules(d, t);
        resultsList.add(new double[] {d, t, n, getBlinkingRate(manager.getTraces())});
        if (debugMode) {
          summarise(IJ::log, manager.getTraces(), manager.getTotalFiltered(), d, t);
        }
      }
    }

    if (debugMode) {
      IJ.log("-=-=-=-");
    }

    IJ.showStatus("");
    IJ.showProgress(1.0);

    return resultsList;
  }

  private static double getBlinkingRate(Trace[] traces) {
    final SummaryStatistics stats = new SummaryStatistics();
    for (final Trace trace : traces) {
      stats.addValue(trace.getBlinks());
    }
    return stats.getMean();
  }

  private static double[] getIntervals(double min, double max, int optimiserSteps) {
    if (max < min) {
      final double tmp = max;
      max = min;
      min = tmp;
    }
    final double range = max - min;
    if (range == 0) {
      return new double[] {min};
    }

    final double[] values = new double[optimiserSteps + ((min != 0) ? 1 : 0)];
    int index = 0;
    if (min != 0) {
      values[index++] = min;
    }

    // Build a set of steps from min to max

    // Calculate a factor so that:
    // f^steps = range + 1
    // => factor^n is in bounds [1:range+1] when n <= steps
    final double f = Math.pow(range + 1, 1.0 / optimiserSteps);

    // Set the first increment, i.e. f^1
    double x = f;

    for (int i = 0; i < optimiserSteps; i++) {
      // Set the value starting from min.
      // This is equivalent to: values[i] = min + Math.pow(f, i+1) - 1
      // Note that the bounds is [1:range+1] and so 1 is subtracted
      values[index++] = min + x - 1;
      x *= f;
    }

    return values;
  }

  private static int[] convert(double[] intervals) {
    final IntOpenHashSet set = new IntOpenHashSet(intervals.length);
    for (final double d : intervals) {
      set.add((int) Math.round(d));
    }

    set.remove(0); // Do not allow zero

    final int[] values = set.toIntArray();
    Arrays.sort(values);
    return values;
  }

  /**
   * Find the contour that intersects zero on the fractional difference plot. Find the point on the
   * contour nearest the origin.
   *
   * @param results the results
   * @return true, if successful
   */
  private boolean findOptimalParameters(List<double[]> results) {
    // This method only works if there are many results and if the results
    // cover enough of the search space to go from above zero (i.e. not enough traces)
    // to below zero (i.e. too many traces)

    final int maxx = timeThresholds.length;
    final int maxy = ddistanceThresholds.length;

    // --------
    // Find zero crossings using linear interpolation
    zeroCrossingPoints = new ArrayList<>();
    // --------

    // Pass across all time points
    boolean noZeroCrossingAtT0 = false;
    boolean noZeroCrossingAtTn = false;
    for (int x = 0; x < maxx; x++) {
      // Find zero crossings on distance points
      final double[] data = new double[maxy];
      for (int y = 0; y < maxy; y++) {
        final int i = y * maxx + x;
        final double[] result = results.get(i);
        data[y] = result[2];
      }
      final double zeroCrossing = findZeroCrossing(data, ddistanceThresholds);
      if (zeroCrossing > 0) {
        zeroCrossingPoints.add(new double[] {timeThresholds[x], zeroCrossing});
      } else if (x == 0) {
        noZeroCrossingAtT0 = true;
      } else if (x == maxx - 1) {
        noZeroCrossingAtTn = true;
      }
    }

    // If there were not enough zero crossings then the ranges are wrong
    if (zeroCrossingPoints.size() < 3) {
      IJ.log(String.format("Very few zero crossings (%d). Increase the optimisation space",
          zeroCrossingPoints.size()));
      return false;
    }

    // --------
    // Use relative distances to find the zero crossing with the smallest distance from origin
    // and set this as the optimal parameters
    // --------
    double minD = Double.MAX_VALUE;
    final double maxTimeThresholdInFrames = settings.getMaxTimeThreshold();
    // The optimiser works using frames so convert back to the correct units
    final TypeConverter<TimeUnit> convert = UnitConverterUtils.createConverter(TimeUnit.FRAME,
        settings.getTimeUnit(), getExposureTimeInMilliSeconds());

    for (final double[] point : zeroCrossingPoints) {
      final double dx = point[0] / maxTimeThresholdInFrames;
      final double dy = point[1] / settings.getMaxDistanceThreshold();
      final double d = dx * dx + dy * dy;
      if (d < minD) {
        minD = d;
        settings.setDistanceThreshold(point[1]);
        settings.setTimeThreshold(convert.convert(point[0]));
      }
    }

    // --------
    // Add more points to make the plotted line look better when showing the plot.
    // --------

    // Pass across all distance points
    boolean noZeroCrossingAtD0 = false;
    boolean noZeroCrossingAtDn = false;
    final double[] tThresholdsD = SimpleArrayUtils.toDouble(timeThresholds);
    for (int y = 0; y < maxy; y++) {
      // Find zero crossings on time points
      final double[] data = new double[maxx];
      for (int x = 0; x < maxx; x++) {
        final int i = y * maxx + x;
        final double[] result = results.get(i);
        data[x] = result[2];
      }
      final double zeroCrossing = findZeroCrossing(data, tThresholdsD);
      if (zeroCrossing > 0) {
        zeroCrossingPoints.add(new double[] {zeroCrossing, ddistanceThresholds[y]});
      } else if (y == 0) {
        noZeroCrossingAtD0 = true;
      } else if (y == maxy - 1) {
        noZeroCrossingAtDn = true;
      }
    }

    sortPoints();

    // --------
    // Output a message suggesting if the limits should be updated.
    // --------
    final StringBuilder sb = new StringBuilder();
    boolean reduceTime = false;
    boolean reduceDistance = false;
    if (noZeroCrossingAtDn && settings.getMinTimeThreshold() > 1) {
      sb.append(" * No zero crossing at max distance\n");
      reduceTime = true;
    }
    if (noZeroCrossingAtTn && settings.getMinDistanceThreshold() > 0) {
      sb.append(" * No zero crossing at max time\n");
      reduceDistance = true;
    }
    if (!noZeroCrossingAtD0 && settings.getMinDistanceThreshold() > 0) {
      sb.append(" * Zero crossing at min distance\n");
      reduceDistance = true;
    }
    if (!noZeroCrossingAtT0 && settings.getMinTimeThreshold() > 1) {
      sb.append(" * Zero crossing at min time\n");
      reduceTime = true;
    }
    if (reduceTime) {
      sb.append(" => Reduce the min time threshold\n");
    }
    if (reduceDistance) {
      sb.append(" => Reduce the min distance threshold\n");
    }
    if (sb.length() > 0) {
      sb.insert(0, "\nWarning:\n");
      sb.append("\n");
      IJ.log(sb.toString());
    }

    // interpolateZeroCrossingPoints();

    return true;
  }

  private static double findZeroCrossing(double[] data, double[] axis) {
    if (data[0] < 0) {
      return -1;
    }
    for (int i = 1; i < data.length; i++) {
      if (data[i] < 0) {
        final double fraction = data[i - 1] / (data[i - 1] - data[i]);
        return fraction * axis[i] + (1 - fraction) * axis[i - 1];
      }
    }

    return -1;
  }

  private void sortPoints() {
    // Sort by x coord, then y
    Collections.sort(zeroCrossingPoints, (o1, o2) -> {
      if (o1[0] < o2[0]) {
        return -1;
      }
      if (o1[0] > o2[0]) {
        return 1;
      }
      if (o1[1] < o2[1]) {
        return -1;
      }
      if (o1[1] > o2[1]) {
        return 1;
      }
      return 0;
    });
  }

  @SuppressWarnings("unused")
  private void interpolateZeroCrossingPoints() {
    final double[] x = new double[zeroCrossingPoints.size()];
    final double[] y = new double[zeroCrossingPoints.size()];
    for (int i = 0; i < x.length; i++) {
      final double[] point = zeroCrossingPoints.get(i);
      x[i] = point[0];
      y[i] = point[1];
    }
    final PolynomialSplineFunction fx = new SplineInterpolator().interpolate(x, y);
    double minX = x[0];
    final double maxX = x[x.length - 1];
    final double xinc = (maxX - minX) / 50;
    for (minX = minX + xinc; minX < maxX; minX += xinc) {
      zeroCrossingPoints.add(new double[] {minX, fx.value(minX)});
    }
    sortPoints();
  }

  /**
   * Build an image using the values within the results to set X,Y and value.
   *
   * @param results the results
   */
  private void createPlotResults(List<double[]> results) {
    final int width = 400;
    final int height = 400;
    switch (OptimiserPlot.get(settings.getOptimiserPlot())) {
      case NONE:
        return;
      case BILINEAR:
        fp = createBilinearPlot(results, width, height);
        break;
      default:
        fp = createNnPlot(results, width, height);
    }

    // Create a calibration to map the pixel position back to distance/time
    cal = new Calibration();
    final double xRange =
        getRange(settings.getMaxTimeThreshold(), settings.getMinTimeThreshold(), origX, width);
    final double yRange = getRange(settings.getMaxDistanceThreshold(),
        settings.getMinDistanceThreshold(), origY, height);
    cal.pixelWidth = xRange / width;
    cal.pixelHeight = yRange / height;
    cal.xOrigin = origX - settings.getMinTimeThreshold() / cal.pixelWidth;
    cal.yOrigin = origY - settings.getMinDistanceThreshold() / cal.pixelHeight;
    cal.setXUnit("frame");
    cal.setYUnit("pixel");

    showPlot();
  }

  /**
   * Shows the plot.
   */
  private void showPlot() {
    if (OptimiserPlot.get(settings.getOptimiserPlot()) == OptimiserPlot.NONE) {
      return;
    }

    // Display the image
    final String title = pluginTitle + ": | N - N_actual | / N_actual";
    ImagePlus imp = WindowManager.getImage(title);
    if (imp != null) {
      fp.setColorModel(imp.getProcessor().getColorModel());
      imp.setProcessor(fp);
    } else {
      imp = new ImagePlus(title, fp);
      imp.show();
      WindowManager.setTempCurrentImage(imp);
      final LutLoader lut = new LutLoader();
      lut.run("fire");
      WindowManager.setTempCurrentImage(null);
    }

    imp.setCalibration(cal);
    addZeroCrossingPoints(imp);
    imp.updateAndDraw();
  }

  private void addZeroCrossingPoints(ImagePlus imp) {
    PolygonRoi roi = null;
    imp.setRoi(roi);
    if (zeroCrossingPoints == null || zeroCrossingPoints.isEmpty()) {
      return;
    }
    final Calibration impCal = imp.getCalibration();
    final int npoints = zeroCrossingPoints.size();
    final float[] xpoints = new float[npoints];
    final float[] ypoints = new float[npoints];
    for (int i = 0; i < npoints; i++) {
      final double[] point = zeroCrossingPoints.get(i);
      // Convert to pixel coordinates.
      xpoints[i] = (float) (impCal.xOrigin + (point[0] / impCal.pixelWidth));
      ypoints[i] = (float) (impCal.yOrigin + (point[1] / impCal.pixelHeight));
    }
    roi = new PolygonRoi(xpoints, ypoints, npoints, Roi.POLYLINE);
    imp.setRoi(roi);
  }

  private FloatProcessor createNnPlot(List<double[]> results, int width, int height) {
    final FloatProcessor ip = new FloatProcessor(width, height);

    // Create lookup table that map the tested threshold values to a position in the image
    final int[] xLookup = createLookup(timeThresholds, settings.getMinTimeThreshold(), width);
    final int[] yLookup =
        createLookup(ddistanceThresholds, settings.getMinDistanceThreshold(), height);
    origX = (settings.getMinTimeThreshold() != 0) ? xLookup[1] : 0;
    origY = (settings.getMinDistanceThreshold() != 0) ? yLookup[1] : 0;

    final int gridWidth = timeThresholds.length;
    final int gridHeight = ddistanceThresholds.length;
    for (int y = 0, i = 0; y < gridHeight; y++) {
      for (int x = 0; x < gridWidth; x++, i++) {
        final int x1 = xLookup[x];
        final int x2 = xLookup[x + 1];
        final int y1 = yLookup[y];
        final int y2 = yLookup[y + 1];
        final double[] result = results.get(i);
        ip.setValue(Math.abs(result[2]));
        ip.setRoi(x1, y1, x2 - x1, y2 - y1);
        ip.fill();
      }
    }
    return ip;
  }

  private static int[] createLookup(int[] values, int min, int scale) {
    final double[] newValues = SimpleArrayUtils.toDouble(values);
    return createLookup(newValues, min, scale);
  }

  private static int[] createLookup(double[] values, double min, int scale) {
    // To allow the lowest result to be plotted, add space at the edge
    // equal to the next interval
    if (min != 0 && values.length > 1) {
      min -= values[1] - values[0];
    }

    final int[] lookup = new int[values.length + 1];
    final double range = values[values.length - 1] - min;
    final double scaleFactor = scale / range;
    for (int i = 1; i < values.length; i++) {
      lookup[i] = (int) Math.round(scaleFactor * (values[i - 1] - min));
    }
    lookup[values.length] = scale;
    return lookup;
  }

  private FloatProcessor createBilinearPlot(List<double[]> results, int width, int height) {
    final FloatProcessor ip = new FloatProcessor(width, height);

    // Create lookup table that map the tested threshold values to a position in the image
    final int[] xLookup = createLookup(timeThresholds, settings.getMinTimeThreshold(), width);
    final int[] yLookup =
        createLookup(ddistanceThresholds, settings.getMinDistanceThreshold(), height);
    origX = (settings.getMinTimeThreshold() != 0) ? xLookup[1] : 0;
    origY = (settings.getMinDistanceThreshold() != 0) ? yLookup[1] : 0;

    final int gridWidth = timeThresholds.length;
    final int gridHeight = ddistanceThresholds.length;
    for (int y = 0, prevY = 0; y < gridHeight; y++) {
      for (int x = 0, prevX = 0; x < gridWidth; x++) {
        // Get the 4 flanking values
        final double x1y1 = results.get(prevY * gridWidth + prevX)[2];
        final double x1y2 = results.get(y * gridWidth + prevX)[2];
        final double x2y1 = results.get(prevY * gridWidth + x)[2];
        final double x2y2 = results.get(y * gridWidth + x)[2];

        // Pixel range
        final int x1 = xLookup[x];
        final int x2 = xLookup[x + 1];
        final int y1 = yLookup[y];
        final int y2 = yLookup[y + 1];

        final double xRange = x2 - x1;
        final double yRange = y2 - y1;

        for (int yy = y1; yy < y2; yy++) {
          final double yFraction = (yy - y1) / yRange;
          for (int xx = x1; xx < x2; xx++) {
            // Interpolate
            final double xFraction = (xx - x1) / xRange;
            final double v1 = x1y1 * (1 - xFraction) + x2y1 * xFraction;
            final double v2 = x1y2 * (1 - xFraction) + x2y2 * xFraction;
            final double value = v1 * (1 - yFraction) + v2 * yFraction;
            ip.setf(xx, yy, (float) value);
          }
        }

        prevX = x;
      }
      prevY = y;
    }

    // Convert to absolute for easier visualisation
    final float[] data = (float[]) ip.getPixels();
    for (int i = 0; i < data.length; i++) {
      data[i] = Math.abs(data[i]);
    }

    return ip;
  }

  private static double getRange(double max, double min, int orig, int width) {
    double range = max - min;
    if (range <= 0) {
      range = 1;
    }
    return range * width / (width - orig);
  }
}
