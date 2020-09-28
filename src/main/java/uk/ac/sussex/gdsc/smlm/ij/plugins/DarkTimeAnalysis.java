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
import ij.Prefs;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;
import uk.ac.sussex.gdsc.core.clustering.Cluster;
import uk.ac.sussex.gdsc.core.clustering.ClusteringAlgorithm;
import uk.ac.sussex.gdsc.core.clustering.ClusteringEngine;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.SimpleImageJTrackProgress;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.logging.TrackProgress;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.StoredData;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.Trace;
import uk.ac.sussex.gdsc.smlm.results.TraceManager;

/**
 * Computes a graph of the dark time and estimates the time threshold for the specified point in the
 * cumulative histogram.
 */
public class DarkTimeAnalysis implements PlugIn {
  private static final String TITLE = "Dark-time Analysis";

  private static final String[] METHOD;
  private static final ClusteringAlgorithm[] algorithms =
      new ClusteringAlgorithm[] {ClusteringAlgorithm.CENTROID_LINKAGE_TIME_PRIORITY,
          ClusteringAlgorithm.CENTROID_LINKAGE_DISTANCE_PRIORITY,
          ClusteringAlgorithm.PARTICLE_CENTROID_LINKAGE_TIME_PRIORITY,
          ClusteringAlgorithm.PARTICLE_CENTROID_LINKAGE_DISTANCE_PRIORITY};

  static {
    final ArrayList<String> methods = new ArrayList<>();
    methods.add("Tracing");
    for (final ClusteringAlgorithm c : algorithms) {
      methods.add("Clustering (" + c.toString() + ")");
    }
    METHOD = methods.toArray(new String[methods.size()]);
  }

  private double msPerFrame;

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    String inputOption;
    int method;
    double searchDistance;
    double maxDarkTime;
    double percentile;
    int histogramBins;

    Settings() {
      // Set defaults
      inputOption = "";
      searchDistance = 100;
      percentile = 99;
    }

    Settings(Settings source) {
      inputOption = source.inputOption;
      method = source.method;
      searchDistance = source.searchDistance;
      maxDarkTime = source.maxDarkTime;
      percentile = source.percentile;
      histogramBins = source.histogramBins;
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

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    // Require some fit results and selected regions
    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "There are no fitting results in memory");
      return;
    }

    if (!showDialog()) {
      return;
    }

    // Assume pixels for now
    final MemoryPeakResults results =
        ResultsManager.loadInputResults(settings.inputOption, true, DistanceUnit.PIXEL);
    IJ.showStatus("");
    if (MemoryPeakResults.isEmpty(results)) {
      IJ.error(TITLE, "No results could be loaded");
      return;
    }
    if (!results.hasCalibration()) {
      IJ.error(TITLE, "Results are not calibrated");
      return;
    }
    msPerFrame = results.getCalibrationReader().getExposureTime();
    if (!(msPerFrame > 0)) {
      IJ.error(TITLE, "ms/frame must be strictly positive: " + msPerFrame);
      return;
    }
    ImageJUtils.log("%s: %d localisations", TITLE, results.size());

    analyse(results);
  }

  private boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("dark-time-analysis"));

    settings = Settings.load();

    gd.addMessage("Compute the cumulative dark-time histogram");
    ResultsManager.addInput(gd, settings.inputOption, InputSource.MEMORY);

    gd.addChoice("Method", METHOD, METHOD[settings.method]);
    gd.addSlider("Search_distance (nm)", 5, 150, settings.searchDistance);
    gd.addNumericField("Max_dark_time (seconds)", settings.maxDarkTime, 2);
    gd.addSlider("Percentile", 0, 100, settings.percentile);
    gd.addSlider("Histogram_bins", -1, 100, settings.histogramBins);
    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.inputOption = gd.getNextChoice();
    settings.method = gd.getNextChoiceIndex();
    settings.searchDistance = gd.getNextNumber();
    settings.maxDarkTime = gd.getNextNumber();
    settings.percentile = gd.getNextNumber();
    settings.histogramBins = (int) gd.getNextNumber();
    settings.save();

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Search distance", settings.searchDistance);
      ParameterUtils.isPositive("Percentile", settings.percentile);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return true;
  }

  private void analyse(MemoryPeakResults results) {
    // Find min and max time frames
    results.sort();
    final int min = results.getFirstFrame();
    final int max = results.getLastFrame();

    // Trace results:
    // TODO - The search distance could have units to avoid assuming the results are in pixels
    final double d = settings.searchDistance / results.getCalibrationReader().getNmPerPixel();
    int range = max - min + 1;
    if (settings.maxDarkTime > 0) {
      range = Math.max(1, (int) Math.round(settings.maxDarkTime * 1000 / msPerFrame));
    }

    final TrackProgress tracker = SimpleImageJTrackProgress.getInstance();
    tracker.status("Analysing ...");
    tracker.log("Analysing (d=%s nm (%s px) t=%s s (%d frames)) ...",
        MathUtils.rounded(settings.searchDistance), MathUtils.rounded(d),
        MathUtils.rounded(range * msPerFrame / 1000.0), range);

    Trace[] traces;
    if (settings.method == 0) {
      final TraceManager tm = new TraceManager(results);
      tm.setTracker(tracker);
      tm.traceMolecules(d, range);
      traces = tm.getTraces();
    } else {
      final ClusteringEngine engine =
          new ClusteringEngine(Prefs.getThreads(), algorithms[settings.method - 1], tracker);
      final List<Cluster> clusters =
          engine.findClusters(TraceMolecules.convertToClusterPoints(results), d, range);
      traces = TraceMolecules.convertToTraces(results, clusters);
    }

    tracker.status("Computing histogram ...");

    // Build dark-time histogram
    final int[] times = new int[range];
    final StoredData stats = new StoredData();
    for (final Trace trace : traces) {
      if (trace.getBlinks() > 1) {
        for (final int t : trace.getOffTimes()) {
          times[t]++;
        }
        stats.add(trace.getOffTimes());
      }
    }

    plotDarkTimeHistogram(stats);

    // Cumulative histogram
    for (int i = 1; i < times.length; i++) {
      times[i] += times[i - 1];
    }
    final int total = times[times.length - 1];

    // Plot dark-time up to 100%
    double[] x = new double[range];
    double[] y = new double[range];
    int truncate = 0;
    for (int i = 0; i < x.length; i++) {
      x[i] = i * msPerFrame;
      if (times[i] == total) {
        // Final value at 100%
        y[i] = 100.0;
        truncate = i + 1;
        break;
      }
      y[i] = (100.0 * times[i]) / total;
    }
    if (truncate > 0) {
      x = Arrays.copyOf(x, truncate);
      y = Arrays.copyOf(y, truncate);
    }

    final String title = "Cumulative Dark-time";
    final Plot plot = new Plot(title, "Time (ms)", "Percentile");
    plot.addPoints(x, y, Plot.LINE);
    ImageJUtils.display(title, plot);

    // Report percentile
    for (int i = 0; i < y.length; i++) {
      if (y[i] >= settings.percentile) {
        ImageJUtils.log("Dark-time Percentile %.1f @ %s ms = %s s", settings.percentile,
            MathUtils.rounded(x[i]), MathUtils.rounded(x[i] / 1000));
        break;
      }
    }

    tracker.status("");
  }

  private void plotDarkTimeHistogram(StoredData stats) {
    if (settings.histogramBins >= 0) {
      // Convert the X-axis to milliseconds
      final double[] xValues = stats.getValues();
      for (int i = 0; i < xValues.length; i++) {
        xValues[i] *= msPerFrame;
      }

      // Ensure the bin width is never less than 1
      new HistogramPlotBuilder("Dark-time", StoredDataStatistics.create(xValues), "Time (ms)")
          .setIntegerBins(true).setNumberOfBins(settings.histogramBins).show();
    }
  }
}
