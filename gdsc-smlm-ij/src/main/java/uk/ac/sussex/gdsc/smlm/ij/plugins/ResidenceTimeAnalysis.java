/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Package
 *
 * Software for single molecule localisation microscopy (SMLM) in ImageJ
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.text.TextWindow;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Consumer;
import java.util.function.DoubleSupplier;
import java.util.function.DoubleUnaryOperator;
import java.util.function.Function;
import java.util.function.IntToDoubleFunction;
import java.util.function.Predicate;
import java.util.function.ToDoubleFunction;
import java.util.function.ToIntFunction;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.numbers.rootfinder.BrentSolver;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.AliasMethodDiscreteSampler;
import org.apache.commons.rng.sampling.distribution.DiscreteSampler;
import org.apache.commons.rng.sampling.distribution.ZigguratSampler;
import org.apache.commons.rng.simple.RandomSource;
import uk.ac.sussex.gdsc.core.ij.ImageJPluginLoggerHelper;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.MultiDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.Seed;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.ValidationUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.fitting.ResidenceTimeFitting;
import uk.ac.sussex.gdsc.smlm.fitting.ResidenceTimeFitting.FitResult;
import uk.ac.sussex.gdsc.smlm.fitting.ResidenceTimeFitting.Model;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.Trace;
import uk.ac.sussex.gdsc.smlm.results.count.FrameCounter;
import uk.ac.sussex.gdsc.smlm.results.sort.IdFramePeakResultComparator;

/**
 * Perform analysis of residence times for stationary particles to identify the dissociation rate.
 */
public class ResidenceTimeAnalysis implements PlugIn {
  private static final String TITLE = "Residence Time Analysis";
  private static final String HELP_KEY = "residence-time-analysis";

  private static final AtomicReference<TextWindow> SUMMARY_TABLE_REF = new AtomicReference<>();

  // Used for the multiMode option
  private static final AtomicReference<List<String>> SELECTED_REF = new AtomicReference<>();

  /**
   * The ROI colors.
   *
   * <p>This is Paul Tol's Vibrant colour scheme (color-blind safe). It has been reordered for
   * contrast between the initial colors. See https://personal.sron.nl/~pault/.
   */
  private static final Color[] COLORS = {new Color(0, 119, 187), // Blue
      new Color(238, 119, 51), // Orange
      new Color(0, 153, 136), // Teal
      new Color(204, 51, 17), // Red
      new Color(51, 187, 238), // Cyan
      new Color(238, 51, 119) // Magenta
  };

  /** The plugin settings. */
  Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> INSTANCE = new AtomicReference<>(new Settings());

    int samples;
    double k0;
    double k1;
    double f0;
    double exposureTime;
    Seed seed;
    int firstFrame;
    int lastFrame;
    int minSize;
    double meanDistance;
    double maxDistance;
    double maxTraceLength;
    int minBinCount;
    double apparentDissociationRate;
    int bootstrapRepeats;
    Seed bootstrapSeed;
    int maxPopulations;
    boolean showQQplot;

    Settings() {
      // Set defaults
      samples = 5000;
      k0 = 0.2;
      k1 = 5;
      f0 = 0.15;
      exposureTime = 50;
      seed = Seed.from("1c0f4f6003fc7f01");
      firstFrame = 1;
      lastFrame = -1;
      minSize = 2;
      meanDistance = 100;
      maxDistance = 400;
      bootstrapRepeats = 100;
      bootstrapSeed = Seed.from("fdfcf7c90ee01f9e");
      maxPopulations = 2;
    }

    Settings(Settings source) {
      samples = source.samples;
      k0 = source.k0;
      k1 = source.k1;
      f0 = source.f0;
      exposureTime = source.exposureTime;
      seed = source.seed;
      firstFrame = source.firstFrame;
      lastFrame = source.lastFrame;
      minSize = source.minSize;
      meanDistance = source.meanDistance;
      maxDistance = source.maxDistance;
      maxTraceLength = source.maxTraceLength;
      minBinCount = source.minBinCount;
      apparentDissociationRate = source.apparentDissociationRate;
      bootstrapRepeats = source.bootstrapRepeats;
      bootstrapSeed = source.bootstrapSeed;
      showQQplot = source.showQQplot;
      maxPopulations = source.maxPopulations;
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
    }
  }

  /**
   * Interface for extracting a double value from an object and an integer.
   *
   * @param <T> the type of object
   */
  private interface ObjIntToDoubleFunction<T> {
    /**
     * Apply the function.
     *
     * @param o the object
     * @param i the integer
     * @return the double
     */
    double apply(T o, int i);
  }

  /**
   * Custom statistics class to collect mean, standard deviation and max.
   */
  private static class MyStatistics {
    private final Statistics stats = new Statistics();
    private double max;

    void add(double value) {
      stats.add(value);
      max = max > value ? max : value;
    }

    double mean() {
      return stats.getMean();
    }

    double sd() {
      return stats.getStandardDeviation();
    }

    double max() {
      return max;
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (ImageJUtils.isExtraOptions()) {
      runSimulation();
      return;
    }
    if (!MemoryPeakResults.isAnyInMemory(ResultsManager::hasId)) {
      IJ.error(TITLE, "No traced localisations in memory");
      return;
    }

    settings = Settings.load();
    // Saved by reference so just save now
    settings.save();

    final List<MemoryPeakResults> allResults = new ArrayList<>();

    // Option to pick multiple input datasets together using a list box.
    if (!showMultiDialog(allResults) || !showDialog()) {
      return;
    }

    // Here results are all calibrated with the same exposure time and pixel pitch.
    // Convert results to residence times with options to exclude non-stationary objects.
    int[] counts = extractResidenceTimes(allResults);
    final double exposureTime = allResults.get(0).getCalibrationReader().getExposureTime();
    final int originalLength = counts.length;
    counts = filterCounts(counts, exposureTime);

    String title = allResults.get(0).getName();
    if (allResults.size() > 1) {
      title += " + " + TextUtils.pleural(allResults.size() - 1, "other");
    }

    runAnalysis(title, exposureTime, counts, counts.length < originalLength);
  }

  private static boolean showMultiDialog(List<MemoryPeakResults> allResults) {
    // Show a list box containing all the results.
    // This should remember the last set of chosen items.
    final MultiDialog md = ResultsManager.createMultiDialog(TITLE, ResultsManager::hasId);
    md.setSelected(SELECTED_REF.get());
    md.setHelpUrl(HelpUrls.getUrl(HELP_KEY));

    md.showDialog();

    if (md.wasCancelled()) {
      return false;
    }

    final List<String> selected = md.getSelectedResults();
    if (selected.isEmpty()) {
      IJ.error(TITLE, "No results were selected");
      return false;
    }
    SELECTED_REF.set(selected);

    for (final String name : selected) {
      final MemoryPeakResults r = MemoryPeakResults.getResults(name);
      if (r != null) {
        allResults.add(r);
      }
    }

    // Check calibration exists for the first set of results
    if (allResults.isEmpty() || !checkCalibration(allResults.get(0))) {
      return false;
    }

    // Check the calibration is the same for the rest
    final CalibrationReader cal = allResults.get(0).getCalibrationReader();
    final double nmPerPixel = cal.getNmPerPixel();
    final double exposureTime = cal.getExposureTime();
    final DistanceUnit distanceUnit = cal.getDistanceUnit();
    for (int i = 1; i < allResults.size(); i++) {
      final MemoryPeakResults results = allResults.get(i);

      if (!results.hasCalibration()
          || results.getCalibrationReader().getExposureTime() != exposureTime
          || results.getNmPerPixel() != nmPerPixel || results.getDistanceUnit() != distanceUnit) {
        IJ.error(TITLE,
            "The exposure time, pixel pitch and distance unit must match across all the results");
        return false;
      }
    }

    return true;
  }

  /**
   * Check the results have a calibrated exposure time and pixel pitch. If not then show a dialog to
   * collect the calibration.
   *
   * @param results the results
   * @return True if calibrated
   */
  private static boolean checkCalibration(MemoryPeakResults results) {
    if (results.getCalibration() == null || !results.getCalibrationReader().hasExposureTime()
        || !results.getCalibrationReader().hasNmPerPixel()) {
      final CalibrationWriter cal = results.getCalibrationWriterSafe();

      final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
      gd.addMessage("Uncalibrated results! Please enter the calibration:");
      gd.addNumericField("Exposure_time (ms)", cal.getExposureTime(), 2);
      gd.addNumericField("Pixel_pitch (nm)", cal.getNmPerPixel(), 2);
      gd.showDialog();
      if (gd.wasCanceled()) {
        return false;
      }
      cal.setExposureTime(gd.getNextNumber());
      cal.setNmPerPixel(gd.getNextNumber());
      if (gd.invalidNumber() || cal.getExposureTime() <= 0 || cal.getNmPerPixel() <= 0) {
        return false;
      }
      results.setCalibration(cal.getCalibration());
    }
    return true;
  }

  /**
   * Show a dialog to collect settings to filter traces to stationary molecules.
   *
   * @return true if successful
   */
  private boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl(HELP_KEY));

    gd.addMessage("Filter traces to stationary molecules");
    gd.addNumericField("First_frame", settings.firstFrame, 0);
    gd.addNumericField("Last_frame", settings.lastFrame, 0);
    gd.addSlider("Min_size", 2, 20, settings.minSize);
    gd.addNumericField("Mean_distance", settings.meanDistance, 2, 6, "nm");
    gd.addNumericField("Max_distance", settings.maxDistance, 2, 6, "nm");
    gd.addSlider("Max_populations", 1, 3, settings.maxPopulations);
    gd.addNumericField("Max_trace_length", settings.maxTraceLength, 2, 6, "sec");
    gd.addNumericField("Min_bin_count", settings.minBinCount, 0);
    gd.addNumericField("Apparent_dissociation_rate", settings.apparentDissociationRate, -3, 6,
        "sec^-1");
    gd.addNumericField("Boostrap_repeats", settings.bootstrapRepeats, 0);
    gd.addHexField("Bootstrap_seed", settings.bootstrapSeed.toBytes());
    gd.addCheckbox("Show_QQ_plot", settings.showQQplot);

    gd.showDialog();
    if (gd.wasCanceled() || gd.invalidNumber()) {
      return false;
    }
    settings.firstFrame = (int) gd.getNextNumber();
    settings.lastFrame = (int) gd.getNextNumber();
    settings.minSize = (int) gd.getNextNumber();
    settings.meanDistance = gd.getNextNumber();
    settings.maxDistance = gd.getNextNumber();
    settings.maxPopulations = (int) gd.getNextNumber();
    settings.maxTraceLength = (int) gd.getNextNumber();
    settings.minBinCount = (int) gd.getNextNumber();
    settings.apparentDissociationRate = gd.getNextNumber();
    settings.bootstrapRepeats = (int) gd.getNextNumber();
    settings.bootstrapSeed = Seed.from(gd.getNextHexBytes());
    settings.showQQplot = gd.getNextBoolean();

    return !gd.invalidNumber();
  }

  /**
   * Extract the residence times. This function only uses the start frame of the peak results.
   *
   * @param allResults the results
   * @return the histogram of residence times
   */
  private int[] extractResidenceTimes(List<MemoryPeakResults> allResults) {
    final IntArrayList times = new IntArrayList();
    final Function<MemoryPeakResults, List<Trace>> traceFunction = createTraceFunction();
    final CalibrationReader cr = allResults.get(0).getCalibrationReader();
    final double nmPerPixel = cr.getNmPerPixel();
    final double exposureTime = cr.getExposureTime() * 1e-3;
    final ToIntFunction<Trace> framesFunction = createFramesFunction(nmPerPixel);

    final Logger logger = ImageJPluginLoggerHelper.getLogger(getClass());
    for (final MemoryPeakResults r : allResults) {
      final List<Trace> traces = traceFunction.apply(r);
      final int before = times.size();
      final MyStatistics frames = new MyStatistics();
      final MyStatistics maxJumps = new MyStatistics();
      final MyStatistics meanJumps = new MyStatistics();
      traces.forEach(trace -> {
        // Collect the largest jump distance and check for duplicate frames
        PeakResult p = trace.getHead();
        double max = 0;
        for (int i = 1; i < trace.size(); i++) {
          final PeakResult result = trace.get(i);
          if (p.getFrame() == result.getFrame()) {
            logger.warning(() -> String.format(
                "Dataset %s : Ignoring trace ID %d (multiple localisations in the same frame)",
                r.getName(), trace.getId()));
            return;
          }
          max = Math.max(max, p.distance2(result));
          p = result;
        }
        frames.add(trace.getTail().getFrame() - trace.getHead().getFrame() + 1);
        maxJumps.add(Math.sqrt(max));
        meanJumps.add(trace.getMeanDistance());

        final int t = framesFunction.applyAsInt(trace);
        if (t >= 0) {
          times.add(t);
        }
      });
      logger.info(() -> String.format("Dataset %s : %s : %s : %s sec : %s nm : %s nm", r.getName(),
          TextUtils.pleural(traces.size(), "trace"),
          TextUtils.pleural(times.size() - before, "time"), format(frames, exposureTime, "Time"),
          format(maxJumps, nmPerPixel, "Max jump"), format(meanJumps, nmPerPixel, "Mean jump")));
    }

    // Convert to histogram
    final int[] h = new int[times.intStream().max().orElse(0) + 1];
    times.intStream().forEach(i -> h[i]++);
    return h;
  }

  /**
   * Creates the trace function to extract the molecule traces. Note that the traces may be length 1
   * but should not be empty.
   *
   * <p>Any traces that start or end outside the allowed limits are discarded.
   *
   * @return the function
   */
  private Function<MemoryPeakResults, List<Trace>> createTraceFunction() {
    Consumer<MemoryPeakResults> checkCalibration;
    if (settings.maxDistance > 0 || settings.meanDistance > 0) {
      checkCalibration = r -> {
        if (r.getDistanceUnit() != DistanceUnit.PIXEL) {
          // The native units for results is pixels. This is a check to prevent
          // analysis errors should this be changed in the future.
          throw new IllegalStateException(String.format(
              "Stationary analysis requires results calibrated in pixels. Results '%s' in %s",
              r.getName(), r.getDistanceUnit()));
        }
      };
    } else {
      checkCalibration = r -> {
        /* do nothing */ };
    }
    final int firstFrame = settings.firstFrame;
    final int lastFrame = settings.lastFrame;
    return r -> {
      checkCalibration.accept(r);
      // Only use the results with a trace Id
      final LocalList<PeakResult> results = new LocalList<>(r.size());
      final int[] max = {Integer.MIN_VALUE};
      r.forEach((PeakResult p) -> {
        if (p.getId() > 0) {
          results.add(p);
          max[0] = Math.max(p.getEndFrame(), max[0]);
        }
      });
      final LocalList<Trace> traces = new LocalList<>();
      if (!results.isEmpty()) {
        results.sort(IdFramePeakResultComparator.INSTANCE);
        final FrameCounter counter = new FrameCounter(results.unsafeGet(0).getId());
        final Trace[] t = {new Trace()};
        t[0].setId(counter.currentFrame());
        traces.add(t[0]);
        results.forEach(p -> {
          if (counter.advance(p.getId())) {
            t[0] = new Trace();
            t[0].setId(counter.currentFrame());
            traces.add(t[0]);
          }
          t[0].add(p);
        });
        Predicate<Trace> filter = createFrameFilter(firstFrame, lastFrame, max[0]);
        if (filter != null) {
          // Debug
          final int[] count = {0};
          final Logger logger = ImageJPluginLoggerHelper.getLogger(ResidenceTimeAnalysis.class);
          filter = filter.and(tr -> {
            count[0]++;
            logger.fine(() -> String.format("Dataset %s : Filter trace %d [%d-%d]", r.getName(),
                tr.getId(), tr.getHead().getFrame(), tr.getTail().getEndFrame()));
            return true;
          });
          traces.removeIf(filter);
          logger.info(() -> String.format("Dataset %s : Filtered %s using first/last frame",
              r.getName(), TextUtils.pleural(count[0], "trace")));
        }
      }
      return traces;
    };
  }

  /**
   * Creates the frame filter.
   *
   * @param firstFrame the first frame
   * @param lastFrame the last frame
   * @param max the maximum frame in the data
   * @return the predicate
   */
  private static Predicate<Trace> createFrameFilter(int firstFrame, int lastFrame, int max) {
    Predicate<Trace> filter = null;
    // Note: It makes no sense to have a first frame below 1 so this disables the filter
    if (firstFrame > 0) {
      filter = t -> t.getHead().getFrame() <= firstFrame;
    }
    if (lastFrame != 0) {
      final int end = lastFrame < 0 ? lastFrame + 1 + max : lastFrame;
      final Predicate<Trace> filter2 = t -> t.getTail().getEndFrame() >= end;
      filter = filter == null ? filter2 : filter.or(filter2);
    }
    return filter;
  }

  /**
   * Creates the function to convert the trace to the length in frames. The length is the duration
   * in frames that the molecule was observed to be stationary. A value of zero indicates that the
   * trace ended immediately (i.e. it was not present in a subsequent frame).
   *
   * <p>The duration is offset by the minimum trace size. For example if all traces must be size 5
   * then a trace with 5 localisations will return 0 (or greater if there were gaps in the trace).
   *
   * <p>This function will return -1 if the trace does not pass the criteria for a stationary
   * molecule.
   *
   * @param nmPerPixel the nm per pixel
   * @return the length (or -1)
   */
  private ToIntFunction<Trace> createFramesFunction(double nmPerPixel) {
    Predicate<Trace> filter = new Predicate<Trace>() {
      @Override
      public boolean test(Trace t) {
        return true;
      }

      @Override
      public Predicate<Trace> and(Predicate<? super Trace> other) {
        return other::test;
      }
    };
    if (settings.maxDistance > 0) {
      final double maxDistance2 = Math.pow(settings.maxDistance / nmPerPixel, 2);
      filter = t -> {
        PeakResult p = t.getHead();
        for (int i = 1; i < t.size(); i++) {
          final PeakResult result = t.get(i);
          if (p.distance2(result) > maxDistance2) {
            return false;
          }
          p = result;
        }
        return true;
      };
    }
    if (settings.meanDistance > 0) {
      final double meanDistance = settings.meanDistance / nmPerPixel;
      filter = filter.and(t -> t.getMeanDistance() <= meanDistance);
    }
    final int min = Math.max(1, settings.minSize) - 1;
    final Predicate<Trace> traceFilter = filter;
    return t -> {
      if (t.size() > min && traceFilter.test(t)) {
        // frame length = end - start - offset
        // (The shortest allowed trace should return 0)
        return t.getTail().getFrame() - t.getHead().getFrame() - min;
      }
      return -1;
    };
  }

  /**
   * Format the statistics for display.
   *
   * @param stats the stats
   * @param scale the scale
   * @param name the name
   * @return the string
   */
  private static String format(MyStatistics stats, double scale, String name) {
    return new StringBuilder(name).append(" = ").append(MathUtils.round(stats.mean() * scale))
        .append(" +/- ").append(MathUtils.round(stats.sd() * scale)).append(" (")
        .append(MathUtils.round(stats.max() * scale)).append(")").toString();
  }

  /**
   * Do a simulation of bound molecules and then run analysis on the data.
   */
  private void runSimulation() {
    if (!showSimulationDialog()) {
      return;
    }
    final StringBuilder sb = new StringBuilder(256);
    int[] counts = createSimulatation(sb);
    final int originalLength = counts.length;
    counts = filterCounts(counts, settings.exposureTime);
    runAnalysis(sb.toString(), settings.exposureTime, counts, counts.length < originalLength);
  }

  private boolean showSimulationDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl(HELP_KEY));

    settings = Settings.load();
    // Saved by reference so just save now
    settings.save();

    gd.addHexField("Seed", settings.seed.toBytes());
    gd.addNumericField("Samples", settings.samples, 0);
    gd.addNumericField("k0", settings.k0, -2, 6, "/sec");
    gd.addNumericField("k1", settings.k1, -2, 6, "/sec");
    gd.addNumericField("f0", settings.f0, 3, 6, "");
    gd.addNumericField("Exposure_time", settings.exposureTime, 2, 6, "msec");
    gd.addSlider("Max_populations", 1, 3, settings.maxPopulations);
    gd.addNumericField("Max_trace_length", settings.maxTraceLength, 2, 6, "sec");
    gd.addNumericField("Min_bin_count", settings.minBinCount, 0);
    gd.addNumericField("Boostrap_repeats", settings.bootstrapRepeats, 0);
    gd.addHexField("Bootstrap_seed", settings.bootstrapSeed.toBytes());
    gd.addCheckbox("Show_QQ_plot", settings.showQQplot);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.seed = Seed.from(gd.getNextHexBytes());
    settings.samples = (int) Math.abs(gd.getNextNumber());
    settings.k0 = gd.getNextNumber();
    settings.k1 = gd.getNextNumber();
    settings.f0 = gd.getNextNumber();
    settings.exposureTime = gd.getNextNumber();
    settings.maxPopulations = (int) gd.getNextNumber();
    settings.maxTraceLength = (int) gd.getNextNumber();
    settings.minBinCount = (int) gd.getNextNumber();
    settings.bootstrapRepeats = (int) gd.getNextNumber();
    settings.bootstrapSeed = Seed.from(gd.getNextHexBytes());
    settings.showQQplot = gd.getNextBoolean();

    return !gd.invalidNumber();
  }

  /**
   * Creates the simulation.
   *
   * @param sb the output summary of the simulation
   * @return the residence time counts
   */
  private int[] createSimulatation(StringBuilder sb) {
    final double m0 = 1 / settings.k0;
    ValidationUtils.checkStrictlyPositive(m0, "Mean residence time 0");
    final UniformRandomProvider rng = createRng(settings.seed.toBytes());
    final Logger logger = ImageJPluginLoggerHelper.getLogger(getClass());
    DoubleSupplier gen;
    final double t = settings.exposureTime / 1000;
    if (settings.f0 < 1) {
      ValidationUtils.checkArgument(0 < settings.f0, "Simulation fraction must be in (0, 1): %s",
          settings.f0);
      final ZigguratSampler.Exponential exp = ZigguratSampler.Exponential.of(rng);
      final double f0 = settings.f0;
      final double m1 = 1 / settings.k1;
      ValidationUtils.checkStrictlyPositive(m1, "Mean residence time 1");
      sb.append(
          String.format("Simulation: seed=%s, n=%d, k0=%s, k1=%s, f0=%s, time=%s (m0=%s, m1=%s)",
              settings.seed, settings.samples, settings.k0, settings.k1, f0, t, m0, m1));
      logger.info(sb::toString);
      gen = () -> exp.sample() * (rng.nextDouble() < f0 ? m0 : m1);
    } else {
      sb.append(String.format("Simulation: seed=%s, n=%d, k0=%s, time=%s (m0=%s)", settings.seed,
          settings.samples, settings.k0, t, m0));
      logger.info(sb::toString);
      final ZigguratSampler.Exponential exp = ZigguratSampler.Exponential.of(rng, m0);
      gen = exp::sample;
    }
    final int[] counts =
        DoubleStream.generate(gen).limit(settings.samples).mapToInt(x -> (int) (x / t)).toArray();
    // Convert to a histogram
    final int max = Arrays.stream(counts).max().getAsInt();
    final int[] h = new int[max + 1];
    Arrays.stream(counts).forEach(i -> h[i]++);
    return h;
  }

  /**
   * Creates the RNG from the hex seed.
   *
   * @param seed the seed
   * @return the RNG
   */
  private static UniformRandomProvider createRng(byte[] seed) {
    return RandomSource.XO_RO_SHI_RO_128_PP.create(seed);
  }

  /**
   * Run analysis on the data.
   *
   * <p>If the histogram has been clipped, the fitting uses truncated exponential distributions.
   *
   * @param title the title
   * @param exposureTime the exposure time (in milliseconds)
   * @param counts the counts
   * @param clipped true if the counts were clipped at the upper bound
   */
  private void runAnalysis(String title, double exposureTime, int[] counts, boolean clipped) {
    final Plot plot = new Plot("Residence Time", "Time (ms)", "Frequency");
    final float[] x = SimpleArrayUtils.newArray(counts.length, 0f, (float) exposureTime);
    plot.addPoints(x, SimpleArrayUtils.toFloat(counts), Plot.BAR);
    final WindowOrganiser wo = new WindowOrganiser();

    final ResidenceTimeFitting rt = ResidenceTimeFitting.of(exposureTime / 1000, counts, clipped);
    final Logger logger = ImageJPluginLoggerHelper.getLogger(getClass());
    final Consumer<String> summary = createSummaryTable();
    final int samples = Arrays.stream(counts).sum();
    final StringBuilder legend = new StringBuilder("Data");

    // Bootstrap should fit the same samples
    final int[][] bootstrapSamples = createBootstrapCounts(counts, samples, logger);

    // Divide each bar of the histogram into steps when plotting the SF curve.
    // Use a power of 2 for fast masking to an index inside a rolling table.
    final int steps = 32;
    final int mask = steps - 1;
    final float[] t =
        SimpleArrayUtils.newArray((counts.length - 1) * steps, 0, (float) (exposureTime / steps));

    final DoubleUnaryOperator residenceTime = createResidenceTimeFunction();
    final Predicate<Pair<FitResult, Model>> fitAccept =
        createFitAcceptPredicate(exposureTime, counts);

    final int max = MathUtils.clip(1, 3, settings.maxPopulations);
    Model best = null;
    double bestAic = Double.POSITIVE_INFINITY;
    for (int n = 1; n <= max; n++) {
      final Pair<FitResult, Model> f = rt.fit(n, logger);
      if (!fitAccept.test(f)) {
        continue;
      }
      final FitResult r = f.getLeft();
      final Model m = f.getRight();

      // Add to fitted SF curve to plot:
      // pmf(t) = sf(t) - sf(t+dt)
      // Use a rolling table of values to cover 1 bar of the histogram.
      final double[] table =
          IntStream.range(0, steps).mapToDouble(i -> m.sf(t[i] / 1000)).toArray();
      final double[] yy = IntStream.range(0, t.length).mapToDouble(i -> {
        final double sf = m.sf((t[i] + exposureTime) / 1000);
        final int j = i & mask;
        final double p = (table[j] - sf) * samples;
        table[j] = sf;
        return p;
      }).toArray();
      plot.setColor(COLORS[n - 1]);
      plot.addPoints(t, SimpleArrayUtils.toFloat(yy), Plot.LINE);
      legend.append("\tn=").append(n);

      // Bootstrap the data to create confidence intervals for model parameters.
      // Do this in parallel. This could be made an option.
      final long start = System.nanoTime();
      final Ticker ticker = ImageJUtils.createTicker(bootstrapSamples.length,
          Runtime.getRuntime().availableProcessors(), "Fitting Bootstrap samples");
      final List<Model> models = Arrays.stream(bootstrapSamples).parallel().map(h -> {
        final Pair<FitResult, Model> fit =
            ResidenceTimeFitting.of(exposureTime / 1000, h, clipped).fit(m.getSize());
        ticker.tick();
        return fit;
      }).filter(fitAccept).map(Pair::getRight).collect(Collectors.toList());
      logger.log(getTimingLogLevel(),
          () -> String.format("Fitted %d/%d Bootstrap samples (n=%d) in %s", models.size(),
              bootstrapSamples.length, m.getSize(),
              TextUtils.nanosToString(System.nanoTime() - start)));
      IJ.showStatus("");

      final double aic = MathUtils.getAkaikeInformationCriterion(r.getLogLikelihood(), r.getP());
      if (aic < bestAic) {
        bestAic = aic;
        best = m;
      }

      // Report results
      final StringBuilder sb = new StringBuilder(256);
      sb.append(title).append('\t').append(samples).append('\t').append(exposureTime).append('\t')
          .append(n);
      append(sb, n, m::getRate);
      appendCI(sb, n, models, Model::getRate);
      append(sb, n, m::getFraction);
      appendCI(sb, n, models, Model::getFraction);
      sb.append('\t').append(MathUtils.rounded(Math.max(settings.apparentDissociationRate, 0)));
      append(sb, n, i -> residenceTime.applyAsDouble(m.getRate(i)));
      appendCI(sb, n, models, (mm, i) -> residenceTime.applyAsDouble(mm.getRate(i)));
      sb.append('\t').append(r.getLogLikelihood()).append('\t').append(aic).append('\t').append(
          MathUtils.getBayesianInformationCriterion(r.getLogLikelihood(), r.getN(), r.getP()));
      summary.accept(sb.toString());

      if (settings.showQQplot) {
        showQqPlot(exposureTime, counts, m, wo);
      }
    }
    // Warn about poor fitting on truncated data
    if (best != null) {
      final int n = best.getSize();
      final double[] times = new double[n];
      for (int i = 0; i < n; i++) {
        // Time in seconds
        times[i] = residenceTime.applyAsDouble(best.getRate(i));
      }
      if (times[0] * 4e3 > t[t.length - 1]) {
        logger.warning(() -> "Best fit has a mean residence time above 1/4 of the input data range;"
            + " truncated distribution fitting may not be accurate.");
      }
      plot.setColor(COLORS[n - 1]);
      plot.addLabel(0, 0, Arrays.stream(times).mapToObj(MathUtils::rounded)
          .collect(Collectors.joining(", ", "Residence time: ", "")));
    }
    plot.setColor(Color.black);
    plot.addLegend(legend.toString(), "top-right");
    ImageJUtils.display(plot.getTitle(), plot, wo);
    wo.tile();
  }

  private static Consumer<String> createSummaryTable() {
    return ImageJUtils.refresh(SUMMARY_TABLE_REF,
        () -> new TextWindow(TITLE + " Summary",
            "Title\tSamples\tTime (ms)\tn\tk (sec^-1)\tCI\tfraction\tCI\t"
                + "ka\tResidence time (sec)\tCI\tLL\tAIC\tBIC",
            "", 1200, 300))::append;
  }

  /**
   * Creates the bootstrap counts. The original counts are used to create a discrete probability
   * distribution which is sampled for each bootstrap repeat.
   *
   * @param counts the original counts
   * @param samples the total number of samples
   * @param logger the logger
   * @return the counts
   */
  private int[][] createBootstrapCounts(int[] counts, final int samples, Logger logger) {
    if (settings.bootstrapRepeats > 1) {
      logger.info(() -> String.format("Bootstrapping: seed=%s, repeats=%d", settings.bootstrapSeed,
          settings.bootstrapRepeats));
      final long start = System.nanoTime();
      final UniformRandomProvider rng = createRng(settings.bootstrapSeed.toBytes());
      final double[] probabilities =
          Arrays.stream(counts).mapToDouble(i -> (double) i / samples).toArray();
      final DiscreteSampler s = AliasMethodDiscreteSampler.of(rng, probabilities);
      final int[][] bootstrapSamples = new int[settings.bootstrapRepeats][counts.length];
      for (final int[] h : bootstrapSamples) {
        IntStream.generate(s::sample).limit(samples).forEach(j -> h[j]++);
      }
      logger.log(getTimingLogLevel(),
          () -> String.format("Created %d Bootstrap samples of size %s in %s",
              bootstrapSamples.length, samples,
              TextUtils.nanosToString(System.nanoTime() - start)));
      return bootstrapSamples;
    }
    return new int[0][0];
  }

  /**
   * Creates the residence time function using the observed dissociation rate and the apparent
   * dissociation rate.
   *
   * <pre>
   * rho = 1 / (k_obs - k_a)
   * </pre>
   *
   * @return the function
   */
  private DoubleUnaryOperator createResidenceTimeFunction() {
    final double kb = settings.apparentDissociationRate;
    return kb > 0 ? k -> 1 / Math.max(0, k - kb) : k -> 1 / k;
  }

  /**
   * Create a predicate to test whether to accept a fit.
   *
   * @param exposureTime the exposure time (in milliseconds)
   * @param counts the counts
   * @return the predicate
   */
  private static Predicate<Pair<FitResult, Model>> createFitAcceptPredicate(double exposureTime,
      int[] counts) {
    final double max = counts.length * exposureTime / 1000;
    // If the entire distribution is within a single bin of the histogram then the fit
    // is invalid. Set a threshold based on 95% of the CDF of an exponential in the exposure
    // time. exponential(mean=1) inverse cdf(0.95) = 2.9957
    final double min = exposureTime / 1000 / 3;
    return fit -> {
      if (fit == null) {
        return false;
      }
      final Model m = fit.getRight();
      final int n = m.getSize();
      // Assume the model returns the residence times in descending order
      return m.getResidenceTime(0) < max && m.getResidenceTime(n - 1) > min;
    };
  }

  /**
   * Appends the value to the StringBuilder. The text is prefixed with a tab character and values
   * are delimited by commas.
   *
   * @param sb the StringBuilder
   * @param n the number of values
   * @param value the value mapping function
   */
  private static void append(StringBuilder sb, int n, IntToDoubleFunction value) {
    sb.append('\t');
    for (int i = 0; i < n; i++) {
      if (i != 0) {
        sb.append(", ");
      }
      sb.append(MathUtils.rounded(value.applyAsDouble(i)));
    }
  }

  /**
   * Appends the value to the StringBuilder. The text is prefixed with a tab character and values
   * are delimited by commas.
   *
   * @param sb the StringBuilder
   * @param n the number of values
   * @param value the value mapping function
   */
  private static void appendCI(StringBuilder sb, int n, List<Model> models,
      ObjIntToDoubleFunction<Model> value) {
    if (models.isEmpty()) {
      sb.append('\t');
      return;
    }
    append(sb, n, i -> getConfidenceInterval(models, mm -> value.apply(mm, i)));
  }

  /**
   * Gets the 95% confidence interval for the value from the model.
   *
   * @return the confidence interval
   */
  private static double getConfidenceInterval(List<Model> models, ToDoubleFunction<Model> value) {
    final Statistics ss = new Statistics();
    models.forEach(m -> ss.add(value.applyAsDouble(m)));
    return ss.getConfidenceInterval(0.95);
  }

  /**
   * Gets the log level for timing messages. ImageJ debug mode sets this to INFO, otherwise FINE.
   *
   * @return the timing log level
   */
  private static Level getTimingLogLevel() {
    return IJ.debugMode ? Level.INFO : Level.FINE;
  }

  /**
   * Filter the histogram of counts. Uses the max trace length, and/or the min bin count to set the
   * last valid bin, to truncate the histogram.
   *
   * @param counts the counts
   * @param exposureTime the exposure time (in milliseconds)
   * @return the new counts
   */
  private int[] filterCounts(int[] counts, double exposureTime) {
    int limit = counts.length;
    int minBinCount = settings.minBinCount;
    if (settings.maxTraceLength > 0) {
      // trace length (sec) / exposure time (ms)
      limit = Math.min(limit, (int) Math.round(settings.maxTraceLength * 1e3 / exposureTime));
      // Always clip to the last bin containing a count of 1.
      minBinCount = Math.max(1, minBinCount);
    }
    if (minBinCount > 0) {
      while (limit > 0 && counts[limit - 1] < minBinCount) {
        limit--;
      }
    }
    if (counts.length > limit) {
      counts = Arrays.copyOf(counts, limit);
    }
    return counts;
  }

  /**
   * Show a QQ plot of the quantiles.
   *
   * @param exposureTime the exposure time (in milliseconds)
   * @param counts the counts
   * @param m the model
   * @param wo the window organiser
   */
  private static void showQqPlot(double exposureTime, int[] counts, Model m, WindowOrganiser wo) {
    // Get an approximation for the search for the quantiles
    // cdf[i] = time (i+1) * exposureTime
    final double et = exposureTime / 1000;
    // QQ plot using the quantile definition (k - 1/3) / (n + 1/3).
    // See https://en.wikipedia.org/wiki/Q%E2%80%93Q_plot#Heuristics
    final int n = Arrays.stream(counts).sum();
    final double a = 1.0 / 3;
    final double na = n + a;
    int k = 0;
    final DoubleArrayList x = new DoubleArrayList();
    final DoubleArrayList y = new DoubleArrayList();
    final DoubleArrayList xp = new DoubleArrayList();
    // Moving bracket for the quantile
    double upper = et;
    final BrentSolver solver = new BrentSolver(1e-4, 1e-6, 1e-6);
    for (int i = 0; i < counts.length; i++) {
      if (counts[i] == 0) {
        continue;
      }
      y.add((i + 1) * et);
      k += counts[i];
      // For the survival function we subtract from 1
      final double p = (k - a) / na;
      final double q = 1 - p;
      // Raise bracket
      while (m.sf(upper) > q) {
        upper += et;
      }
      // Find
      final double z = solver.findRoot(xx -> m.sf(xx) - q, upper - et, upper);
      x.add(z);
      xp.add(p);
    }
    final Plot plot = new Plot("Residence Time QQ n=" + m.getSize(), "Model (sec)", "Data (sec)");
    plot.setColor(COLORS[0]);
    plot.addPoints(x.toDoubleArray(), y.toDoubleArray(), Plot.CIRCLE);
    plot.setColor(COLORS[1]);
    final double min = Math.min(et, y.getDouble(0)) / 2;
    final double maxx = x.getDouble(x.size() - 1);
    final double maxy = y.getDouble(y.size() - 1);
    final double max = Math.max(maxx, maxy);
    plot.drawLine(min, min, max * 2, max * 2);
    // For a log scale we must set the limits manually.
    // Currently disabled but setting the limits allows it to be set using More... > Set Range...
    // plot.setLogScaleX();
    // plot.setLogScaleY();
    plot.setLimits(min, maxx * 1.025, min, maxy * 1.025);
    // Add known model quantiles
    final double[] quantiles = {0.5, 0.75, 0.90, 0.95, 0.99};
    final double[] cdf = xp.toDoubleArray();
    final DoubleArrayList points = new DoubleArrayList();
    for (final double p : quantiles) {
      // Estimate the bracket
      int index = Arrays.binarySearch(cdf, p);
      if (index < 0) {
        index = ~index;
      }
      double lo;
      double hi;
      if (index == 0) {
        lo = 0;
        hi = x.getDouble(index);
      } else if (index == x.size()) {
        lo = x.getDouble(index);
        // The last point should correspond to ~99% of the distribution. For an exponential
        // if we double the value it should bracket the quantile.
        hi = lo * 2;
      } else {
        // Valid bracket
        lo = x.getDouble(index - 1);
        hi = x.getDouble(index);
      }
      final double q = 1 - p;
      final DoubleUnaryOperator func = xx -> m.sf(xx) - q;
      // Ensure the bracket is wide enough
      while (Double.compare(func.applyAsDouble(lo) * func.applyAsDouble(hi), 0.0) >= 0) {
        lo *= 0.9375;
        hi *= 1.25;
      }
      points.add(solver.findRoot(func, lo, hi));
    }
    final float[] qp = SimpleArrayUtils.toFloat(points.toDoubleArray());
    plot.setColor(COLORS[2]);
    plot.addPoints(qp, qp, Plot.DIAMOND);
    plot.setColor(Color.BLACK);
    plot.addLabel(0, 0, Arrays.stream(quantiles).mapToObj(String::valueOf)
        .collect(Collectors.joining(", ", "Quantiles: ", "")));
    ImageJUtils.display(plot.getTitle(), plot, wo);
  }
}
