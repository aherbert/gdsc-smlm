/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Package
 *
 * Software for single molecule localisation microscopy (SMLM) in ImageJ
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.process.LUT;
import ij.text.TextWindow;
import it.unimi.dsi.fastutil.floats.FloatArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import java.awt.Color;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.exception.TooManyIterationsException;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.numbers.core.DD;
import org.apache.commons.rng.JumpableUniformRandomProvider;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import org.apache.commons.rng.simple.RandomSource;
import uk.ac.sussex.gdsc.core.ij.ImageJPluginLoggerHelper;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.MultiDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper.LutColour;
import uk.ac.sussex.gdsc.core.logging.LoggerUtils;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;
import uk.ac.sussex.gdsc.core.utils.rng.UniformRandomProviders;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.fitting.DiffusionAnalysis;
import uk.ac.sussex.gdsc.smlm.function.ChiSquaredDistributionTable;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.MemoryResultsList;
import uk.ac.sussex.gdsc.smlm.math3.optim.nonlinear.scalar.noderiv.CustomPowellOptimizer;
import uk.ac.sussex.gdsc.smlm.results.IdPeakResult;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyrResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.sort.IdFramePeakResultComparator;

/**
 * Fit the diffusion distance of a set of molecule tracks using different time delays.
 *
 * <p>Fitting accounts for the expected diffusion of molecules out of the depth-of-field after a
 * given time delay. This is based on the Spot-On model described in the methods section of the
 * paper:
 *
 * <p>Hansen, A.S., Woringer, M., Grimm, J.B., Lavis, L.D., Tjian, R., and Darzacq, X. (2018) Robust
 * model-based analysis of single-particle tracking experiments with Spot-On. eLife 7, e33125.
 * doi:10.7554/eLife.33125.
 */
public class TrackDiffusionAnalysis implements PlugIn {
  private static final String TITLE = "Track Diffusion Analysis";

  private static final AtomicReference<TextWindow> TABLE_REF = new AtomicReference<>();

  /** The exposure time in seconds. */
  private double exposureTime;
  private Logger logger;

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> INSTANCE = new AtomicReference<>(new Settings());

    List<String> selected;

    double depthOfField;
    int gap;
    int maxT;
    int offsets;

    double binWidth;

    double a;
    double b;

    int groupSize;
    int repeats;
    double minD;
    double maxD;

    double exposureTime;
    int numberOfMolecules;
    int simT;
    double precision;
    double f1;
    double f2;
    double d1;
    double d2;
    double d3;

    boolean fitPrecision;
    boolean fitThreeState;
    boolean fitMle;
    double significanceLevel;

    Settings() {
      selected = new LocalList<>();

      depthOfField = 750;
      gap = 1;
      maxT = 7;
      offsets = 0;

      binWidth = 0.01;

      groupSize = 5;
      repeats = 3;
      minD = 1;
      maxD = 5;

      exposureTime = 10;
      numberOfMolecules = 50000;
      simT = 10;
      precision = 30;
      f1 = 0.2;
      d1 = 0.01;
      d2 = 1.5;
      d3 = 5;
      significanceLevel = 0.05;
    }

    Settings(Settings source) {
      selected = new LocalList<>(source.selected);

      depthOfField = source.depthOfField;
      gap = source.gap;
      maxT = source.maxT;
      offsets = source.offsets;

      binWidth = source.binWidth;

      a = source.a;
      b = source.b;

      groupSize = source.groupSize;
      repeats = source.repeats;
      minD = source.minD;
      maxD = source.maxD;

      exposureTime = source.exposureTime;
      numberOfMolecules = source.numberOfMolecules;
      simT = source.simT;
      precision = source.precision;
      f1 = source.f1;
      f2 = source.f2;
      d1 = source.d1;
      d2 = source.d2;
      d3 = source.d3;

      fitPrecision = source.fitPrecision;
      fitThreeState = source.fitThreeState;
      fitMle = source.fitMle;
      significanceLevel = source.significanceLevel;
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
     * Save the settings. This can be called only once as it saves via a reference.
     */
    void save() {
      INSTANCE.set(this);
    }
  }

  /**
   * Store the localisation position.
   */
  static class Position {
    int t;
    float x;
    float y;

    Position(int t, float x, float y) {
      this.t = t;
      this.x = x;
      this.y = y;
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (ImageJUtils.isExtraOptions()) {
      simulateTracks();
      return;
    }

    final MemoryResultsList items = new MemoryResultsList(MemoryPeakResults::hasId);

    if (items.isEmpty()) {
      IJ.error(TITLE, "No traced localisations in memory");
      return;
    }

    settings = Settings.load();
    settings.save();

    final List<MemoryPeakResults> allResults = new LocalList<>();

    if (!showDialog() || !showMultiDialog(allResults, items)) {
      return;
    }

    logger = ImageJPluginLoggerHelper.getLogger(getClass());

    final float[][] distances = getDistances(allResults);
    final float[][] pdf = createPdf(distances);

    final int threadCount = Prefs.getThreads();
    final ExecutorService executor = Executors.newFixedThreadPool(threadCount);

    try {
      PointValuePair result;
      // Update the fit method to accept distances or pdf and fit MLE or SS
      if (settings.fitMle) {
        final float[][] gDistances = groupDistances(distances);
        result = fitDistances(gDistances, executor, true);
      } else {
        result = fitDistances(pdf, executor, false);
      }
      plotDistances(distances, pdf, result, executor);
    } finally {
      executor.shutdown();
    }

    ImageJUtils.finished(TITLE + " done");
  }

  private void simulateTracks() {
    if (!showSimulationDialog()) {
      return;
    }

    final int threadCount = Prefs.getThreads();
    final ExecutorService executor = Executors.newFixedThreadPool(threadCount);

    try {
      IJ.showStatus("Simulating tracks...");

      final double halfDz = settings.depthOfField / 2000;
      final double dt = settings.exposureTime / 1000;
      final int g = settings.gap;
      final int maxT = settings.simT;

      // Sample populations
      final double p1 = settings.f1;
      final double p2 = settings.f2 > 0 ? p1 + settings.f2 : 1;

      // Precision in um (used for the units of the tracks)
      final double precision = settings.precision / 1000;

      // Simulate tracks across the depth of field.
      // Each track is scaled using the diffusion coefficient and the molecule tested if
      // it remains in the depth-of-field.

      final double[] sampleD = {settings.d1, settings.d2, settings.d3};
      // Create the Gaussian step for each diffusion coefficient
      final double[] step = Arrays.stream(sampleD).map(
          // Std.dev of Gaussian for the step size: sqrt(2D * dt)
          d -> Math.sqrt(2 * d * dt)).toArray();

      final List<Future<MemoryPeakResults>> futures = new LinkedList<>();

      final int total = settings.numberOfMolecules;
      final Ticker ticker = ImageJUtils.createTicker(total, threadCount);

      final AtomicInteger position = new AtomicInteger(total);
      final AtomicInteger nextId = new AtomicInteger();
      // Require 4 dimension equi-distributed sampler
      final JumpableUniformRandomProvider parentRng =
          (JumpableUniformRandomProvider) RandomSource.XO_RO_SHI_RO_1024_PP.create();
      for (int n = 0; n < threadCount; n++) {
        final UniformRandomProvider rng = parentRng.jump();
        futures.add(executor.submit(() -> {
          final NormalizedGaussianSampler sampler =
              SamplerUtils.createNormalizedGaussianSampler(rng);
          final MemoryPeakResults observed = new MemoryPeakResults(total);

          for (;;) {
            final int pos = position.decrementAndGet();
            if (pos < 0) {
              break;
            }
            // Molecule type
            final double p = rng.nextDouble();
            double d;
            if (p < p1) {
              d = step[0];
            } else if (p < p2) {
              d = step[1];
            } else {
              d = step[2];
            }
            // Simulate the track from the origin within the depth of field
            int id = nextId.incrementAndGet();
            double x = 0;
            double y = 0;
            double z = rng.nextDouble(-halfDz, halfDz);
            int last = 0;
            for (int i = 1; i <= maxT; i++) {
              x += sampler.sample() * d;
              y += sampler.sample() * d;
              z += sampler.sample() * d;
              if (Math.abs(z) < halfDz) {
                // Record molecule
                // Check the gap from the last time it was in the depth-of-field
                if (i - last > g) {
                  // Start a new track
                  id = nextId.incrementAndGet();
                }
                last = i;
                // Add localisation precision
                observed.add(new IdPeakResult(i, (float) (x + sampler.sample() * precision),
                    (float) (y + sampler.sample() * precision), 1f, id));
              }
            }
            ticker.tick();
          }
          return observed;
        }));
      }

      final MemoryPeakResults results = new MemoryPeakResults(total);
      final CalibrationWriter cw = new CalibrationWriter();
      cw.setExposureTime(settings.exposureTime);
      cw.setNmPerPixel(100);
      cw.setDistanceUnit(DistanceUnit.UM);
      results.setCalibration(cw.getCalibration());
      results.setName(TITLE);
      MemoryPeakResults.addResults(results);

      // Uncomment to allow rendering the dataset (which requires pixel units)
      // results.convertToPreferredUnits();

      // Pre-select the simulation
      if (settings.selected.isEmpty()) {
        settings.selected.add(TITLE);
        // These are fitted values for the simulation defaults
        settings.a = 0.19;
        settings.b = 0.06;
      }

      // Build final results in memory
      for (final Future<MemoryPeakResults> f : futures) {
        try {
          results.add(f.get());
        } catch (InterruptedException | ExecutionException e) {
          throw new RuntimeException(e);
        }
      }

    } finally {
      executor.shutdown();
    }

    ImageJUtils.finished();
  }

  private boolean showSimulationDialog() {
    settings = Settings.load();
    settings.save();

    final GenericDialog gd = new GenericDialog(TITLE);

    gd.addNumericField("Depth_of_field", settings.depthOfField, 1, 6, "nm");
    gd.addNumericField("Exposure_time", settings.exposureTime, 1, 6, "ms");
    gd.addSlider("Gap", 1, 5, settings.gap);

    gd.addNumericField("Number_of_molecules", settings.numberOfMolecules);
    gd.addSlider("Max_t", 5, 15, settings.simT);

    gd.addNumericField("Precision", settings.precision, 1, 6, "nm");
    gd.addNumericField("F1", settings.f1, 3);
    gd.addNumericField("F2", settings.f2, 3);
    gd.addNumericField("D1", settings.d1, 3, 6, "um^2/s");
    gd.addNumericField("D2", settings.d2, 3, 6, "um^2/s");
    gd.addNumericField("D3", settings.d3, 3, 6, "um^2/s");

    gd.addHelp(HelpUrls.getUrl("diffusion-depth-of-field"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    settings.depthOfField = gd.getNextNumber();
    settings.exposureTime = gd.getNextNumber();
    settings.gap = (int) gd.getNextNumber();

    settings.numberOfMolecules = (int) gd.getNextNumber();
    settings.simT = (int) gd.getNextNumber();

    settings.precision = gd.getNextNumber();
    settings.f1 = gd.getNextNumber();
    settings.f2 = gd.getNextNumber();
    settings.d1 = gd.getNextNumber();
    settings.d2 = gd.getNextNumber();
    settings.d3 = gd.getNextNumber();

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Depth of Field", settings.depthOfField);
      ParameterUtils.isAboveZero("Exposure time", settings.exposureTime);
      ParameterUtils.isAboveZero("Gap", settings.gap);

      ParameterUtils.isAboveZero("Number of molecules", settings.numberOfMolecules);
      ParameterUtils.isAboveZero("Max T", settings.simT);

      ParameterUtils.isPositive("Precision", settings.precision);
      ParameterUtils.isAboveZero("F1", settings.f1);
      ParameterUtils.isPositive("F2", settings.f2);
      ParameterUtils.isPositive("D1", settings.d1);
      ParameterUtils.isPositive("D2", settings.d2);
      ParameterUtils.isPositive("D3", settings.d3);

      ParameterUtils.isBelow("F1 + F2", settings.f1 + settings.f2, 1);

    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    if (gd.invalidNumber()) {
      return false;
    }

    return true;
  }

  private boolean showDialog() {
    final GenericDialog gd = new GenericDialog(TITLE);

    gd.addNumericField("Depth_of_field", settings.depthOfField, 1, 6, "nm");
    gd.addSlider("Max_t", 3, 10, settings.maxT);
    gd.addSlider("Offsets", 0, 5, settings.offsets);
    gd.addNumericField("Precision", settings.precision, 1, 6, "nm");

    gd.addNumericField("Bin_width", settings.binWidth, -3);

    gd.addNumericField("A", settings.a, 4);
    gd.addNumericField("B", settings.b, 4);

    gd.addNumericField("Group_size", settings.groupSize, 0);
    gd.addSlider("Repeats", 1, 5, settings.repeats);
    gd.addNumericField("Min_D", settings.minD, 3, 6, "um^2/s");
    gd.addNumericField("Max_D", settings.maxD, 3, 6, "um^2/s");

    gd.addCheckbox("Fit_precision", settings.fitPrecision);
    gd.addCheckbox("Fit_three_state", settings.fitThreeState);
    gd.addCheckbox("Fit_MLE", settings.fitMle);
    gd.addNumericField("significanceLevel", settings.significanceLevel, -3);

    gd.addHelp(HelpUrls.getUrl("track-diffusion-analysis"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    settings.depthOfField = gd.getNextNumber();
    settings.maxT = (int) gd.getNextNumber();
    settings.offsets = (int) gd.getNextNumber();
    settings.precision = gd.getNextNumber();

    settings.binWidth = gd.getNextNumber();

    settings.a = gd.getNextNumber();
    settings.b = gd.getNextNumber();

    settings.groupSize = (int) gd.getNextNumber();
    settings.repeats = (int) gd.getNextNumber();
    settings.minD = gd.getNextNumber();
    settings.maxD = gd.getNextNumber();

    settings.fitPrecision = gd.getNextBoolean();
    settings.fitThreeState = gd.getNextBoolean();
    settings.fitMle = gd.getNextBoolean();
    settings.significanceLevel = gd.getNextNumber();

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Depth of Field", settings.depthOfField);
      ParameterUtils.isAboveZero("Max T", settings.maxT);
      ParameterUtils.isPositive("Offsets", settings.offsets);
      ParameterUtils.isPositive("Precision", settings.precision);

      ParameterUtils.isAboveZero("Bin width", settings.binWidth);

      ParameterUtils.isPositive("A", settings.a);
      ParameterUtils.isPositive("B", settings.b);

      ParameterUtils.isPositive("Group size", settings.groupSize);
      ParameterUtils.isAboveZero("Repeats", settings.repeats);
      ParameterUtils.isAboveZero("Min D", settings.minD);
      ParameterUtils.isAboveZero("Max D", settings.maxD);
      ParameterUtils.isAboveZero("Significance level", settings.significanceLevel);

      ParameterUtils.isEqualOrAbove("Max D", settings.maxD, settings.minD);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    if (gd.invalidNumber()) {
      return false;
    }

    return true;
  }

  private boolean showMultiDialog(List<MemoryPeakResults> allResults, MemoryResultsList items) {
    // Show a list box containing all the results. This should remember the last set of chosen
    // items.
    final MultiDialog md = new MultiDialog(TITLE, items);
    md.setDisplayConverter(items.getDisplayConverter());
    md.setSelected(settings.selected);
    md.setHelpUrl(HelpUrls.getUrl("track-diffusion-analysis"));

    md.showDialog();

    if (md.wasCancelled()) {
      return false;
    }

    final List<String> selected = md.getSelectedResults();
    if (selected.isEmpty()) {
      IJ.error(TITLE, "No results were selected");
      return false;
    }
    settings.selected = selected;

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
    this.exposureTime = exposureTime * 1e-3;
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
      if (gd.wasCanceled() || gd.invalidNumber()) {
        return false;
      }
      cal.setExposureTime(gd.getNextNumber());
      cal.setNmPerPixel(gd.getNextNumber());
      if (cal.getExposureTime() <= 0 || cal.getNmPerPixel() <= 0) {
        return false;
      }
      results.setCalibration(cal.getCalibration());
    }
    return true;
  }

  /**
   * Get the jump distances for each time delay.
   *
   * @param allResults the results
   * @return distances
   */
  private float[][] getDistances(List<MemoryPeakResults> allResults) {
    final int maxT = settings.maxT;
    final FloatArrayList[] distances =
        IntStream.rangeClosed(0, maxT).mapToObj(FloatArrayList::new).toArray(FloatArrayList[]::new);
    final List<Position> track = new LocalList<>();
    for (final MemoryPeakResults results : allResults) {
      results.sort(IdFramePeakResultComparator.INSTANCE);
      final int[] id = {-1};
      results.forEach(DistanceUnit.UM, (XyrResultProcedure) (x, y, r) -> {
        if (r.getId() != id[0]) {
          id[0] = r.getId();
          // Process track
          final int maxStart = Math.min(track.size() - 1, settings.offsets + 1);
          for (int i = 0; i < maxStart; i++) {
            final Position origin = track.get(i);
            for (int j = i + 1; j < track.size(); j++) {
              final Position position = track.get(j);
              final int g = position.t - origin.t;
              if (g > maxT) {
                break;
              }
              distances[g].add(MathUtils.distance(origin.x, origin.y, position.x, position.y));
            }
          }
          track.clear();
        }
        track.add(new Position(r.getFrame(), x, y));
      });
      // Process final track
      for (int i = 0; i < track.size(); i++) {
        final Position origin = track.get(i);
        for (int j = i + 1; j < track.size(); j++) {
          final Position position = track.get(i);
          final int g = position.t - origin.t;
          if (g > maxT) {
            break;
          }
          distances[g].add(MathUtils.distance(origin.x, origin.y, position.x, position.y));
        }
      }
    }

    // Convert to sorted arrays
    final float[][] sortedDistances = IntStream.rangeClosed(1, maxT).mapToObj(i -> {
      float[] data = distances[i].toFloatArray();
      Arrays.sort(data);
      // Remove zero distances. These have a log-likelihood of -infinity.
      // Since they have zero probability these will be artifacts from the data
      // such as an average localisation position over time.
      int index = 0;
      while (index < data.length && data[index] == 0) {
        index++;
      }
      if (index > 0) {
        LoggerUtils.log(logger, Level.INFO, "Removing zero distances at time delay=[%d] %d", i,
            index);
        data = Arrays.copyOfRange(data, index, data.length);
      }
      return data;
    }).toArray(float[][]::new);

    final int[] sizes = Arrays.stream(sortedDistances).mapToInt(x -> x.length).toArray();
    LoggerUtils.log(logger, Level.INFO, "Distance counts: %s", Arrays.toString(sizes));

    // Avoid fitting when there is no data
    final int time = SimpleArrayUtils.findIndex(sizes, i -> i == 0);
    if (time >= 0) {
      throw new IllegalStateException("No sizes for time delay: " + (time + 1));
    }

    return sortedDistances;
  }

  /**
   * Group the distances into average distances in each group. This can be used to reduce the number
   * of observations for fitting.
   *
   * @param distances the distances (must be sorted)
   * @return the grouped distances
   */
  private float[][] groupDistances(float[][] distances) {
    final int groupSize = settings.groupSize;
    if (groupSize < 2) {
      return distances;
    }
    final List<FloatArrayList> groupDistances = new LocalList<>();
    for (final float[] d : distances) {
      final FloatArrayList data = new FloatArrayList();
      for (int i = 0; i < d.length; i += groupSize) {
        double sum = 0;
        final int end = Math.min(i + groupSize, d.length);
        for (int j = i; j < end; j++) {
          sum += d[j];
        }
        data.add((float) (sum / (end - i)));
      }
      groupDistances.add(data);
    }

    final int[] sizes = groupDistances.stream().mapToInt(FloatArrayList::size).toArray();
    LoggerUtils.log(logger, Level.INFO, "Grouped distance counts: %s", Arrays.toString(sizes));

    return groupDistances.stream().map(FloatArrayList::toFloatArray).toArray(float[][]::new);
  }

  private PointValuePair fitDistances(float[][] distances, ExecutorService executor, boolean mle) {
    final PointValuePair result2 = fitTwoStateDistances(distances, executor, mle);
    if (result2 == null) {
      return null;
    }
    addToResultTable(result2);
    final PointValuePair result3 = fitThreeStateDistances(distances, executor, mle);
    if (result3 != null) {
      addToResultTable(result3);
      double ll2 = result2.getValue();
      double ll3 = result3.getValue();

      if (!mle) {
        // Choose best model for SS using approximate log-likelihood
        final int n = distances.length * distances[0].length;
        ll2 = MathUtils.getLogLikelihood(ll2, n);
        ll3 = MathUtils.getLogLikelihood(ll3, n);
      }

      // Select best fit using log-likelihood ratio test
      final double llr = 2 * (ll3 - ll2);
      if (llr >= 0) {
        // The difference in the number of fitted parameters will be 2:
        // i.e. 1 extra diffusion coefficient and population fraction
        final double q = ChiSquaredDistributionTable.computeQValue(llr, 2);
        final boolean reject = q > settings.significanceLevel;
        LoggerUtils.log(logger, Level.INFO,
            "Two-state -> three-state : LL = %s -> %s, LLR = %s, q-value = %s (Reject=%b)",
            MathUtils.rounded(ll2, 4), MathUtils.rounded(ll3, 4),
            MathUtils.rounded(llr, 4), MathUtils.rounded(q, 4), reject);
        if (!reject) {
          return result3;
        }
      }
    }
    return result2;
  }

  private PointValuePair fitTwoStateDistances(float[][] distances, ExecutorService executor,
      boolean mle) {
    IJ.showStatus("Fitting two-state model...");

    final double dt = exposureTime;
    final double dz = settings.depthOfField / 1000;
    final double precision = settings.precision / 1000;
    final UniformRandomProvider rng = UniformRandomProviders.create();

    final MaxEval maxEval = new MaxEval(2000);
    final CustomPowellOptimizer powellOptimizer = createCustomPowellOptimizer();

    double best = Double.NEGATIVE_INFINITY;
    PointValuePair result = null;

    // fit: f1, d1, f2, d2, sigma
    // Note that f2 = 1 - f1 so the number of fitted parameters is 4.
    // However for convenience f2 is also a fitted parameter for the optimiser.
    final ObjectiveFunction fun;
    final GoalType goalType;
    if (mle) {
      fun = new ObjectiveFunction(new TwoStateModelFunction(dt, dz, settings.a, settings.b,
          precision, distances, executor));
      goalType = GoalType.MAXIMIZE;
    } else {
      fun = new ObjectiveFunction(new TwoStateSSFunction(dt, dz, settings.a, settings.b, precision,
          settings.binWidth, distances, executor));
      goalType = GoalType.MINIMIZE;
      best = Double.POSITIVE_INFINITY;
    }
    final SimpleBounds bounds = new SimpleBounds(addPrecision(new double[] {0, 0.001, 0, 0.1}, 0),
        addPrecision(new double[] {1, 0.1, 1, Double.POSITIVE_INFINITY}, 0.1));
    final CustomPowellOptimizer.BasisStep step = new CustomPowellOptimizer.BasisStep(
        addPrecision(new double[] {0.1, 0.01, 0.1, 0.5}, 0.001));

    final Ticker ticker = ImageJUtils.createTicker(settings.repeats, 1);
    for (int n = 1; n <= settings.repeats; n++) {
      try {
        // Guess in range
        final double[] start = addPrecision(new double[] {
            // Bound fraction
            rng.nextDouble(0.1, 0.9), rng.nextDouble(0.1),
            // Free fraction
            0, rng.nextDouble(settings.minD, settings.maxD)},
            // precision
            rng.nextDouble(0.01, 0.03));
        // Set f2
        start[2] = 1 - start[0];

        final PointValuePair solution =
            powellOptimizer.optimize(maxEval, fun, bounds, step, goalType, new InitialGuess(start));

        if (mle) {
          LoggerUtils.log(logger, Level.INFO,
              "Two-state fit [%d]: MLE = %s, Akaike IC = %s (%d evaluations)", n,
              solution.getValue(), MathUtils.getAkaikeInformationCriterion(solution.getValue(), 4),
              powellOptimizer.getEvaluations());
          if (solution.getValue() > best) {
            best = solution.getValue();
            result = solution;
          }
        } else {
          LoggerUtils.log(logger, Level.INFO, "Two-state fit [%d]: SS = %s (%d evaluations)", n,
              solution.getValue(), powellOptimizer.getEvaluations());
          if (solution.getValue() < best) {
            best = solution.getValue();
            result = solution;
          }
        }
      } catch (final TooManyEvaluationsException ex) {
        LoggerUtils.log(logger, Level.INFO,
            "Powell optimiser [%d] failed : Too many evaluation (%d)", n,
            powellOptimizer.getEvaluations());
      } catch (final TooManyIterationsException ex) {
        LoggerUtils.log(logger, Level.INFO,
            "Powell optimiser [%d] failed : Too many iterations (%d)", n,
            powellOptimizer.getIterations());
      } catch (final ConvergenceException ex) {
        LoggerUtils.log(logger, Level.INFO, "Powell optimiser [%d] failed to fit : %s", n,
            ex.getMessage());
      }
      ticker.tick();
    }

    if (result != null) {
      // Normalise fractions
      sort2state(result.getPointRef());
    }

    return result;
  }

  private static void sort2state(final double[] r) {
    final double sum = r[0] + r[2];
    r[0] /= sum;
    r[2] /= sum;
    // Sort (in event that the free fraction is slower than the bound)
    if (r[3] < r[1]) {
      final double f = r[0];
      final double d = r[1];
      r[0] = r[2];
      r[1] = r[3];
      r[2] = f;
      r[3] = d;
    }
  }

  private PointValuePair fitThreeStateDistances(float[][] distances, ExecutorService executor,
      boolean mle) {
    if (!settings.fitThreeState) {
      return null;
    }

    IJ.showStatus("Fitting three-state model...");

    final double dt = exposureTime;
    final double dz = settings.depthOfField / 1000;
    final double precision = settings.precision / 1000;
    final UniformRandomProvider rng = UniformRandomProviders.create();

    final MaxEval maxEval = new MaxEval(2000);
    final CustomPowellOptimizer powellOptimizer = createCustomPowellOptimizer();

    double best = Double.NEGATIVE_INFINITY;
    PointValuePair result = null;

    // fit: f1, d1, f2, d2, f3, d3, sigma
    // Note that f3 = 1 - f1 - f2 so the number of fitted parameters is 6.
    // However the optimiser does not support compound constraints (f1+f2<=1)
    // so for convenience f3 is also a fitted parameter for the optimiser.
    final ObjectiveFunction fun;
    final GoalType goalType;
    if (mle) {
      fun = new ObjectiveFunction(new ThreeStateModelFunction(dt, dz, settings.a, settings.b,
          precision, distances, executor));
      goalType = GoalType.MAXIMIZE;
    } else {
      fun = new ObjectiveFunction(new ThreeStateSSFunction(dt, dz, settings.a, settings.b,
          precision, settings.binWidth, distances, executor));
      goalType = GoalType.MINIMIZE;
      best = Double.POSITIVE_INFINITY;
    }
    final SimpleBounds bounds =
        new SimpleBounds(addPrecision(new double[] {0, 0.001, 0, 0.1, 0, 0.1}, 0), addPrecision(
            new double[] {1, 0.1, 1, Double.POSITIVE_INFINITY, 1, Double.POSITIVE_INFINITY}, 0.1));
    final CustomPowellOptimizer.BasisStep step = new CustomPowellOptimizer.BasisStep(
        addPrecision(new double[] {0.1, 0.01, 0.1, 0.5, 0.1, 0.5}, 0.001));

    final Ticker ticker = ImageJUtils.createTicker(settings.repeats, 1);
    final double range = settings.maxD - settings.minD;
    for (int n = 1; n <= settings.repeats; n++) {
      try {
        // Guess in range
        final double[] start = addPrecision(new double[] {
            // Bound fraction
            rng.nextDouble(0.1, 0.4), rng.nextDouble(0.1),
            // Slow fraction
            rng.nextDouble(0.1, 0.4), rng.nextDouble(settings.minD, settings.minD + range / 3),
            // Fast fraction
            0, rng.nextDouble(settings.maxD - range / 3, settings.maxD)},
            // precision
            precision * rng.nextDouble(0.9, 1.1));
        // Set f3
        start[4] = 1 - start[0] - start[2];

        final PointValuePair solution =
            powellOptimizer.optimize(maxEval, fun, bounds, step, goalType, new InitialGuess(start));

        if (mle) {
          LoggerUtils.log(logger, Level.INFO,
              "Three-state fit [%d]: MLE = %s, Akaike IC = %s (%d evaluations)", n,
              solution.getValue(), MathUtils.getAkaikeInformationCriterion(solution.getValue(), 4),
              powellOptimizer.getEvaluations());
          if (solution.getValue() > best) {
            best = solution.getValue();
            result = solution;
          }
        } else {
          LoggerUtils.log(logger, Level.INFO, "Three-state fit [%d]: SS = %s (%d evaluations)", n,
              solution.getValue(), powellOptimizer.getEvaluations());
          if (solution.getValue() < best) {
            best = solution.getValue();
            result = solution;
          }
        }
      } catch (final TooManyEvaluationsException ex) {
        LoggerUtils.log(logger, Level.INFO,
            "Powell optimiser [%d] failed : Too many evaluation (%d)", n,
            powellOptimizer.getEvaluations());
      } catch (final TooManyIterationsException ex) {
        LoggerUtils.log(logger, Level.INFO,
            "Powell optimiser [%d] failed : Too many iterations (%d)", n,
            powellOptimizer.getIterations());
      } catch (final ConvergenceException ex) {
        LoggerUtils.log(logger, Level.INFO, "Powell optimiser [%d] failed to fit : %s", n,
            ex.getMessage());
      }
      ticker.tick();
    }

    if (result != null) {
      sort3state(result.getPointRef());
    }

    return result;
  }

  private static void sort3state(final double[] r) {
    final double sum = r[0] + r[2] + r[4];
    r[0] /= sum;
    r[2] /= sum;
    r[4] /= sum;
    // Sort (in event that the free fraction is slower than the bound)
    double f;
    double d;
    if (r[3] < r[1]) {
      f = r[0];
      d = r[1];
      r[0] = r[2];
      r[1] = r[3];
      r[2] = f;
      r[3] = d;
    }
    if (r[5] < r[3]) {
      f = r[2];
      d = r[3];
      r[2] = r[4];
      r[3] = r[5];
      r[4] = f;
      r[5] = d;
      if (r[3] < r[1]) {
        f = r[0];
        d = r[1];
        r[0] = r[2];
        r[1] = r[3];
        r[2] = f;
        r[3] = d;
      }
    }
  }

  private static CustomPowellOptimizer createCustomPowellOptimizer() {
    final double rel = 1e-8;
    final double abs = 1e-100;
    final boolean basisConvergence = false;
    return new CustomPowellOptimizer(rel, abs, null, basisConvergence);
  }

  private double[] addPrecision(double[] a, double p) {
    if (settings.fitPrecision) {
      final double[] b = Arrays.copyOf(a, a.length + 1);
      b[a.length] = p;
      return b;
    }
    return a;
  }

  private void addToResultTable(PointValuePair result) {
    createTable().append(addResult(result));
  }

  private TextWindow createTable() {
    return ImageJUtils.refresh(TABLE_REF, () -> {
      return new TextWindow(TITLE + " Analysis", createHeader(), "", 800, 300);
    });
  }

  private String createHeader() {
    return Arrays.stream(new String[] {"dz (nm)", "dt (ms)", "offsets", "precision (nm)",
        "fit precision", "a", "b", "group size", "repeats", "min D (um^2/s)", "max D (um^2/s)", "F",
        "D", "sigma (nm)", "Value"}).collect(Collectors.joining("\t"));
  }

  private String addResult(PointValuePair result) {
    final StringBuilder sb = new StringBuilder();
    //@formatter:off
    sb.append(MathUtils.rounded(settings.depthOfField)).append('\t')
      .append(MathUtils.rounded(exposureTime * 1e3)).append('\t')
      .append(settings.offsets).append('\t')
      .append(MathUtils.rounded(settings.precision)).append('\t')
      .append(settings.fitPrecision).append('\t')
      .append(MathUtils.rounded(settings.a)).append('\t')
      .append(MathUtils.rounded(settings.b)).append('\t')
      .append(settings.groupSize).append('\t')
      .append(settings.repeats).append('\t')
      .append(MathUtils.rounded(settings.minD)).append('\t')
      .append(MathUtils.rounded(settings.maxD))
      ;
    //@formatter:on
    if (result != null) {
      final double[] fit = result.getPointRef();
      final int n = fit.length / 2;
      final String f = IntStream.range(0, n).mapToObj(i -> MathUtils.rounded(fit[2 * i]))
          .collect(Collectors.joining(", "));
      final String d = IntStream.range(0, n).mapToObj(i -> MathUtils.rounded(fit[2 * i + 1]))
          .collect(Collectors.joining(", "));
      final double precision = (fit.length & 1) == 1 ? fit[2 * n] * 1e3 : settings.precision;
      sb.append('\t').append(f).append('\t').append(d).append('\t')
          .append(MathUtils.rounded(precision)).append('\t')
          .append(MathUtils.rounded(result.getValue()));
    }
    return sb.toString();
  }

  private float[][] createPdf(float[][] distances) {
    IJ.showStatus("Creating PDF...");

    final float binWidth = (float) settings.binWidth;
    final float[][] pdf = new float[distances.length][];
    int n = 0;
    for (int i = 0; i < distances.length; i++) {
      final float[] x = distances[i];
      pdf[i] = createPdf(binWidth, x);
      n = Math.max(n, pdf[i].length);
    }

    // Zero pad to max length
    for (int i = 0; i < distances.length; i++) {
      pdf[i] = Arrays.copyOf(pdf[i], n);
    }

    return pdf;
  }

  /**
   * Creates the probabilty density function.
   *
   * @param binWidth the bin width
   * @param distances the distances (must be sorted)
   * @return the pdf
   */
  private static float[] createPdf(float binWidth, float[] distances) {
    final IntArrayList counts = new IntArrayList();
    double max = binWidth;
    int bin = 1;
    int count = 0;
    for (final float d : distances) {
      if (d >= max) {
        counts.add(count);
        bin++;
        max = bin * binWidth;
        while (d >= max) {
          counts.add(0);
          bin++;
          max = bin * binWidth;
        }
        count = 1;
      } else {
        count++;
      }
    }
    counts.add(count);

    final double scale = 1.0 / distances.length;
    final float[] pdf = new float[counts.size()];
    for (int i = 0; i < pdf.length; i++) {
      pdf[i] = (float) (counts.getInt(i) * scale);
    }
    return pdf;
  }

  private void plotDistances(float[][] distances, float[][] pdf, PointValuePair result,
      ExecutorService executor) {
    IJ.showStatus("Plotting results...");
    String title = TITLE + " distance CDF";
    Plot plot = new Plot(title, "Distance (um)", "Probability");
    double maxD = 0.01;
    final LUT lut = LutHelper.createLut(LutColour.RED_MAGENTA_BLUE, false);
    for (int i = 0; i < distances.length; i++) {
      final float[] x = distances[i];
      plot.setColor(LutHelper.getColour(lut, distances.length - i, 1, distances.length));
      final float[] y = SimpleArrayUtils.newArray(x.length, 1f, 1f);
      final double scale = 1.0 / x.length;
      SimpleArrayUtils.apply(y, f -> (float) (f * scale));
      maxD = Math.max(maxD, x[x.length - 1]);
      plot.addPoints(x, y, Plot.LINE);
    }
    plot.setColor(Color.black);
    plot.setLimits(0, maxD, 0, 1.05);

    final WindowOrganiser wo = new WindowOrganiser();
    ImageJUtils.display(title, plot, 0, wo);

    title = TITLE + " distance PDF";
    plot = new Plot(title, "Distance (um)", "Probability");
    final float binWidth = (float) settings.binWidth;
    final float[] r = SimpleArrayUtils.newArray(pdf[0].length, binWidth / 2, binWidth);
    float maxP = 0;
    for (int i = 0; i < distances.length; i++) {
      final float[] y = pdf[i];
      plot.setColor(LutHelper.getColour(lut, distances.length - i, 1, distances.length));
      maxP = MathUtils.maxDefault(maxP, y);
      plot.addPoints(r, y, Plot.BAR);
    }
    plot.setColor(Color.black);
    plot.setLimits(0, maxD, 0, maxP * 1.05);

    ImageJUtils.display(title, plot, 0, wo);
    wo.tile();

    if (result == null) {
      return;
    }

    // Compute fitted PDF.
    // Evaluate the PDF using Simpson integration: use [a, (a+b)/2, b] for each bin [a, b].
    final float[][] evalDistances = new float[distances.length][];
    Arrays.fill(evalDistances, SimpleArrayUtils.newArray(r.length * 2 + 1, 0, binWidth * 0.5f));

    final double dt = exposureTime;
    final double dz = settings.depthOfField / 1000;
    final double precision = settings.precision / 1000;
    double[][] p;
    final double[] fit = result.getPointRef();
    if (fit.length < 6) {
      p = new TwoStateModelFunction(dt, dz, settings.a, settings.b, precision, evalDistances,
          executor).evaluate(fit);
    } else {
      p = new ThreeStateModelFunction(dt, dz, settings.a, settings.b, precision, evalDistances,
          executor).evaluate(fit);
    }

    for (int i = 0; i < evalDistances.length; i++) {
      final double[] pi = p[i];
      final int n = pi.length / 2;
      plot.setColor(LutHelper.getColour(lut, distances.length - i, 1, distances.length));
      final float[] y = new float[n];
      for (int j = 0; j < y.length; j++) {
        final int index = j * 2;
        y[j] = (float) ((pi[index] + 4 * pi[index + 1] + pi[index + 2]) / 6);
      }
      // The computed probability will be lower for increasing time-delay
      // as the diffusing molecules are lost. Normalise to sum to 1 so
      // the curve matches the empirical PDF (which sums to 1).
      final double s = 1.0 / MathUtils.sum(y);
      SimpleArrayUtils.apply(y, f -> (float) (f * s));

      plot.addPoints(r, y, Plot.CIRCLE);
    }
  }

  /**
   * Creates the weights for each time delay. These are used to evenly balance the different number
   * of observations for each time delay. The weight for is maxN / N where N is the number of
   * observations for the time delay.
   *
   * @param distances the distances
   * @return the weights
   */
  private static double[] createWeights(float[][] distances) {
    final int[] sizes = Arrays.stream(distances).mapToInt(x -> x.length).toArray();
    final double maxSize = Arrays.stream(sizes).max().getAsInt();
    return Arrays.stream(sizes).mapToDouble(size -> maxSize / size).toArray();
  }

  // Implementation Note:
  //
  // The functions assumes d1 is the bound fraction.
  // If the fractions are sorted before evaluation then this
  // creates a non-smooth function at the transition which can have
  // poor fitting behaviour for an optimiser that is
  // following a optimisation direction.
  // It is assumed that if d values are initialised sufficiently
  // far apart with d1 << d2 << d3 then the populations will
  // be correctly fit. An option is to also compute p_remaining
  // for the bound fraction. In practice this is always close to 1
  // when d1 is bound. It does avoid issues if the d values overlap.

  /**
   * Function to evaluate the log-likelihood of the observed distances for the two-state diffusion
   * model.
   *
   * <p>This function is slow and can use an {@link ExecutorService} for parallel evaluation.
   */
  static class TwoStateModelFunction implements MultivariateFunction {
    double dt;
    double dz;
    double a;
    double b;
    double precision;
    float[][] distances;
    ExecutorService executor;
    double[] weights;

    TwoStateModelFunction(double dt, double dz, double a, double b, double precision,
        float[][] distances, ExecutorService executor) {
      this.dt = dt;
      this.dz = dz;
      this.a = a;
      this.b = b;
      this.precision = precision;
      this.distances = distances;
      this.executor = executor;
      weights = createWeights(distances);
    }

    @Override
    public double value(double[] point) {
      final double sum = point[0] + point[2];
      final double f1 = point[0] / sum;
      final double d1 = point[1];
      final double f2 = point[2] / sum;
      final double d2 = point[3];
      final double sigma = point.length > 4 ? point[4] : precision;

      DD ll = DD.ZERO;

      if (executor == null) {
        for (int time = 0; time < distances.length; time++) {
          ll = ll.add(compute(f1, d1, f2, d2, sigma, time));
        }
        return ll.doubleValue();
      }

      final List<Future<DD>> futures = new LinkedList<>();

      for (int n = distances.length; --n >= 0;) {
        final int time = n;
        futures.add(executor.submit(() -> compute(f1, d1, f2, d2, sigma, time)));
      }

      // Finish processing data
      for (final Future<DD> f : futures) {
        try {
          ll = ll.add(f.get());
        } catch (InterruptedException | ExecutionException e) {
          throw new RuntimeException(e);
        }
      }
      return ll.doubleValue();
    }

    private DD compute(double f1, double d1, double f2, double d2, double sigma, int time) {
      DD ll = DD.ZERO;
      final float[] d = distances[time];
      final double deltaT = dt * (time + 1);
      // Use corrected depth-of-field
      final double p2 = DiffusionAnalysis.remaining(deltaT, dz + a * Math.sqrt(d2) + b, d2);
      final double denom1 = 1.0 / (4 * (d1 * deltaT + sigma * sigma));
      final double denom2 = 1.0 / (4 * (d2 * deltaT + sigma * sigma));
      for (final float r : d) {
        final double p = f1 * (r * 2 * denom1) * Math.exp(-r * r * denom1)
            + p2 * f2 * (r * 2 * denom2) * Math.exp(-r * r * denom2);
        ll = ll.add(Math.log(p));
      }
      return ll.multiply(weights[time]);
    }

    double[][] evaluate(double[] point) {
      final double sum = point[0] + point[2];
      final double f1 = point[0] / sum;
      final double d1 = point[1];
      final double f2 = point[2] / sum;
      final double d2 = point[3];
      final double sigma = point.length > 4 ? point[4] : precision;

      final double[][] p = new double[distances.length][];

      if (executor == null) {
        for (int time = 0; time < distances.length; time++) {
          p[time] = evaluate(f1, d1, f2, d2, sigma, time);
        }
        return p;
      }

      final List<Future<double[]>> futures = new LinkedList<>();

      for (int n = 0; n < distances.length; n++) {
        final int time = n;
        futures.add(executor.submit(() -> evaluate(f1, d1, f2, d2, sigma, time)));
      }

      // Finish processing data
      int n = 0;
      for (final Future<double[]> f : futures) {
        try {
          p[n++] = f.get();
        } catch (InterruptedException | ExecutionException e) {
          throw new RuntimeException(e);
        }
      }
      return p;
    }

    private double[] evaluate(double f1, double d1, double f2, double d2, double sigma, int time) {
      final float[] d = distances[time];
      final double deltaT = dt * (time + 1);
      // Use corrected depth-of-field
      final double p2 = DiffusionAnalysis.remaining(deltaT, dz + a * Math.sqrt(d2) + b, d2);
      final double denom1 = 1.0 / (4 * (d1 * deltaT + sigma * sigma));
      final double denom2 = 1.0 / (4 * (d2 * deltaT + sigma * sigma));
      final double[] p = new double[d.length];
      for (int i = 0; i < d.length; i++) {
        final double r = d[i];
        p[i] = f1 * (r * 2 * denom1) * Math.exp(-r * r * denom1)
            + p2 * f2 * (r * 2 * denom2) * Math.exp(-r * r * denom2);
      }
      return p;
    }
  }

  /**
   * Function to evaluate the log-likelihood of the observed distances for the three-state diffusion
   * model.
   *
   * <p>This function is slow and can use an {@link ExecutorService} for parallel evaluation.
   */
  static class ThreeStateModelFunction implements MultivariateFunction {
    double dt;
    double dz;
    double a;
    double b;
    double precision;

    float[][] distances;
    ExecutorService executor;
    double[] weights;

    ThreeStateModelFunction(double dt, double dz, double a, double b, double precision,
        float[][] distances, ExecutorService executor) {
      this.dt = dt;
      this.dz = dz;
      this.a = a;
      this.b = b;
      this.precision = precision;
      this.distances = distances;
      this.executor = executor;
      weights = createWeights(distances);
    }

    @Override
    public double value(double[] point) {
      final double sum = point[0] + point[2] + point[4];
      final double f1 = point[0] / sum;
      final double d1 = point[1];
      final double f2 = point[2] / sum;
      final double d2 = point[3];
      final double f3 = point[4] / sum;
      final double d3 = point[5];
      final double sigma = point.length > 6 ? point[6] : precision;

      DD ll = DD.ZERO;

      if (executor == null) {
        for (int time = 0; time < distances.length; time++) {
          ll = ll.add(compute(f1, d1, f2, d2, f3, d3, sigma, time));
        }
        return ll.doubleValue();
      }

      final List<Future<DD>> futures = new LinkedList<>();

      for (int n = distances.length; --n >= 0;) {
        final int time = n;
        futures.add(executor.submit(() -> compute(f1, d1, f2, d2, f3, d3, sigma, time)));
      }

      // Finish processing data
      for (final Future<DD> f : futures) {
        try {
          ll = ll.add(f.get());
        } catch (InterruptedException | ExecutionException e) {
          throw new RuntimeException(e);
        }
      }
      return ll.doubleValue();
    }

    private DD compute(double f1, double d1, double f2, double d2, double f3, double d3,
        double sigma, int time) {
      DD ll = DD.ZERO;
      final float[] d = distances[time];
      final double deltaT = dt * (time + 1);
      // Use corrected depth-of-field
      final double p2 = DiffusionAnalysis.remaining(deltaT, dz + a * Math.sqrt(d2) + b, d2);
      final double p3 = DiffusionAnalysis.remaining(deltaT, dz + a * Math.sqrt(d3) + b, d3);
      final double denom1 = 1.0 / (4 * (d1 * deltaT + sigma * sigma));
      final double denom2 = 1.0 / (4 * (d2 * deltaT + sigma * sigma));
      final double denom3 = 1.0 / (4 * (d3 * deltaT + sigma * sigma));
      for (final float r : d) {
        final double p = f1 * (r * 2 * denom1) * Math.exp(-r * r * denom1)
            + p2 * f2 * (r * 2 * denom2) * Math.exp(-r * r * denom2)
            + p3 * f3 * (r * 2 * denom3) * Math.exp(-r * r * denom3);
        ll = ll.add(Math.log(p));
      }
      return ll.multiply(weights[time]);
    }

    double[][] evaluate(double[] point) {
      final double sum = point[0] + point[2] + point[4];
      final double f1 = point[0] / sum;
      final double d1 = point[1];
      final double f2 = point[2] / sum;
      final double d2 = point[3];
      final double f3 = point[4] / sum;
      final double d3 = point[5];
      final double sigma = point.length > 6 ? point[6] : precision;

      final double[][] p = new double[distances.length][];

      if (executor == null) {
        for (int time = 0; time < distances.length; time++) {
          p[time] = evaluate(f1, d1, f2, d2, f3, d3, sigma, time);
        }
        return p;
      }

      final List<Future<double[]>> futures = new LinkedList<>();

      for (int n = 0; n < distances.length; n++) {
        final int time = n;
        futures.add(executor.submit(() -> evaluate(f1, d1, f2, d2, f3, d3, sigma, time)));
      }

      // Finish processing data
      int n = 0;
      for (final Future<double[]> f : futures) {
        try {
          p[n++] = f.get();
        } catch (InterruptedException | ExecutionException e) {
          throw new RuntimeException(e);
        }
      }
      return p;
    }

    private double[] evaluate(double f1, double d1, double f2, double d2, double f3, double d3,
        double sigma, int time) {
      final float[] d = distances[time];
      final double deltaT = dt * (time + 1);
      // Use corrected depth-of-field
      final double p2 = DiffusionAnalysis.remaining(deltaT, dz + a * Math.sqrt(d2) + b, d2);
      final double p3 = DiffusionAnalysis.remaining(deltaT, dz + a * Math.sqrt(d3) + b, d3);
      final double denom1 = 1.0 / (4 * (d1 * deltaT + sigma * sigma));
      final double denom2 = 1.0 / (4 * (d2 * deltaT + sigma * sigma));
      final double denom3 = 1.0 / (4 * (d3 * deltaT + sigma * sigma));
      final double[] p = new double[d.length];
      for (int i = 0; i < d.length; i++) {
        final double r = d[i];
        p[i] = f1 * (r * 2 * denom1) * Math.exp(-r * r * denom1)
            + p2 * f2 * (r * 2 * denom2) * Math.exp(-r * r * denom2)
            + p3 * f3 * (r * 2 * denom3) * Math.exp(-r * r * denom3);
      }
      return p;
    }
  }

  /**
   * Function to evaluate the sum-of-squares of the observed PDF for the two-state diffusion model.
   *
   * <p>This function is slow and can use an {@link ExecutorService} for parallel evaluation.
   */
  static class TwoStateSSFunction implements MultivariateFunction {
    double dt;
    double dz;
    double a;
    double b;
    double precision;
    double[] distances;
    float[][] pdf;
    ExecutorService executor;

    TwoStateSSFunction(double dt, double dz, double a, double b, double precision, double dr,
        float[][] pdf, ExecutorService executor) {
      this.dt = dt;
      this.dz = dz;
      this.a = a;
      this.b = b;
      this.precision = precision;
      this.distances = SimpleArrayUtils.newArray(pdf[0].length * 2 + 1, 0, dr * 0.5);
      this.pdf = pdf;
      this.executor = executor;
    }

    @Override
    public double value(double[] point) {
      final double sum = point[0] + point[2];
      final double f1 = point[0] / sum;
      final double d1 = point[1];
      final double f2 = point[2] / sum;
      final double d2 = point[3];
      final double sigma = point.length > 4 ? point[4] : precision;

      double ss = 0;

      if (executor == null) {
        for (int time = 0; time < pdf.length; time++) {
          ss += compute(f1, d1, f2, d2, sigma, time);
        }
        return ss;
      }

      final List<Future<Double>> futures = new LinkedList<>();

      for (int n = pdf.length; --n >= 0;) {
        final int time = n;
        futures.add(executor.submit(() -> compute(f1, d1, f2, d2, sigma, time)));
      }

      // Finish processing data
      for (final Future<Double> f : futures) {
        try {
          ss += f.get();
        } catch (InterruptedException | ExecutionException e) {
          throw new RuntimeException(e);
        }
      }
      return ss;
    }

    private double compute(double f1, double d1, double f2, double d2, double sigma, int time) {
      double ss = 0;
      final double deltaT = dt * (time + 1);
      // Use corrected depth-of-field
      final double p2 = DiffusionAnalysis.remaining(deltaT, dz + a * Math.sqrt(d2) + b, d2);
      final double denom1 = 1.0 / (4 * (d1 * deltaT + sigma * sigma));
      final double denom2 = 1.0 / (4 * (d2 * deltaT + sigma * sigma));
      final double[] d = distances;
      final double[] p = new double[d.length];
      for (int i = 0; i < p.length; i++) {
        final double r = d[i];
        p[i] = f1 * (r * 2 * denom1) * Math.exp(-r * r * denom1)
            + p2 * f2 * (r * 2 * denom2) * Math.exp(-r * r * denom2);
      }
      normalise(p);
      final float[] obs = pdf[time];
      for (int i = 0; i < obs.length; i++) {
        final int index = i << 1;
        // Simpson integration using points [a, (a+b)/2, b] for bin [a, b]
        final double dp = (p[index] + 4 * p[index + 1] + p[index + 2]) / 6 - obs[i];
        ss += dp * dp;
      }
      return ss;
    }

    /**
     * Normalise to sum to 1 so the curve matches the empirical PDF (which sums to 1).
     *
     * @param p the pdf
     */
    static void normalise(double[] p) {
      // Use composite Simpson's 1/3 rule integration
      double s4 = 0;
      double s2 = 0;
      for (int i = p.length - 2; i >= 1; i -= 2) {
        s4 += p[i];
      }
      for (int i = p.length - 3; i >= 2; i -= 2) {
        s2 += p[i];
      }
      final double sum = p[0] + 4 * s4 + 2 * s2 + p[p.length - 1];
      final double inverseSum = 6 / sum;
      for (int i = 0; i < p.length; i++) {
        p[i] *= inverseSum;
      }
    }
  }

  /**
   * Function to evaluate the sum-of-squares of the observed PDF for the two-state diffusion model.
   *
   * <p>This function is slow and can use an {@link ExecutorService} for parallel evaluation.
   */
  static class ThreeStateSSFunction implements MultivariateFunction {
    double dt;
    double dz;
    double a;
    double b;
    double precision;
    double[] distances;
    float[][] pdf;
    ExecutorService executor;

    ThreeStateSSFunction(double dt, double dz, double a, double b, double precision, double dr,
        float[][] pdf, ExecutorService executor) {
      this.dt = dt;
      this.dz = dz;
      this.a = a;
      this.b = b;
      this.precision = precision;
      this.distances = SimpleArrayUtils.newArray(pdf[0].length * 2 + 1, 0, dr * 0.5);
      this.pdf = pdf;
      this.executor = executor;
    }

    @Override
    public double value(double[] point) {
      final double sum = point[0] + point[2] + point[4];
      final double f1 = point[0] / sum;
      final double d1 = point[1];
      final double f2 = point[2] / sum;
      final double d2 = point[3];
      final double f3 = point[4] / sum;
      final double d3 = point[5];
      final double sigma = point.length > 6 ? point[6] : precision;

      double ss = 0;

      if (executor == null) {
        for (int time = 0; time < pdf.length; time++) {
          ss += compute(f1, d1, f2, d2, f3, d3, sigma, time);
        }
        return ss;
      }

      final List<Future<Double>> futures = new LinkedList<>();

      for (int n = pdf.length; --n >= 0;) {
        final int time = n;
        futures.add(executor.submit(() -> compute(f1, d1, f2, d2, f3, d3, sigma, time)));
      }

      // Finish processing data
      for (final Future<Double> f : futures) {
        try {
          ss += f.get();
        } catch (InterruptedException | ExecutionException e) {
          throw new RuntimeException(e);
        }
      }
      return ss;
    }

    private double compute(double f1, double d1, double f2, double d2, double f3, double d3,
        double sigma, int time) {
      double ss = 0;
      final double deltaT = dt * (time + 1);
      // Use corrected depth-of-field
      final double p2 = DiffusionAnalysis.remaining(deltaT, dz + a * Math.sqrt(d2) + b, d2);
      final double p3 = DiffusionAnalysis.remaining(deltaT, dz + a * Math.sqrt(d3) + b, d3);
      final double denom1 = 1.0 / (4 * (d1 * deltaT + sigma * sigma));
      final double denom2 = 1.0 / (4 * (d2 * deltaT + sigma * sigma));
      final double denom3 = 1.0 / (4 * (d3 * deltaT + sigma * sigma));
      final double[] d = distances;
      final double[] p = new double[d.length];
      for (int i = 0; i < p.length; i++) {
        final double r = d[i];
        p[i] = f1 * (r * 2 * denom1) * Math.exp(-r * r * denom1)
            + p2 * f2 * (r * 2 * denom2) * Math.exp(-r * r * denom2)
            + p3 * f3 * (r * 2 * denom3) * Math.exp(-r * r * denom3);
      }
      TwoStateSSFunction.normalise(p);
      final float[] obs = pdf[time];
      for (int i = 0; i < obs.length; i++) {
        final int index = i << 1;
        // Simpson integration using points [a, (a+b)/2, b] for bin [a, b]
        final double dp = (p[index] + 4 * p[index + 1] + p[index + 2]) / 6 - obs[i];
        ss += dp * dp;
      }
      return ss;
    }
  }
}
