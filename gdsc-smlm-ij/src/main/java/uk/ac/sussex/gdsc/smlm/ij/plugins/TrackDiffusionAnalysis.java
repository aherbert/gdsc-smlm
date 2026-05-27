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
import java.util.stream.IntStream;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.rng.JumpableUniformRandomProvider;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import org.apache.commons.rng.simple.RandomSource;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.MultiDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper.LutColour;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.fitting.DiffusionAnalysis;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.MemoryResultsList;
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

  private double exposureTime;

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

    Settings() {
      selected = new LocalList<>();

      depthOfField = 750;
      gap = 1;
      maxT = 7;
      offsets = 0;

      binWidth = 0.01;

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
    final int threadCount = Prefs.getThreads();
    final ExecutorService executor = Executors.newFixedThreadPool(threadCount);

    try {

      final float[][] distances = getDistances(allResults);

      plotDistances(distances);

      // final double[][] probability = simulateRemaining(threadCount, executor);
      //
      // // Fit the observed probability
      // final double[] ab = fitRemaining(probability, threadCount, executor);
      //
      // addToResultTable(ab);
      //
      // final double[][] fittedProbability = createProbability(ab);
      //
      // plotRemaining(probability, fittedProbability);

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
      // Uncomment to allow rendering the dataset (which requires pixel units)
      // results.convertToPreferredUnits();
      MemoryPeakResults.addResults(results);

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

    gd.addNumericField("Bin_width", settings.binWidth, -3);

    gd.addNumericField("A", settings.a, 4);
    gd.addNumericField("B", settings.b, 4);

    gd.addSlider("Repeats", 0, 5, settings.repeats);
    gd.addNumericField("Min_D", settings.minD, 3, 6, "um^2/s");
    gd.addNumericField("Max_D", settings.maxD, 3, 6, "um^2/s");

    gd.addHelp(HelpUrls.getUrl("track-diffusion-analysis"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    settings.depthOfField = gd.getNextNumber();
    settings.maxT = (int) gd.getNextNumber();
    settings.offsets = (int) gd.getNextNumber();

    settings.binWidth = gd.getNextNumber();

    settings.a = gd.getNextNumber();
    settings.b = gd.getNextNumber();

    settings.repeats = (int) gd.getNextNumber();
    settings.minD = gd.getNextNumber();
    settings.maxD = gd.getNextNumber();

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Depth of Field", settings.depthOfField);
      ParameterUtils.isAboveZero("Max T", settings.maxT);
      ParameterUtils.isPositive("Offsets", settings.offsets);

      ParameterUtils.isAboveZero("Bin width", settings.binWidth);

      ParameterUtils.isPositive("A", settings.a);
      ParameterUtils.isPositive("B", settings.b);

      ParameterUtils.isPositive("Repeats", settings.repeats);
      ParameterUtils.isAboveZero("Min D", settings.minD);
      ParameterUtils.isAboveZero("Max D", settings.maxD);

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

    final int[] sizes = IntStream.rangeClosed(1, maxT).map(i -> distances[i].size()).toArray();
    ImageJUtils.log("Distance counts: %s", Arrays.toString(sizes));

    // Convert to sorted arrays
    return IntStream.rangeClosed(1, maxT).mapToObj(i -> {
      final float[] data = distances[i].toFloatArray();
      Arrays.sort(data);
      return data;
    }).toArray(float[][]::new);
  }

  private void plotDistances(float[][] distances) {
    String title = TITLE + " distance CDF";
    Plot plot = new Plot(title, "Distance (um)", "Probability");
    double maxD = 0.01;
    final LUT lut = LutHelper.createLut(LutColour.RED_MAGENTA_BLUE, false);
    for (int i = 0; i < distances.length; i++) {
      final float[] x = distances[i];
      if (x.length == 0) {
        continue;
      }
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
    float maxP = 0;
    for (int i = 0; i < distances.length; i++) {
      final float[] x = distances[i];
      if (x.length == 0) {
        continue;
      }
      plot.setColor(LutHelper.getColour(lut, distances.length - i, 1, distances.length));
      final float[] y = createPdf(binWidth, x);
      final float[] r = SimpleArrayUtils.newArray(x.length, binWidth / 2, binWidth);
      maxP = MathUtils.maxDefault(maxP, y);
      plot.addPoints(r, y, Plot.BAR);
    }
    plot.setColor(Color.black);
    plot.setLimits(0, maxD, 0, maxP * 1.05);

    ImageJUtils.display(title, plot, 0, wo);
    wo.tile();
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

  /**
   * Function to evaluate the squared distance between the simulated probability and the computed
   * probability using P(dt, dz_corr, D) where z_corr is dz + a * sqrt(D) + b.
   *
   * <p>This function is slow and can use an {@link ExecutorService} for parallel evaluation.
   */
  static class DepthOfFieldFunction implements MultivariateFunction {
    double dt;
    double dz;
    double[] diffusionCoefficients;
    double[][] probability;
    int threadCount;
    ExecutorService executor;

    DepthOfFieldFunction(double dt, double dz, double[] diffusionCoefficients,
        double[][] probability, int threadCount, ExecutorService executor) {
      this.dt = dt;
      this.dz = dz;
      this.diffusionCoefficients = diffusionCoefficients;
      this.probability = probability;
      this.threadCount = threadCount;
      this.executor = executor;
    }

    @Override
    public double value(double[] point) {
      final double a = point[0];
      final double b = point[1];
      if (threadCount < 2 || executor == null) {
        double ss = 0;
        for (int i = 0; i < probability.length; i++) {
          final double d = diffusionCoefficients[i];
          for (int j = 0; j < probability[i].length; j++) {
            final double dz_corr = dz + a * Math.sqrt(d) + b;
            final double p = DiffusionAnalysis.remaining(dt * (j + 1), dz_corr, d);
            final double dp = probability[i][j] - p;
            ss += dp * dp;
          }
        }
        return ss;
      }

      final List<Future<Double>> futures = new LinkedList<>();
      final int maxj = probability[0].length;
      final AtomicInteger position = new AtomicInteger(probability.length * maxj);

      for (int n = 0; n < threadCount; n++) {
        futures.add(executor.submit(() -> {
          double ss = 0;
          for (;;) {
            final int p = position.decrementAndGet();
            if (p < 0) {
              break;
            }
            final int i = p / maxj;
            final int j = p % maxj;

            final double d = diffusionCoefficients[i];
            final double dz_corr = dz + a * Math.sqrt(d) + b;
            final double dp =
                probability[i][j] - DiffusionAnalysis.remaining(dt * (j + 1), dz_corr, d);
            ss += dp * dp;
          }
          return ss;
        }));
      }

      // Finish processing data
      double ss = 0;
      for (final Future<Double> f : futures) {
        try {
          ss += f.get();
        } catch (InterruptedException | ExecutionException e) {
          throw new RuntimeException(e);
        }
      }
      return ss;
    }
  }
}
