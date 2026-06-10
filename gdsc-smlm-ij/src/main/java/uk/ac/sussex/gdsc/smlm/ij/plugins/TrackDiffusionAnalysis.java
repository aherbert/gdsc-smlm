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
import ij.gui.PlotWindow;
import ij.plugin.PlugIn;
import ij.process.LUT;
import ij.text.TextWindow;
import it.unimi.dsi.fastutil.floats.FloatArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import java.awt.Color;
import java.awt.Font;
import java.awt.TextField;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.exception.TooManyIterationsException;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.rng.JumpableUniformRandomProvider;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.NormalizedGaussianSampler;
import org.apache.commons.rng.simple.RandomSource;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
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
import uk.ac.sussex.gdsc.core.utils.StoredData;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.rng.SamplerUtils;
import uk.ac.sussex.gdsc.core.utils.rng.UniformRandomProviders;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.fitting.DiffusionAnalysis;
import uk.ac.sussex.gdsc.smlm.function.ChiSquaredDistributionTable;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.MemoryResultsList;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.TrackDiffusionAnalysisSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.math3.optim.nonlinear.scalar.noderiv.CustomPowellOptimizer;
import uk.ac.sussex.gdsc.smlm.results.IdPeakResult;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.count.FrameCounter;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedure;
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
  private TrackDiffusionAnalysisSettings.Builder settings;

  /**
   * Store the localisation position.
   */
  private static class Position {
    int t;
    float x;
    float y;

    Position(int t, float x, float y) {
      this.t = t;
      this.x = x;
      this.y = y;
    }
  }

  /**
   * Store the localisation with the z position.
   */
  private static class XyzPeakResult extends IdPeakResult {
    private static final long serialVersionUID = 20260605;

    XyzPeakResult(int t, double x, double y, double z, int id) {
      // Fixed intensity of 1
      super(t, (float) x, (float) y, 1f, id);
      setZPosition((float) z);
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

    settings = SettingsManager.readTrackDiffusionAnalysisSettings(0).toBuilder();

    final List<MemoryPeakResults> allResults = new LocalList<>();

    if (!showDialog() || !showMultiDialog(allResults, items)) {
      return;
    }

    logger = ImageJPluginLoggerHelper.getLogger(getClass());

    final float[][] distances = getDistances(allResults);
    final int[][] counts = createCounts(distances);
    final float[][] pdf = createPdf(counts);

    final int threadCount = Prefs.getThreads();
    final ExecutorService executor = Executors.newFixedThreadPool(threadCount);

    try {
      final PointValuePair result = fitDistances(counts, pdf, executor, settings.getFitMle());
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

      final double halfDz = settings.getDepthOfField() / 2000;
      final double dt = settings.getExposureTime() / 1000;
      final int g = settings.getGap();
      final int maxT = settings.getSimT();

      // Sample populations
      final double p1 = settings.getF1();
      final double p2 = settings.getF2() > 0 ? p1 + settings.getF2() : 1;

      // Precision in um (used for the units of the tracks)
      final double precision = settings.getPrecision() / 1000;

      // Simulate tracks across the depth of field.
      // Each track is scaled using the diffusion coefficient and the molecule tested if
      // it remains in the depth-of-field.

      final double[] sampleD = {settings.getD1(), settings.getD2(), settings.getD3()};
      // Create the Gaussian step for each diffusion coefficient
      final double[] step = Arrays.stream(sampleD).map(
          // Std.dev of Gaussian for the step size: sqrt(2D * dt)
          d -> Math.sqrt(2 * d * dt)).toArray();

      final List<Future<MemoryPeakResults>> futures = new LocalList<>(threadCount);

      final int total = settings.getNumberOfMolecules();
      final Ticker ticker = ImageJUtils.createTicker(total, threadCount);

      final AtomicInteger position = new AtomicInteger();
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
            final int pos = position.incrementAndGet();
            if (pos > total) {
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
            final int id = nextId.incrementAndGet();
            double x = 0;
            double y = 0;
            double z = rng.nextDouble(-halfDz, halfDz);
            int last = 0;
            // Ensure all tracks have unique frames to allow tracing the dataset
            // Each track should be separated by maxT frames so tracing even with
            // small frame gaps will not join stationary molecules at the origin.
            final int frame = pos * (2 * maxT);
            // Record the origin
            observed.add(new XyzPeakResult(frame, x + sampler.sample() * precision,
                y + sampler.sample() * precision, z, id));
            for (int i = 1; i <= maxT; i++) {
              x += sampler.sample() * d;
              y += sampler.sample() * d;
              z += sampler.sample() * d;
              if (Math.abs(z) < halfDz) {
                // Record molecule
                // Check the gap from the last time it was in the depth-of-field
                if (i - last > g) {
                  break;
                  // Start a new track
                  // id = nextId.incrementAndGet();
                }
                last = i;
                // Add localisation precision
                observed.add(new XyzPeakResult(frame + i, x + sampler.sample() * precision,
                    y + sampler.sample() * precision, z, id));
              }
            }
            ticker.tick();
          }
          return observed;
        }));
      }

      final MemoryPeakResults results = new MemoryPeakResults(total);
      final CalibrationWriter cw = new CalibrationWriter();
      cw.setExposureTime(settings.getExposureTime());
      cw.setNmPerPixel(100);
      cw.setDistanceUnit(DistanceUnit.UM);
      cw.setIntensityUnit(IntensityUnit.COUNT);
      results.setCalibration(cw.getCalibration());
      results.setName(TITLE);
      MemoryPeakResults.addResults(results);

      // Uncomment to allow rendering the dataset (which requires pixel units)
      // results.convertToPreferredUnits();

      // Build final results in memory
      for (final Future<MemoryPeakResults> f : futures) {
        try {
          results.add(f.get());
        } catch (InterruptedException | ExecutionException e) {
          throw new RuntimeException(e);
        }
      }

      final StoredData data = new StoredData();
      // Already sorted by ID. Store each start z.
      final FrameCounter c = new FrameCounter();
      results.forEach((PeakResultProcedure) peak -> {
        if (c.advance(peak.getId())) {
          data.add(peak.getZPosition());
        }
      });
      new HistogramPlotBuilder("Track Origin", data, "z (um)").show();

    } finally {
      executor.shutdown();
    }

    ImageJUtils.finished();
  }

  private boolean showSimulationDialog() {
    settings = SettingsManager.readTrackDiffusionAnalysisSettings(0).toBuilder();

    final GenericDialog gd = new GenericDialog(TITLE);

    gd.addNumericField("Depth_of_field", settings.getDepthOfField(), 1, 6, "nm");
    gd.addNumericField("Exposure_time", settings.getExposureTime(), 1, 6, "ms");
    gd.addSlider("Gap", 1, 5, settings.getGap());

    gd.addNumericField("Number_of_molecules", settings.getNumberOfMolecules());
    gd.addSlider("Max_t", 5, 15, settings.getSimT());

    gd.addNumericField("Precision", settings.getPrecision(), 1, 6, "nm");
    gd.addNumericField("F1", settings.getF1(), 3);
    gd.addNumericField("F2", settings.getF2(), 3);
    gd.addNumericField("D1", settings.getD1(), 3, 6, "um^2/s");
    gd.addNumericField("D2", settings.getD2(), 3, 6, "um^2/s");
    gd.addNumericField("D3", settings.getD3(), 3, 6, "um^2/s");

    gd.addHelp(HelpUrls.getUrl("diffusion-depth-of-field"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    settings.setDepthOfField(gd.getNextNumber());
    settings.setExposureTime(gd.getNextNumber());
    settings.setGap((int) gd.getNextNumber());

    settings.setNumberOfMolecules((int) gd.getNextNumber());
    settings.setSimT((int) gd.getNextNumber());

    settings.setPrecision(gd.getNextNumber());
    settings.setF1(gd.getNextNumber());
    settings.setF2(gd.getNextNumber());
    settings.setD1(gd.getNextNumber());
    settings.setD2(gd.getNextNumber());
    settings.setD3(gd.getNextNumber());

    SettingsManager.writeSettings(settings.build());

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Depth of Field", settings.getDepthOfField());
      ParameterUtils.isAboveZero("Exposure time", settings.getExposureTime());
      ParameterUtils.isAboveZero("Gap", settings.getGap());

      ParameterUtils.isAboveZero("Number of molecules", settings.getNumberOfMolecules());
      ParameterUtils.isAboveZero("Max T", settings.getSimT());

      ParameterUtils.isPositive("Precision", settings.getPrecision());
      ParameterUtils.isAboveZero("F1", settings.getF1());
      ParameterUtils.isPositive("F2", settings.getF2());
      ParameterUtils.isPositive("D1", settings.getD1());
      ParameterUtils.isPositive("D2", settings.getD2());
      ParameterUtils.isPositive("D3", settings.getD3());

      ParameterUtils.isBelow("F1 + F2", settings.getF1() + settings.getF2(), 1);

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
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

    final TextField tfZ =
        gd.addAndGetNumericField("Depth_of_field", settings.getDepthOfField(), 1, 6, "nm");
    gd.addSlider("Max_t", 3, 10, settings.getMaxT());
    gd.addSlider("Offsets", 0, 5, settings.getOffsets());
    gd.addNumericField("Precision", settings.getPrecision(), 1, 6, "nm");

    gd.addNumericField("Bin_width", settings.getBinWidth(), -3);

    gd.addMessage("z_corr = z + a * sqrt(D) + b");
    final TextField tfA = gd.addAndGetNumericField("A", settings.getA(), 4);
    final TextField tfB = gd.addAndGetNumericField("B", settings.getB(), 4);
    gd.addButton("Calibrate", e -> {
      // Run on non GUI thread and while running disable this plugin dialog
      gd.setEnabled(false);
      ForkJoinPool.commonPool().execute(() -> {
        try {
          final double dz = Double.parseDouble(tfZ.getText());
          final double[] ab = new DiffusionDepthOfField().run(dz);
          if (ab != null) {
            tfA.setText(MathUtils.rounded(ab[0]));
            tfB.setText(MathUtils.rounded(ab[1]));
          }
        } finally {
          gd.setEnabled(true);
        }
      });
    });

    gd.addSlider("Repeats", 1, 5, settings.getRepeats());
    gd.addNumericField("Min_D", settings.getMinD(), 3, 6, "um^2/s");
    gd.addNumericField("Max_D", settings.getMaxD(), 3, 6, "um^2/s");

    gd.addCheckbox("Fit_precision", settings.getFitPrecision());
    gd.addCheckbox("Fit_three_state", settings.getFitThreeState());
    gd.addCheckbox("Fit_MLE", settings.getFitMle());
    gd.addNumericField("significanceLevel", settings.getSignificanceLevel(), -3);

    gd.addCheckbox("Show_CDF", settings.getShowCdf());
    gd.addCheckbox("Seperate_plots", settings.getShowSeparatePlots());

    gd.addHelp(HelpUrls.getUrl("track-diffusion-analysis"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    settings.setDepthOfField(gd.getNextNumber());
    settings.setMaxT((int) gd.getNextNumber());
    settings.setOffsets((int) gd.getNextNumber());
    settings.setPrecision(gd.getNextNumber());

    settings.setBinWidth(gd.getNextNumber());

    settings.setA(gd.getNextNumber());
    settings.setB(gd.getNextNumber());

    settings.setRepeats((int) gd.getNextNumber());
    settings.setMinD(gd.getNextNumber());
    settings.setMaxD(gd.getNextNumber());

    settings.setFitPrecision(gd.getNextBoolean());
    settings.setFitThreeState(gd.getNextBoolean());
    settings.setFitMle(gd.getNextBoolean());
    settings.setSignificanceLevel(gd.getNextNumber());

    settings.setShowCdf(gd.getNextBoolean());
    settings.setShowSeparatePlots(gd.getNextBoolean());

    SettingsManager.writeSettings(settings.build());

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Depth of Field", settings.getDepthOfField());
      ParameterUtils.isAboveZero("Max T", settings.getMaxT());
      ParameterUtils.isPositive("Offsets", settings.getOffsets());
      ParameterUtils.isPositive("Precision", settings.getPrecision());

      ParameterUtils.isAboveZero("Bin width", settings.getBinWidth());

      ParameterUtils.isPositive("A", settings.getA());
      ParameterUtils.isPositive("B", settings.getB());

      ParameterUtils.isAboveZero("Repeats", settings.getRepeats());
      ParameterUtils.isAboveZero("Min D", settings.getMinD());
      ParameterUtils.isAboveZero("Max D", settings.getMaxD());
      ParameterUtils.isAboveZero("Significance level", settings.getSignificanceLevel());

      ParameterUtils.isEqualOrAbove("Max D", settings.getMaxD(), settings.getMinD());
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
    md.setSelected(settings.getSelectedList());
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
    settings.clearSelected();
    selected.forEach(settings::addSelected);
    SettingsManager.writeSettings(settings.build());

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
    final int maxT = settings.getMaxT();
    final FloatArrayList[] distances =
        IntStream.rangeClosed(0, maxT).mapToObj(FloatArrayList::new).toArray(FloatArrayList[]::new);
    final LocalList<Position> track = new LocalList<>();

    for (final MemoryPeakResults results : allResults) {
      results.sort(IdFramePeakResultComparator.INSTANCE);
      final int[] id = {-1};
      results.forEach(DistanceUnit.UM, (XyrResultProcedure) (x, y, r) -> {
        if (r.getId() != id[0]) {
          id[0] = r.getId();
          // Process track
          final int maxStart = Math.min(track.size() - 1, settings.getOffsets() + 1);
          for (int i = 0; i < maxStart; i++) {
            final Position origin = track.unsafeGet(i);
            for (int j = i + 1; j < track.size(); j++) {
              final Position position = track.unsafeGet(j);
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
        final Position origin = track.unsafeGet(i);
        for (int j = i + 1; j < track.size(); j++) {
          final Position position = track.unsafeGet(i);
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

  private PointValuePair fitDistances(int[][] counts, float[][] pdf, ExecutorService executor,
      boolean mle) {
    final PointValuePair result2 = fitTwoStateDistances(counts, pdf, executor, mle);
    if (result2 == null) {
      return null;
    }
    addToResultTable(result2);
    final PointValuePair result3 = fitThreeStateDistances(counts, pdf, executor, mle);
    if (result3 != null) {
      addToResultTable(result3);

      if (!mle) {
        // Select best fit using BIC
        final int numberOfPoints = pdf[0].length * pdf.length;
        final double bic2 = MathUtils.getDeltaBayesianInformationCriterion(result2.getValue(),
            numberOfPoints, result2.getPointRef().length - 1);
        final double bic3 = MathUtils.getDeltaBayesianInformationCriterion(result3.getValue(),
            numberOfPoints, result3.getPointRef().length - 1);
        final boolean reject = bic3 >= bic2;
        LoggerUtils.log(logger, Level.INFO, "Two-state -> three-state : BIC = %s -> %s (Reject=%b)",
            bic2, bic3, reject);
        if (!reject) {
          return result3;
        }
      } else {
        // Select best fit using log-likelihood ratio test
        final double ll2 = result2.getValue();
        final double ll3 = result3.getValue();
        final double llr = 2 * (ll3 - ll2);
        if (llr >= 0) {
          // The difference in the number of fitted parameters will be 2:
          // i.e. 1 extra diffusion coefficient and population fraction
          final double q = ChiSquaredDistributionTable.computeQValue(llr, 2);
          final boolean reject = q > settings.getSignificanceLevel();
          LoggerUtils.log(logger, Level.INFO,
              "Two-state -> three-state : LL = %s -> %s, LLR = %s, q-value = %s (Reject=%b)", ll2,
              ll3, MathUtils.rounded(llr, 4), MathUtils.rounded(q, 4), reject);
          if (!reject) {
            return result3;
          }
        } else {
          LoggerUtils.log(logger, Level.INFO,
              "Two-state -> three-state : LL = %s -> %s, LLR = %s : No improvement", ll2, ll3,
              MathUtils.rounded(llr, 4));
        }
      }
    }
    return result2;
  }

  private PointValuePair fitTwoStateDistances(int[][] counts, float[][] pdf,
      ExecutorService executor, boolean mle) {
    IJ.showStatus("Fitting two-state model...");

    final double dt = exposureTime;
    final double dz = settings.getDepthOfField() / 1000;
    final double precision = settings.getPrecision() / 1000;
    final UniformRandomProvider rng = UniformRandomProviders.create();

    final MaxEval maxEval = new MaxEval(2000);
    final CustomPowellOptimizer powellOptimizer = createCustomPowellOptimizer();

    double best = Double.NEGATIVE_INFINITY;
    PointValuePair result = null;

    // fit: f1, d1, d2, sigma
    // Note that f2 = 1 - f1
    final ObjectiveFunction fun;
    final GoalType goalType;
    final int numberOfPoints;
    if (mle) {
      fun = new ObjectiveFunction(new TwoStateFunction(dt, dz, settings.getA(), settings.getB(),
          precision, settings.getBinWidth(), counts, executor)::logLikelihood);
      goalType = GoalType.MAXIMIZE;
      numberOfPoints = (int) (counts[0].length * MathUtils.sum(counts[0]));
    } else {
      fun = new ObjectiveFunction(new TwoStateFunction(dt, dz, settings.getA(), settings.getB(),
          precision, settings.getBinWidth(), pdf, executor)::sumOfSquares);
      goalType = GoalType.MINIMIZE;
      best = Double.POSITIVE_INFINITY;
      numberOfPoints = pdf[0].length * pdf.length;
    }
    final double minD = Math.min(0.1, settings.getMinD());
    final SimpleBounds bounds = new SimpleBounds(addPrecision(new double[] {0, 0.0, minD}, 0),
        addPrecision(new double[] {1, 0.1, Double.POSITIVE_INFINITY}, 0.1));
    final CustomPowellOptimizer.BasisStep step =
        new CustomPowellOptimizer.BasisStep(addPrecision(new double[] {0.1, 0.01, 0.5}, 0.001));

    final Ticker ticker = ImageJUtils.createTicker(settings.getRepeats(), 1);
    for (int n = 1; n <= settings.getRepeats(); n++) {
      try {
        // Guess in range
        final double[] start = addPrecision(new double[] {
            // Bound fraction
            rng.nextDouble(0.1, 0.9), rng.nextDouble(0.1),
            // Free fraction
            rng.nextDouble(settings.getMinD(), settings.getMaxD())},
            // precision
            rng.nextDouble(0.01, 0.03));

        final PointValuePair solution =
            powellOptimizer.optimize(maxEval, fun, bounds, step, goalType, new InitialGuess(start));

        if (mle) {
          LoggerUtils.log(logger, Level.INFO,
              "Two-state fit [%d]: MLE = %s, BIC = %s (%d evaluations)", n, solution.getValue(),
              MathUtils.getBayesianInformationCriterion(solution.getValue(), start.length,
                  numberOfPoints),
              powellOptimizer.getEvaluations());
          if (solution.getValue() > best) {
            best = solution.getValue();
            result = solution;
          }
        } else {
          LoggerUtils.log(logger, Level.INFO,
              "Two-state fit [%d]: SS = %s, delta BIC = %s (%d evaluations)", n,
              solution.getValue(),
              MathUtils.getDeltaBayesianInformationCriterion(solution.getValue(), numberOfPoints,
                  start.length),
              powellOptimizer.getEvaluations());
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
      // Computed as [f1, d1, d2 [, sigma]]
      // Return as [f1, d1, f2, d2 [, sigma]]
      final double[] a = result.getPointRef();
      final double[] p = new double[a.length + 1];
      System.arraycopy(a, 0, p, 0, 2);
      p[2] = 1 - p[0];
      System.arraycopy(a, 2, p, 3, a.length - 2);
      sort2state(p);
      result = new PointValuePair(p, result.getValue());
    }

    return result;
  }

  private static void sort2state(final double[] r) {
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

  private PointValuePair fitThreeStateDistances(int[][] counts, float[][] pdf,
      ExecutorService executor, boolean mle) {
    if (!settings.getFitThreeState()) {
      return null;
    }

    IJ.showStatus("Fitting three-state model...");

    final double dt = exposureTime;
    final double dz = settings.getDepthOfField() / 1000;
    final double precision = settings.getPrecision() / 1000;
    final UniformRandomProvider rng = UniformRandomProviders.create();

    final MaxEval maxEval = new MaxEval(2000);
    final CustomPowellOptimizer powellOptimizer = createCustomPowellOptimizer();

    double best = Double.NEGATIVE_INFINITY;
    PointValuePair result = null;

    // TODO: Update to not fit f3

    // fit: f1, d1, f2, d2, f3, d3, sigma
    // Note that f3 = 1 - f1 - f2 so the number of fitted parameters is 6.
    // However the optimiser does not support compound constraints (f1+f2<=1)
    // so for convenience f3 is also a fitted parameter for the optimiser.
    final ObjectiveFunction fun;
    final GoalType goalType;
    final int numberOfPoints;
    if (mle) {
      fun = new ObjectiveFunction(new ThreeStateFunction(dt, dz, settings.getA(), settings.getB(),
          precision, settings.getBinWidth(), counts, executor)::logLikelihood);
      goalType = GoalType.MAXIMIZE;
      numberOfPoints = (int) (counts[0].length * MathUtils.sum(counts[0]));
    } else {
      fun = new ObjectiveFunction(new ThreeStateFunction(dt, dz, settings.getA(), settings.getB(),
          precision, settings.getBinWidth(), pdf, executor)::sumOfSquares);
      goalType = GoalType.MINIMIZE;
      best = Double.POSITIVE_INFINITY;
      numberOfPoints = pdf[0].length * pdf.length;
    }
    final double minD = Math.min(0.1, settings.getMinD());
    final SimpleBounds bounds =
        new SimpleBounds(addPrecision(new double[] {0, 0.0, 0, minD, 0, minD}, 0), addPrecision(
            new double[] {1, 0.1, 1, Double.POSITIVE_INFINITY, 1, Double.POSITIVE_INFINITY}, 0.1));
    final CustomPowellOptimizer.BasisStep step = new CustomPowellOptimizer.BasisStep(
        addPrecision(new double[] {0.1, 0.01, 0.1, 0.5, 0.1, 0.5}, 0.001));

    final Ticker ticker = ImageJUtils.createTicker(settings.getRepeats(), 1);
    final double range = settings.getMaxD() - settings.getMinD();
    for (int n = 1; n <= settings.getRepeats(); n++) {
      try {
        // Guess in range
        final double[] start = addPrecision(new double[] {
            // Bound fraction
            rng.nextDouble(0.1, 0.4), rng.nextDouble(0.1),
            // Slow fraction
            rng.nextDouble(0.1, 0.4),
            rng.nextDouble(settings.getMinD(), settings.getMinD() + range / 3),
            // Fast fraction
            0, rng.nextDouble(settings.getMaxD() - range / 3, settings.getMaxD())},
            // precision
            precision * rng.nextDouble(0.9, 1.1));
        // Set f3
        start[4] = 1 - start[0] - start[2];

        final PointValuePair solution =
            powellOptimizer.optimize(maxEval, fun, bounds, step, goalType, new InitialGuess(start));

        if (mle) {
          LoggerUtils.log(logger, Level.INFO,
              "Three-state fit [%d]: MLE = %s, BIC = %s (%d evaluations)", n, solution.getValue(),
              MathUtils.getBayesianInformationCriterion(solution.getValue(), start.length - 1,
                  numberOfPoints),
              powellOptimizer.getEvaluations());
          if (solution.getValue() > best) {
            best = solution.getValue();
            result = solution;
          }
        } else {
          LoggerUtils.log(logger, Level.INFO,
              "Three-state fit [%d]: SS = %s, delta BIC = %s (%d evaluations)", n,
              solution.getValue(),
              MathUtils.getDeltaBayesianInformationCriterion(solution.getValue(), numberOfPoints,
                  start.length - 1),
              powellOptimizer.getEvaluations());
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
    if (settings.getFitPrecision()) {
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
      return new TextWindow(TITLE + " Analysis", createHeader(), "", 1500, 300);
    });
  }

  private String createHeader() {
    return Arrays
        .stream(new String[] {"Dataset", "dz (nm)", "dt (ms)", "max t", "offsets", "a", "b",
            "repeats", "min D", "max D", "F", "D (um^2/s)", "sigma (nm)", "Value"})
        .collect(Collectors.joining("\t"));
  }

  private String addResult(PointValuePair result) {
    final StringBuilder sb = new StringBuilder();
    sb.append(settings.getSelected(0));
    if (settings.getSelectedCount() > 1) {
      sb.append(" + ").append(TextUtils.pleural(settings.getSelectedCount() - 1, "other"));
    }
    sb.append('\t');
    //@formatter:off
    sb.append(MathUtils.rounded(settings.getDepthOfField())).append('\t')
      .append(MathUtils.rounded(exposureTime * 1e3)).append('\t')
      .append(settings.getMaxT()).append('\t')
      .append(settings.getOffsets()).append('\t')
      .append(MathUtils.rounded(settings.getA())).append('\t')
      .append(MathUtils.rounded(settings.getB())).append('\t')
      .append(settings.getRepeats()).append('\t')
      .append(MathUtils.rounded(settings.getMinD())).append('\t')
      .append(MathUtils.rounded(settings.getMaxD()))
      ;
    //@formatter:on
    if (result != null) {
      final double[] fit = result.getPointRef();
      final int n = fit.length / 2;
      final String f = IntStream.range(0, n).mapToObj(i -> MathUtils.rounded(fit[2 * i]))
          .collect(Collectors.joining(", "));
      final String d = IntStream.range(0, n).mapToObj(i -> MathUtils.rounded(fit[2 * i + 1]))
          .collect(Collectors.joining(", "));
      final double precision = (fit.length & 1) == 1 ? fit[2 * n] * 1e3 : settings.getPrecision();
      sb.append('\t').append(f).append('\t').append(d).append('\t')
          .append(MathUtils.rounded(precision)).append('\t').append(result.getValue());
    }
    return sb.toString();
  }

  private int[][] createCounts(float[][] distances) {
    IJ.showStatus("Creating histograms...");

    final float binWidth = (float) settings.getBinWidth();
    final int[][] counts = new int[distances.length][];
    for (int i = 0; i < distances.length; i++) {
      final float[] x = distances[i];
      counts[i] = createCounts(binWidth, x);
    }

    return counts;
  }

  /**
   * Creates the histogram of the distances.
   *
   * @param binWidth the bin width
   * @param distances the distances (must be sorted)
   * @return the counts
   */
  private static int[] createCounts(float binWidth, float[] distances) {
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
    return counts.toIntArray();
  }

  private float[][] createPdf(int[][] counts) {
    IJ.showStatus("Creating PDF...");

    // Zero pad to max length
    final int n = Arrays.stream(counts).mapToInt(x -> x.length).max().getAsInt();
    final float[][] pdf = new float[counts.length][];
    for (int i = 0; i < counts.length; i++) {
      final int[] x = counts[i];
      final double w = 1.0 / MathUtils.sum(x);
      final float[] p = new float[n];
      for (int j = 0; j < x.length; j++) {
        p[j] = (float) (x[j] * w);
      }
      pdf[i] = p;
    }

    // Check the t=1 PDF is not truncated
    final double maxP = MathUtils.max(pdf[0]);
    final float endP = pdf[0][counts[0].length - 1];
    if (endP / maxP > 1e-3) {
      logger.warning(() -> String.format(
          "Possible truncation of observed distances: pdf(r=%s, dt=%s ms)=%s (fraction of max %s)",
          MathUtils.rounded(settings.getBinWidth() * counts[0].length),
          MathUtils.rounded(exposureTime * 1e3), MathUtils.rounded(endP),
          MathUtils.rounded(endP / maxP)));
    }

    return pdf;
  }

  private void plotDistances(float[][] distances, float[][] pdf, PointValuePair result,
      ExecutorService executor) {
    IJ.showStatus("Plotting results...");

    final double toMillies = exposureTime * 1e3;
    final String labels = IntStream.rangeClosed(1, distances.length)
        .mapToObj(i -> MathUtils.rounded(i * toMillies) + "ms").collect(Collectors.joining("\n"));
    final WindowOrganiser wo = new WindowOrganiser();
    final LUT lut = LutHelper.createLut(LutColour.RED_MAGENTA_BLUE, false);

    Plot[] cdfPlot = null;
    if (settings.getShowCdf()) {
      final String[] title = new String[distances.length];
      cdfPlot = new Plot[distances.length];
      if (settings.getShowSeparatePlots()) {
        for (int i = 0; i < distances.length; i++) {
          title[i] = TITLE + " distance CDF " + MathUtils.rounded((i + 1) * toMillies) + "ms";
          cdfPlot[i] = new Plot(title[i], "Distance (um)", "Probability");
        }
      } else {
        title[0] = TITLE + " distance CDF";
        cdfPlot[0] = new Plot(title[0], "Distance (um)", "Probability");
        Arrays.fill(cdfPlot, cdfPlot[0]);
      }
      double maxD = 0.01;
      for (int i = 0; i < distances.length; i++) {
        final float[] x = distances[i];
        cdfPlot[i].setColor(LutHelper.getColour(lut, distances.length - i, 1, distances.length));
        final float[] y = SimpleArrayUtils.newArray(x.length, 1f, 1f);
        final double scale = 1.0 / x.length;
        SimpleArrayUtils.apply(y, f -> (float) (f * scale));
        maxD = Math.max(maxD, x[x.length - 1]);
        cdfPlot[i].addPoints(x, y, Plot.LINE);
        cdfPlot[i].setColor(Color.black);
      }
      for (int i = 0; i < distances.length; i++) {
        cdfPlot[i].setLimits(0, maxD, 0, 1.05);
      }

      if (settings.getShowSeparatePlots()) {
        for (int i = 0; i < distances.length; i++) {
          cdfPlot[i].setFrameSize(PlotWindow.plotWidth, PlotWindow.plotHeight / 2);
          ImageJUtils.display(title[i], cdfPlot[i], 0, wo);
        }
      } else {
        final Font font = cdfPlot[0].getCurrentFont();
        cdfPlot[0].setFontSize(9);
        cdfPlot[0].addLegend(labels, "bottom-right");
        cdfPlot[0].setFont(font);
        ImageJUtils.display(title[0], cdfPlot[0], 0, wo);
      }
    }

    final String[] title = new String[distances.length];
    final Plot[] pdfPlot = new Plot[distances.length];
    if (settings.getShowSeparatePlots()) {
      for (int i = 0; i < distances.length; i++) {
        title[i] = TITLE + " distance PDF " + MathUtils.rounded((i + 1) * toMillies) + "ms";
        pdfPlot[i] = new Plot(title[i], "Distance (um)", "Probability");
      }
    } else {
      title[0] = TITLE + " distance PDF";
      pdfPlot[0] = new Plot(title[0], "Distance (um)", "Probability");
      Arrays.fill(pdfPlot, pdfPlot[0]);
    }
    float binWidth = (float) settings.getBinWidth();
    float[] r = SimpleArrayUtils.newArray(pdf[0].length, binWidth / 2, binWidth);
    float maxP = 0;
    for (int i = 0; i < distances.length; i++) {
      final float[] y = pdf[i];
      pdfPlot[i].setColor(LutHelper.getColour(lut, distances.length - i, 1, distances.length));
      maxP = MathUtils.maxDefault(maxP, y);
      pdfPlot[i].addPoints(r, y, Plot.BAR);
      pdfPlot[i].setColor(Color.black);
    }
    for (int i = 0; i < distances.length; i++) {
      pdfPlot[i].setLimits(0, r[r.length - 1], 0, maxP * 1.05);
    }

    if (settings.getShowSeparatePlots()) {
      for (int i = 0; i < distances.length; i++) {
        pdfPlot[i].setFrameSize(PlotWindow.plotWidth, PlotWindow.plotHeight / 2);
        ImageJUtils.display(title[i], pdfPlot[i], 0, wo);
      }
    } else {
      final Font font = pdfPlot[0].getCurrentFont();
      pdfPlot[0].setFontSize(9);
      pdfPlot[0].addLegend(labels, "top-right");
      pdfPlot[0].setFont(font);
      ImageJUtils.display(title[0], pdfPlot[0], 0, wo);
    }

    wo.tile();

    if (result == null) {
      return;
    }

    // Compute fitted PDF using a smaller bin width.
    // Scale using a power of 2.
    final int scale = 1 << 2;
    binWidth /= scale;
    r = SimpleArrayUtils.newArray(pdf[0].length * scale, binWidth / 2, binWidth);

    final double dt = exposureTime;
    final double dz = settings.getDepthOfField() / 1000;
    final double precision = settings.getPrecision() / 1000;
    double[][] p;
    final double[] fit = result.getPointRef();
    if (fit.length < 6) {
      final double f1 = fit[0];
      final double d1 = fit[1];
      final double f2 = fit[2];
      final double d2 = fit[3];
      final double sigma = fit.length > 4 ? fit[4] : precision;
      p = new TwoStateFunction(dt, dz, settings.getA(), settings.getB(), precision, binWidth,
          r.length, executor).pdf(distances.length, f1, d1, f2, d2, sigma);
    } else {
      p = new ThreeStateFunction(dt, dz, settings.getA(), settings.getB(), precision, binWidth,
          r.length, executor).pdf(distances.length, fit);
    }

    for (int i = 0; i < distances.length; i++) {
      final double[] pi = p[i];
      // pdfPlot[i].setLineWidth(2f);
      pdfPlot[i].setColor(LutHelper.getColour(lut, distances.length - i, 1, distances.length));

      // The PDF is normalised to sum to 1 so increase the height to compensate the
      // reduced bin width so the plots align.
      float[] yy = SimpleArrayUtils.toFloat(pi);
      SimpleArrayUtils.apply(yy, v -> v * scale);
      pdfPlot[i].addPoints(r, yy, Plot.LINE);

      // Cumulative sum
      if (cdfPlot == null) {
        continue;
      }
      // cdfPlot[i].setLineWidth(2f);
      cdfPlot[i].setColor(LutHelper.getColour(lut, distances.length - i, 1, distances.length));
      double sum = 0;
      yy = yy.clone();
      for (int j = 0; j < yy.length; j++) {
        sum += pi[j];
        yy[j] = (float) sum;
      }
      cdfPlot[i].addPoints(r, yy, Plot.DOT);
    }
  }

  /**
   * Creates the weights for each time delay. These are used to evenly balance the different number
   * of observations for each time delay. The weight for is maxN / N where N is the number of
   * observations for the time delay.
   *
   * @param counts the observed counts
   * @return the weights
   */
  private static double[] createWeights(int[][] counts) {
    final int[] sizes = Arrays.stream(counts).mapToInt(x -> (int) MathUtils.sum(x)).toArray();
    final double maxSize = MathUtils.max(sizes);
    return Arrays.stream(sizes).mapToDouble(size -> maxSize / size).toArray();
  }

  // Implementation Note:
  //
  // Note that the Spot-On PDF function is the probability density of the observed distribution.
  // It is not a true probability density function that sums to 1. This is because the mixture
  // of probabilities for each fraction discards some of the probability for the moving
  // fractions (out-of-focus correction). As the fraction of moving particles increases the sum of
  // the PDF will decrease. Thus the PDF cannot be used directly in maximum likelihood estimation.
  // In order to use the PDF as a likelihood function it is integrated over the range
  // of the function where PDF(r) significantly contributes to the sum. Then the values
  // are normalised so the sum is 1.
  //
  // This corrected PDF is either fit to the observed PDF using least squares fitting, or
  // used as a likelihood function for maximum likelihood estimation (MLE). For performance
  // MLE uses binned counts of the observations with the entire bin using the same likelihood.
  //
  // Note that least squares fitting weights all observations the same as the empirical PDF
  // for each time delay is the same size. MLE fitting sums n * log(likelihood) where n
  // is the count of observations at each histogram bin. Since shorter time delays
  // have more observations this will bias MLE. To compensate the time delays are weighted
  // using their total number of observations so that the weighted number of observations
  // for each time delay are the same (i.e. sum w * n * log(likelihood)).
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
   * Base function for the probability density.
   */
  private static abstract class BaseFunction {
    double dt;
    double dz;
    double a;
    double b;
    double precision;
    double dr;
    int n;
    float[][] pdf;
    int[][] counts;
    double[] weights;
    ExecutorService executor;

    BaseFunction(double dt, double dz, double a, double b, double precision, double dr,
        float[][] pdf, ExecutorService executor) {
      this(dt, dz, a, b, precision, dr, Arrays.stream(pdf).mapToInt(x -> x.length).max().getAsInt(),
          executor);
      this.pdf = pdf;
    }

    BaseFunction(double dt, double dz, double a, double b, double precision, double dr,
        int[][] counts, ExecutorService executor) {
      this(dt, dz, a, b, precision, dr,
          Arrays.stream(counts).mapToInt(x -> x.length).max().getAsInt(), executor);
      this.counts = counts;
      weights = createWeights(counts);
    }

    BaseFunction(double dt, double dz, double a, double b, double precision, double dr, int n,
        ExecutorService executor) {
      this.dt = dt;
      this.dz = dz;
      this.a = a;
      this.b = b;
      this.precision = precision;
      this.dr = 0.5 * dr;
      this.n = n;
      this.executor = executor;
    }
  }
  /**
   * Function to evaluate the two-state diffusion model. Creates a binned PDF for the model.
   * Evaluate the sum-of-squares of the observed PDF, or the log-likelihood of the observed counts.
   *
   * <p>This function is slow and can use an {@link ExecutorService} for parallel evaluation.
   */
  static class TwoStateFunction extends BaseFunction {

    TwoStateFunction(double dt, double dz, double a, double b, double precision, double dr,
        float[][] pdf, ExecutorService executor) {
      super(dt, dz, a, b, precision, dr, pdf, executor);
    }

    TwoStateFunction(double dt, double dz, double a, double b, double precision, double dr,
        int[][] counts, ExecutorService executor) {
      super(dt, dz, a, b, precision, dr, counts, executor);
    }

    TwoStateFunction(double dt, double dz, double a, double b, double precision, double dr, int n,
        ExecutorService executor) {
      super(dt, dz, a, b, precision, dr, n, executor);
    }

    double[][] pdf(int maxT, double f1, double d1, double f2, double d2, double sigma) {
      final List<Future<double[]>> futures = new LocalList<>(maxT);
      for (int n = 0; n < maxT; n++) {
        final int time = n;
        futures.add(executor.submit(() -> computePdf(f1, d1, f2, d2, sigma, time)));
      }

      // Finish processing data
      final double[][] p = new double[maxT][];
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

    double sumOfSquares(double[] point) {
      final double f1 = point[0];
      final double d1 = point[1];
      final double d2 = point[2];
      final double sigma = point.length > 3 ? point[3] : precision;
      final double f2 = 1 - f1;

      final List<Future<Double>> futures = new LocalList<>(pdf.length);
      for (int n = pdf.length; --n >= 0;) {
        final int time = n;
        futures.add(executor.submit(() -> computeSS(f1, d1, f2, d2, sigma, time)));
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

    private double computeSS(double f1, double d1, double f2, double d2, double sigma, int time) {
      final double[] p = computePdf(f1, d1, f2, d2, sigma, time);
      final float[] obs = pdf[time];
      double ss = 0;
      for (int i = 0; i < obs.length; i++) {
        final double dp = p[i] - obs[i];
        ss += dp * dp;
      }
      return ss;
    }

    double logLikelihood(double[] point) {
      final double f1 = point[0];
      final double d1 = point[1];
      final double d2 = point[2];
      final double sigma = point.length > 3 ? point[3] : precision;
      final double f2 = 1 - f1;

      final List<Future<Double>> futures = new LocalList<>(counts.length);
      for (int n = counts.length; --n >= 0;) {
        final int time = n;
        futures.add(executor.submit(() -> computeLL(f1, d1, f2, d2, sigma, time)));
      }

      // Finish processing data
      double ll = 0;
      for (final Future<Double> f : futures) {
        try {
          ll += f.get();
        } catch (InterruptedException | ExecutionException e) {
          throw new RuntimeException(e);
        }
      }
      return ll;
    }

    private double computeLL(double f1, double d1, double f2, double d2, double sigma, int time) {
      final double[] p = computePdf(f1, d1, f2, d2, sigma, time);
      final int[] obs = counts[time];
      double ll = 0;
      for (int i = 0; i < obs.length; i++) {
        ll += obs[i] * Math.log(p[i]);
      }
      return ll * weights[time];
    }

    private double[] computePdf(double f1, double d1, double f2, double d2, double sigma,
        int time) {
      final double deltaT = dt * (time + 1);
      // Use corrected depth-of-field
      final double p2 = DiffusionAnalysis.remaining(deltaT, dz + a * Math.sqrt(d2) + b, d2);
      final double denom1 = 1.0 / (4 * (d1 * deltaT + sigma * sigma));
      final double denom2 = 1.0 / (4 * (d2 * deltaT + sigma * sigma));
      final double[] p = new double[n];
      // We have 2n+1 distances r for n observations.
      // Integrate using Simpson's rule: f(a) + 4*f((a+b)/2) + f(b)
      double sum = 0;
      // i=0 : pdf(r=0) = 0
      double last = 0;
      double r = 0;
      int j = 0;
      for (int i = 0; i < p.length; i++) {
        r = ++j * dr;
        final double x = f1 * (r * 2 * denom1) * Math.exp(-r * r * denom1)
            + p2 * f2 * (r * 2 * denom2) * Math.exp(-r * r * denom2);
        r = ++j * dr;
        final double y = f1 * (r * 2 * denom1) * Math.exp(-r * r * denom1)
            + p2 * f2 * (r * 2 * denom2) * Math.exp(-r * r * denom2);
        p[i] = last + 4 * x + y;
        sum += p[i];
        last = y;
      }
      // Continue the integration until terms are insignificant
      while (last / sum > 1e-5) {
        r = ++j * dr;
        final double x = f1 * (r * 2 * denom1) * Math.exp(-r * r * denom1)
            + p2 * f2 * (r * 2 * denom2) * Math.exp(-r * r * denom2);
        r = ++j * dr;
        final double y = f1 * (r * 2 * denom1) * Math.exp(-r * r * denom1)
            + p2 * f2 * (r * 2 * denom2) * Math.exp(-r * r * denom2);
        sum += last + 4 * x + y;
        last = y;
      }
      final double inverseSum = 1 / sum;
      for (int i = 0; i < p.length; i++) {
        p[i] *= inverseSum;
      }
      return p;
    }
  }

  /**
   * Function to evaluate the two-state diffusion model. Creates a binned PDF for the model.
   * Evaluate the sum-of-squares of the observed PDF, or the log-likelihood of the observed counts.
   *
   * <p>This function is slow and can use an {@link ExecutorService} for parallel evaluation.
   */
  static class ThreeStateFunction extends BaseFunction {

    ThreeStateFunction(double dt, double dz, double a, double b, double precision, double dr,
        float[][] pdf, ExecutorService executor) {
      super(dt, dz, a, b, precision, dr, pdf, executor);
    }

    ThreeStateFunction(double dt, double dz, double a, double b, double precision, double dr,
        int[][] counts, ExecutorService executor) {
      super(dt, dz, a, b, precision, dr, counts, executor);
    }

    ThreeStateFunction(double dt, double dz, double a, double b, double precision, double dr, int n,
        ExecutorService executor) {
      super(dt, dz, a, b, precision, dr, n, executor);
    }

    double[][] pdf(int maxT, double[] point) {
      final double sum = point[0] + point[2] + point[4];
      final double f1 = point[0] / sum;
      final double d1 = point[1];
      final double f2 = point[2] / sum;
      final double d2 = point[3];
      final double f3 = point[4] / sum;
      final double d3 = point[5];
      final double sigma = point.length > 6 ? point[6] : precision;

      final List<Future<double[]>> futures = new LocalList<>(maxT);
      for (int n = 0; n < maxT; n++) {
        final int time = n;
        futures.add(executor.submit(() -> computePdf(f1, d1, f2, d2, f3, d3, sigma, time)));
      }

      // Finish processing data
      final double[][] p = new double[maxT][];
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

    double sumOfSquares(double[] point) {
      final double sum = point[0] + point[2] + point[4];
      final double f1 = point[0] / sum;
      final double d1 = point[1];
      final double f2 = point[2] / sum;
      final double d2 = point[3];
      final double f3 = point[4] / sum;
      final double d3 = point[5];
      final double sigma = point.length > 6 ? point[6] : precision;

      final List<Future<Double>> futures = new LocalList<>(pdf.length);
      for (int n = pdf.length; --n >= 0;) {
        final int time = n;
        futures.add(executor.submit(() -> computeSS(f1, d1, f2, d2, f3, d3, sigma, time)));
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

    private double computeSS(double f1, double d1, double f2, double d2, double f3, double d3,
        double sigma, int time) {
      final double[] p = computePdf(f1, d1, f2, d2, f3, d3, sigma, time);
      final float[] obs = pdf[time];
      double ss = 0;
      for (int i = 0; i < obs.length; i++) {
        final double dp = p[i] - obs[i];
        ss += dp * dp;
      }
      return ss;
    }

    double logLikelihood(double[] point) {
      final double sum = point[0] + point[2] + point[4];
      final double f1 = point[0] / sum;
      final double d1 = point[1];
      final double f2 = point[2] / sum;
      final double d2 = point[3];
      final double f3 = point[4] / sum;
      final double d3 = point[5];
      final double sigma = point.length > 6 ? point[6] : precision;

      final List<Future<Double>> futures = new LocalList<>(counts.length);
      for (int n = counts.length; --n >= 0;) {
        final int time = n;
        futures.add(executor.submit(() -> computeLL(f1, d1, f2, d2, f3, d3, sigma, time)));
      }

      // Finish processing data
      double ll = 0;
      for (final Future<Double> f : futures) {
        try {
          ll += f.get();
        } catch (InterruptedException | ExecutionException e) {
          throw new RuntimeException(e);
        }
      }
      return ll;
    }

    private double computeLL(double f1, double d1, double f2, double d2, double f3, double d3,
        double sigma, int time) {
      final double[] p = computePdf(f1, d1, f2, d2, f3, d3, sigma, time);
      final int[] obs = counts[time];
      double ll = 0;
      for (int i = 0; i < obs.length; i++) {
        ll += obs[i] * Math.log(p[i]);
      }
      return ll * weights[time];
    }

    private double[] computePdf(double f1, double d1, double f2, double d2, double f3, double d3,
        double sigma, int time) {
      final double deltaT = dt * (time + 1);
      // Use corrected depth-of-field
      final double p2 = DiffusionAnalysis.remaining(deltaT, dz + a * Math.sqrt(d2) + b, d2);
      final double p3 = DiffusionAnalysis.remaining(deltaT, dz + a * Math.sqrt(d3) + b, d3);
      final double denom1 = 1.0 / (4 * (d1 * deltaT + sigma * sigma));
      final double denom2 = 1.0 / (4 * (d2 * deltaT + sigma * sigma));
      final double denom3 = 1.0 / (4 * (d3 * deltaT + sigma * sigma));
      final double[] p = new double[n];
      // We have 2n+1 distances r for n observations.
      // Integrate using Simpson's rule: f(a) + 4*f((a+b)/2) + f(b)
      double sum = 0;
      // i=0 : pdf(r=0) = 0
      double last = 0;
      double r = 0;
      int j = 0;
      for (int i = 0; i < p.length; i++) {
        r = ++j * dr;
        final double x = f1 * (r * 2 * denom1) * Math.exp(-r * r * denom1)
            + p2 * f2 * (r * 2 * denom2) * Math.exp(-r * r * denom2)
            + p3 * f3 * (r * 2 * denom3) * Math.exp(-r * r * denom3);
        r = ++j * dr;
        final double y = f1 * (r * 2 * denom1) * Math.exp(-r * r * denom1)
            + p2 * f2 * (r * 2 * denom2) * Math.exp(-r * r * denom2)
            + p3 * f3 * (r * 2 * denom3) * Math.exp(-r * r * denom3);
        p[i] = last + 4 * x + y;
        sum += p[i];
        last = y;
      }
      // Continue the integration until terms are insignificant
      while (last / sum > 1e-5) {
        r = ++j * dr;
        final double x = f1 * (r * 2 * denom1) * Math.exp(-r * r * denom1)
            + p2 * f2 * (r * 2 * denom2) * Math.exp(-r * r * denom2)
            + p3 * f3 * (r * 2 * denom3) * Math.exp(-r * r * denom3);
        r = ++j * dr;
        final double y = f1 * (r * 2 * denom1) * Math.exp(-r * r * denom1)
            + p2 * f2 * (r * 2 * denom2) * Math.exp(-r * r * denom2)
            + p3 * f3 * (r * 2 * denom3) * Math.exp(-r * r * denom3);
        sum += last + 4 * x + y;
        last = y;
      }
      final double inverseSum = 1 / sum;
      for (int i = 0; i < p.length; i++) {
        p[i] *= inverseSum;
      }
      return p;
    }
  }
}
