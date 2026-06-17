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
import java.util.function.DoublePredicate;
import java.util.function.DoubleUnaryOperator;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.exception.TooManyIterationsException;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.SimpleValueChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.rng.JumpableUniformRandomProvider;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.ContinuousSampler;
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
import uk.ac.sussex.gdsc.core.utils.rng.RandomGeneratorAdapter;
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
  private static final int MODE_CDF = 1;
  private static final int MODE_PDF_MLE = 2;
  private static final String[] OPTIMISER_MODES = {"Powell", "CMA-ES", "BOBYQA"};
  private static final int OPT_MODE_CMAES = 1;
  private static final int OPT_MODE_BOBYQA = 2;

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

    final IntArrayList gapCounts = new IntArrayList();
    final float[][] distances = getDistances(allResults, gapCounts);
    // Create counts and PDF for plotting/fitting; and optionally CDF for fitting only
    final int[][] counts = createCounts(distances, (float) settings.getBinWidth());
    final float[][] pdf = createPdf(counts);
    float[][] df = pdf;
    float[][] cdf = null;
    // Note: The CDF is binned for plot display as showing all distances is slow and
    // the plot resolution cannot show all distances.
    if (settings.getShowCdf() || settings.getFitMode() == MODE_CDF) {
      cdf = createCdf(createCounts(distances, (float) settings.getCdfBinWidth()));
      // Fit the CDF as the distribution function
      if (settings.getFitMode() == MODE_CDF) {
        df = cdf;
      }
    }

    final int threadCount = Prefs.getThreads();
    final ExecutorService executor = Executors.newFixedThreadPool(threadCount);

    try {
      final PointValuePair result = fitDistances(counts, df, executor, settings.getFitMode());
      if (result != null) {
        final double[] fit = result.getPointRef();
        double d;;
        double sigma;
        if (fit.length < 6) {
          d = fit[3];
          sigma = fit.length > 4 ? fit[4] : settings.getPrecision();
        } else {
          d = fit[5];
          sigma = fit.length > 6 ? fit[6] : settings.getPrecision();
        }
        checkTraceSettings(distances[0][distances[0].length - 1], d, sigma, gapCounts);
      }
      plotDistances(cdf, pdf, result, executor);
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
      final boolean allowRestarts = settings.getAllowRestarts();
      final double halfLife = settings.getHalfLife();

      // Create the depth of field detector curve
      final DoubleUnaryOperator detectorCurve = settings.getUseDetector()
          ? DiffusionDepthOfField.loadDetectorCurve(settings.getDetectorCurveFilename())
          : null;

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
        // If no restarts then there is no requirement for a half-life.
        // The half-life parameter can be zero and the track will reach max T.
        final ContinuousSampler lifeSampler = halfLife <= 0 ? () -> maxT
            : SamplerUtils.createExponentialSampler(parentRng.jump(), halfLife);
        // Use an optional detector curve
        final DoublePredicate dof;
        if (settings.getUseDetector()) {
          dof = DiffusionDepthOfField.createDetector(detectorCurve, parentRng.jump());
        } else {
          dof = z -> Math.abs(z) <= halfDz;
        }
        futures.add(executor.submit(() -> {
          final NormalizedGaussianSampler sampler =
              SamplerUtils.createNormalizedGaussianSampler(rng);
          final MemoryPeakResults observed = new MemoryPeakResults(total);
          for (;;) {
            final int pos = position.incrementAndGet();
            if (pos > total) {
              break;
            }
            // Create a lifetime within the depth-of-field.
            final int lifetime = (int) Math.min(lifeSampler.sample(), 20 * maxT);
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
            int id = -1;
            double x = 0;
            double y = 0;
            double z = rng.nextDouble(-halfDz, halfDz);
            // Last frame the molecule was inside the DoF
            int lastT = 0;
            // Ensure all tracks have unique frames to allow tracing the dataset
            // Each track should be separated by maxT frames so tracing even with
            // small frame gaps will not join stationary molecules at the origin.
            final int frame = pos * (21 * maxT);
            // Diffuse until detected
            int t = 0;
            while (!dof.test(z) && t < lifetime) {
              t++;
              // Note: Updates to (x,y) do not matter here
              z += sampler.sample() * d;
            }
            if (t == lifetime) {
              // Never detected
              break;
            }
            // Record the origin
            id = nextId.incrementAndGet();
            lastT = t;
            observed.add(new XyzPeakResult(frame + t, x + sampler.sample() * precision,
                y + sampler.sample() * precision, z, id));
            // Diffuse until lifetime expires
            while (t < lifetime) {
              t++;
              x += sampler.sample() * d;
              y += sampler.sample() * d;
              z += sampler.sample() * d;
              if (dof.test(z)) {
                // Record molecule
                // Check the gap from the last time it was in the depth-of-field
                if (t - lastT > g) {
                  if (allowRestarts) {
                    // Optionally allow restarting a new track
                    id = nextId.incrementAndGet();
                  } else {
                    break;
                  }
                }
                lastT = t;
                // Add localisation precision
                observed.add(new XyzPeakResult(frame + t, x + sampler.sample() * precision,
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

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

    gd.addNumericField("Depth_of_field", settings.getDepthOfField(), 1, 6, "nm");
    gd.addNumericField("Exposure_time", settings.getExposureTime(), 1, 6, "ms");
    gd.addSlider("Gap", 1, 5, settings.getGap());

    gd.addNumericField("Number_of_molecules", settings.getNumberOfMolecules());
    gd.addSlider("Max_t", 5, 15, settings.getSimT());
    gd.addCheckbox("Use_dectector", settings.getUseDetector());
    gd.addFilenameField("Detector_curve", settings.getDetectorCurveFilename());
    gd.addCheckbox("Allow_restarts", settings.getAllowRestarts());
    gd.addNumericField("Half_life", settings.getHalfLife());

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
    settings.setUseDetector(gd.getNextBoolean());
    settings.setDetectorCurveFilename(gd.getNextString());
    settings.setAllowRestarts(gd.getNextBoolean());
    settings.setHalfLife(gd.getNextNumber());

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

      if (settings.getAllowRestarts()) {
        ParameterUtils.isAboveZero("Half-life with restarts", settings.getHalfLife());
      }

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

    gd.addChoice("Fit_mode", new String[] {"PDF", "CDF", "PDF_MLE"}, settings.getFitMode());
    gd.addNumericField("Bin_width", settings.getBinWidth(), -3);
    gd.addNumericField("CDF_bin_width", settings.getCdfBinWidth(), -3);

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

    gd.addChoice("Optimiser_mode", OPTIMISER_MODES, settings.getOptimiserMode());
    gd.addSlider("Repeats", 1, 5, settings.getRepeats());
    gd.addNumericField("Min_D", settings.getMinD(), 3, 6, "um^2/s");
    gd.addNumericField("Max_D", settings.getMaxD(), 3, 6, "um^2/s");

    gd.addCheckbox("Fit_precision", settings.getFitPrecision());
    gd.addCheckbox("Fit_three_state", settings.getFitThreeState());
    gd.addNumericField("Significance_level", settings.getSignificanceLevel(), -3);

    gd.addCheckbox("Show_CDF", settings.getShowCdf());
    gd.addCheckbox("Seperate_plots", settings.getShowSeparatePlots());
    gd.addNumericField("Plot_max_r", settings.getPlotMaxR());

    gd.addHelp(HelpUrls.getUrl("track-diffusion-analysis"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    settings.setDepthOfField(gd.getNextNumber());
    settings.setMaxT((int) gd.getNextNumber());
    settings.setOffsets((int) gd.getNextNumber());
    settings.setPrecision(gd.getNextNumber());

    settings.setFitMode(gd.getNextChoiceIndex());
    settings.setBinWidth(gd.getNextNumber());
    settings.setCdfBinWidth(gd.getNextNumber());

    settings.setA(gd.getNextNumber());
    settings.setB(gd.getNextNumber());

    settings.setOptimiserMode(gd.getNextChoiceIndex());
    settings.setRepeats((int) gd.getNextNumber());
    settings.setMinD(gd.getNextNumber());
    settings.setMaxD(gd.getNextNumber());

    settings.setFitPrecision(gd.getNextBoolean());
    settings.setFitThreeState(gd.getNextBoolean());
    settings.setSignificanceLevel(gd.getNextNumber());

    settings.setShowCdf(gd.getNextBoolean());
    settings.setShowSeparatePlots(gd.getNextBoolean());
    settings.setPlotMaxR(gd.getNextNumber());

    SettingsManager.writeSettings(settings.build());

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Depth of Field", settings.getDepthOfField());
      ParameterUtils.isAboveZero("Max T", settings.getMaxT());
      ParameterUtils.isPositive("Offsets", settings.getOffsets());
      ParameterUtils.isPositive("Precision", settings.getPrecision());

      ParameterUtils.isAboveZero("Bin width", settings.getBinWidth());
      ParameterUtils.isAboveZero("CDF bin width", settings.getCdfBinWidth());

      ParameterUtils.isPositive("A", settings.getA());
      ParameterUtils.isPositive("B", settings.getB());

      ParameterUtils.isAboveZero("Repeats", settings.getRepeats());
      ParameterUtils.isAboveZero("Min D", settings.getMinD());
      ParameterUtils.isAboveZero("Max D", settings.getMaxD());
      ParameterUtils.isAboveZero("Significance level", settings.getSignificanceLevel());

      ParameterUtils.isEqualOrAbove("Max D", settings.getMaxD(), settings.getMinD());
      if (settings.getShowCdf() || settings.getFitMode() == MODE_CDF) {
        ParameterUtils.isEqualOrBelow("CDF bin width", settings.getCdfBinWidth(),
            settings.getBinWidth());
      }
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
   * Get the jump distances for each time delay. Also count the frame gaps in the used part of the
   * traces. In conjunction with the maximum observed distance, this allows analysis of the trace
   * settings.
   *
   * @param allResults the results
   * @param gapCounts the gap counts
   * @return distances
   */
  private float[][] getDistances(List<MemoryPeakResults> allResults, IntArrayList gapCounts) {
    final int maxT = settings.getMaxT();
    final FloatArrayList[] distances =
        IntStream.rangeClosed(0, maxT).mapToObj(FloatArrayList::new).toArray(FloatArrayList[]::new);
    final LocalList<Position> track = new LocalList<>();
    // Count frame gaps in the traces
    final int[] gap = new int[11];

    for (final MemoryPeakResults results : allResults) {
      results.sort(IdFramePeakResultComparator.INSTANCE);
      final FrameCounter id = new FrameCounter(-1);
      results.forEach(DistanceUnit.UM, (XyrResultProcedure) (x, y, r) -> {
        if (id.advance(r.getId())) {
          if (!track.isEmpty()) {
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
            // Count the gaps in the trace. This is used to determine if trace settings
            // are suitable for fitted diffusion coefficients
            final int limit = Math.min(track.size(), maxStart + maxT + 1);
            int t = track.unsafeGet(0).t;
            for (int i = 1; i < limit; i++) {
              final int tt = track.unsafeGet(i).t;
              if (tt - t < gap.length) {
                gap[tt - t]++;
              }
              t = tt;
            }
            track.clear();
          }
        }
        track.add(new Position(r.getFrame(), x, y));
      });
      // Process final track
      if (!track.isEmpty()) {
        final int maxStart = Math.min(track.size() - 1, settings.getOffsets() + 1);
        for (int i = 0; i < maxStart; i++) {
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
        final int limit = Math.min(track.size(), maxStart + maxT + 1);
        int t = track.unsafeGet(0).t;
        for (int i = 1; i < limit; i++) {
          final int tt = track.unsafeGet(i).t;
          if (tt - t < gap.length) {
            gap[tt - t]++;
          }
          t = tt;
        }
        track.clear();
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


    // Record the gap lengths
    final int end = SimpleArrayUtils.findLastIndex(gap, v -> v != 0);
    for (int i = 1; i <= end; i++) {
      gapCounts.add(gap[i]);
    }
    LoggerUtils.log(logger, Level.INFO, "Gap counts: %s", Arrays.toString(gapCounts.toIntArray()));

    // Avoid fitting when there is no data
    final int time = SimpleArrayUtils.findIndex(sizes, i -> i == 0);
    if (time >= 0) {
      throw new IllegalStateException("No sizes for time delay: " + (time + 1));
    }

    return sortedDistances;
  }

  private PointValuePair fitDistances(int[][] counts, float[][] df, ExecutorService executor,
      int mode) {
    // TODO: Try a CMA-ES optimiser
    final PointValuePair result2 = fitTwoStateDistances(counts, df, executor, mode);
    if (result2 == null) {
      return null;
    }
    addToResultTable(result2);
    final PointValuePair result3 = fitThreeStateDistances(counts, df, executor, mode);
    if (result3 != null) {
      addToResultTable(result3);

      if (mode != MODE_PDF_MLE) {
        // Select best fit using BIC
        final int numberOfPoints = df[0].length * df.length;
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

  private PointValuePair fitTwoStateDistances(int[][] counts, float[][] df,
      ExecutorService executor, int mode) {
    IJ.showStatus("Fitting two-state model...");

    final double dt = exposureTime;
    final double dz = settings.getDepthOfField() / 1000;
    final double precision = settings.getPrecision() / 1000;
    final UniformRandomProvider rng = UniformRandomProviders.create();

    double best = Double.NEGATIVE_INFINITY;
    PointValuePair result = null;

    final LocalList<OptimizationData> args = new LocalList<>();
    // CMA-ES can take 30 N to 300 N^2 evaluations for N parameters
    args.add(new MaxEval(10000));

    // fit: f1, d1, d2, sigma
    // Note that f2 = 1 - f1
    final int numberOfPoints;
    if (mode == MODE_PDF_MLE) {
      args.add(new ObjectiveFunction(new TwoStateFunction(dt, dz, settings.getA(), settings.getB(),
          precision, settings.getBinWidth(), counts, executor)::logLikelihood));
      args.add(GoalType.MAXIMIZE);
      numberOfPoints = (int) (counts[0].length * MathUtils.sum(counts[0]));
    } else {
      final boolean isCdf = mode == MODE_CDF;
      final double binWidth = mode == MODE_CDF ? settings.getCdfBinWidth() : settings.getBinWidth();
      args.add(new ObjectiveFunction(new TwoStateFunction(dt, dz, settings.getA(), settings.getB(),
          precision, binWidth, df, isCdf, executor)::sumOfSquares));
      args.add(GoalType.MINIMIZE);
      best = Double.POSITIVE_INFINITY;
      numberOfPoints = df[0].length * df.length;
    }
    final double minD = Math.min(0.1, settings.getMinD());
    args.add(new SimpleBounds(addPrecision(new double[] {0, 0.0, minD}, 0),
        addPrecision(new double[] {1, 0.1, Double.POSITIVE_INFINITY}, 0.1)));

    MultivariateOptimizer optimizer;
    if (settings.getOptimiserMode() == OPT_MODE_CMAES) {
      optimizer = createCmaesOptimizer(rng);
      // The sigma determines the search range for the variables.
      // It should be 1/3 of the initial search region.
      args.add(new CMAESOptimizer.Sigma(addPrecision(new double[] {0.1, 0.01, 0.3}, 0.001)));
      args.add(new CMAESOptimizer.PopulationSize((int) (4 + Math.floor(3 * Math.log(4)))));
    } else if (settings.getOptimiserMode() == OPT_MODE_BOBYQA) {
      optimizer = createBobyqaOptimizer(4);
    } else {
      optimizer = createCustomPowellOptimizer();
      args.add(
          new CustomPowellOptimizer.BasisStep(addPrecision(new double[] {0.1, 0.01, 0.3}, 0.001)));
    }

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

        args.push(new InitialGuess(start));
        final PointValuePair solution =
            optimizer.optimize(args.toArray(new OptimizationData[args.size()]));
        args.pop();

        if (mode == MODE_PDF_MLE) {
          LoggerUtils.log(logger, Level.INFO,
              "Two-state fit [%d]: MLE = %s, BIC = %s (%d evaluations)", n, solution.getValue(),
              MathUtils.getBayesianInformationCriterion(solution.getValue(), start.length,
                  numberOfPoints),
              optimizer.getEvaluations());
          if (solution.getValue() > best) {
            best = solution.getValue();
            result = solution;
          }
        } else {
          LoggerUtils.log(logger, Level.INFO,
              "Two-state fit [%d]: SS = %s, delta BIC = %s (%d evaluations)", n,
              solution.getValue(), MathUtils.getDeltaBayesianInformationCriterion(
                  solution.getValue(), numberOfPoints, start.length),
              optimizer.getEvaluations());
          if (solution.getValue() < best) {
            best = solution.getValue();
            result = solution;
          }
        }
      } catch (final TooManyEvaluationsException ex) {
        LoggerUtils.log(logger, Level.INFO, "Optimiser [%d] failed : Too many evaluation (%d)", n,
            optimizer.getEvaluations());
      } catch (final TooManyIterationsException ex) {
        LoggerUtils.log(logger, Level.INFO, "Optimiser [%d] failed : Too many iterations (%d)", n,
            optimizer.getIterations());
      } catch (final ConvergenceException ex) {
        LoggerUtils.log(logger, Level.INFO, "Optimiser [%d] failed to fit : %s", n,
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

  private PointValuePair fitThreeStateDistances(int[][] counts, float[][] df,
      ExecutorService executor, int mode) {
    if (!settings.getFitThreeState()) {
      return null;
    }

    IJ.showStatus("Fitting three-state model...");

    final double dt = exposureTime;
    final double dz = settings.getDepthOfField() / 1000;
    final double precision = settings.getPrecision() / 1000;
    final UniformRandomProvider rng = UniformRandomProviders.create();

    double best = Double.NEGATIVE_INFINITY;
    PointValuePair result = null;

    final LocalList<OptimizationData> args = new LocalList<>();
    // CMA-ES can take 30 N to 30 N^2 evaluations for N parameters
    args.add(new MaxEval(20000));

    // TODO: Update to not fit f3

    // fit: f1, d1, f2, d2, f3, d3, sigma
    // Note that f3 = 1 - f1 - f2 so the number of fitted parameters is 6.
    // However the optimiser does not support compound constraints (f1+f2<=1)
    // so for convenience f3 is also a fitted parameter for the optimiser.
    final int numberOfPoints;
    if (mode == MODE_PDF_MLE) {
      args.add(new ObjectiveFunction(new ThreeStateFunction(dt, dz, settings.getA(),
          settings.getB(), precision, settings.getBinWidth(), counts, executor)::logLikelihood));
      args.add(GoalType.MAXIMIZE);
      numberOfPoints = (int) (counts[0].length * MathUtils.sum(counts[0]));
    } else {
      final boolean isCdf = mode == MODE_CDF;
      final double binWidth = mode == MODE_CDF ? settings.getCdfBinWidth() : settings.getBinWidth();
      args.add(new ObjectiveFunction(new ThreeStateFunction(dt, dz, settings.getA(),
          settings.getB(), precision, binWidth, df, isCdf, executor)::sumOfSquares));
      args.add(GoalType.MINIMIZE);
      best = Double.POSITIVE_INFINITY;
      numberOfPoints = df[0].length * df.length;
    }
    final double minD = Math.min(0.1, settings.getMinD());
    args.add(
        new SimpleBounds(addPrecision(new double[] {0, 0.0, 0, minD, 0, minD}, 0), addPrecision(
            new double[] {1, 0.1, 1, Double.POSITIVE_INFINITY, 1, Double.POSITIVE_INFINITY}, 0.1)));

    MultivariateOptimizer optimizer;
    if (settings.getOptimiserMode() == OPT_MODE_CMAES) {
      optimizer = createCmaesOptimizer(rng);
      // The sigma determines the search range for the variables.
      // It should be 1/3 of the initial search region.
      args.add(new CMAESOptimizer.Sigma(
          addPrecision(new double[] {0.1, 0.01, 0.1, 0.3, 0.1, 0.3}, 0.001)));
      args.add(new CMAESOptimizer.PopulationSize((int) (4 + Math.floor(3 * Math.log(7)))));
    } else if (settings.getOptimiserMode() == OPT_MODE_BOBYQA) {
      optimizer = createBobyqaOptimizer(7);
    } else {
      optimizer = createCustomPowellOptimizer();
      args.add(new CustomPowellOptimizer.BasisStep(
          addPrecision(new double[] {0.1, 0.01, 0.1, 0.3, 0.1, 0.3}, 0.001)));
    }

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

        args.push(new InitialGuess(start));
        final PointValuePair solution =
            optimizer.optimize(args.toArray(new OptimizationData[args.size()]));
        args.pop();

        if (mode == MODE_PDF_MLE) {
          LoggerUtils.log(logger, Level.INFO,
              "Three-state fit [%d]: MLE = %s, BIC = %s (%d evaluations)", n, solution.getValue(),
              MathUtils.getBayesianInformationCriterion(solution.getValue(), start.length - 1,
                  numberOfPoints),
              optimizer.getEvaluations());
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
              optimizer.getEvaluations());
          if (solution.getValue() < best) {
            best = solution.getValue();
            result = solution;
          }
        }
      } catch (final TooManyEvaluationsException ex) {
        LoggerUtils.log(logger, Level.INFO, "Optimiser[%d] failed : Too many evaluation (%d)", n,
            optimizer.getEvaluations());
      } catch (final TooManyIterationsException ex) {
        LoggerUtils.log(logger, Level.INFO, "Optimiser[%d] failed : Too many iterations (%d)", n,
            optimizer.getIterations());
      } catch (final ConvergenceException ex) {
        LoggerUtils.log(logger, Level.INFO, "Optimiser[%d] failed to fit : %s", n, ex.getMessage());
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

  private static CMAESOptimizer createCmaesOptimizer(UniformRandomProvider rng) {
    final double rel = 1e-8;
    final double abs = 1e-100;
    final int maxIterations = 2000;
    final double stopFitness = 0;
    final boolean isActiveCma = true;
    final int diagonalOnly = 20;
    final int checkFeasableCount = 1;
    final RandomGenerator random = new RandomGeneratorAdapter(rng);
    final boolean generateStatistics = false;
    final ConvergenceChecker<PointValuePair> checker = new SimpleValueChecker(rel, abs);
    return new CMAESOptimizer(maxIterations, stopFitness, isActiveCma, diagonalOnly,
        checkFeasableCount, random, generateStatistics, checker);
  }

  private static BOBYQAOptimizer createBobyqaOptimizer(int n) {
    final int numberOfInterpolationpoints = n + 2;
    // Optional parameters: initialTrustRegionRadius; stoppingTrustRegionRadius.
    // It is unclear what the initial radius should be. The default is 10.
    // When the optimiser starts it sets the value to the minimum bound difference
    // [upper - lower] / 3. So when using a tight bounds this settings does not
    // matter unless we require it to be smaller.
    // The default stopping trust radius is 1e-8.
    return new BOBYQAOptimizer(numberOfInterpolationpoints);
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
            "optimiser", "repeats", "min D", "max D", "F", "D (um^2/s)", "sigma (nm)", "Value"})
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
      .append(getOptimiserMode(settings.getOptimiserMode())).append('\t')
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

  private static String getOptimiserMode(int mode) {
    if (mode > 0 && mode < OPTIMISER_MODES.length) {
      return OPTIMISER_MODES[mode];
    }
    return OPTIMISER_MODES[0];
  }

  private static int[][] createCounts(float[][] distances, float binWidth) {
    IJ.showStatus("Creating histograms...");

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

  private float[][] createCdf(int[][] counts) {
    IJ.showStatus("Creating CDF...");

    // Zero pad to max length
    final int n = Arrays.stream(counts).mapToInt(x -> x.length).max().getAsInt();
    final float[][] cdf = new float[counts.length][];
    for (int i = 0; i < counts.length; i++) {
      final int[] x = counts[i];
      final double w = 1.0 / MathUtils.sum(x);
      final float[] p = new float[n];
      long c = 0;
      for (int j = 0; j < x.length; j++) {
        c += x[j];
        p[j] = (float) (c * w);
      }
      for (int j = x.length; j < n; j++) {
        p[j] = 1;
      }
      cdf[i] = p;
    }

    return cdf;
  }

  /**
   * Check the tracing settings. This checks the distance used for tracing covers enough of the
   * cumulative mean-squared distance distribution.
   *
   * @param r the maximum observed jump distance for 1 frame (um)
   * @param d the diffusion coefficient (um^2/s)
   * @param sigma the localisation precision
   * @param gapCounts the gap counts
   */
  private void checkTraceSettings(double r, double d, double sigma, IntArrayList gapCounts) {
    final long sum = MathUtils.sum(gapCounts.toIntArray());
    for (int i = 0; i < gapCounts.size(); i++) {
      final double t = (i + 1) * exposureTime;
      // cdf(r^2) = 1 - exp(-r^2 / 4dt+s^2)
      final double p = -Math.expm1(-r * r / (4 * d * t + sigma * sigma));
      Level level = Level.INFO;
      String msg = "Observed";
      if (p < 0.95) {
        level = Level.WARNING;
        msg = "Possible truncation of observed";
      }
      LoggerUtils.log(logger, level,
          "%s distances (gap=%d; %s) for max D: cdf(r=%s, D=%s, s=%s, t=%s)=%s", msg, i + 1,
          MathUtils.rounded((double) gapCounts.getInt(i) / sum), MathUtils.rounded(r),
          MathUtils.rounded(d), MathUtils.rounded(sigma), MathUtils.rounded(t),
          MathUtils.rounded(p));
    }
  }

  private void plotDistances(float[][] cdf, float[][] pdf, PointValuePair result,
      ExecutorService executor) {
    IJ.showStatus("Plotting results...");

    final double toMillies = exposureTime * 1e3;
    final String labels = IntStream.rangeClosed(1, pdf.length)
        .mapToObj(i -> MathUtils.rounded(i * toMillies) + "ms").collect(Collectors.joining("\n"));
    final WindowOrganiser wo = new WindowOrganiser();
    final LUT lut = LutHelper.createLut(LutColour.RED_MAGENTA_BLUE, false);

    Plot[] cdfPlot = null;
    if (settings.getShowCdf()) {
      final String[] title = new String[cdf.length];
      cdfPlot = new Plot[cdf.length];
      if (settings.getShowSeparatePlots()) {
        for (int i = 0; i < cdf.length; i++) {
          title[i] = TITLE + " distance CDF " + MathUtils.rounded((i + 1) * toMillies) + "ms";
          cdfPlot[i] = new Plot(title[i], "Distance (um)", "Probability");
        }
      } else {
        title[0] = TITLE + " distance CDF";
        cdfPlot[0] = new Plot(title[0], "Distance (um)", "Probability");
        Arrays.fill(cdfPlot, cdfPlot[0]);
      }
      final float binWidth = (float) settings.getCdfBinWidth();
      final float[] r = SimpleArrayUtils.newArray(cdf[0].length + 1, 0, binWidth);
      for (int i = 0; i < cdf.length; i++) {
        // Add a zero at r=0
        final float[] y = new float[r.length];
        System.arraycopy(cdf[i], 0, y, 1, y.length - 1);
        cdfPlot[i].setColor(LutHelper.getColour(lut, cdf.length - i, 1, cdf.length));
        cdfPlot[i].addPoints(r, y, Plot.LINE);
        cdfPlot[i].setColor(Color.black);
      }
      double maxR = r[r.length - 1];
      if (settings.getPlotMaxR() > 0) {
        maxR = MathUtils.clip(0.01, maxR, settings.getPlotMaxR());
      }
      for (int i = 0; i < cdf.length; i++) {
        cdfPlot[i].setLimits(0, maxR, 0, 1.05);
      }

      if (settings.getShowSeparatePlots()) {
        for (int i = 0; i < cdf.length; i++) {
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

    final String[] title = new String[pdf.length];
    final Plot[] pdfPlot = new Plot[pdf.length];
    if (settings.getShowSeparatePlots()) {
      for (int i = 0; i < pdf.length; i++) {
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
    for (int i = 0; i < pdf.length; i++) {
      final float[] y = pdf[i];
      pdfPlot[i].setColor(LutHelper.getColour(lut, pdf.length - i, 1, pdf.length));
      maxP = MathUtils.maxDefault(maxP, y);
      pdfPlot[i].addPoints(r, y, Plot.BAR);
      pdfPlot[i].setColor(Color.black);
    }
    double maxR = r[r.length - 1];
    if (settings.getPlotMaxR() > 0) {
      maxR = MathUtils.clip(0.01, maxR, settings.getPlotMaxR());
    }
    for (int i = 0; i < pdf.length; i++) {
      pdfPlot[i].setLimits(0, maxR, 0, maxP * 1.05);
    }

    if (settings.getShowSeparatePlots()) {
      for (int i = 0; i < pdf.length; i++) {
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

    // Compute fitted PDF using a smaller bin width for a more detailed plot.
    // If fitting the CDF, then use the CDF bin width if it increases resolution.
    // Scale using a power of 2.
    binWidth /= 8;
    if (cdfPlot != null && settings.getCdfBinWidth() < binWidth) {
      binWidth = (float) settings.getCdfBinWidth();
    }
    final float scale = (float) (settings.getBinWidth() / binWidth);
    r = SimpleArrayUtils.newArray((int) Math.ceil(scale * pdf[0].length), binWidth / 2, binWidth);

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
          r.length, executor).pdf(pdf.length, f1, d1, f2, d2, sigma);
    } else {
      p = new ThreeStateFunction(dt, dz, settings.getA(), settings.getB(), precision, binWidth,
          r.length, executor).pdf(pdf.length, fit);
    }

    for (int i = 0; i < pdf.length; i++) {
      final double[] pi = p[i];
      pdfPlot[i].setColor(LutHelper.getColour(lut, pdf.length - i, 1, pdf.length));

      // The PDF is normalised to sum to 1 so increase the height to compensate the
      // reduced bin width so the plots align.
      float[] yy = SimpleArrayUtils.toFloat(pi);
      SimpleArrayUtils.apply(yy, v -> v * scale);
      pdfPlot[i].addPoints(r, yy, Plot.LINE);

      // Cumulative sum
      if (cdfPlot == null) {
        continue;
      }
      cdfPlot[i].setColor(LutHelper.getColour(lut, cdf.length - i, 1, cdf.length));
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
  // This is avoided by using bounds on the fitted parameters.

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
    float[][] df;
    boolean isCdf;
    int[][] counts;
    double[] weights;
    ExecutorService executor;
    /**
     * The fraction of the current sum to terminate integration. This effectively controls how far
     * into the distance r to proceed with the integration.
     *
     * <p>Initial trials used a value of 1e-5 to ensure most of the PDF was computed. This has been
     * disabled using a term error of 1.
     *
     * <p>Note that the observed PDF/CDF may be a truncation of the full density as larger distances
     * are not observed. If the full PDF is computed then this will result in a penalty applied to
     * all distances above the final observed distance. This affects the sum-of-squares (SS) for the
     * PDF by adding a small error term equal to the computed PDF for all r until integration is
     * terminated. The effect on the SS for the CDF is larger. The observed CDF will be 1 at the
     * last observed distance. The computed CDF at this point is initially far from 1, and increases
     * to 1 over the additional distance terms computed until integration terminates. The penalty is
     * far greater for the CDF due to summing the square of these larger differences to the observed
     * density (compared to the PDF). Fitting will avoid this penalty by reducing the truncation of
     * the computed PDF by using smaller diffusion coefficients.
     *
     * <p>The best fit to the observed distribution is made by computing the PDF using the maximum
     * distance in the observed PDF and normalising to that sum. This effectively fits the curve of
     * the observed PDF.
     *
     * <p>The maximum distance is used across all time delays so each curve has the same bin width.
     * An alternative is to use a fixed number of bins and use a variable bin width for each. A
     * compromise is to only integrate the computed PDF over the max distance observed for each time
     * delay. This will alter the number of observed points.
     */
    double termError = 1;

    BaseFunction(double dt, double dz, double a, double b, double precision, double dr,
        float[][] df, boolean isCdf, ExecutorService executor) {
      this(dt, dz, a, b, precision, dr, Arrays.stream(df).mapToInt(x -> x.length).max().getAsInt(),
          executor);
      this.df = df;
      this.isCdf = isCdf;
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
   * Evaluate the sum-of-squares of the observed PDF, the observed CDF, or the log-likelihood of the
   * observed counts.
   *
   * <p>This function is slow and can use an {@link ExecutorService} for parallel evaluation.
   */
  static class TwoStateFunction extends BaseFunction {

    TwoStateFunction(double dt, double dz, double a, double b, double precision, double dr,
        float[][] df, boolean isCdf, ExecutorService executor) {
      super(dt, dz, a, b, precision, dr, df, isCdf, executor);
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

      final List<Future<Double>> futures = new LocalList<>(df.length);
      for (int n = df.length; --n >= 0;) {
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
      final float[] obs = df[time];
      double ss = 0;
      if (isCdf) {
        double s = 0;
        for (int i = 0; i < obs.length; i++) {
          s += p[i];
          final double dp = s - obs[i];
          ss += dp * dp;
        }
      } else {
        for (int i = 0; i < obs.length; i++) {
          final double dp = p[i] - obs[i];
          ss += dp * dp;
        }
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
      final double error = this.termError;
      while (last / sum > error) {
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
   * Evaluate the sum-of-squares of the observed PDF, the observed CDF, or the log-likelihood of the
   * observed counts.
   *
   * <p>This function is slow and can use an {@link ExecutorService} for parallel evaluation.
   */
  static class ThreeStateFunction extends BaseFunction {

    ThreeStateFunction(double dt, double dz, double a, double b, double precision, double dr,
        float[][] df, boolean isCdf, ExecutorService executor) {
      super(dt, dz, a, b, precision, dr, df, isCdf, executor);
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

      final List<Future<Double>> futures = new LocalList<>(df.length);
      for (int n = df.length; --n >= 0;) {
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
      final float[] obs = df[time];
      double ss = 0;
      if (isCdf) {
        double s = 0;
        for (int i = 0; i < obs.length; i++) {
          s += p[i];
          final double dp = s - obs[i];
          ss += dp * dp;
        }
      } else {
        for (int i = 0; i < obs.length; i++) {
          final double dp = p[i] - obs[i];
          ss += dp * dp;
        }
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
      final double error = this.termError;
      while (last / sum > error) {
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
