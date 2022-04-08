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
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.text.TextWindow;
import java.awt.Checkbox;
import java.awt.Color;
import java.awt.TextField;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyIterationsException;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer.Optimum;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.RealVector;
import uk.ac.sussex.gdsc.core.data.utils.Converter;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
import uk.ac.sussex.gdsc.core.ij.ImageJPluginLoggerHelper;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.SimpleImageJTrackProgress;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog.OptionListener;
import uk.ac.sussex.gdsc.core.ij.gui.MultiDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.utils.FileUtils;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationHelper;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.ConfigurationException;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.fitting.JumpDistanceAnalysis;
import uk.ac.sussex.gdsc.smlm.fitting.JumpDistanceAnalysis.CurveLogger;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.GUIProtos.ClusteringSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.ArrayPeakResultStore;
import uk.ac.sussex.gdsc.smlm.results.DynamicMultipleTargetTracing;
import uk.ac.sussex.gdsc.smlm.results.DynamicMultipleTargetTracing.DmttConfiguration;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResultStoreList;
import uk.ac.sussex.gdsc.smlm.results.Trace;
import uk.ac.sussex.gdsc.smlm.results.TraceManager;
import uk.ac.sussex.gdsc.smlm.utils.StdMath;

/**
 * Run a tracing algorithm on the peak results to trace molecules across the frames.
 */
public class TraceDiffusion implements PlugIn, CurveLogger {
  private static final String TITLE = "Trace Diffusion";

  private static AtomicReference<TextWindow> summaryTableRef = new AtomicReference<>();

  // Used for the macro extensions
  private static AtomicReference<double[][]> jumpDistanceParametersRef = new AtomicReference<>();

  // Used for the multiMode option
  private static AtomicReference<List<String>> selectedRef = new AtomicReference<>();

  /** The names for the different tracing modes. */
  private static final String[] TRACE_MODE =
      {"Nearest neighbour", "Dynamic Multiple Target Tracing"};

  private boolean directoryChosen;
  private ClusteringSettings.Builder clusteringSettings;
  private MemoryPeakResults results;
  private boolean extraOptions;
  private boolean multiMode;
  private int myMinN = 1;

  // The number of additional datasets
  private int additionalDatasets;

  // Store exposure time in seconds
  private double exposureTime;
  private double precision;
  private double beta;
  private double fitValue = Double.NaN;

  // Used to tile new plot windows
  private final WindowOrganiser windowOrganiser = new WindowOrganiser();

  private String jdTitle = TITLE + " Jump Distance";
  private Plot jdPlot;

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    static final String[] NAMES = {"Total Signal", "Signal/Frame", "t-On (s)"};
    static final boolean[] ROUNDED = {false, false, true};
    static final int TOTAL_SIGNAL = 0;
    static final int SIGNAL_PER_FRAME = 1;
    static final int T_ON = 2;

    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    String inputOption;
    String header;

    boolean[] displayHistograms;
    boolean[] alwaysRemoveOutliers;

    boolean displayMsdHistogram;
    boolean displayDHistogram;
    boolean displayTraceLength;
    boolean displayTraceSize;

    boolean saveTraceDistances;
    boolean saveRawData;
    String rawDataDirectory;
    String distancesFilename;
    double significanceLevel;
    double minFraction;
    double minDifference;
    int minN;
    int maxN;
    boolean debugFitting;
    String tracesFilename;
    String title;

    Settings() {
      // Set defaults
      inputOption = "";
      displayHistograms = new boolean[NAMES.length];
      alwaysRemoveOutliers = new boolean[NAMES.length];
      alwaysRemoveOutliers[TOTAL_SIGNAL] = false;
      displayMsdHistogram = true;
      displayDHistogram = true;
      rawDataDirectory = "";
      distancesFilename = "";
      significanceLevel = 0.05;
      minFraction = 0.1;
      minDifference = 2;
      minN = 1;
      maxN = 5;
      tracesFilename = "";
      title = "";
    }

    Settings(Settings source) {
      inputOption = source.inputOption;
      header = source.header;
      displayHistograms = source.displayHistograms.clone();
      alwaysRemoveOutliers = source.alwaysRemoveOutliers.clone();
      displayMsdHistogram = source.displayMsdHistogram;
      displayDHistogram = source.displayDHistogram;
      displayTraceLength = source.displayTraceLength;
      displayTraceSize = source.displayTraceSize;
      saveTraceDistances = source.saveTraceDistances;
      saveRawData = source.saveRawData;
      rawDataDirectory = source.rawDataDirectory;
      distancesFilename = source.distancesFilename;
      significanceLevel = source.significanceLevel;
      minFraction = source.minFraction;
      minDifference = source.minDifference;
      minN = source.minN;
      maxN = source.maxN;
      debugFitting = source.debugFitting;
      tracesFilename = source.tracesFilename;
      title = source.title;
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

    jumpDistanceParametersRef.set(null);

    extraOptions = ImageJUtils.isExtraOptions();
    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "No localisations in memory");
      return;
    }

    settings = Settings.load();
    // Saved by reference so just save now
    settings.save();

    final ArrayList<MemoryPeakResults> allResults = new ArrayList<>();

    // Option to pick multiple input datasets together using a list box.
    if ("multi".equals(arg) && !showMultiDialog(allResults)) {
      return;
    }

    // This shows the dialog for selecting trace options
    if (!showTraceDialog(allResults)) {
      return;
    }

    if (allResults.isEmpty()) {
      return;
    }

    ImageJUtils.log(TITLE + "...");

    // - Trace each single dataset (and store in memory)
    // - Combine trace results held in memory
    final Trace[] traces = getTraces(allResults);

    // -=-=-
    // Analyse the traces
    // -=-=-

    // Only show the second dialog if we have traces.
    // This still allows a zero entry in the results table.
    if (traces.length > 0 && !showDialog()) {
      return;
    }

    final int count = traces.length;
    double[] fitMsdResult = null;
    int numberOfDataPoints = 0;
    double[][] jdParams = null;
    if (count > 0) {
      calculatePrecision(traces, allResults.size() > 1);

      // --- MSD Analysis ---

      // Conversion constants
      final double px2ToUm2 = MathUtils.pow2(results.getCalibrationReader().getNmPerPixel()) / 1e6;
      final double px2ToUm2PerSecond = px2ToUm2 / exposureTime;

      // Get the maximum trace length
      int length = clusteringSettings.getMinimumTraceLength();
      if (!clusteringSettings.getTruncate()) {
        for (final Trace trace : traces) {
          if (length < trace.size()) {
            length = trace.size();
          }
        }
      }

      // Get the localisation error (4s^2) in um^2
      final double error =
          (clusteringSettings.getPrecisionCorrection()) ? 4 * precision * precision / 1e6 : 0;
      // Pre-calculate MSD correction factors. This accounts for the fact that the distance moved
      // in the start/end frames is reduced due to the averaging of the particle location over the
      // entire frame into a single point. The true MSD may be restored by applying a factor.
      // Note: These are used for the calculation of the diffusion coefficients per molecule and
      // the MSD passed to the Jump Distance analysis. However the error is not included in the
      // jump distance analysis so will be subtracted from the fitted D coefficients later.
      final double[] factors;
      if (clusteringSettings.getMsdCorrection()) {
        factors = new double[length];
        for (int t = 1; t < length; t++) {
          factors[t] = JumpDistanceAnalysis.getConversionfactor(t);
        }
      } else {
        factors = SimpleArrayUtils.newArray(length, 0.0, 1.0);
      }

      // Extract the mean-squared distance statistics
      final Statistics[] stats = new Statistics[length];
      for (int i = 0; i < stats.length; i++) {
        stats[i] = new Statistics();
      }

      final ArrayList<double[]> distances =
          (settings.saveTraceDistances || settings.displayTraceLength)
              ? new ArrayList<>(traces.length)
              : null;

      // Store all the jump distances at the specified interval
      final StoredDataStatistics jumpDistances = new StoredDataStatistics();
      final int jumpDistanceInterval = clusteringSettings.getJumpDistance();

      // Compute squared distances
      final StoredDataStatistics msdPerMoleculeAllVsAll = new StoredDataStatistics();
      final StoredDataStatistics msdPerMoleculeAdjacent = new StoredDataStatistics();
      for (final Trace trace : traces) {
        final PeakResultStoreList results = trace.getPoints();
        // Sum the MSD and the time
        final int traceLength =
            (clusteringSettings.getTruncate()) ? clusteringSettings.getMinimumTraceLength()
                : trace.size();

        // Get the mean for each time separation
        final double[] sumDistance = new double[traceLength + 1];
        final double[] sumTime = new double[sumDistance.length];

        // Do the distances to the origin (saving if necessary)
        final float x0 = results.get(0).getXPosition();
        final float y0 = results.get(0).getYPosition();
        if (distances != null) {
          final double[] msd = new double[traceLength - 1];
          for (int j = 1; j < traceLength; j++) {
            final int t = j;
            final double d = distance2(x0, y0, results.get(j));
            msd[j - 1] = px2ToUm2 * d;
            if (t == jumpDistanceInterval) {
              jumpDistances.add(msd[j - 1]);
            }
            sumDistance[t] += d;
            sumTime[t] += t;
          }
          distances.add(msd);
        } else {
          for (int j = 1; j < traceLength; j++) {
            final int t = j;
            final double d = distance2(x0, y0, results.get(j));
            if (t == jumpDistanceInterval) {
              jumpDistances.add(px2ToUm2 * d);
            }
            sumDistance[t] += d;
            sumTime[t] += t;
          }
        }

        if (clusteringSettings.getInternalDistances()) {
          // Do the internal distances
          for (int i = 1; i < traceLength; i++) {
            final float x = results.get(i).getXPosition();
            final float y = results.get(i).getYPosition();
            for (int j = i + 1; j < traceLength; j++) {
              final int t = j - i;
              final double d = distance2(x, y, results.get(j));
              if (t == jumpDistanceInterval) {
                jumpDistances.add(px2ToUm2 * d);
              }
              sumDistance[t] += d;
              sumTime[t] += t;
            }
          }

          // Add the average distance per time separation to the population
          for (int t = 1; t < traceLength; t++) {
            // Note: (traceLength - t) == count
            stats[t].add(sumDistance[t] / (traceLength - t));
          }
        } else {
          // Add the distance per time separation to the population
          for (int t = 1; t < traceLength; t++) {
            stats[t].add(sumDistance[t]);
          }
        }

        // Fix this for the precision and MSD adjustment.
        // It may be necessary to:
        // - sum the raw distances for each time interval (this is sumDistance[t])
        // - subtract the precision error
        // - apply correction factor for the n-frames to get actual MSD
        // - sum the actual MSD

        double sumD = 0;
        final double sumD_adjacent = Math.max(0, sumDistance[1] - error) * factors[1];
        double sumT = 0;
        final double sumT_adjacent = sumTime[1];
        for (int t = 1; t < traceLength; t++) {
          sumD += Math.max(0, sumDistance[t] - error) * factors[t];
          sumT += sumTime[t];
        }

        // Calculate the average displacement for the trace (do not simply use the largest
        // time separation since this will miss moving molecules that end up at the origin)

        msdPerMoleculeAllVsAll.add(px2ToUm2PerSecond * sumD / sumT);
        msdPerMoleculeAdjacent.add(px2ToUm2PerSecond * sumD_adjacent / sumT_adjacent);
      }

      StoredDataStatistics dperMoleculeAllVsAll = null;
      StoredDataStatistics dperMoleculeAdjacent = null;
      if (settings.saveTraceDistances
          || (clusteringSettings.getShowHistograms() && settings.displayDHistogram)) {
        dperMoleculeAllVsAll = calculateDiffusionCoefficient(msdPerMoleculeAllVsAll);
        dperMoleculeAdjacent = calculateDiffusionCoefficient(msdPerMoleculeAdjacent);
      }

      if (settings.saveTraceDistances) {
        saveTraceDistances(traces.length, distances, msdPerMoleculeAllVsAll, msdPerMoleculeAdjacent,
            dperMoleculeAllVsAll, dperMoleculeAdjacent);
      }

      if (settings.displayTraceLength) {
        final StoredDataStatistics lengths = calculateTraceLengths(distances);
        showHistogram(lengths, "Trace length (um)");
      }

      if (settings.displayTraceSize) {
        final StoredDataStatistics sizes = calculateTraceSizes(traces);
        showHistogram(sizes, "Trace size", true);
      }

      // Plot the per-trace histogram of MSD and D
      if (clusteringSettings.getShowHistograms()) {
        if (settings.displayMsdHistogram) {
          showHistogram(msdPerMoleculeAllVsAll, "MSD/Molecule (all-vs-all)");
          showHistogram(msdPerMoleculeAdjacent, "MSD/Molecule (adjacent)");
        }
        if (settings.displayDHistogram) {
          showHistogram(dperMoleculeAllVsAll, "D/Molecule (all-vs-all)");
          showHistogram(dperMoleculeAdjacent, "D/Molecule (adjacent)");
        }
      }

      // Calculate the mean squared distance (MSD)
      final double[] x = new double[stats.length];
      final double[] y = new double[x.length];
      final double[] sd = new double[x.length];
      // Intercept is the 4s^2 (in um^2)
      y[0] = 4 * precision * precision / 1e6;
      for (int i = 1; i < stats.length; i++) {
        x[i] = i * exposureTime;
        y[i] = stats[i].getMean() * px2ToUm2;
        // sd[i] = stats[i].getStandardDeviation() * px2ToUm2;
        sd[i] = stats[i].getStandardError() * px2ToUm2;
      }

      final String title = TITLE + " MSD";
      final Plot plot = plotMsd(x, y, sd, title);

      // Fit the MSD using a linear fit
      fitMsdResult = fitMsd(x, y, title, plot);

      // Jump Distance analysis
      if (settings.saveRawData) {
        saveStatistics(jumpDistances, "Jump Distance", "Distance (um^2)", false);
      }

      // Calculate the cumulative jump-distance histogram
      final double[][] jdHistogram =
          JumpDistanceAnalysis.cumulativeHistogram(jumpDistances.getValues());

      // Always show the jump distance histogram
      jdTitle = TITLE + " Jump Distance";
      jdPlot = new Plot(jdTitle, "Distance (um^2)", "Cumulative Probability");
      jdPlot.addPoints(jdHistogram[0], jdHistogram[1], Plot.LINE);
      display(jdTitle, jdPlot);

      // Fit Jump Distance cumulative probability
      numberOfDataPoints = jumpDistances.getN();
      jdParams = fitJumpDistance(jumpDistances, jdHistogram);
      jumpDistanceParametersRef.set(jdParams);
    }

    summarise(traces, fitMsdResult, numberOfDataPoints, jdParams);
  }

  /**
   * Calculate trace lengths.
   *
   * @param distances the distances for each trace
   * @return the trace lengths
   */
  private static StoredDataStatistics calculateTraceLengths(ArrayList<double[]> distances) {
    final StoredDataStatistics lengths = new StoredDataStatistics();
    for (final double[] trace : distances) {
      double sum = 0;
      for (final double d : trace) {
        sum += Math.sqrt(d);
      }
      lengths.add(sum);
    }
    return lengths;
  }

  private static StoredDataStatistics calculateTraceSizes(Trace[] traces) {
    final StoredDataStatistics sizes = new StoredDataStatistics();
    for (final Trace trace : traces) {
      sizes.add(trace.size());
    }
    return sizes;
  }

  private void display(String title, Plot plot) {
    ImageJUtils.display(title, plot, windowOrganiser);
  }

  private void showHistogram(StoredDataStatistics stats, String title, boolean integerData) {
    showHistogram(stats, title, false, false, integerData);
  }

  private void showHistogram(StoredDataStatistics stats, String title) {
    showHistogram(stats, title, false, false, false);
  }

  private void showHistogram(StoredDataStatistics stats, String title, boolean alwaysRemoveOutliers,
      boolean rounded, boolean integerData) {
    if (settings.saveRawData) {
      saveStatistics(stats, title, title, rounded);
    }

    new HistogramPlotBuilder(TITLE, stats, title).setIntegerBins(integerData)
        .setRemoveOutliersOption(
            (clusteringSettings.getRemoveOutliers() || alwaysRemoveOutliers) ? 1 : 0)
        .setNumberOfBins(clusteringSettings.getHistogramBins()).show(windowOrganiser);
  }

  /**
   * Calculate the average precision of localisation in the traces.
   *
   * @param traces the traces
   * @param multi true if the traces were from multiple input results
   */
  private void calculatePrecision(Trace[] traces, boolean multi) {
    // Check the diffusion simulation for a precision
    if (DiffusionRateTest.isSimulated(results.getName()) && !multi) {
      precision = DiffusionRateTest.getLastSimulationPrecision();
    } else {
      precision = 999;
      try {
        final Gaussian2DPeakResultCalculator calculator = Gaussian2DPeakResultHelper.create(
            results.getPsf(), results.getCalibration(), Gaussian2DPeakResultHelper.LSE_PRECISION);
        // Get the average precision of the localisations
        precision = 0;
        int count = 0;
        for (final Trace trace : traces) {
          for (int k = 0; k < trace.size(); k++) {
            final PeakResult r = trace.get(k);
            precision += calculator.getLsePrecision(r.getParameters(), r.getNoise());
          }
          count += trace.size();
        }
        precision /= count;
      } catch (final ConfigurationException ex) {
        // Ignore this and we will ask the user for the precision
      }
    }

    if (precision > 100) {
      final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
      gd.addMessage("The average precision of the traced results is "
          + MathUtils.rounded(precision, 4) + " nm.\nPlease verify the precision.");
      gd.addSlider("Precision (nm)", 5, 100, precision);
      gd.showDialog();
      if (!(gd.wasCanceled() || gd.invalidNumber())) {
        precision = Math.abs(gd.getNextNumber());
      }
    }
  }

  /**
   * Calculate the diffusion coefficient (D) of the molecule. This is done by using the mean-squared
   * deviation between frames divided by the time interval (delta) between frames. This is divided
   * by 4 to produce the diffusion coefficient from two-dimensional distance analysis.
   *
   * <p>See Uphoff, et al, 2013. Single-molecule DNA repair in live bacteria, PNAS 110, 8063-8068
   *
   * @param msdPerMoleculeAdjacent the MSD per molecule adjacent
   * @return The D per molecule
   */
  private static StoredDataStatistics
      calculateDiffusionCoefficient(StoredDataStatistics msdPerMoleculeAdjacent) {
    final StoredDataStatistics dPerMolecule = new StoredDataStatistics();
    final double diffusionCoefficientConversion = 1.0 / 4.0;
    for (final double msd : msdPerMoleculeAdjacent.getValues()) {
      dPerMolecule.add(msd * diffusionCoefficientConversion);
    }
    return dPerMolecule;
  }

  private void saveTraceDistances(int traceCount, ArrayList<double[]> distances,
      StoredDataStatistics msdPerMolecule, StoredDataStatistics msdPerMoleculeAdjacent,
      StoredDataStatistics dstarPerMolecule, StoredDataStatistics dstarPerMoleculeAdjacent) {
    settings.distancesFilename =
        ImageJUtils.getFilename("Trace_Distances_File", settings.distancesFilename);
    if (settings.distancesFilename != null) {
      settings.distancesFilename = FileUtils.replaceExtension(settings.distancesFilename, "xls");

      try (BufferedWriter out = Files.newBufferedWriter(Paths.get(settings.distancesFilename))) {
        final double[] msd = msdPerMolecule.getValues();
        final double[] msd2 = msdPerMoleculeAdjacent.getValues();
        final double[] dStar = dstarPerMolecule.getValues();
        final double[] dStar2 = dstarPerMoleculeAdjacent.getValues();
        out.write(String.format("#%d traces : Precision = %s nm : Exposure time = %s s", traceCount,
            MathUtils.rounded(precision, 4), MathUtils.rounded(exposureTime, 4)));
        out.newLine();
        out.write(String.format(
            "#TraceId\tMSD all-vs-all (um^2/s)\tMSD adjacent (um^2/s)\tD all-vs-all(um^2/s)\t"
                + "D adjacent(um^2/s)\tDistances (um^2) per %ss ... ",
            MathUtils.rounded(exposureTime, 4)));
        out.newLine();
        for (int i = 0; i < msd.length; i++) {
          out.write(Integer.toString(i + 1));
          out.write('\t');
          out.write(MathUtils.rounded(msd[i], 4));
          out.write('\t');
          out.write(MathUtils.rounded(msd2[i], 4));
          out.write('\t');
          out.write(MathUtils.rounded(dStar[i], 4));
          out.write('\t');
          out.write(MathUtils.rounded(dStar2[i], 4));
          for (final double d : distances.get(i)) {
            out.write('\t');
            out.write(MathUtils.rounded(d, 4));
          }
          out.newLine();
        }
      } catch (final IOException ex) {
        Logger.getLogger(getClass().getName()).log(Level.WARNING,
            "Failed to save trace distances: " + settings.distancesFilename, ex);
      }
    }
  }

  private void saveMsd(double[] x, double[] y, double[] se) {
    if (!directoryChosen) {
      settings.rawDataDirectory =
          ImageJUtils.getDirectory("Data_directory", settings.rawDataDirectory);
    }
    directoryChosen = true;
    if (settings.rawDataDirectory == null) {
      return;
    }
    final Path filename = Paths.get(settings.rawDataDirectory, "MSD.txt");

    try (BufferedWriter out = Files.newBufferedWriter(filename)) {
      out.write("Time (s)\tDistance (um^2)\tS.E.");
      out.newLine();
      for (int i = 0; i < x.length; i++) {
        out.write(MathUtils.rounded(x[i]));
        out.write('\t');
        out.write(Double.toString(y[i]));
        out.write('\t');
        out.write(Double.toString(se[i]));
        out.newLine();
      }
    } catch (final IOException ex) {
      Logger.getLogger(getClass().getName()).log(Level.WARNING,
          "Failed to save MSD file: " + filename.toString(), ex);
    }
  }

  private void saveStatistics(StoredDataStatistics stats, String title, String label,
      boolean rounded) {
    if (!directoryChosen) {
      settings.rawDataDirectory =
          ImageJUtils.getDirectory("Data_directory", settings.rawDataDirectory);
    }
    directoryChosen = true;
    if (settings.rawDataDirectory == null) {
      return;
    }
    final String filename =
        settings.rawDataDirectory + title.replace("/", " per ").replace("*", "star") + ".txt";

    try (BufferedWriter out = Files.newBufferedWriter(Paths.get(filename))) {
      out.write(label);
      out.newLine();
      final double[] data = stats.getValues();
      Arrays.sort(data);
      if (rounded) {
        for (final double d : data) {
          out.write(MathUtils.rounded(d, 4));
          out.newLine();
        }
      } else {
        for (final double d : data) {
          out.write(Double.toString(d));
          out.newLine();
        }
      }
    } catch (final IOException ex) {
      Logger.getLogger(getClass().getName()).log(Level.WARNING,
          "Failed to save statistics: " + filename, ex);
    }
  }

  private void saveFit(double[] x, double[] y, String title) {
    if (!directoryChosen) {
      settings.rawDataDirectory =
          ImageJUtils.getDirectory("Data_directory", settings.rawDataDirectory);
    }
    directoryChosen = true;
    if (settings.rawDataDirectory == null) {
      return;
    }
    final String filename = settings.rawDataDirectory + "Fit." + title + ".txt";

    try (BufferedWriter out = Files.newBufferedWriter(Paths.get(filename))) {
      out.write("JumpDistance\tCumulativeP");
      out.newLine();
      for (int i = 0; i < x.length; i++) {
        out.write(Double.toString(x[i]));
        out.write('\t');
        out.write(Double.toString(y[i]));
        out.newLine();
      }
    } catch (final IOException ex) {
      Logger.getLogger(getClass().getName()).log(Level.WARNING, "Failed to save fit: " + filename,
          ex);
    }
  }

  private static double distance2(final float x, final float y, PeakResult r2) {
    final double dx = x - r2.getXPosition();
    final double dy = y - r2.getYPosition();
    return dx * dx + dy * dy;
  }

  /**
   * Split traces to contiguous traces and filter traces that are not the minimum length. Re-assigns
   * the ID for the output traces.
   *
   * @param name the name
   * @param traces the traces
   * @param minimumTraceLength the minimum trace length
   * @param ignoreEnds the ignore ends
   * @return The new traces
   */
  private static Trace[] filterTraces(String name, Trace[] traces, int minimumTraceLength,
      boolean ignoreEnds) {
    final LocalList<Trace> list = new LocalList<>(traces.length);

    final int minLength = (ignoreEnds) ? minimumTraceLength + 2 : minimumTraceLength;
    final Consumer<Trace> action = (ignoreEnds) ? t -> {
      if (t.size() >= minLength) {
        t.removeEnds();
        list.add(t);
        t.setId(list.size());
      }
    } : t -> {
      if (t.size() >= minLength) {
        list.add(t);
        t.setId(list.size());
      }
    };
    final Consumer<ArrayPeakResultStore> action2 = (ignoreEnds) ? r -> {
      if (r.size() >= minLength) {
        r.remove(0);
        r.remove(r.size() - 1);
        final Trace t = new Trace(r);
        list.add(t);
        t.setId(list.size());
      }
    } : r -> {
      if (r.size() >= minLength) {
        final Trace t = new Trace(r);
        list.add(t);
        t.setId(list.size());
      }
    };

    final ArrayPeakResultStore results = new ArrayPeakResultStore(11);
    for (final Trace trace : traces) {
      if (trace.size() < minLength) {
        // Too short
        continue;
      }
      if (trace.size() == trace.getTail().getFrame() - trace.getHead().getFrame() + 1) {
        // Contiguous
        action.accept(trace);
      } else {
        // Split the trace
        int t1 = trace.getHead().getFrame();
        for (int i = 0; i < trace.size(); i++) {
          final PeakResult peak = trace.get(i);
          final int t2 = peak.getFrame();
          if (t2 - t1 > 1) {
            // Non-contiguous
            action2.accept(results);
            results.clear();
          }
          t1 = t2;
          results.add(peak);
        }
        // Final trace
        action2.accept(results);
        results.clear();
      }
    }

    ImageJUtils.log(
        "Filtered results '%s' : %s split and filtered to %d using "
            + "minimum length %d (Ignore ends = %b)",
        name, TextUtils.pleural(traces.length, "trace"), list.size(), minimumTraceLength,
        ignoreEnds);
    return list.toArray(new Trace[0]);
  }

  private String createSettingsComment() {
    return String.format("Molecule tracing : distance-threshold = %f nm",
        clusteringSettings.getDistanceThreshold());
  }

  private void summarise(Trace[] traces, double[] fitMsdResult, int n, double[][] jdParams) {
    IJ.showStatus("Calculating summary ...");

    final Statistics[] stats = new Statistics[Settings.NAMES.length];
    for (int i = 0; i < stats.length; i++) {
      stats[i] =
          (clusteringSettings.getShowHistograms()) ? new StoredDataStatistics() : new Statistics();
    }
    for (final Trace trace : traces) {
      stats[Settings.T_ON].add(trace.getOnTime() * exposureTime);
      final double signal = trace.getSignal() / results.getGain();
      stats[Settings.TOTAL_SIGNAL].add(signal);
      stats[Settings.SIGNAL_PER_FRAME].add(signal / trace.size());
    }

    // Add to the summary table
    final StringBuilder sb = new StringBuilder(settings.title);
    sb.append('\t').append(createCombinedName());
    sb.append('\t');
    sb.append(MathUtils.rounded(exposureTime * 1000, 3)).append('\t');
    appendClusteringSettings(sb).append('\t');
    sb.append(clusteringSettings.getMinimumTraceLength()).append('\t');
    sb.append(clusteringSettings.getIgnoreEnds()).append('\t');
    sb.append(clusteringSettings.getTruncate()).append('\t');
    sb.append(clusteringSettings.getInternalDistances()).append('\t');
    sb.append(clusteringSettings.getFitLength()).append('\t');
    sb.append(clusteringSettings.getMsdCorrection()).append('\t');
    sb.append(clusteringSettings.getPrecisionCorrection()).append('\t');
    sb.append(clusteringSettings.getMle()).append('\t');
    sb.append(traces.length).append('\t');
    sb.append(MathUtils.rounded(precision, 4)).append('\t');
    double diffCoeff = 0; // D
    double precision = 0;
    if (fitMsdResult != null) {
      diffCoeff = fitMsdResult[0];
      precision = fitMsdResult[1];
    }
    sb.append(MathUtils.rounded(diffCoeff, 4)).append('\t');
    sb.append(MathUtils.rounded(precision * 1000, 4)).append('\t');
    sb.append(MathUtils.rounded(clusteringSettings.getJumpDistance() * exposureTime)).append('\t');
    sb.append(n).append('\t');
    sb.append(MathUtils.rounded(beta, 4)).append('\t');
    if (jdParams == null) {
      sb.append("\t\t\t");
    } else {
      sb.append(format(jdParams[0])).append('\t');
      sb.append(format(jdParams[1])).append('\t');
      sb.append(MathUtils.rounded(fitValue)).append('\t');
    }

    for (int i = 0; i < stats.length; i++) {
      sb.append(MathUtils.rounded(stats[i].getMean(), 3)).append('\t');
    }
    createSummaryTable().accept(sb.toString());

    if (java.awt.GraphicsEnvironment.isHeadless()) {
      return;
    }

    if (clusteringSettings.getShowHistograms()) {
      IJ.showStatus("Calculating histograms ...");

      for (int i = 0; i < Settings.NAMES.length; i++) {
        if (settings.displayHistograms[i]) {
          showHistogram((StoredDataStatistics) stats[i], Settings.NAMES[i],
              settings.alwaysRemoveOutliers[i], Settings.ROUNDED[i], false);
        }
      }
    }

    windowOrganiser.tile();

    IJ.showStatus("Finished " + TITLE);
  }

  private StringBuilder appendClusteringSettings(StringBuilder sb) {
    if (clusteringSettings.getTraceDiffusionMode() == 1) {
      sb.append("Dynamic MTT : D=")
          .append(MathUtils.rounded(clusteringSettings.getDiffusionCoefficentMaximum(), 3))
          .append("um^2/s");
      if (!clusteringSettings.getDisableLocalDiffusionModel()
          || !clusteringSettings.getDisableIntensityModel()) {
        sb.append(", window=").append(clusteringSettings.getTemporalWindow());
      }
      if (!clusteringSettings.getDisableLocalDiffusionModel()) {
        sb.append(", wLocal=")
            .append(MathUtils.rounded(clusteringSettings.getLocalDiffusionWeight(), 2));
      }
      if (!clusteringSettings.getDisableIntensityModel()) {
        sb.append(", wOn=").append(MathUtils.rounded(clusteringSettings.getOnIntensityWeight(), 2));
      }
      sb.append(", decay=")
          .append(MathUtils.rounded(clusteringSettings.getDisappearanceDecayFactor(), 2))
          .append(", disappear=").append(clusteringSettings.getDisappearanceThreshold());
    } else {
      sb.append("Nearest Neighbour : ")
          .append(MathUtils.rounded(clusteringSettings.getDistanceThreshold(), 3)).append("nm");
      if (clusteringSettings.getDistanceExclusion() > clusteringSettings.getDistanceThreshold()) {
        sb.append(" (ex. ").append(MathUtils.rounded(clusteringSettings.getDistanceExclusion(), 3))
            .append("nm)");
      }
    }
    return sb;
  }

  private static String format(double[] jumpD) {
    if (jumpD == null || jumpD.length == 0) {
      return "";
    }
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < jumpD.length; i++) {
      if (i != 0) {
        sb.append(", ");
      }
      sb.append(MathUtils.rounded(jumpD[i], 4));
    }
    return sb.toString();
  }

  private Consumer<String> createSummaryTable() {
    if (java.awt.GraphicsEnvironment.isHeadless()) {
      if (settings.header == null) {
        settings.header = createHeader();
        IJ.log(settings.header);
      }
      return IJ::log;
    }
    return ImageJUtils.refresh(summaryTableRef,
        () -> new TextWindow(TITLE + " Data Summary", createHeader(), "", 800, 300))::append;
  }

  private static String createHeader() {
    final StringBuilder sb =
        new StringBuilder("Title\tDataset\tExposure time (ms)\tTrace settings\t");
    sb.append("Min.Length\tIgnoreEnds\tTruncate\tInternal\tFit Length");
    sb.append("\tMSD corr.\ts corr.\tMLE\tTraces\ts (nm)\tD (um^2/s)\tfit s (nm)");
    sb.append("\tJump Distance (s)\tN\tBeta\tJump D (um^2/s)\tFractions\tFit Score");
    for (int i = 0; i < Settings.NAMES.length; i++) {
      sb.append('\t').append(Settings.NAMES[i]);
    }
    return sb.toString();
  }

  private boolean showTraceDialog(ArrayList<MemoryPeakResults> allResults) {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("trace-diffusion"));

    if (!multiMode) {
      ResultsManager.addInput(gd, settings.inputOption, InputSource.MEMORY);
    }

    clusteringSettings = SettingsManager.readClusteringSettings(0).toBuilder();

    gd.addChoice("Mode", TRACE_MODE, clusteringSettings.getTraceDiffusionMode(),
        new OptionListener<Integer>() {
          @Override
          public boolean collectOptions(Integer value) {
            clusteringSettings.setTraceDiffusionMode(value);
            return collectOptions(false);
          }

          @Override
          public boolean collectOptions() {
            return collectOptions(true);
          }

          private boolean collectOptions(boolean silent) {
            final ExtendedGenericDialog egd =
                new ExtendedGenericDialog("Trace diffusion options", null);
            // Only 2 modes
            if (clusteringSettings.getTraceDiffusionMode() == 1) {
              // Dynamic Multiple Target Tracing
              final TextField tfD = egd.addAndGetNumericField("Diffusion_coefficient",
                  clusteringSettings.getDiffusionCoefficentMaximum(), 3, 6, "um^2/s");
              final TextField tfW = egd.addAndGetNumericField("Temporal_window",
                  clusteringSettings.getTemporalWindow(), 0, 6, "frames");
              final TextField tfLdw = egd.addAndGetNumericField("Local_diffusion_weight",
                  clusteringSettings.getLocalDiffusionWeight(), 2);
              final TextField tfOiw = egd.addAndGetNumericField("On_intensity_weight",
                  clusteringSettings.getOnIntensityWeight(), 2);
              final TextField tfDdf = egd.addAndGetNumericField("Disappearance_decay_factor",
                  clusteringSettings.getDisappearanceDecayFactor(), 0, 6, "frames");
              final TextField tfDt = egd.addAndGetNumericField("Disappearance_threshold",
                  clusteringSettings.getDisappearanceThreshold(), 0, 6, "frames");
              final Checkbox cbDld = egd.addAndGetCheckbox("Disable_local_diffusion_model",
                  clusteringSettings.getDisableLocalDiffusionModel());
              final Checkbox cbDim = egd.addAndGetCheckbox("Disable_intensity_model",
                  clusteringSettings.getDisableIntensityModel());

              // Allow reset to defaults
              egd.addAndGetButton("Defaults", e -> {
                final DmttConfiguration config = DmttConfiguration.newBuilder(1).build();
                tfD.setText(String.valueOf(clusteringSettings.getDiffusionCoefficentMaximum()));
                tfW.setText(String.valueOf(config.getTemporalWindow()));
                tfLdw.setText(String.valueOf(config.getLocalDiffusionWeight()));
                tfOiw.setText(String.valueOf(config.getOnIntensityWeight()));
                tfDdf.setText(String.valueOf(config.getDisappearanceDecayFactor()));
                tfDt.setText(String.valueOf(config.getDisappearanceThreshold()));
                cbDld.setState(config.isDisableLocalDiffusionModel());
                cbDim.setState(config.isDisableIntensityModel());
              });
            } else {
              // Nearest Neighbour
              egd.addNumericField("Distance_Threshold (nm)",
                  clusteringSettings.getDistanceThreshold(), 0);
              egd.addNumericField("Distance_Exclusion (nm)",
                  clusteringSettings.getDistanceExclusion(), 0);
            }

            egd.setSilent(silent);
            egd.showDialog(true, gd);
            if (egd.wasCanceled()) {
              return false;
            }

            if (clusteringSettings.getTraceDiffusionMode() == 1) {
              // Dynamic Multiple Target Tracing
              clusteringSettings.setDiffusionCoefficentMaximum(egd.getNextNumber());
              clusteringSettings.setTemporalWindow((int) egd.getNextNumber());
              clusteringSettings.setLocalDiffusionWeight(egd.getNextNumber());
              clusteringSettings.setOnIntensityWeight(egd.getNextNumber());
              clusteringSettings.setDisappearanceDecayFactor(egd.getNextNumber());
              clusteringSettings.setDisappearanceThreshold((int) egd.getNextNumber());
              clusteringSettings.setDisableLocalDiffusionModel(egd.getNextBoolean());
              clusteringSettings.setDisableIntensityModel(egd.getNextBoolean());
            } else {
              // Nearest Neighbour
              clusteringSettings.setDistanceThreshold(egd.getNextNumber());
              clusteringSettings.setDistanceExclusion(Math.abs(egd.getNextNumber()));
            }
            return true;
          }
        });

    gd.addSlider("Min_trace_length", 2, 20, clusteringSettings.getMinimumTraceLength());
    gd.addCheckbox("Ignore_ends", clusteringSettings.getIgnoreEnds());
    gd.addCheckbox("Save_traces", clusteringSettings.getSaveTraces());

    gd.showDialog();

    if (gd.wasCanceled() || !readTraceDialog(gd)) {
      return false;
    }

    // Update the settings
    SettingsManager.writeSettings(clusteringSettings.build());

    // Load the results
    if (!multiMode) {
      final MemoryPeakResults results =
          ResultsManager.loadInputResults(settings.inputOption, true, null, null);
      if (MemoryPeakResults.isEmpty(results)) {
        IJ.error(TITLE, "No results could be loaded");
        IJ.showStatus("");
        return false;
      }

      if (!checkCalibration(results)) {
        return false;
      }

      allResults.add(results);
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

  private boolean readTraceDialog(ExtendedGenericDialog gd) {
    if (!multiMode) {
      settings.inputOption = ResultsManager.getInputSource(gd);
    }
    clusteringSettings.setTraceDiffusionMode(gd.getNextChoiceIndex());
    clusteringSettings.setMinimumTraceLength((int) Math.abs(gd.getNextNumber()));
    clusteringSettings.setIgnoreEnds(gd.getNextBoolean());
    clusteringSettings.setSaveTraces(gd.getNextBoolean());

    gd.collectOptions();

    if (gd.invalidNumber()) {
      return false;
    }

    // Check arguments
    try {
      if (clusteringSettings.getTraceDiffusionMode() == 1) {
        // Validate using the Builder
        createDmttConfiguration();
      } else {
        ParameterUtils.isAboveZero("Distance threshold", clusteringSettings.getDistanceThreshold());
      }
      ParameterUtils.isAbove("Min trace length", clusteringSettings.getMinimumTraceLength(), 1);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
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
    return DmttConfiguration.newBuilder(clusteringSettings.getDiffusionCoefficentMaximum())
        .setTemporalWindow(clusteringSettings.getTemporalWindow())
        .setLocalDiffusionWeight(clusteringSettings.getLocalDiffusionWeight())
        .setOnIntensityWeight(clusteringSettings.getOnIntensityWeight())
        .setDisappearanceDecayFactor(clusteringSettings.getDisappearanceDecayFactor())
        .setDisappearanceThreshold(clusteringSettings.getDisappearanceThreshold())
        .setDisableLocalDiffusionModel(clusteringSettings.getDisableLocalDiffusionModel())
        .setDisableIntensityModel(clusteringSettings.getDisableIntensityModel()).build();
  }

  private Trace[] getTraces(ArrayList<MemoryPeakResults> allResults) {
    this.results = allResults.get(0);

    // Results should be checked for calibration by this point
    exposureTime = results.getCalibrationReader().getExposureTime() / 1000;

    final Function<MemoryPeakResults, Trace[]> traceFunction = createTraceFunction();

    final ArrayList<Trace> allTraces = new ArrayList<>();
    additionalDatasets = -1;
    for (final MemoryPeakResults r : allResults) {
      additionalDatasets++;

      Trace[] traces = traceFunction.apply(r);

      traces = filterTraces(r.getName(), traces, clusteringSettings.getMinimumTraceLength(),
          clusteringSettings.getIgnoreEnds());
      allTraces.addAll(Arrays.asList(traces));

      // --- Save results ---
      if (traces.length != 0) {
        // Save the traces to memory
        TraceMolecules.saveResults(r, traces, "Tracks");

        if (clusteringSettings.getSaveTraces()) {
          // Sort traces by time to assist the results source in extracting frames sequentially.
          // Do this before saving to assist in debugging using the saved traces file.
          TraceMolecules.sortByTime(traces);
          final String newFilename = TraceMolecules.saveTraces(r, traces, createSettingsComment(),
              settings.tracesFilename, additionalDatasets);
          // Only keep the main filename in memory
          if (additionalDatasets == 0) {
            settings.tracesFilename = newFilename;
          }
        }
      }
    }

    final Trace[] all = allTraces.toArray(new Trace[0]);

    if (additionalDatasets != 0) {
      ImageJUtils.log("Multiple inputs provide %d traces", allTraces.size());

      final MemoryPeakResults tracedResults =
          TraceManager.toPeakResults(all, results.getCalibration(), true);
      tracedResults.copySettings(results);
      tracedResults.setName(createCombinedName() + " Tracks");
      MemoryPeakResults.addResults(tracedResults);
    }

    return all;
  }

  /**
   * Creates the trace function for the configured trace diffusion mode.
   *
   * @return the function
   */
  private Function<MemoryPeakResults, Trace[]> createTraceFunction() {
    if (clusteringSettings.getTraceDiffusionMode() == 1) {
      final DmttConfiguration config = createDmttConfiguration();
      return r -> new DynamicMultipleTargetTracing(r).traceMolecules(config).toArray(new Trace[0]);
    }

    // Nearest neighbour

    // Convert from NM to the native units of the results
    final Converter c =
        CalibrationHelper.getDistanceConverter(results.getCalibration(), DistanceUnit.NM);
    final double distanceThreshold = c.convertBack(clusteringSettings.getDistanceThreshold());
    final double distanceExclusion = c.convertBack(clusteringSettings.getDistanceExclusion());

    return r -> {
      final TraceManager manager = new TraceManager(r);

      // Run the tracing
      manager.setTracker(SimpleImageJTrackProgress.getInstance());
      manager.setDistanceExclusion(distanceExclusion);
      manager.traceMolecules(distanceThreshold, 1);
      return manager.getTraces();
    };
  }

  private String createCombinedName() {
    if (additionalDatasets > 0) {
      return results.getName() + " + " + TextUtils.pleural(additionalDatasets, "other");
    }
    return results.getName();
  }

  private boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("trace-diffusion"));

    clusteringSettings = SettingsManager.readClusteringSettings(0).toBuilder();

    gd.addCheckbox("Truncate_traces", clusteringSettings.getTruncate());
    gd.addCheckbox("Internal_distances", clusteringSettings.getInternalDistances());
    // gd.addCheckbox("Sub-sample_distances", settings.subSampledDistances);
    gd.addSlider("Fit_length", 2, 20, clusteringSettings.getFitLength());
    gd.addCheckbox("MSD_correction", clusteringSettings.getMsdCorrection());
    gd.addCheckbox("Precision_correction", clusteringSettings.getPrecisionCorrection());
    gd.addCheckbox("Maximum_likelihood", clusteringSettings.getMle());
    gd.addSlider("MLE_significance_level", 0, 0.5, settings.significanceLevel);
    gd.addSlider("Fit_restarts", 0, 10, clusteringSettings.getFitRestarts());
    gd.addSlider("Jump_distance", 1, 20, clusteringSettings.getJumpDistance());
    gd.addSlider("Minimum_difference", 0, 10, settings.minDifference);
    gd.addSlider("Minimum_fraction", 0, 1, settings.minFraction);
    if (extraOptions) {
      gd.addSlider("Minimum_N", 1, 10, settings.minN);
    }
    gd.addSlider("Maximum_N", 2, 10, settings.maxN);
    gd.addCheckbox("Debug_fitting", settings.debugFitting);
    gd.addCheckbox("Save_trace_distances", settings.saveTraceDistances);
    gd.addCheckbox("Save_raw_data", settings.saveRawData);
    gd.addCheckbox("Show_histograms", clusteringSettings.getShowHistograms());
    gd.addStringField("Title", settings.title);

    gd.showDialog();

    if (gd.wasCanceled() || !readDialog(gd)) {
      return false;
    }

    // Update the settings
    SettingsManager.writeSettings(clusteringSettings.build());

    return true;
  }

  private boolean readDialog(ExtendedGenericDialog gd) {
    clusteringSettings.setTruncate(gd.getNextBoolean());
    clusteringSettings.setInternalDistances(gd.getNextBoolean());
    // settings.subSampledDistances = gd.getNextBoolean();
    clusteringSettings.setFitLength((int) Math.abs(gd.getNextNumber()));
    clusteringSettings.setMsdCorrection(gd.getNextBoolean());
    clusteringSettings.setPrecisionCorrection(gd.getNextBoolean());
    clusteringSettings.setMle(gd.getNextBoolean());
    settings.significanceLevel = Math.abs(gd.getNextNumber());
    clusteringSettings.setFitRestarts((int) Math.abs(gd.getNextNumber()));
    clusteringSettings.setJumpDistance((int) Math.abs(gd.getNextNumber()));
    settings.minDifference = Math.abs(gd.getNextNumber());
    settings.minFraction = Math.abs(gd.getNextNumber());
    if (extraOptions) {
      myMinN = settings.minN = (int) Math.abs(gd.getNextNumber());
    }
    settings.maxN = (int) Math.abs(gd.getNextNumber());
    settings.debugFitting = gd.getNextBoolean();
    settings.saveTraceDistances = gd.getNextBoolean();
    settings.saveRawData = gd.getNextBoolean();
    clusteringSettings.setShowHistograms(gd.getNextBoolean());
    settings.title = gd.getNextString();

    if (gd.invalidNumber()) {
      return false;
    }

    if (clusteringSettings.getShowHistograms()) {
      gd = new ExtendedGenericDialog(TITLE);
      gd.addMessage("Select the histograms to display");
      gd.addCheckbox("Remove_outliers", clusteringSettings.getRemoveOutliers());
      gd.addNumericField("Histogram_bins", clusteringSettings.getHistogramBins(), 0);
      for (int i = 0; i < settings.displayHistograms.length; i++) {
        gd.addCheckbox(Settings.NAMES[i].replace(' ', '_'), settings.displayHistograms[i]);
      }
      gd.addCheckbox("MSD/Molecule", settings.displayMsdHistogram);
      gd.addCheckbox("D/Molecule", settings.displayDHistogram);
      gd.addCheckbox("Trace_length", settings.displayTraceLength);
      gd.addCheckbox("Trace_size", settings.displayTraceSize);
      gd.showDialog();
      if (gd.wasCanceled()) {
        return false;
      }
      clusteringSettings.setRemoveOutliers(gd.getNextBoolean());
      clusteringSettings.setHistogramBins((int) Math.abs(gd.getNextNumber()));
      for (int i = 0; i < settings.displayHistograms.length; i++) {
        settings.displayHistograms[i] = gd.getNextBoolean();
      }
      settings.displayMsdHistogram = gd.getNextBoolean();
      settings.displayDHistogram = gd.getNextBoolean();
      settings.displayTraceLength = gd.getNextBoolean();
      settings.displayTraceSize = gd.getNextBoolean();
    }

    // Check arguments
    try {
      // Parameters.isAboveZero("Histogram bins", settings.getHistogramBins());
      ParameterUtils.isAbove("Fit length", clusteringSettings.getFitLength(), 1);
      ParameterUtils.isAboveZero("Jump distance", clusteringSettings.getJumpDistance());
      ParameterUtils.isEqualOrAbove("Maximum N", settings.maxN, myMinN);
      if (clusteringSettings.getMle()) {
        ParameterUtils.isAboveZero("Significance level", settings.significanceLevel);
      }
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return true;
  }

  private Plot plotMsd(double[] x, double[] y, double[] sd, String title) {
    if (settings.saveRawData) {
      saveMsd(x, y, sd);
    }

    final Plot plot = new Plot(title, "Time (s)", "Distance (um^2)");
    plot.addPoints(x, y, Plot.LINE);
    // Set limits before any plotting
    double max = 0;
    for (int i = 1; i < x.length; i++) {
      final double value = y[i] + sd[i];
      max = Math.max(max, value);
    }
    plot.setLimits(0, x[x.length - 1] + exposureTime * 0.5, 0, max);
    plot.setColor(Color.blue);
    for (int i = 1; i < x.length; i++) {
      plot.drawLine(x[i], y[i] - sd[i], x[i], y[i] + sd[i]);
    }
    plot.setColor(Color.red);
    display(title, plot);
    return plot;
  }

  /**
   * Fit the MSD using a linear fit that must pass through 0,0.
   *
   * <p>Update the plot by adding the fit line.
   *
   * @param x the x
   * @param y the y
   * @param title the title
   * @param plot the plot
   * @return [D, precision]
   */
  private double[] fitMsd(double[] x, double[] y, String title, Plot plot) {
    // The Weimann paper (Plos One e64287) fits:
    // MSD(n dt) = 4D n dt + 4s^2
    // n = number of jumps
    // dt = time difference between frames
    // s = localisation precision
    // Thus we should fit an intercept as well.

    // From the fit D = gradient / (4*exposureTime)

    double diffCoeff = 0; // D
    double intercept = 0;
    double precision = 0;

    final LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();
    Optimum lvmSolution;
    double ic = Double.NaN;

    // Fit with no intercept
    try {
      final LinearFunction function = new LinearFunction(x, y, clusteringSettings.getFitLength());
      final double[] parameters = new double[] {function.guess()};

      //@formatter:off
      final LeastSquaresProblem problem = new LeastSquaresBuilder()
          .maxEvaluations(Integer.MAX_VALUE)
          .maxIterations(3000)
          .start(parameters)
          .target(function.getY())
          .weight(new DiagonalMatrix(function.getWeights()))
          .model(function, function::jacobian)
          .build();
      //@formatter:on

      lvmSolution = optimizer.optimize(problem);

      final double ss = lvmSolution.getResiduals().dotProduct(lvmSolution.getResiduals());
      // double ss = 0;
      // double[] obs = function.getY();
      // double[] exp = lvmSolution.getValue();
      // for (int i = 0; i < obs.length; i++)
      // ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);

      ic = getAkaikeInformationCriterionFromResiduals(ss, function.getY().length, 1);

      final double gradient = lvmSolution.getPoint().getEntry(0);
      diffCoeff = gradient / 4;

      ImageJUtils.log(
          "Linear fit (%d points) : Gradient = %s, D = %s um^2/s, SS = %s, "
              + "IC = %s (%d evaluations)",
          function.getY().length, MathUtils.rounded(gradient, 4), MathUtils.rounded(diffCoeff, 4),
          MathUtils.rounded(ss), MathUtils.rounded(ic), lvmSolution.getEvaluations());
    } catch (final TooManyIterationsException ex) {
      ImageJUtils.log("Failed to fit : Too many iterations (%s)", ex.getMessage());
    } catch (final ConvergenceException ex) {
      ImageJUtils.log("Failed to fit : %s", ex.getMessage());
    }

    // Fit with intercept.
    // Optionally include the intercept (which is the estimated precision).
    final boolean fitIntercept = true;
    try {
      final LinearFunctionWithIntercept function =
          new LinearFunctionWithIntercept(x, y, clusteringSettings.getFitLength(), fitIntercept);

      //@formatter:off
      final LeastSquaresProblem problem = new LeastSquaresBuilder()
          .maxEvaluations(Integer.MAX_VALUE)
          .maxIterations(3000)
          .start(function.guess())
          .target(function.getY())
          .weight(new DiagonalMatrix(function.getWeights()))
          .model(function, function::jacobian)
          .build();
      //@formatter:on

      lvmSolution = optimizer.optimize(problem);

      final RealVector residuals = lvmSolution.getResiduals();
      final double ss = residuals.dotProduct(residuals);
      // double ss = 0;
      // double[] obs = function.getY();
      // double[] exp = lvmSolution.getValue();
      // for (int i = 0; i < obs.length; i++)
      // ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);

      final double ic2 = getAkaikeInformationCriterionFromResiduals(ss, function.getY().length, 2);
      final double gradient = lvmSolution.getPoint().getEntry(0);
      final double s = lvmSolution.getPoint().getEntry(1);
      final double intercept2 = 4 * s * s;

      if (ic2 < ic || Double.isNaN(ic)) {
        if (settings.debugFitting) {
          // Convert fitted precision in um to nm
          ImageJUtils.log(
              "Linear fit with intercept (%d points) : Gradient = %s, Intercept = %s, "
                  + "D = %s um^2/s, precision = %s nm, SS = %s, IC = %s (%d evaluations)",
              function.getY().length, MathUtils.rounded(gradient, 4),
              MathUtils.rounded(intercept2, 4), MathUtils.rounded(gradient / 4, 4),
              MathUtils.rounded(s * 1000, 4), MathUtils.rounded(ss), MathUtils.rounded(ic2),
              lvmSolution.getEvaluations());
        }

        intercept = intercept2;
        diffCoeff = gradient / 4;
        precision = s;
      }
    } catch (final TooManyIterationsException ex) {
      ImageJUtils.log("Failed to fit with intercept : Too many iterations (%s)", ex.getMessage());
    } catch (final ConvergenceException ex) {
      ImageJUtils.log("Failed to fit with intercept : %s", ex.getMessage());
    }

    if (clusteringSettings.getMsdCorrection()) {
      // Fit with intercept including the MSD correction in the intercept.
      // For the MSD correction we fit including the correction factor (n-1/3)/n:
      // MSD = 4Dt n * (n - 1/3)/n + 4 s^2
      // MSD = 4Dt n - (4Dt) / 3 + 4 s^2
      // i.e. the intercept is allowed to be a small negative.
      try {
        // This function fits the jump distance (n) not the time (nt) so update x
        final double[] x2 = new double[x.length];
        for (int i = 0; i < x2.length; i++) {
          x2[i] = x[i] / exposureTime;
        }

        final LinearFunctionWithMsdCorrectedIntercept function =
            new LinearFunctionWithMsdCorrectedIntercept(x2, y, clusteringSettings.getFitLength(),
                fitIntercept);

        //@formatter:off
        final LeastSquaresProblem problem = new LeastSquaresBuilder()
          .maxEvaluations(Integer.MAX_VALUE)
          .maxIterations(3000)
          .start(function.guess())
          .target(function.getY())
          .weight(new DiagonalMatrix(function.getWeights()))
          .model(function, function::jacobian)
          .build();
        // @formatter:on

        lvmSolution = optimizer.optimize(problem);

        final RealVector residuals = lvmSolution.getResiduals();
        final double ss = residuals.dotProduct(residuals);
        // double ss = 0;
        // double[] obs = function.getY();
        // double[] exp = lvmSolution.getValue();
        // for (int i = 0; i < obs.length; i++)
        // ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);

        final double ic2 =
            getAkaikeInformationCriterionFromResiduals(ss, function.getY().length, 2);
        double gradient = lvmSolution.getPoint().getEntry(0);
        final double s = lvmSolution.getPoint().getEntry(1);
        final double intercept2 = 4 * s * s - gradient / 3;

        // Q. Is this working?
        // Try fixed precision fitting. Is the gradient correct?
        // Revisit all the equations to see if they are wrong.
        // Try adding the x[0] datapoint using the precision.
        // Change the formula to not be linear at x[0] and to just fit the precision, i.e. the
        // intercept2 = 4 * s * s - gradient / 3 is wrong as the
        // equation is not linear below n=1.

        // Incorporate the exposure time into the gradient to allow comparison to other fits
        gradient /= exposureTime;

        if (ic2 < ic || Double.isNaN(ic)) {
          if (settings.debugFitting) {
            // Convert fitted precision in um to nm
            ImageJUtils.log(
                "Linear fit with MSD corrected intercept (%d points) : Gradient = %s, "
                    + "Intercept = %s, D = %s um^2/s, precision = %s nm, SS = %s, "
                    + "IC = %s (%d evaluations)",
                function.getY().length, MathUtils.rounded(gradient, 4),
                MathUtils.rounded(intercept2, 4), MathUtils.rounded(gradient / 4, 4),
                MathUtils.rounded(s * 1000, 4), MathUtils.rounded(ss), MathUtils.rounded(ic2),
                lvmSolution.getEvaluations());
          }

          intercept = intercept2;
          diffCoeff = gradient / 4;
          precision = s;
        }
      } catch (final TooManyIterationsException ex) {
        ImageJUtils.log("Failed to fit with intercept : Too many iterations (%s)", ex.getMessage());
      } catch (final ConvergenceException ex) {
        ImageJUtils.log("Failed to fit with intercept : %s", ex.getMessage());
      }
    }

    // Add the fit to the plot
    if (diffCoeff > 0) {
      plot.setColor(Color.magenta);
      plot.drawLine(0, intercept, x[x.length - 1], 4 * diffCoeff * x[x.length - 1] + intercept);
      display(title, plot);

      checkTraceSettings(diffCoeff);
    }

    return new double[] {diffCoeff, precision};
  }

  /**
   * Get the Akaike Information Criterion (AIC) for a least squares estimate. This assumes that the
   * residuals are distributed according to independent identical normal distributions (with zero
   * mean).
   *
   * @param sumOfSquaredResiduals the sum of squared residuals from the least-squares fit
   * @param numberOfPoints The number of data points
   * @param numberOfParameters The number of fitted parameters
   * @return The corrected Akaike Information Criterion
   * @see <a
   *      href="https://en.wikipedia.org/wiki/Akaike_information_criterion#Comparison_with_least_squares">https://
   *      en.wikipedia.org/wiki/Akaike_information_criterion#Comparison_with_least_squares</a>
   */
  private static double getAkaikeInformationCriterionFromResiduals(double sumOfSquaredResiduals,
      int numberOfPoints, int numberOfParameters) {
    return MathUtils.getAkaikeInformationCriterion(
        MathUtils.getLogLikelihood(sumOfSquaredResiduals, numberOfPoints), numberOfParameters);
  }

  /**
   * Check the tracing settings. For nearest neighbour tracing this checks the distance used for
   * tracing covers enough of the cumulative mean-squared distance distribution.
   *
   * @param d the diffusion coefficient (um^2/s)
   */
  private void checkTraceSettings(double d) {
    if (clusteringSettings.getTraceDiffusionMode() == 1) {
      final double dMax = clusteringSettings.getDiffusionCoefficentMaximum();
      ImageJUtils.log("Checking trace settings: Dmax = %s um^s/s, D = %s um^2/s",
          MathUtils.rounded(dMax), MathUtils.rounded(d));
      if (d > dMax * 2) {
        ImageJUtils.log("WARNING *** The diffusion coefficient may not be large enough! ***");
      }
    } else {
      final double t = exposureTime;
      // Cumul P(r^2) = 1 - exp(-r^2 / 4dt)
      final double r = clusteringSettings.getDistanceThreshold() / 1000;
      final double msd = 4 * d * t;
      final double p = 1 - StdMath.exp(-r * r / msd);
      ImageJUtils.log("Checking trace distance: r = %s nm, D = %s um^2/s, Cumul p(r^2|frame) = %s",
          clusteringSettings.getDistanceThreshold(), MathUtils.rounded(d), MathUtils.rounded(p));
      if (p < 0.95) {
        ImageJUtils.log("WARNING *** The tracing distance may not be large enough! ***");
      }
    }
  }

  private static class LinearFunction implements MultivariateVectorFunction {
    double[] x;
    double[] y;
    double[][] jacobian;

    public LinearFunction(double[] x, double[] y, int length) {
      final int to = Math.min(x.length, 1 + length);
      this.x = Arrays.copyOfRange(x, 1, to);
      this.y = Arrays.copyOfRange(y, 1, to);
      jacobian = calculateJacobian();
    }

    // Adapted from http://commons.apache.org/proper/commons-math/userguide/optimization.html

    /**
     * Guess the gradient.
     *
     * @return An estimate for the linear gradient.
     */
    private double guess() {
      return y[y.length - 1] / x[x.length - 1];
    }

    private double[] getWeights() {
      final double[] w = new double[x.length];
      Arrays.fill(w, 1);
      return w;
    }

    private double[] getY() {
      return y;
    }

    @Override
    public double[] value(double[] variables) {
      final double[] values = new double[x.length];
      for (int i = 0; i < values.length; i++) {
        values[i] = x[i] * variables[0];
      }
      return values;
    }

    /**
     * Get the Jacobian.
     *
     * @param variables The variables (ignored)
     * @return the Jacobian
     */
    double[][] jacobian(double[] variables) {
      return jacobian;
    }

    double[][] calculateJacobian() {
      // Compute the gradients using calculus differentiation:
      // y = ax + c
      // dyDa = x
      final double[][] jacobian = new double[x.length][1];

      for (int i = 0; i < jacobian.length; ++i) {
        jacobian[i][0] = x[i];
      }

      return jacobian;
    }
  }

  private static class LinearFunctionWithIntercept implements MultivariateVectorFunction {
    final double[] x;
    final double[] y;
    final boolean fitIntercept;

    public LinearFunctionWithIntercept(double[] x, double[] y, int length, boolean fitIntercept) {
      this.fitIntercept = fitIntercept;
      final int to = Math.min(x.length, 1 + length);
      // Optionally include the intercept
      final int from = (fitIntercept) ? 0 : 1;
      this.x = Arrays.copyOfRange(x, from, to);
      this.y = Arrays.copyOfRange(y, from, to);
    }

    /**
     * Guess the gradient and intercept.
     *
     * @return An estimate for the linear gradient and intercept.
     */
    public double[] guess() {
      final int n1 = (fitIntercept) ? 1 : 0;

      if (y.length == n1 + 1) {
        return new double[] {y[n1] / x[n1], 0};
      }

      final double a = (y[y.length - 1] - y[n1]) / (x[x.length - 1] - x[n1]);
      // y = ax + 4c^2
      // y = ax + intercept
      // intercept = y - ax
      // = 4c^2
      final double intercept = y[y.length - 1] - a * x[x.length - 1];
      final double c = (intercept < 0) ? 0 : Math.sqrt(intercept / 4);

      return new double[] {a, c};
    }

    public double[] getWeights() {
      final double[] w = new double[x.length];
      Arrays.fill(w, 1);
      return w;
    }

    public double[] getY() {
      return y;
    }

    @Override
    public double[] value(double[] variables) {
      // y = ax + 4c^2
      final double[] values = new double[x.length];
      final double a = variables[0];
      final double intercept = 4 * variables[1] * variables[1];
      for (int i = 0; i < values.length; i++) {
        values[i] = a * x[i] + intercept;
      }
      return values;
    }

    double[][] jacobian(double[] variables) {
      // Compute the gradients using calculus differentiation:
      // y = ax + 4c^2
      // dyDa = x
      // dy_dc = 8c
      final double[][] jacobian = new double[x.length][2];
      final double dy_dc = 8 * variables[1];

      for (int i = 0; i < jacobian.length; ++i) {
        jacobian[i][0] = x[i];
        jacobian[i][1] = dy_dc;
      }

      return jacobian;
    }
  }

  private static class LinearFunctionWithMsdCorrectedIntercept
      implements MultivariateVectorFunction {
    static final double THIRD = 1 / 3.0;
    final double[] x;
    final double[] y;
    final boolean fitIntercept;

    public LinearFunctionWithMsdCorrectedIntercept(double[] x, double[] y, int length,
        boolean fitIntercept) {
      this.fitIntercept = fitIntercept;
      final int to = Math.min(x.length, 1 + length);
      // Optionally include the intercept
      final int from = (fitIntercept) ? 0 : 1;
      this.x = Arrays.copyOfRange(x, from, to);
      this.y = Arrays.copyOfRange(y, from, to);
    }

    /**
     * Guess the gradient and intercept.
     *
     * @return An estimate for the linear gradient and intercept.
     */
    public double[] guess() {
      final int n1 = (fitIntercept) ? 1 : 0;

      if (y.length == n1 + 1) {
        return new double[] {y[n1] / x[n1], 0};
      }

      final double a = (y[y.length - 1] - y[n1]) / (x[x.length - 1] - x[n1]);
      // y = ax - a/3 + 4c^2
      // y = ax + intercept
      // intercept = y - ax
      // = 4c^2 - a/3
      // 4c^2 = intercept + a/3
      double intercept = y[y.length - 1] - a * x[x.length - 1];
      intercept += a * THIRD;
      final double c = (intercept < 0) ? 0 : Math.sqrt(intercept / 4);

      return new double[] {a, c};
    }

    public double[] getWeights() {
      final double[] w = new double[x.length];
      Arrays.fill(w, 1);
      return w;
    }

    public double[] getY() {
      return y;
    }

    @Override
    public double[] value(double[] variables) {
      // When x>=1:
      // y = ax - a/3 + 4c^2
      // When x==0:
      // y = 4c^2
      final double[] values = new double[x.length];
      final double a = variables[0];
      final double intercept = 4 * variables[1] * variables[1];
      final double error = intercept - a * THIRD;
      int index = 0;
      // Special case for fitting the intercept since the line is not linear below n=1
      if (fitIntercept) {
        values[index++] = intercept;
      }
      for (; index < values.length; index++) {
        values[index] = a * x[index] + error;
      }
      return values;
    }

    double[][] jacobian(double[] variables) {
      // Compute the gradients using calculus differentiation:
      // y = ax - a/3 + 4c^2
      // dyDa = x - 1/3
      // dy_dc = 8c
      final double[][] jacobian = new double[x.length][2];
      final double dy_dc = 8 * variables[1];

      for (int i = 0; i < jacobian.length; ++i) {
        jacobian[i][0] = x[i] - THIRD;
        jacobian[i][1] = dy_dc;
      }

      return jacobian;
    }
  }

  /**
   * Fit the jump distance histogram.
   *
   * <p>Update the plot by adding the fit line(s).
   *
   * @param jumpDistances (in um^2)
   * @param jdHistogram the jump distance histogram
   * @return The fitted coefficients and fractions
   */
  private double[][] fitJumpDistance(StoredDataStatistics jumpDistances, double[][] jdHistogram) {
    final double msd = jumpDistances.getMean();
    final double meanDistance = Math.sqrt(msd) * 1e3;
    // TODO:
    // Q. Should the beta be expressed using the mean-distance or MSD?
    // Q. Should it be normalised to the frame length. If not then the beta will be invariant on
    // jump distance length
    beta = meanDistance / precision;
    // Set the minimum diffusion coefficient using the precision:
    // Note: 4D = MSD
    // 4D = precision^2
    // D = precision^2 / 4
    // D = (precision/2)^2
    final double minD = MathUtils.pow2(precision / 2000.0); // Extra 1000 factor to convert nm to um
    ImageJUtils.log(
        "Jump Distance analysis : N = %d, Time = %d frames (%s seconds). "
            + "MSD = %s um^2/jump, Mean Distance = %s nm/jump, Precision = %s nm, "
            + "Beta = %s, minD = %s um^2/jump",
        jumpDistances.getN(), clusteringSettings.getJumpDistance(),
        MathUtils.rounded(clusteringSettings.getJumpDistance() * exposureTime, 4),
        MathUtils.rounded(msd, 4), MathUtils.rounded(meanDistance, 4),
        MathUtils.rounded(precision, 4), MathUtils.rounded(beta, 4), MathUtils.rounded(minD, 4));

    final Logger logger = ImageJPluginLoggerHelper.getLogger(getClass());
    if (settings.debugFitting) {
      logger.setLevel(Level.FINE);
    }
    final JumpDistanceAnalysis jd = new JumpDistanceAnalysis(logger);
    jd.setFitRestarts(clusteringSettings.getFitRestarts());
    jd.setMinFraction(settings.minFraction);
    jd.setMinDifference(settings.minDifference);
    jd.setMinN(myMinN);
    jd.setMaxN(settings.maxN);
    jd.setMinD(minD);
    jd.setSignificanceLevel(settings.significanceLevel);
    // Update the plot with the fit
    jd.setCurveLogger(this);

    // Set the calibration
    jd.setN(clusteringSettings.getJumpDistance());
    jd.setDeltaT(exposureTime);
    if (clusteringSettings.getPrecisionCorrection()) {
      jd.setError(precision, true);
    }
    jd.setMsdCorrection(clusteringSettings.getMsdCorrection());

    double[][] fit;
    if (clusteringSettings.getMle()) {
      fit = jd.fitJumpDistancesMle(jumpDistances.getValues(), jdHistogram);
    } else {
      fit = jd.fitJumpDistanceHistogram(jumpDistances.getMean(), jdHistogram);
    }

    // Get the raw fitted D and convert it to a calibrated D*
    if (fit != null) {
      fit[0] = jd.calculateApparentDiffusionCoefficient(fit[0]);
      // Check the largest D
      checkTraceSettings(fit[0][0]);
      fitValue = jd.getLastFitValue();
    }

    return fit;
  }

  @Override
  public int getNumberOfCurvePoints() {
    return 300;
  }

  @Override
  public void saveSinglePopulationCurve(double[][] curve) {
    final double[] x = curve[0];
    final double[] y = curve[1];
    if (settings.saveRawData) {
      saveFit(x, y, "Single");
    }
    addToJumpDistancePlot(x, y, Color.magenta);
  }

  @Override
  public void saveMixedPopulationCurve(double[][] curve) {
    final double[] x = curve[0];
    final double[] y = curve[1];
    if (settings.saveRawData) {
      saveFit(x, y, "Mixed");
    }
    addToJumpDistancePlot(x, y, Color.yellow);
  }

  private void addToJumpDistancePlot(double[] x, double[] y, Color color) {
    jdPlot.setColor(color);
    jdPlot.addPoints(x, y, Plot.LINE);
    display(jdTitle, jdPlot);
  }

  /**
   * Macro extension function.
   *
   * <p>Get the number of fitted species from the last call to fit the jump distances.
   *
   * @param args An array of output variables.
   *
   *        <ul>
   *
   *        <li>[0]: Double[1] - output: the number of species
   *
   *        </ul>
   * @return Empty string
   */
  public static String getNumberOfSpecies(Object[] args) {
    int species = 0;
    final double[][] jumpDistanceParameters = jumpDistanceParametersRef.get();
    if (jumpDistanceParameters != null) {
      species = jumpDistanceParameters[0].length;
    }
    final Double[] array = (Double[]) args[0];
    array[0] = Double.valueOf(species);
    return "";
  }

  /**
   * Macro extension function.
   *
   * <p>Get the diffusion coefficient for the requested species from the last call to fit the jump
   * distances.
   *
   * @param args An array of input and output variables.
   *
   *        <ul>
   *
   *        <li>[0]: Double[1] - input: the index of the species
   *
   *        <li>[1]: Double[1] - output: the coefficient
   *
   *        </ul>
   * @return Empty string
   */
  public static String getD(Object[] args) {
    double value = 0;
    final double[][] jumpDistanceParameters = jumpDistanceParametersRef.get();
    if (jumpDistanceParameters != null) {
      final int index = ((Double) args[0]).intValue();
      if (index >= 0 && index < jumpDistanceParameters[0].length) {
        value = jumpDistanceParameters[0][index];
      }
    }
    ((Double[]) args[1])[0] = Double.valueOf(value);
    return "";
  }

  /**
   * Macro extension function.
   *
   * <p>Get the population fraction for the requested species from the last call to fit the jump
   * distances.
   *
   * @param args An array of input and output variables.
   *
   *        <ul>
   *
   *        <li>[0]: Double[1] - input: the index of the species
   *
   *        <li>[1]: Double[1] - output: the population fraction
   *
   *        </ul>
   * @return Empty string
   */
  public static String getF(Object[] args) {
    double value = 0;
    final double[][] jumpDistanceParameters = jumpDistanceParametersRef.get();
    if (jumpDistanceParameters != null) {
      final int index = ((Double) args[0]).intValue();
      if (index >= 0 && index < jumpDistanceParameters[1].length) {
        value = jumpDistanceParameters[1][index];
      }
    }
    ((Double[]) args[1])[0] = Double.valueOf(value);
    return "";
  }

  /**
   * Macro extension function.
   *
   * <p>Get the diffusion coefficient and population fraction for the requested species from the
   * last call to fit the jump distances.
   *
   * @param args An array of input and output variables.
   *
   *        <ul>
   *
   *        <li>[0]: Double[1] - input: the index of the species
   *
   *        <li>[1]: Double[1] - output: the coefficient
   *
   *        <li>[2]: Double[1] - output: the population fraction
   *
   *        </ul>
   * @return Empty string
   */
  public static String getSpecies(Object[] args) {
    double value = 0;
    double value2 = 0;
    final double[][] jumpDistanceParameters = jumpDistanceParametersRef.get();
    if (jumpDistanceParameters != null) {
      final int i = ((Double) args[0]).intValue();
      if (i >= 0 && i < jumpDistanceParameters[0].length) {
        value = jumpDistanceParameters[0][i];
        value2 = jumpDistanceParameters[1][i];
      }
    }
    ((Double[]) args[1])[0] = Double.valueOf(value);
    ((Double[]) args[2])[0] = Double.valueOf(value2);
    return "";
  }

  private boolean showMultiDialog(ArrayList<MemoryPeakResults> allResults) {
    multiMode = true;

    // Show a list box containing all the results. This should remember the last set of chosen
    // items.
    final MultiDialog md = ResultsManager.createMultiDialog(TITLE);
    md.setSelected(selectedRef.get());
    md.setHelpUrl(HelpUrls.getUrl("trace-diffusion-multi"));

    md.showDialog();

    if (md.wasCancelled()) {
      return false;
    }

    final List<String> selected = md.getSelectedResults();
    if (selected.isEmpty()) {
      IJ.error(TITLE, "No results were selected");
      return false;
    }
    selectedRef.set(selected);

    for (final String name : selected) {
      final MemoryPeakResults r = MemoryPeakResults.getResults(name);
      if (r != null) {
        allResults.add(r);
      }
    }

    if (allResults.isEmpty()) {
      return false;
    }

    // Check calibration exists for the first set of results
    if (!checkCalibration(allResults.get(0))) {
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
}
