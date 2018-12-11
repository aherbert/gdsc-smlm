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

import uk.ac.sussex.gdsc.core.data.utils.Converter;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.HistogramPlotBuilder;
import uk.ac.sussex.gdsc.core.ij.ImageJPluginLoggerHelper;
import uk.ac.sussex.gdsc.core.ij.ImageJTrackProgress;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.Plot2;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.StoredDataStatistics;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationHelper;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.ConfigurationException;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.ClusteringSettings;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.fitting.JumpDistanceAnalysis;
import uk.ac.sussex.gdsc.smlm.fitting.JumpDistanceAnalysis.CurveLogger;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultCalculator;
import uk.ac.sussex.gdsc.smlm.results.Gaussian2DPeakResultHelper;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResultStoreList;
import uk.ac.sussex.gdsc.smlm.results.Trace;
import uk.ac.sussex.gdsc.smlm.results.TraceManager;

import ij.IJ;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.text.TextWindow;

import org.apache.commons.math3.analysis.MultivariateMatrixFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyIterationsException;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer.Optimum;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.util.FastMath;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Run a tracing algorithm on the peak results to trace molecules across the frames.
 */
public class TraceDiffusion implements PlugIn, CurveLogger {
  private static final String TITLE = "Trace Diffusion";
  private static String inputOption = "";
  private static String header = null;
  private static TextWindow summaryTable = null;

  private static final String[] NAMES = new String[] {"Total Signal", "Signal/Frame", "t-On (s)"};
  private static final boolean[] ROUNDED = new boolean[] {false, false, true};
  private static boolean[] displayHistograms = new boolean[NAMES.length];

  static {
    for (int i = 0; i < displayHistograms.length; i++) {
      displayHistograms[i] = false;
    }
  }
  private static final int TOTAL_SIGNAL = 0;
  private static final int SIGNAL_PER_FRAME = 1;
  private static final int T_ON = 2;

  private static boolean[] alwaysRemoveOutliers;

  static {
    alwaysRemoveOutliers = new boolean[NAMES.length];
    alwaysRemoveOutliers[TOTAL_SIGNAL] = false;
  }

  private static boolean displayMSDHistogram = true;
  private static boolean displayDHistogram = true;
  private static boolean displayTraceLength = false;
  private static boolean displayTraceSize = false;

  private static boolean saveTraceDistances = false;
  private static boolean saveRawData = false;
  private static String rawDataDirectory = "";
  private boolean directoryChosen = false;
  private static String distancesFilename = "";
  private static double significanceLevel = 0.05;
  private static double minFraction = 0.1;
  private static double minDifference = 2;
  private static int minN = 1;
  private static int maxN = 5;
  private static boolean debugFitting = false;
  private static String tracesFilename = "";
  private static String title = "";

  private ClusteringSettings.Builder settings;
  private MemoryPeakResults results;
  private boolean extraOptions;
  private boolean multiMode;
  private int myMinN = 1;

  // The number of additional datasets
  private int additionalDatasets = 0;

  // Store exposure time in seconds
  private double exposureTime = 0;
  private double precision;
  private double beta;
  private double fitValue = Double.NaN;

  // Used to tile new plot windows
  private final WindowOrganiser windowOrganiser = new WindowOrganiser();

  private String jdTitle = TITLE + " Jump Distance";
  private Plot2 jdPlot;

  // Used for the macro extensions
  private static double[][] jumpDistanceParameters = null;

  // Used for the multiMode option
  private static ArrayList<String> selected;

  /** {@inheritDoc} */
  @Override
  public void run(String arg) {
    SMLMUsageTracker.recordPlugin(this.getClass(), arg);

    jumpDistanceParameters = null;

    extraOptions = ImageJUtils.isExtraOptions();
    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "No localisations in memory");
      return;
    }

    final ArrayList<MemoryPeakResults> allResults = new ArrayList<>();

    // Option to pick multiple input datasets together using a list box.
    if ("multi".equals(arg)) {
      if (!showMultiDialog(allResults)) {
        return;
      }
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
    if (traces.length > 0) {
      if (!showDialog()) {
        return;
      }
    }

    final int count = traces.length;
    double[] fitMSDResult = null;
    int n = 0;
    double[][] jdParams = null;
    if (count > 0) {
      calculatePrecision(traces, allResults.size() > 1);

      // --- MSD Analysis ---

      // Conversion constants
      final double px2ToUm2 = MathUtils.pow2(results.getCalibrationReader().getNmPerPixel()) / 1e6;
      final double px2ToUm2PerSecond = px2ToUm2 / exposureTime;

      // Get the maximum trace length
      int length = settings.getMinimumTraceLength();
      if (!settings.getTruncate()) {
        for (final Trace trace : traces) {
          if (length < trace.size()) {
            length = trace.size();
          }
        }
      }

      // Get the localisation error (4s^2) in um^2
      final double error =
          (settings.getPrecisionCorrection()) ? 4 * precision * precision / 1e6 : 0;
      // Pre-calculate MSD correction factors. This accounts for the fact that the distance moved
      // in the start/end frames is reduced due to the averaging of the particle location over the
      // entire frame into a single point. The true MSD may be restored by applying a factor.
      // Note: These are used for the calculation of the diffusion coefficients per molecule and
      // the MSD passed to the Jump Distance analysis. However the error is not included in the
      // jump distance analysis so will be subtracted from the fitted D coefficients later.
      final double[] factors;
      if (settings.getMsdCorrection()) {
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
          (saveTraceDistances || displayTraceLength) ? new ArrayList<>(traces.length) : null;

      // Store all the jump distances at the specified interval
      final StoredDataStatistics jumpDistances = new StoredDataStatistics();
      final int jumpDistanceInterval = settings.getJumpDistance();

      // Compute squared distances
      final StoredDataStatistics msdPerMoleculeAllVsAll = new StoredDataStatistics();
      final StoredDataStatistics msdPerMoleculeAdjacent = new StoredDataStatistics();
      for (final Trace trace : traces) {
        final PeakResultStoreList results = trace.getPoints();
        // Sum the MSD and the time
        final int traceLength =
            (settings.getTruncate()) ? settings.getMinimumTraceLength() : trace.size();

        // Get the mean for each time separation
        final double[] sumDistance = new double[traceLength + 1];
        final double[] sumTime = new double[sumDistance.length];

        // Do the distances to the origin (saving if necessary)
        {
          final float x = results.get(0).getXPosition();
          final float y = results.get(0).getYPosition();
          if (distances != null) {
            final double[] msd = new double[traceLength - 1];
            for (int j = 1; j < traceLength; j++) {
              final int t = j;
              final double d = distance2(x, y, results.get(j));
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
              final double d = distance2(x, y, results.get(j));
              if (t == jumpDistanceInterval) {
                jumpDistances.add(px2ToUm2 * d);
              }
              sumDistance[t] += d;
              sumTime[t] += t;
            }
          }
        }

        if (settings.getInternalDistances()) {
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

      StoredDataStatistics dPerMoleculeAllVsAll = null;
      StoredDataStatistics dPerMoleculeAdjacent = null;
      if (saveTraceDistances || (settings.getShowHistograms() && displayDHistogram)) {
        dPerMoleculeAllVsAll = calculateDiffusionCoefficient(msdPerMoleculeAllVsAll);
        dPerMoleculeAdjacent = calculateDiffusionCoefficient(msdPerMoleculeAdjacent);
      }

      if (saveTraceDistances) {
        saveTraceDistances(traces.length, distances, msdPerMoleculeAllVsAll, msdPerMoleculeAdjacent,
            dPerMoleculeAllVsAll, dPerMoleculeAdjacent);
      }

      if (displayTraceLength) {
        final StoredDataStatistics lengths = calculateTraceLengths(distances);
        showHistogram(lengths, "Trace length (um)");
      }

      if (displayTraceSize) {
        final StoredDataStatistics sizes = calculateTraceSizes(traces);
        showHistogram(sizes, "Trace size", true);
      }

      // Plot the per-trace histogram of MSD and D
      if (settings.getShowHistograms()) {
        if (displayMSDHistogram) {
          showHistogram(msdPerMoleculeAllVsAll, "MSD/Molecule (all-vs-all)");
          showHistogram(msdPerMoleculeAdjacent, "MSD/Molecule (adjacent)");
        }
        if (displayDHistogram) {
          showHistogram(dPerMoleculeAllVsAll, "D/Molecule (all-vs-all)");
          showHistogram(dPerMoleculeAdjacent, "D/Molecule (adjacent)");
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
      final Plot2 plot = plotMSD(x, y, sd, title);

      // Fit the MSD using a linear fit
      fitMSDResult = fitMSD(x, y, title, plot);

      // Jump Distance analysis
      if (saveRawData) {
        saveStatistics(jumpDistances, "Jump Distance", "Distance (um^2)", false);
      }

      // Calculate the cumulative jump-distance histogram
      final double[][] jdHistogram =
          JumpDistanceAnalysis.cumulativeHistogram(jumpDistances.getValues());

      // Always show the jump distance histogram
      jdTitle = TITLE + " Jump Distance";
      jdPlot = new Plot2(jdTitle, "Distance (um^2)", "Cumulative Probability", jdHistogram[0],
          jdHistogram[1]);
      display(jdTitle, jdPlot);

      // Fit Jump Distance cumulative probability
      n = jumpDistances.getN();
      jumpDistanceParameters = jdParams = fitJumpDistance(jumpDistances, jdHistogram);
    }

    summarise(traces, fitMSDResult, n, jdParams);
  }

  /**
   * Calculate trace lengths.
   *
   * @param distances the distances for each trace
   * @return the trace lengths
   */
  public static StoredDataStatistics calculateTraceLengths(ArrayList<double[]> distances) {
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

  private void display(String title, Plot2 plot) {
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
    if (saveRawData) {
      saveStatistics(stats, title, title, rounded);
    }

    new HistogramPlotBuilder(TITLE, stats, title).setIntegerBins(integerData)
        .setRemoveOutliersOption((settings.getRemoveOutliers() || alwaysRemoveOutliers) ? 1 : 0)
        .setNumberOfBins(settings.getHistogramBins()).show(windowOrganiser);
  }

  /**
   * Calculate the average precision of localisation in the traces.
   *
   * @param traces the traces
   * @param multi the multi
   */
  private void calculatePrecision(Trace[] traces, boolean multi) {
    // Check the diffusion simulation for a precision
    if (DiffusionRateTest.isSimulated(results.getName()) && !multi) {
      precision = DiffusionRateTest.lastSimulatedPrecision;
    } else {
      precision = 999;
      try {
        final Gaussian2DPeakResultCalculator calculator = Gaussian2DPeakResultHelper.create(
            results.getPSF(), results.getCalibration(), Gaussian2DPeakResultHelper.LSE_PRECISION);
        // Get the average precision of the localisations
        precision = 0;
        int n = 0;
        for (final Trace trace : traces) {
          for (int k = 0; k < trace.size(); k++) {
            final PeakResult r = trace.get(k);
            precision += calculator.getLSEPrecision(r.getParameters(), r.getNoise());
          }
          n += trace.size();
        }
        precision /= n;
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
  private static StoredDataStatistics calculateDiffusionCoefficient(
      StoredDataStatistics msdPerMoleculeAdjacent) {
    final StoredDataStatistics dPerMolecule = new StoredDataStatistics();
    final double diffusionCoefficientConversion = 1.0 / 4.0;
    for (final double msd : msdPerMoleculeAdjacent.getValues()) {
      dPerMolecule.add(msd * diffusionCoefficientConversion);
    }
    return dPerMolecule;
  }

  private void saveTraceDistances(int nTraces, ArrayList<double[]> distances,
      StoredDataStatistics msdPerMolecule, StoredDataStatistics msdPerMoleculeAdjacent,
      StoredDataStatistics dStarPerMolecule, StoredDataStatistics dStarPerMoleculeAdjacent) {
    distancesFilename = ImageJUtils.getFilename("Trace_Distances_File", distancesFilename);
    if (distancesFilename != null) {
      distancesFilename = ImageJUtils.replaceExtension(distancesFilename, "xls");

      try (BufferedWriter out = new BufferedWriter(
          new OutputStreamWriter(new FileOutputStream(distancesFilename), "UTF-8"))) {
        final double[] msd = msdPerMolecule.getValues();
        final double[] msd2 = msdPerMoleculeAdjacent.getValues();
        final double[] dStar = dStarPerMolecule.getValues();
        final double[] dStar2 = dStarPerMoleculeAdjacent.getValues();
        out.write(String.format("#%d traces : Precision = %s nm : Exposure time = %s s", nTraces,
            MathUtils.rounded(precision, 4), MathUtils.rounded(exposureTime, 4)));
        out.newLine();
        out.write(String.format(
            "#TraceId\tMSD all-vs-all (um^2/s)\tMSD adjacent (um^2/s)\tD all-vs-all(um^2/s)\tD adjacent(um^2/s)\tDistances (um^2) per %ss ... ",
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
      } catch (final Exception ex) {
        // Ignore
      }
    }
  }

  private void saveMSD(double[] x, double[] y, double[] se) {
    if (!directoryChosen) {
      rawDataDirectory = ImageJUtils.getDirectory("Data_directory", rawDataDirectory);
    }
    directoryChosen = true;
    if (rawDataDirectory == null) {
      return;
    }
    final String filename = rawDataDirectory + "MSD.txt";

    try (BufferedWriter out =
        new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename), "UTF-8"))) {
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
    } catch (final Exception ex) {
      // Ignore
    }
  }

  private void saveStatistics(StoredDataStatistics stats, String title, String label,
      boolean rounded) {
    if (!directoryChosen) {
      rawDataDirectory = ImageJUtils.getDirectory("Data_directory", rawDataDirectory);
    }
    directoryChosen = true;
    if (rawDataDirectory == null) {
      return;
    }
    final String filename =
        rawDataDirectory + title.replace("/", " per ").replace("*", "star") + ".txt";

    try (BufferedWriter out =
        new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename), "UTF-8"))) {
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
    } catch (final Exception ex) {
      // Ignore
    }
  }

  private void saveFit(double[] x, double[] y, String title) {
    if (!directoryChosen) {
      rawDataDirectory = ImageJUtils.getDirectory("Data_directory", rawDataDirectory);
    }
    directoryChosen = true;
    if (rawDataDirectory == null) {
      return;
    }
    final String filename = rawDataDirectory + "Fit." + title + ".txt";

    try (BufferedWriter out =
        new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename), "UTF-8"))) {
      out.write("JumpDistance\tCumulativeP");
      out.newLine();
      for (int i = 0; i < x.length; i++) {
        out.write(Double.toString(x[i]));
        out.write('\t');
        out.write(Double.toString(y[i]));
        out.newLine();
      }
    } catch (final Exception ex) {
      // Ignore
    }
  }

  private static double distance2(final float x, final float y, PeakResult r2) {
    final double dx = x - r2.getXPosition();
    final double dy = y - r2.getYPosition();
    return dx * dx + dy * dy;
  }

  /**
   * Filter traces that are not the minimum length.
   *
   * @param name the name
   * @param traces the traces
   * @param minimumTraceLength the minimum trace length
   * @param ignoreEnds the ignore ends
   * @return The new traces
   */
  private static Trace[] filterTraces(String name, Trace[] traces, int minimumTraceLength,
      boolean ignoreEnds) {
    final int minLength = (ignoreEnds) ? minimumTraceLength + 2 : minimumTraceLength;
    int count = 0;
    for (int i = 0; i < traces.length; i++) {
      if (traces[i].size() >= minLength) {
        if (ignoreEnds) {
          traces[i].removeEnds();
        }
        traces[count++] = traces[i];
      }
    }

    ImageJUtils.log(
        "Filtered results '%s' : %s filtered to %d using minimum length %d (Ignore ends = %b)",
        name, TextUtils.pleural(traces.length, "trace"), count, minimumTraceLength, ignoreEnds);
    return Arrays.copyOf(traces, count);
  }

  private String createSettingsComment() {
    return String.format("Molecule tracing : distance-threshold = %f nm",
        settings.getDistanceThreshold());
  }

  private void summarise(Trace[] traces, double[] fitMSDResult, int n, double[][] jdParams) {
    IJ.showStatus("Calculating summary ...");

    // Create summary table
    createSummaryTable();

    final Statistics[] stats = new Statistics[NAMES.length];
    for (int i = 0; i < stats.length; i++) {
      stats[i] = (settings.getShowHistograms()) ? new StoredDataStatistics() : new Statistics();
    }
    for (final Trace trace : traces) {
      stats[T_ON].add(trace.getOnTime() * exposureTime);
      final double signal = trace.getSignal() / results.getGain();
      stats[TOTAL_SIGNAL].add(signal);
      stats[SIGNAL_PER_FRAME].add(signal / trace.size());
    }

    // Add to the summary table
    final StringBuilder sb = new StringBuilder(title);
    sb.append('\t').append(createCombinedName());
    sb.append('\t');
    sb.append(MathUtils.rounded(exposureTime * 1000, 3)).append('\t');
    sb.append(MathUtils.rounded(settings.getDistanceThreshold(), 3)).append('\t');
    sb.append(MathUtils.rounded(settings.getDistanceExclusion(), 3)).append('\t');
    sb.append(settings.getMinimumTraceLength()).append('\t');
    sb.append(settings.getIgnoreEnds()).append('\t');
    sb.append(settings.getTruncate()).append('\t');
    sb.append(settings.getInternalDistances()).append('\t');
    sb.append(settings.getFitLength()).append('\t');
    sb.append(settings.getMsdCorrection()).append('\t');
    sb.append(settings.getPrecisionCorrection()).append('\t');
    sb.append(settings.getMle()).append('\t');
    sb.append(traces.length).append('\t');
    sb.append(MathUtils.rounded(precision, 4)).append('\t');
    double D = 0;
    double s = 0;
    if (fitMSDResult != null) {
      D = fitMSDResult[0];
      s = fitMSDResult[1];
    }
    sb.append(MathUtils.rounded(D, 4)).append('\t');
    sb.append(MathUtils.rounded(s * 1000, 4)).append('\t');
    sb.append(MathUtils.rounded(settings.getJumpDistance() * exposureTime)).append('\t');
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
    if (java.awt.GraphicsEnvironment.isHeadless()) {
      IJ.log(sb.toString());
      return;
    }
    summaryTable.append(sb.toString());

    if (settings.getShowHistograms()) {
      IJ.showStatus("Calculating histograms ...");

      for (int i = 0; i < NAMES.length; i++) {
        if (displayHistograms[i]) {
          showHistogram((StoredDataStatistics) stats[i], NAMES[i], alwaysRemoveOutliers[i],
              ROUNDED[i], false);
        }
      }
    }

    windowOrganiser.tile();

    IJ.showStatus("Finished " + TITLE);
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

  private static void createSummaryTable() {
    if (java.awt.GraphicsEnvironment.isHeadless()) {
      if (header == null) {
        header = createHeader();
        IJ.log(header);
      }
    } else if (summaryTable == null || !summaryTable.isVisible()) {
      summaryTable = new TextWindow(TITLE + " Data Summary", createHeader(), "", 800, 300);
      summaryTable.setVisible(true);
    }
  }

  private static String createHeader() {
    final StringBuilder sb =
        new StringBuilder("Title\tDataset\tExposure time (ms)\tD-threshold (nm)");
    sb.append("\tEx-threshold (nm)\t");
    sb.append("Min.Length\tIgnoreEnds\tTruncate\tInternal\tFit Length");
    sb.append("\tMSD corr.\ts corr.\tMLE\tTraces\ts (nm)\tD (um^2/s)\tfit s (nm)");
    sb.append("\tJump Distance (s)\tN\tBeta\tJump D (um^2/s)\tFractions\tFit Score");
    for (int i = 0; i < NAMES.length; i++) {
      sb.append('\t').append(NAMES[i]);
    }
    return sb.toString();
  }

  private boolean showTraceDialog(ArrayList<MemoryPeakResults> allResults) {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    if (!multiMode) {
      ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);
    }

    settings = SettingsManager.readClusteringSettings(0).toBuilder();

    gd.addNumericField("Distance_Threshold (nm)", settings.getDistanceThreshold(), 0);
    gd.addNumericField("Distance_Exclusion (nm)", settings.getDistanceExclusion(), 0);
    gd.addSlider("Min_trace_length", 2, 20, settings.getMinimumTraceLength());
    gd.addCheckbox("Ignore_ends", settings.getIgnoreEnds());
    gd.addCheckbox("Save_traces", settings.getSaveTraces());

    gd.showDialog();

    if (gd.wasCanceled() || !readTraceDialog(gd)) {
      return false;
    }

    // Update the settings
    SettingsManager.writeSettings(settings.build());

    // Load the results
    if (!multiMode) {
      final MemoryPeakResults results =
          ResultsManager.loadInputResults(inputOption, true, null, null);
      if (results == null || results.size() == 0) {
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
      inputOption = ResultsManager.getInputSource(gd);
    }
    settings.setDistanceThreshold(gd.getNextNumber());
    settings.setDistanceExclusion(Math.abs(gd.getNextNumber()));
    settings.setMinimumTraceLength((int) Math.abs(gd.getNextNumber()));
    settings.setIgnoreEnds(gd.getNextBoolean());
    settings.setSaveTraces(gd.getNextBoolean());

    if (gd.invalidNumber()) {
      return false;
    }

    // Check arguments
    try {
      Parameters.isAboveZero("Distance threshold", settings.getDistanceThreshold());
      Parameters.isAbove("Min trace length", settings.getMinimumTraceLength(), 1);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return true;
  }

  private Trace[] getTraces(ArrayList<MemoryPeakResults> allResults) {
    this.results = allResults.get(0);

    // Results should be checked for calibration by this point
    exposureTime = results.getCalibrationReader().getExposureTime() / 1000;

    // Convert from NM to the native units of the results
    final Converter c =
        CalibrationHelper.getDistanceConverter(results.getCalibration(), DistanceUnit.NM);
    final double distanceThreshold = c.convertBack(settings.getDistanceThreshold());
    final double distanceExclusion = c.convertBack(settings.getDistanceExclusion());

    final ArrayList<Trace> allTraces = new ArrayList<>();
    additionalDatasets = -1;
    for (final MemoryPeakResults r : allResults) {
      additionalDatasets++;

      final TraceManager manager = new TraceManager(r);

      // Run the tracing
      manager.setTracker(new ImageJTrackProgress());
      // convert from
      manager.setDistanceExclusion(distanceExclusion);
      manager.traceMolecules(distanceThreshold, 1);
      Trace[] traces = manager.getTraces();

      traces = filterTraces(r.getName(), traces, settings.getMinimumTraceLength(),
          settings.getIgnoreEnds());
      allTraces.addAll(Arrays.asList(traces));

      // --- Save results ---
      if (traces.length > 0) {
        // Save the traces to memory
        TraceMolecules.saveResults(r, traces, "Tracks");

        if (settings.getSaveTraces()) {
          // Sort traces by time to assist the results source in extracting frames sequentially.
          // Do this before saving to assist in debugging using the saved traces file.
          TraceMolecules.sortByTime(traces);
          final String newFilename = TraceMolecules.saveTraces(r, traces, createSettingsComment(),
              tracesFilename, additionalDatasets);
          // Only keep the main filename in memory
          if (additionalDatasets == 0) {
            tracesFilename = newFilename;
          }
        }
      }
    }

    final Trace[] all = allTraces.toArray(new Trace[allTraces.size()]);

    if (additionalDatasets > 0) {
      ImageJUtils.log("Multiple inputs provide %d traces", allTraces.size());

      final MemoryPeakResults tracedResults =
          TraceManager.toPeakResults(all, results.getCalibration(), true);
      tracedResults.copySettings(results);
      tracedResults.setName(createCombinedName() + " Tracks");
      MemoryPeakResults.addResults(tracedResults);
    }

    return all;
  }

  private String createCombinedName() {
    if (additionalDatasets > 0) {
      return results.getName() + " + " + TextUtils.pleural(additionalDatasets, "other");
    }
    return results.getName();
  }

  private boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    settings = SettingsManager.readClusteringSettings(0).toBuilder();

    gd.addCheckbox("Truncate_traces", settings.getTruncate());
    gd.addCheckbox("Internal_distances", settings.getInternalDistances());
    // gd.addCheckbox("Sub-sample_distances", settings.subSampledDistances);
    gd.addSlider("Fit_length", 2, 20, settings.getFitLength());
    gd.addCheckbox("MSD_correction", settings.getMsdCorrection());
    gd.addCheckbox("Precision_correction", settings.getPrecisionCorrection());
    gd.addCheckbox("Maximum_likelihood", settings.getMle());
    gd.addSlider("MLE_significance_level", 0, 0.5, significanceLevel);
    gd.addSlider("Fit_restarts", 0, 10, settings.getFitRestarts());
    gd.addSlider("Jump_distance", 1, 20, settings.getJumpDistance());
    gd.addSlider("Minimum_difference", 0, 10, minDifference);
    gd.addSlider("Minimum_fraction", 0, 1, minFraction);
    if (extraOptions) {
      gd.addSlider("Minimum_N", 1, 10, minN);
    }
    gd.addSlider("Maximum_N", 2, 10, maxN);
    gd.addCheckbox("Debug_fitting", debugFitting);
    gd.addCheckbox("Save_trace_distances", saveTraceDistances);
    gd.addCheckbox("Save_raw_data", saveRawData);
    gd.addCheckbox("Show_histograms", settings.getShowHistograms());
    gd.addStringField("Title", title);

    gd.showDialog();

    if (gd.wasCanceled() || !readDialog(gd)) {
      return false;
    }

    // Update the settings
    SettingsManager.writeSettings(settings.build());

    return true;
  }

  private boolean readDialog(ExtendedGenericDialog gd) {
    settings.setTruncate(gd.getNextBoolean());
    settings.setInternalDistances(gd.getNextBoolean());
    // settings.subSampledDistances = gd.getNextBoolean();
    settings.setFitLength((int) Math.abs(gd.getNextNumber()));
    settings.setMsdCorrection(gd.getNextBoolean());
    settings.setPrecisionCorrection(gd.getNextBoolean());
    settings.setMle(gd.getNextBoolean());
    significanceLevel = Math.abs(gd.getNextNumber());
    settings.setFitRestarts((int) Math.abs(gd.getNextNumber()));
    settings.setJumpDistance((int) Math.abs(gd.getNextNumber()));
    minDifference = Math.abs(gd.getNextNumber());
    minFraction = Math.abs(gd.getNextNumber());
    if (extraOptions) {
      myMinN = minN = (int) Math.abs(gd.getNextNumber());
    }
    maxN = (int) Math.abs(gd.getNextNumber());
    debugFitting = gd.getNextBoolean();
    saveTraceDistances = gd.getNextBoolean();
    saveRawData = gd.getNextBoolean();
    settings.setShowHistograms(gd.getNextBoolean());
    title = gd.getNextString();

    if (gd.invalidNumber()) {
      return false;
    }

    if (settings.getShowHistograms()) {
      gd = new ExtendedGenericDialog(TITLE);
      gd.addMessage("Select the histograms to display");
      gd.addCheckbox("Remove_outliers", settings.getRemoveOutliers());
      gd.addNumericField("Histogram_bins", settings.getHistogramBins(), 0);
      for (int i = 0; i < displayHistograms.length; i++) {
        gd.addCheckbox(NAMES[i].replace(' ', '_'), displayHistograms[i]);
      }
      gd.addCheckbox("MSD/Molecule", displayMSDHistogram);
      gd.addCheckbox("D/Molecule", displayDHistogram);
      gd.addCheckbox("Trace_length", displayTraceLength);
      gd.addCheckbox("Trace_size", displayTraceSize);
      gd.showDialog();
      if (gd.wasCanceled()) {
        return false;
      }
      settings.setRemoveOutliers(gd.getNextBoolean());
      settings.setHistogramBins((int) Math.abs(gd.getNextNumber()));
      for (int i = 0; i < displayHistograms.length; i++) {
        displayHistograms[i] = gd.getNextBoolean();
      }
      displayMSDHistogram = gd.getNextBoolean();
      displayDHistogram = gd.getNextBoolean();
      displayTraceLength = gd.getNextBoolean();
      displayTraceSize = gd.getNextBoolean();
    }

    // Check arguments
    try {
      // Parameters.isAboveZero("Histogram bins", settings.getHistogramBins());
      Parameters.isAbove("Fit length", settings.getFitLength(), 1);
      Parameters.isAboveZero("Jump distance", settings.getJumpDistance());
      Parameters.isEqualOrAbove("Maximum N", maxN, myMinN);
      if (settings.getMle()) {
        Parameters.isAboveZero("Significance level", significanceLevel);
      }
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return true;
  }

  private Plot2 plotMSD(double[] x, double[] y, double[] sd, String title) {
    if (saveRawData) {
      saveMSD(x, y, sd);
    }

    final Plot2 plot = new Plot2(title, "Time (s)", "Distance (um^2)", x, y);
    // Set limits before any plotting
    double max = 0;
    for (int i = 1; i < x.length; i++) {
      final double value = y[i] + sd[i];
      max = FastMath.max(max, value);
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
  @SuppressWarnings("null")
  private double[] fitMSD(double[] x, double[] y, String title, Plot2 plot) {
    // The Weimann paper (Plos One e64287) fits:
    // MSD(n dt) = 4D n dt + 4s^2
    // n = number of jumps
    // dt = time difference between frames
    // s = localisation precision
    // Thus we should fit an intercept as well.

    // From the fit D = gradient / (4*exposureTime)

    double D = 0;
    double intercept = 0;
    double precision = 0;

    final LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();
    Optimum lvmSolution;
    double ic = 0;

    // Fit with no intercept
    try {
      final LinearFunction function = new LinearFunction(x, y, settings.getFitLength());
      final double[] parameters = new double[] {function.guess()};

      //@formatter:off
      final LeastSquaresProblem problem = new LeastSquaresBuilder()
          .maxEvaluations(Integer.MAX_VALUE)
          .maxIterations(3000)
          .start(parameters)
          .target(function.getY())
          .weight(new DiagonalMatrix(function.getWeights()))
          .model(function, new MultivariateMatrixFunction() {
            @Override
            public double[][] value(double[] point) throws IllegalArgumentException
            {
              return function.jacobian(point);
            }} )
          .build();
      //@formatter:on

      lvmSolution = optimizer.optimize(problem);

      final double ss = lvmSolution.getResiduals().dotProduct(lvmSolution.getResiduals());
      // double ss = 0;
      // double[] obs = function.getY();
      // double[] exp = lvmSolution.getValue();
      // for (int i = 0; i < obs.length; i++)
      // ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);

      ic = MathUtils.getAkaikeInformationCriterionFromResiduals(ss, function.getY().length, 1);

      final double gradient = lvmSolution.getPoint().getEntry(0);
      D = gradient / 4;

      ImageJUtils.log(
          "Linear fit (%d points) : Gradient = %s, D = %s um^2/s, SS = %s, IC = %s (%d evaluations)",
          function.getY().length, MathUtils.rounded(gradient, 4), MathUtils.rounded(D, 4),
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
          new LinearFunctionWithIntercept(x, y, settings.getFitLength(), fitIntercept);

      //@formatter:off
      final LeastSquaresProblem problem = new LeastSquaresBuilder()
          .maxEvaluations(Integer.MAX_VALUE)
          .maxIterations(3000)
          .start(function.guess())
          .target(function.getY())
          .weight(new DiagonalMatrix(function.getWeights()))
          .model(function, new MultivariateMatrixFunction() {
            @Override
            public double[][] value(double[] point) throws IllegalArgumentException
            {
              return function.jacobian(point);
            }} )
          .build();
      //@formatter:on

      lvmSolution = optimizer.optimize(problem);

      final double ss = lvmSolution.getResiduals().dotProduct(lvmSolution.getResiduals());
      // double ss = 0;
      // double[] obs = function.getY();
      // double[] exp = lvmSolution.getValue();
      // for (int i = 0; i < obs.length; i++)
      // ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);

      final double ic2 =
          MathUtils.getAkaikeInformationCriterionFromResiduals(ss, function.getY().length, 2);
      final double gradient = lvmSolution.getPoint().getEntry(0);
      final double s = lvmSolution.getPoint().getEntry(1);
      final double intercept2 = 4 * s * s;

      if (ic2 < ic || debugFitting) {
        // Convert fitted precision in um to nm
        ImageJUtils.log(
            "Linear fit with intercept (%d points) : Gradient = %s, Intercept = %s, D = %s um^2/s, precision = %s nm, SS = %s, IC = %s (%d evaluations)",
            function.getY().length, MathUtils.rounded(gradient, 4),
            MathUtils.rounded(intercept2, 4), MathUtils.rounded(gradient / 4, 4),
            MathUtils.rounded(s * 1000, 4), MathUtils.rounded(ss), MathUtils.rounded(ic2),
            lvmSolution.getEvaluations());
      }

      if (lvmSolution == null || ic2 < ic) {
        intercept = intercept2;
        D = gradient / 4;
        precision = s;
      }
    } catch (final TooManyIterationsException ex) {
      ImageJUtils.log("Failed to fit with intercept : Too many iterations (%s)", ex.getMessage());
    } catch (final ConvergenceException ex) {
      ImageJUtils.log("Failed to fit with intercept : %s", ex.getMessage());
    }

    if (settings.getMsdCorrection()) {
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

        final LinearFunctionWithMSDCorrectedIntercept function =
            new LinearFunctionWithMSDCorrectedIntercept(x2, y, settings.getFitLength(),
                fitIntercept);

        //@formatter:off
final LeastSquaresProblem problem = new LeastSquaresBuilder()
          .maxEvaluations(Integer.MAX_VALUE)
          .maxIterations(3000)
          .start(function.guess())
          .target(function.getY())
          .weight(new DiagonalMatrix(function.getWeights()))
          .model(function, new MultivariateMatrixFunction() {
            @Override
            public double[][] value(double[] point) throws IllegalArgumentException
            {
              return function.jacobian(point);
            }} )
          .build();
//@formatter:on

        lvmSolution = optimizer.optimize(problem);

        final double ss = lvmSolution.getResiduals().dotProduct(lvmSolution.getResiduals());
        // double ss = 0;
        // double[] obs = function.getY();
        // double[] exp = lvmSolution.getValue();
        // for (int i = 0; i < obs.length; i++)
        // ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);

        final double ic2 =
            MathUtils.getAkaikeInformationCriterionFromResiduals(ss, function.getY().length, 2);
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

        if (ic2 < ic || debugFitting) {
          // Convert fitted precision in um to nm
          ImageJUtils.log(
              "Linear fit with MSD corrected intercept (%d points) : Gradient = %s, Intercept = %s, D = %s um^2/s, precision = %s nm, SS = %s, IC = %s (%d evaluations)",
              function.getY().length, MathUtils.rounded(gradient, 4),
              MathUtils.rounded(intercept2, 4), MathUtils.rounded(gradient / 4, 4),
              MathUtils.rounded(s * 1000, 4), MathUtils.rounded(ss), MathUtils.rounded(ic2),
              lvmSolution.getEvaluations());
        }

        if (lvmSolution == null || ic2 < ic) {
          intercept = intercept2;
          D = gradient / 4;
          precision = s;
        }
      } catch (final TooManyIterationsException ex) {
        ImageJUtils.log("Failed to fit with intercept : Too many iterations (%s)", ex.getMessage());
      } catch (final ConvergenceException ex) {
        ImageJUtils.log("Failed to fit with intercept : %s", ex.getMessage());
      }
    }

    // Add the fit to the plot
    if (D > 0) {
      plot.setColor(Color.magenta);
      plot.drawLine(0, intercept, x[x.length - 1], 4 * D * x[x.length - 1] + intercept);
      display(title, plot);

      checkTraceDistance(D);
    }

    return new double[] {D, precision};
  }

  /**
   * Check the distance used for tracing covers enough of the cumulative mean-squared distance
   * distribution.
   *
   * @param d the distance
   */
  private void checkTraceDistance(double d) {
    final double t = exposureTime;
    // Cumul P(r^2) = 1 - exp(-r^2 / 4dt)
    final double r = settings.getDistanceThreshold() / 1000;
    final double msd = 4 * d * t;
    final double p = 1 - FastMath.exp(-r * r / msd);
    ImageJUtils.log("Checking trace distance: r = %s nm, D = %s um^2/s, Cumul p(r^2|frame) = %s",
        settings.getDistanceThreshold(), MathUtils.rounded(d), MathUtils.rounded(p));
    if (p < 0.95) {
      ImageJUtils.log("WARNING *** The tracing distance may not be large enough! ***");
    }
  }

  private class LinearFunction implements MultivariateVectorFunction {
    double[] x;
    double[] y;
    double[][] jacobian;

    public LinearFunction(double[] x, double[] y, int length) {
      final int to = FastMath.min(x.length, 1 + length);
      this.x = Arrays.copyOfRange(x, 1, to);
      this.y = Arrays.copyOfRange(y, 1, to);
      jacobian = calculateJacobian();
    }

    // Adapted from http://commons.apache.org/proper/commons-math/userguide/optimization.html

    /**
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

    /** {@inheritDoc} */
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
      // dy_da = x
      final double[][] jacobian = new double[x.length][1];

      for (int i = 0; i < jacobian.length; ++i) {
        jacobian[i][0] = x[i];
      }

      return jacobian;
    }
  }

  private class LinearFunctionWithIntercept implements MultivariateVectorFunction {
    final double[] x;
    final double[] y;
    final boolean fitIntercept;

    public LinearFunctionWithIntercept(double[] x, double[] y, int length, boolean fitIntercept) {
      this.fitIntercept = fitIntercept;
      final int to = FastMath.min(x.length, 1 + length);
      // Optionally include the intercept
      final int from = (fitIntercept) ? 0 : 1;
      this.x = Arrays.copyOfRange(x, from, to);
      this.y = Arrays.copyOfRange(y, from, to);
    }

    /**
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

    /** {@inheritDoc} */
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
      // dy_da = x
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

  private class LinearFunctionWithMSDCorrectedIntercept implements MultivariateVectorFunction {
    final double THIRD = 1 / 3.0;
    final double[] x;
    final double[] y;
    final boolean fitIntercept;

    public LinearFunctionWithMSDCorrectedIntercept(double[] x, double[] y, int length,
        boolean fitIntercept) {
      this.fitIntercept = fitIntercept;
      final int to = FastMath.min(x.length, 1 + length);
      // Optionally include the intercept
      final int from = (fitIntercept) ? 0 : 1;
      this.x = Arrays.copyOfRange(x, from, to);
      this.y = Arrays.copyOfRange(y, from, to);
    }

    /**
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

    /** {@inheritDoc} */
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
      int i = 0;
      // Special case for fitting the intercept since the line is not linear below n=1
      if (fitIntercept) {
        values[i++] = intercept;
      }
      for (; i < values.length; i++) {
        values[i] = a * x[i] + error;
      }
      return values;
    }

    double[][] jacobian(double[] variables) {
      // Compute the gradients using calculus differentiation:
      // y = ax - a/3 + 4c^2
      // dy_da = x - 1/3
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
        "Jump Distance analysis : N = %d, Time = %d frames (%s seconds). MSD = %s um^2/jump, Mean Distance = %s nm/jump, Precision = %s nm, Beta = %s, minD = %s um^2/jump",
        jumpDistances.getN(), settings.getJumpDistance(),
        MathUtils.rounded(settings.getJumpDistance() * exposureTime, 4), MathUtils.rounded(msd, 4),
        MathUtils.rounded(meanDistance, 4), MathUtils.rounded(precision, 4),
        MathUtils.rounded(beta, 4), MathUtils.rounded(minD, 4));

    final Logger logger = ImageJPluginLoggerHelper.getLogger(getClass());
    if (debugFitting) {
      logger.setLevel(Level.FINE);
    }
    final JumpDistanceAnalysis jd = new JumpDistanceAnalysis(logger);
    jd.setFitRestarts(settings.getFitRestarts());
    jd.setMinFraction(minFraction);
    jd.setMinDifference(minDifference);
    jd.setMinN(myMinN);
    jd.setMaxN(maxN);
    jd.setMinD(minD);
    jd.setSignificanceLevel(significanceLevel);
    // Update the plot with the fit
    jd.setCurveLogger(this);

    // Set the calibration
    jd.setN(settings.getJumpDistance());
    jd.setDeltaT(exposureTime);
    if (settings.getPrecisionCorrection()) {
      jd.setError(precision, true);
    }
    jd.setMsdCorrection(settings.getMsdCorrection());

    double[][] fit;
    if (settings.getMle()) {
      fit = jd.fitJumpDistancesMLE(jumpDistances.getValues(), jdHistogram);
    } else {
      fit = jd.fitJumpDistanceHistogram(jumpDistances.getMean(), jdHistogram);
    }

    // Get the raw fitted D and convert it to a calibrated D*
    if (fit != null) {
      fit[0] = jd.calculateApparentDiffusionCoefficient(fit[0]);
      // Check the largest D
      checkTraceDistance(fit[0][0]);
      fitValue = jd.getFitValue();
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
    if (saveRawData) {
      saveFit(x, y, "Single");
    }
    addToJumpDistancePlot(x, y, Color.magenta);
  }

  @Override
  public void saveMixedPopulationCurve(double[][] curve) {
    final double[] x = curve[0];
    final double[] y = curve[1];
    if (saveRawData) {
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
   * @param args 0: Double[1] - output the number of species
   * @return Empty string
   */
  public static String getNumberOfSpecies(Object[] args) {
    int n = 0;
    if (jumpDistanceParameters != null) {
      n = jumpDistanceParameters[0].length;
    }
    final Double[] array = (Double[]) args[0];
    array[0] = new Double(n);
    return "";
  }

  /**
   * Macro extension function.
   *
   * <p>Get the diffusion coefficient for the requested species from the last call to fit the jump
   * distances.
   *
   * @param args 0: Double[1] - input the index of the species; 1: Double[1] - output the
   *        coefficient
   * @return Empty string
   */
  public static String getD(Object[] args) {
    double value = 0;
    if (jumpDistanceParameters != null) {
      final int i = ((Double) args[0]).intValue();
      if (i >= 0 && i < jumpDistanceParameters[0].length) {
        value = jumpDistanceParameters[0][i];
      }
    }
    ((Double[]) args[1])[0] = new Double(value);
    return "";
  }

  /**
   * Macro extension function.
   *
   * <p>Get the population fraction for the requested species from the last call to fit the jump
   * distances.
   *
   * @param args 0: Double[1] - input the index of the species; 1: Double[1] - output the population
   *        fraction
   * @return Empty string
   */
  public static String getF(Object[] args) {
    double value = 0;
    if (jumpDistanceParameters != null) {
      final int i = ((Double) args[0]).intValue();
      if (i >= 0 && i < jumpDistanceParameters[1].length) {
        value = jumpDistanceParameters[1][i];
      }
    }
    ((Double[]) args[1])[0] = new Double(value);
    return "";
  }

  /**
   * Macro extension function.
   *
   * <p>Get the diffusion coefficient and population fraction for the requested species from the
   * last call to fit the jump distances.
   *
   * @param args 0: Double[1] - input the index of the species; 1: Double[1] - output the
   *        coefficient; 1: Double[1] - output the population fraction
   * @return Empty string
   */
  public static String getSpecies(Object[] args) {
    double value = 0;
    double value2 = 0;
    if (jumpDistanceParameters != null) {
      final int i = ((Double) args[0]).intValue();
      if (i >= 0 && i < jumpDistanceParameters[0].length) {
        value = jumpDistanceParameters[0][i];
        value2 = jumpDistanceParameters[1][i];
      }
    }
    ((Double[]) args[1])[0] = new Double(value);
    ((Double[]) args[2])[0] = new Double(value2);
    return "";
  }

  private boolean showMultiDialog(ArrayList<MemoryPeakResults> allResults) {
    multiMode = true;

    // Show a list box containing all the results. This should remember the last set of chosen
    // items.
    final MultiDialog md = new MultiDialog(TITLE, new MultiDialog.MemoryResultsItems());
    md.addSelected(selected);

    md.showDialog();

    if (md.wasCancelled()) {
      return false;
    }

    selected = md.getSelectedResults();
    if (selected.isEmpty()) {
      IJ.error(TITLE, "No results were selected");
      return false;
    }

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
      final MemoryPeakResults results = allResults.get(1);

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
