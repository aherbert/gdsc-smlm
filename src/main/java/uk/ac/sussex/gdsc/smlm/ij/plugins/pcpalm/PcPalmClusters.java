/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.ij.plugins.pcpalm;

import uk.ac.sussex.gdsc.core.clustering.Cluster;
import uk.ac.sussex.gdsc.core.clustering.ClusterPoint;
import uk.ac.sussex.gdsc.core.clustering.ClusteringAlgorithm;
import uk.ac.sussex.gdsc.core.clustering.ClusteringEngine;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot;
import uk.ac.sussex.gdsc.core.ij.ImageJPluginLoggerHelper;
import uk.ac.sussex.gdsc.core.ij.ImageJTrackProgress;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.Plot2;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationHelper;
import uk.ac.sussex.gdsc.smlm.fitting.BinomialFitter;
import uk.ac.sussex.gdsc.smlm.ij.plugins.About;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ParameterUtils;
import uk.ac.sussex.gdsc.smlm.ij.plugins.SmlmUsageTracker;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.ExtendedPeakResult;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;

import ij.IJ;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.plugin.frame.Recorder;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.util.FastMath;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.InputMismatchException;
import java.util.List;
import java.util.Locale;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.regex.Pattern;

/**
 * Find clusters of molecules using a partial centroid-linkage hierarchical clustering algorithm.
 *
 *
 * <p>Points are added to the nearest cluster if they are below the distance threshold to the
 * cluster centroid. The cluster centroid is updated. All points above the cluster distance
 * threshold remain as single molecules.
 *
 * <p>The purpose is to join colocalising molecules into clusters.
 *
 * <p>See Puchnar, et al (2013). Counting molecules in single organelles with superresolution
 * microscopy allows tracking of the endosome maturation trajectory. PNAS.
 * doi:10.1073/pnas.1309676110
 */
public class PcPalmClusters implements PlugIn {
  /** The title. */
  static final String TITLE = "PC-PALM Clusters";

  private static int runMode;
  private static double distance = 50;
  private static ClusteringAlgorithm sClusteringAlgorithm =
      ClusteringAlgorithm.PARTICLE_CENTROID_LINKAGE;
  private static int minN = 1;
  private static int maxN;
  private static boolean maximumLikelihood;
  private static boolean showCumulativeHistogram;
  private static boolean multiThread = true;
  private static boolean sWeightedClustering;
  private static boolean saveHistogram;
  private static String histogramFile = "";
  private static String noiseFile = "";
  private static boolean sAutoSave = true;
  private boolean autoSave = true;
  private static boolean calibrateHistogram;
  private static int frames = 1;
  private static double area = 1;
  private static final String[] UNITS = {"pixels^2", "um^2"};
  private static String units = UNITS[0];

  private ClusteringAlgorithm clusteringAlgorithm = ClusteringAlgorithm.PARTICLE_CENTROID_LINKAGE;
  private boolean weightedClustering;

  private boolean fileInput;

  private int numberOfMolecules;
  private double count;

  private static class HistogramData {
    float[][] histogram;
    int frames;
    double area;
    String units;
    String filename = "";

    HistogramData(float[][] histogram, int frames, double area, String units) {
      this.histogram = histogram;
      this.frames = frames;
      this.area = area;
      this.units = units;
    }

    HistogramData(float[][] histogram) {
      this.histogram = histogram;
      units = "";
    }

    boolean isCalibrated() {
      return frames > 0 && area > 0;
    }
  }

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (!showDialog()) {
      return;
    }

    PcPalmMolecules.logSpacer();
    ImageJUtils.log(TITLE);
    PcPalmMolecules.logSpacer();
    final long start = System.currentTimeMillis();

    HistogramData histogramData;
    if (fileInput) {
      histogramData = loadHistogram(histogramFile);
    } else {
      histogramData = doClustering();
    }

    if (histogramData == null) {
      return;
    }
    final float[][] hist = histogramData.histogram;

    // Create a histogram of the cluster sizes
    String title = TITLE + " Molecules/cluster";
    final String xTitle = "Molecules/cluster";
    final String yTitle = "Frequency";

    // Create the data required for fitting and plotting
    float[] xValues = HistogramPlot.createHistogramAxis(hist[0]);
    float[] yValues = HistogramPlot.createHistogramValues(hist[1]);

    // Plot the histogram
    float yMax = MathUtils.max(yValues);
    Plot2 plot = new Plot2(title, xTitle, yTitle, xValues, yValues);
    if (xValues.length > 0) {
      final double xPadding = 0.05 * (xValues[xValues.length - 1] - xValues[0]);
      plot.setLimits(xValues[0] - xPadding, xValues[xValues.length - 1] + xPadding, 0, yMax * 1.05);
    }
    ImageJUtils.display(title, plot);

    final HistogramData noiseData = loadNoiseHistogram(histogramData);
    if (noiseData != null) {
      if (subtractNoise(histogramData, noiseData)) {
        // Update the histogram
        title += " (noise subtracted)";
        xValues = HistogramPlot.createHistogramAxis(hist[0]);
        yValues = HistogramPlot.createHistogramValues(hist[1]);
        yMax = MathUtils.max(yValues);
        plot = new Plot2(title, xTitle, yTitle, xValues, yValues);
        if (xValues.length > 0) {
          final double xPadding = 0.05 * (xValues[xValues.length - 1] - xValues[0]);
          plot.setLimits(xValues[0] - xPadding, xValues[xValues.length - 1] + xPadding, 0,
              yMax * 1.05);
        }
        ImageJUtils.display(title, plot);

        // Automatically save
        if (autoSave) {
          final String newFilename =
              ImageJUtils.replaceExtension(histogramData.filename, ".noise.tsv");
          if (saveHistogram(histogramData, newFilename)) {
            ImageJUtils.log("Saved noise-subtracted histogram to " + newFilename);
          }
        }
      }
    }

    // Fit the histogram
    final double[] fitParameters = fitBinomial(histogramData);
    if (fitParameters != null) {
      // Add the binomial to the histogram
      final int n = (int) fitParameters[0];
      final double p = fitParameters[1];

      ImageJUtils.log("Optimal fit : N=%d, p=%s", n, MathUtils.rounded(p));

      final BinomialDistribution dist = new BinomialDistribution(n, p);

      // A zero-truncated binomial was fitted.
      // pi is the adjustment factor for the probability density.
      final double pi = 1 / (1 - dist.probability(0));

      if (!fileInput) {
        // Calculate the estimated number of clusters from the observed molecules:
        // Actual = (Observed / p-value) / N
        final double actual = (numberOfMolecules / p) / n;
        ImageJUtils.log("Estimated number of clusters : (%d / %s) / %d = %s", numberOfMolecules,
            MathUtils.rounded(p), n, MathUtils.rounded(actual));
      }

      final double[] x = new double[n + 2];
      final double[] y = new double[n + 2];

      // Scale the values to match those on the histogram
      final double normalisingFactor = count * pi;
      for (int i = 0; i <= n; i++) {
        x[i] = i + 0.5;
        y[i] = dist.probability(i) * normalisingFactor;
      }
      x[n + 1] = n + 1.5;
      y[n + 1] = 0;

      // Redraw the plot since the limits may have changed
      plot = new Plot2(title, xTitle, yTitle, xValues, yValues);
      final double xPadding = 0.05 * (xValues[xValues.length - 1] - xValues[0]);
      plot.setLimits(xValues[0] - xPadding, xValues[xValues.length - 1] + xPadding, 0,
          MathUtils.maxDefault(yMax, y) * 1.05);
      plot.setColor(Color.magenta);
      plot.addPoints(x, y, Plot.LINE);
      plot.addPoints(x, y, Plot.CIRCLE);
      plot.setColor(Color.black);
      ImageJUtils.display(title, plot);
    }

    final double seconds = (System.currentTimeMillis() - start) / 1000.0;
    final String msg = TITLE + " complete : " + seconds + "s";
    IJ.showStatus(msg);
    ImageJUtils.log(msg);
    return;
  }

  /**
   * Extract the results from the PCPALM molecules using the area ROI and then do clustering to
   * obtain the histogram of molecules per cluster.
   *
   * @return the histogram data
   */
  private HistogramData doClustering() {
    // Perform clustering analysis to generate the histogram of cluster sizes
    final PcPalmAnalysis analysis = new PcPalmAnalysis();
    final ArrayList<Molecule> molecules = analysis.cropToRoi(WindowManager.getCurrentImage());

    if (molecules.size() < 2) {
      error("No results within the crop region");
      return null;
    }

    ImageJUtils.log("Using %d molecules (Density = %s um^-2) @ %s nm", molecules.size(),
        MathUtils.rounded(molecules.size() / analysis.croppedArea), MathUtils.rounded(distance));

    final long s1 = System.nanoTime();
    final ClusteringEngine engine =
        new ClusteringEngine(1, clusteringAlgorithm, new ImageJTrackProgress());
    if (multiThread) {
      engine.setThreadCount(Prefs.getThreads());
    }
    engine.setTracker(new ImageJTrackProgress());
    IJ.showStatus("Clustering ...");
    final List<Cluster> clusters = engine.findClusters(convertToPoint(molecules), distance);
    IJ.showStatus("");

    if (clusters == null) {
      ImageJUtils.log("Aborted");
      return null;
    }
    numberOfMolecules = molecules.size();
    ImageJUtils.log("Finished : %d total clusters (%s ms)", clusters.size(),
        MathUtils.rounded((System.nanoTime() - s1) / 1e6));

    // Save cluster centroids to a results set in memory. Then they can be plotted.
    final MemoryPeakResults results = new MemoryPeakResults(clusters.size());
    results.setName(TITLE);
    // Set an arbitrary calibration so that the lifetime of the results is stored in the exposure
    // time
    // The results will be handled as a single mega-frame containing all localisation.
    results.setCalibration(CalibrationHelper.create(100, 1, PcPalmMolecules.seconds * 1000));
    int id = 0;
    for (final Cluster c : clusters) {
      results.add(new ExtendedPeakResult((float) c.getX(), (float) c.getY(), c.getSize(), ++id));
    }
    MemoryPeakResults.addResults(results);

    // Get the data for fitting
    final float[] values = new float[clusters.size()];
    for (int i = 0; i < values.length; i++) {
      values[i] = clusters.get(i).getSize();
    }
    final float yMax = (int) Math.ceil(MathUtils.max(values));
    final int nBins = (int) (yMax + 1);
    final float[][] hist = HistogramPlot.calcHistogram(values, 0, yMax, nBins);

    final HistogramData histogramData =
        (calibrateHistogram) ? new HistogramData(hist, frames, area, units)
            : new HistogramData(hist);

    saveHistogram(histogramData);

    return histogramData;
  }

  /**
   * Convert molecules for clustering.
   *
   * @param molecules the molecules
   * @return the list of cluster points
   */
  private List<ClusterPoint> convertToPoint(ArrayList<Molecule> molecules) {
    final ArrayList<ClusterPoint> points = new ArrayList<>(molecules.size());
    int id = 0;
    for (final Molecule m : molecules) {
      points
          .add(ClusterPoint.newClusterPoint(id++, m.x, m.y, (weightedClustering) ? m.photons : 1));
    }
    return points;
  }

  /**
   * Saves the histogram to the user selected file if the save histogram option is enabled.
   *
   * @param histogramData the histogram data
   * @return true, if successful
   */
  private static boolean saveHistogram(HistogramData histogramData) {
    if (!saveHistogram) {
      return false;
    }
    histogramFile = ImageJUtils.getFilename("Histogram_file", histogramFile);
    return saveHistogram(histogramData, histogramFile);
  }

  /**
   * Saves the histogram to the selected file. Updates the filename property of the histogram
   * object.
   *
   * @param histogramData the histogram data
   * @param filename the filename
   * @return true, if successful
   */
  private static boolean saveHistogram(HistogramData histogramData, String filename) {
    if (filename == null) {
      return false;
    }

    final float[][] hist = histogramData.histogram;

    filename = ImageJUtils.replaceExtension(filename, "tsv");
    try (BufferedWriter output = Files.newBufferedWriter(Paths.get(filename))) {
      if (histogramData.isCalibrated()) {
        output.write(String.format("Frames  %d", histogramData.frames));
        output.newLine();
        output.write(String.format("Area    %f", histogramData.area));
        output.newLine();
        output.write(String.format("Units   %s", histogramData.units));
        output.newLine();
      }
      output.write("Size\tFrequency");
      output.newLine();
      for (int i = 0; i < hist[0].length; i++) {
        output.write(String.format("%d\t%s", (int) hist[0][i], MathUtils.rounded(hist[1][i])));
        output.newLine();
      }

      histogramData.filename = filename;
      return true;
    } catch (final Exception ex) {
      ex.printStackTrace();
      IJ.log("Failed to save histogram to file: " + filename);
    }
    return false;
  }

  /**
   * Load the histogram from the file. Assumes the histogram is [int, float] format and creates a
   * contiguous histogram from zero
   *
   * @param filename the filename
   * @return the histogram data
   */
  private static HistogramData loadHistogram(String filename) {
    int count = 0;

    try (BufferedReader input = Files.newBufferedReader(Paths.get(filename))) {
      int frames = 0;
      double area = 0;
      String units = "";

      String line;

      final ArrayList<float[]> data = new ArrayList<>();

      // Read the header and store the calibration if present
      while ((line = input.readLine()) != null) {
        count++;
        if (line.length() == 0) {
          continue;
        }
        if (Character.isDigit(line.charAt(0))) {
          // This is the first record
          break;
        }
        final String[] fields = line.split("[\t, ]+");
        if (fields[0].equalsIgnoreCase("frames")) {
          frames = Integer.parseInt(fields[1]);
        }
        if (fields[0].equalsIgnoreCase("area")) {
          area = Double.parseDouble(fields[1]);
        }
        if (fields[0].equalsIgnoreCase("units")) {
          units = fields[1];
        }
      }

      final Pattern pattern = Pattern.compile("[\t, ]+");
      while (line != null) {
        if (line.length() == 0) {
          continue;
        }
        if (!Character.isDigit(line.charAt(0))) {
          continue;
        }

        // Extract the first 2 fields
        try (Scanner scanner = new Scanner(line)) {
          scanner.useLocale(Locale.US);
          scanner.useDelimiter(pattern);
          final int molecules = scanner.nextInt();
          final float frequency = scanner.nextFloat();

          // Check for duplicates
          for (final float[] d : data) {
            if (d[0] == molecules) {
              error("Duplicate molecules field on line " + count);
              return null;
            }
          }

          data.add(new float[] {molecules, frequency});
        }

        // Get the next line
        line = input.readLine();
        count++;
      }

      if (data.isEmpty()) {
        error("No data in file " + filename);
        return null;
      }

      // Create a contiguous histogram from zero
      int maxN = 0;
      for (final float[] d : data) {
        if (maxN < d[0]) {
          maxN = (int) d[0];
        }
      }

      final float[][] hist = new float[2][maxN + 1];
      for (int n = 0; n <= maxN; n++) {
        hist[0][n] = n;
        for (final float[] d : data) {
          if (n == d[0]) {
            hist[1][n] = d[1];
          }
        }
      }
      final HistogramData histogramData = new HistogramData(hist, frames, area, units);
      histogramData.filename = filename;
      return histogramData;
    } catch (final InputMismatchException ex) {
      error("Incorrect fields on line " + count);
    } catch (final NoSuchElementException ex) {
      error("Incorrect fields on line " + count);
    } catch (final IOException ex) {
      IJ.error(TITLE, "Unable to read from file " + filename);
    }
    return null;
  }

  /**
   * If the histogram is calibrated then ask the user if they wish to subtract a calibrated noise
   * histogram.
   *
   * <p>Loads a noise histogram from a user selected file and check the units match those provided
   *
   * @param histogramData the histogram data
   * @return The histogram (or null)
   */
  private HistogramData loadNoiseHistogram(HistogramData histogramData) {
    if (!histogramData.isCalibrated()) {
      return null;
    }
    final GenericDialog gd = new GenericDialog(TITLE);
    gd.enableYesNoCancel();
    gd.hideCancelButton();
    gd.addMessage("The histogram is calibrated.\n \n"
        + "Do you want to subtract a noise histogram before fitting?");
    final boolean allowSave = new File(histogramData.filename).exists();
    if (allowSave) {
      gd.addCheckbox("Auto_save noise-subtracted histogram", sAutoSave);
    }

    // If this is a macro then the dialog will not have Yes or No pressed.
    // Add a checkbox that can be read from the macro arguments by ImageJ.
    final String macroOption = "subtract";
    if (IJ.isMacro()) {
      gd.addCheckbox(macroOption, true);
    }

    gd.showDialog();
    if (!gd.wasOKed()) {
      return null;
    }
    if (allowSave) {
      autoSave = sAutoSave = gd.getNextBoolean();
    }

    if (IJ.isMacro()) {
      // If the macro option flag is not found then the arguments do not want this to run
      if (!gd.getNextBoolean()) {
        return null;
      }
    } else {
      // Ensure that the 'Yes' result is recorded for macros to detect
      Recorder.recordOption(macroOption);
    }

    noiseFile = ImageJUtils.getFilename("Noise_file", noiseFile);
    if (noiseFile != null) {
      final HistogramData data = loadHistogram(noiseFile);
      // Check the data is calibrated with the same units
      if (data.isCalibrated() && data.units.equalsIgnoreCase(histogramData.units)) {
        return data;
      }
    }
    return null;
  }

  private boolean showDialog() {
    if (PcPalmMolecules.molecules == null || PcPalmMolecules.molecules.size() < 2) {
      ImageJUtils.log(TITLE + " defaulting to File mode");
      fileInput = true;
      // Ensure this gets recorded
      Recorder.recordOption("Method", "File");
    } else {
      final GenericDialog gd = new GenericDialog(TITLE);
      final String[] items = {"Clustering", "File"};

      gd.addMessage("Fit a Binomial distribution to a histogram of cluster sizes.\n \n"
          + "Select the method to generate the histogram:");
      gd.addChoice("Method", items, items[runMode]);
      gd.showDialog();
      if (gd.wasCanceled()) {
        return false;
      }
      runMode = gd.getNextChoiceIndex();
      fileInput = (runMode == 1);
    }

    if (fileInput) {
      if ((histogramFile = ImageJUtils.getFilename("Histogram_file", histogramFile)) == null) {
        return false;
      }
    }

    final GenericDialog gd = new GenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    // Check if the molecules have weights
    boolean haveWeights = false;
    if (!fileInput) {
      haveWeights = checkForWeights();

      gd.addMessage("Find clusters using centroid-linkage clustering.");

      gd.addNumericField("Distance (nm)", distance, 0);
      final String[] clusteringAlgorithmNames =
          SettingsManager.getNames((Object[]) ClusteringAlgorithm.values());
      gd.addChoice("Algorithm", clusteringAlgorithmNames,
          clusteringAlgorithmNames[sClusteringAlgorithm.ordinal()]);
      gd.addCheckbox("Multi_thread", multiThread);
      if (haveWeights) {
        gd.addCheckbox("Weighted_clustering", sWeightedClustering);
      }
    }

    gd.addSlider("Min_N", 1, 10, minN);
    gd.addSlider("Max_N", 0, 10, maxN);
    gd.addCheckbox("Show_cumulative_histogram", showCumulativeHistogram);
    gd.addCheckbox("Maximum_likelihood", maximumLikelihood);

    if (!fileInput) {
      gd.addCheckbox("Save_histogram", saveHistogram);
      gd.addMessage("Histogram calibration (optional)");
      gd.addCheckbox("Calibrate_histogram", calibrateHistogram);
      gd.addNumericField("Frames", frames, 0);
      gd.addNumericField("Area", area, 2);
      gd.addChoice("Units", UNITS, units);
    }

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    if (!fileInput) {
      distance = gd.getNextNumber();
      clusteringAlgorithm =
          sClusteringAlgorithm = ClusteringAlgorithm.values()[gd.getNextChoiceIndex()];
      multiThread = gd.getNextBoolean();
      if (haveWeights) {
        weightedClustering = sWeightedClustering = gd.getNextBoolean();
      }
    }
    minN = (int) Math.abs(gd.getNextNumber());
    maxN = (int) Math.abs(gd.getNextNumber());
    showCumulativeHistogram = gd.getNextBoolean();
    maximumLikelihood = gd.getNextBoolean();
    if (!fileInput) {
      saveHistogram = gd.getNextBoolean();
      calibrateHistogram = gd.getNextBoolean();
      frames = (int) Math.abs(gd.getNextNumber());
      area = Math.abs(gd.getNextNumber());
      units = gd.getNextChoice();
    }

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Min N", minN);
      if (!fileInput) {
        ParameterUtils.isAboveZero("Distance", distance);
        ParameterUtils.isAboveZero("Frames", frames);
        ParameterUtils.isAboveZero("Area", area);
      }
    } catch (final IllegalArgumentException ex) {
      error(ex.getMessage());
      return false;
    }

    return true;
  }

  /**
   * Check if all the molecules have weights (allowing weighted clustering).
   *
   * @return True if all the molecules have weights (allowing weighted clustering).
   */
  private static boolean checkForWeights() {
    for (final Molecule m : PcPalmMolecules.molecules) {
      if (m.photons <= 0) {
        return false;
      }
    }
    return true;
  }

  private static void error(String message) {
    ImageJUtils.log("ERROR : " + message);
    IJ.error(TITLE, message);
  }

  /**
   * Normalise the histograms using the (frames*area). Subtract the noise from the histogram and
   * then rescale.
   *
   * @param histogramData the histogram data
   * @param noiseData the noise data
   * @return true, if successful
   */
  private static boolean subtractNoise(HistogramData histogramData, HistogramData noiseData) {
    final float[] v1 = normalise(histogramData);
    final float[] v2 = normalise(noiseData);
    final int length = v1.length; // FastMath.max(v1.length, v2.length);
    final double factor = (histogramData.frames * histogramData.area);
    for (int i = 0; i < length; i++) {
      histogramData.histogram[1][i] =
          (float) (FastMath.max(0, v1[i] - ((i < v2.length) ? v2[i] : 0)) * factor);
    }
    return true;
  }

  /**
   * Normalise the histogram using the (frames*area).
   *
   * @param data the data
   * @return the normalised data
   */
  private static float[] normalise(HistogramData data) {
    final float[] values = Arrays.copyOf(data.histogram[1], data.histogram[1].length);
    final double normalisingFactor = 1.0 / (data.frames * data.area);
    for (int i = 0; i < values.length; i++) {
      values[i] *= normalisingFactor;
    }
    return values;
  }

  /**
   * Fit a zero-truncated Binomial to the cumulative histogram.
   *
   * @param histogramData the histogram data
   * @return the double[]
   */
  private double[] fitBinomial(HistogramData histogramData) {
    // Get the mean and sum of the input histogram
    double mean;
    double sum = 0;
    count = 0;
    for (int i = 0; i < histogramData.histogram[1].length; i++) {
      count += histogramData.histogram[1][i];
      sum += histogramData.histogram[1][i] * i;
    }
    mean = sum / count;

    final String name = "Zero-truncated Binomial distribution";
    ImageJUtils.log("Mean cluster size = %s", MathUtils.rounded(mean));
    ImageJUtils.log("Fitting cumulative " + name);

    // Convert to a normalised double array for the binomial fitter
    final double[] histogram = new double[histogramData.histogram[1].length];
    for (int i = 0; i < histogramData.histogram[1].length; i++) {
      histogram[i] = histogramData.histogram[1][i] / count;
    }

    // Plot the cumulative histogram
    final String title = TITLE + " Cumulative Distribution";
    Plot2 plot = null;
    if (showCumulativeHistogram) {
      // Create a cumulative histogram for fitting
      final double[] cumulativeHistogram = new double[histogram.length];
      sum = 0;
      for (int i = 0; i < histogram.length; i++) {
        sum += histogram[i];
        cumulativeHistogram[i] = sum;
      }

      final double[] values = SimpleArrayUtils.newArray(histogram.length, 0.0, 1.0);
      plot = new Plot2(title, "N", "Cumulative Probability", values, cumulativeHistogram);
      plot.setLimits(0, histogram.length - 1, 0, 1.05);
      plot.addPoints(values, cumulativeHistogram, Plot.CIRCLE);
      ImageJUtils.display(title, plot);
    }

    // Do fitting for different N
    double bestSs = Double.POSITIVE_INFINITY;
    double[] parameters = null;
    int worse = 0;
    int countN = histogram.length - 1;
    int min = minN;
    final boolean customRange = (minN > 1) || (maxN > 0);
    if (min > countN) {
      min = countN;
    }
    if (maxN > 0 && countN > maxN) {
      countN = maxN;
    }

    ImageJUtils.log("Fitting N from %d to %d%s", min, countN,
        (customRange) ? " (custom-range)" : "");

    // Since varying the N should be done in integer steps do this
    // for n=1,2,3,... until the SS peaks then falls off (is worse then the best
    // score several times in succession)
    final BinomialFitter bf = new BinomialFitter(ImageJPluginLoggerHelper.getLogger(getClass()));
    bf.setMaximumLikelihood(maximumLikelihood);
    for (int n = min; n <= countN; n++) {
      final PointValuePair solution = bf.fitBinomial(histogram, mean, n, true);
      if (solution == null) {
        continue;
      }

      final double p = solution.getPointRef()[0];

      ImageJUtils.log("Fitted %s : N=%d, p=%s. SS=%g", name, n, MathUtils.rounded(p),
          solution.getValue());

      if (bestSs > solution.getValue()) {
        bestSs = solution.getValue();
        parameters = new double[] {n, p};
        worse = 0;
      } else if (bestSs < Double.POSITIVE_INFINITY && ++worse >= 3) {
        break;
      }

      if (showCumulativeHistogram) {
        addToPlot(n, p, title, plot, new Color((float) n / countN, 0, 1f - (float) n / countN));
      }
    }

    // Add best it in magenta
    if (showCumulativeHistogram && parameters != null) {
      addToPlot((int) parameters[0], parameters[1], title, plot, Color.magenta);
    }

    return parameters;
  }

  private static void addToPlot(int trials, double pvalue, String title, Plot2 plot, Color color) {
    final double[] x = new double[trials + 1];
    final double[] y = new double[trials + 1];

    final BinomialDistribution dist = new BinomialDistribution(trials, pvalue);

    final int startIndex = 1;

    // Normalise optionally excluding the x=0 point
    double total = 1;
    if (startIndex > 0) {
      total -= dist.probability(0);
    }

    double cumul = 0;
    for (int i = startIndex; i <= trials; i++) {
      cumul += dist.probability(i) / total;
      x[i] = i;
      y[i] = cumul;
    }

    plot.setColor(color);
    plot.addPoints(x, y, Plot.LINE);
    // plot.addPoints(x, y, Plot2.CIRCLE);
    ImageJUtils.display(title, plot);
  }
}
