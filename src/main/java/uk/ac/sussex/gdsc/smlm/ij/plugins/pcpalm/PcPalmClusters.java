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

package uk.ac.sussex.gdsc.smlm.ij.plugins.pcpalm;

import ij.IJ;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.plugin.frame.Recorder;
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
import java.util.concurrent.atomic.AtomicReference;
import java.util.regex.Pattern;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.util.FastMath;
import uk.ac.sussex.gdsc.core.clustering.Cluster;
import uk.ac.sussex.gdsc.core.clustering.ClusterPoint;
import uk.ac.sussex.gdsc.core.clustering.ClusteringAlgorithm;
import uk.ac.sussex.gdsc.core.clustering.ClusteringEngine;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot;
import uk.ac.sussex.gdsc.core.ij.ImageJPluginLoggerHelper;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.SimpleImageJTrackProgress;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.Plot2;
import uk.ac.sussex.gdsc.core.utils.FileUtils;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationHelper;
import uk.ac.sussex.gdsc.smlm.fitting.BinomialFitter;
import uk.ac.sussex.gdsc.smlm.ij.plugins.HelpUrls;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ParameterUtils;
import uk.ac.sussex.gdsc.smlm.ij.plugins.SmlmUsageTracker;
import uk.ac.sussex.gdsc.smlm.ij.plugins.pcpalm.PcPalmMolecules.MoleculesResults;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.ExtendedPeakResult;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;

/**
 * Find clusters of molecules using a partial centroid-linkage hierarchical clustering algorithm.
 *
 * <p>Points are added to the nearest cluster if they are below the distance threshold to the
 * cluster centroid. The cluster centroid is updated. All points above the cluster distance
 * threshold remain as single molecules.
 *
 * <p>The purpose is to join colocalising molecules into clusters.
 *
 * <p>See Puchner, et al (2013). Counting molecules in single organelles with superresolution
 * microscopy allows tracking of the endosome maturation trajectory. PNAS.
 * doi:10.1073/pnas.1309676110
 */
public class PcPalmClusters implements PlugIn {
  /** The title. */
  static final String TITLE = "PC-PALM Clusters";

  /**
   * The auto save is only set if the histogram data file already exists. The option to save is then
   * a choice of the user. Otherwise default to true.
   */
  private boolean autoSave = true;
  /**
   * The clustering algorithm is only set when not running in file input mode.
   */
  private ClusteringAlgorithm clusteringAlgorithm = ClusteringAlgorithm.PARTICLE_CENTROID_LINKAGE;
  /**
   * The weighted clustering flag is only set if there are weights.
   */
  private boolean weightedClustering;

  private boolean fileInput;

  private int numberOfMolecules;
  private double count;

  /** The plugin settings. */
  private Settings settings;
  /** The results from PC-PALM molecules. */
  private MoleculesResults moleculesResults;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    private static final String[] UNITS = {"pixels^2", "um^2"};

    int runMode;
    double distance;
    ClusteringAlgorithm clusteringAlgorithm;
    int minN;
    int maxN;
    boolean maximumLikelihood;
    boolean showCumulativeHistogram;
    boolean multiThread;
    boolean weightedClustering;
    boolean saveHistogram;
    String histogramFile;
    String noiseFile;
    boolean autoSave;
    boolean calibrateHistogram;
    int frames;
    double area;
    int units;

    Settings() {
      // Set defaults
      distance = 50;
      clusteringAlgorithm = ClusteringAlgorithm.PARTICLE_CENTROID_LINKAGE;
      minN = 1;
      multiThread = true;
      histogramFile = "";
      noiseFile = "";
      autoSave = true;
      frames = 1;
      area = 1;
    }

    Settings(Settings source) {
      runMode = source.runMode;
      distance = source.distance;
      clusteringAlgorithm = source.clusteringAlgorithm;
      minN = source.minN;
      maxN = source.maxN;
      maximumLikelihood = source.maximumLikelihood;
      showCumulativeHistogram = source.showCumulativeHistogram;
      multiThread = source.multiThread;
      weightedClustering = source.weightedClustering;
      saveHistogram = source.saveHistogram;
      histogramFile = source.histogramFile;
      noiseFile = source.noiseFile;
      autoSave = source.autoSave;
      calibrateHistogram = source.calibrateHistogram;
      frames = source.frames;
      area = source.area;
      units = source.units;
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
     * Save the settings. This can be called only once as it saves via a reference.
     */
    void save() {
      lastSettings.set(this);
    }
  }

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
      histogramData = loadHistogram(settings.histogramFile);
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
    if (noiseData != null && subtractNoise(histogramData, noiseData)) {
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
        final String newFilename = FileUtils.replaceExtension(histogramData.filename, ".noise.tsv");
        if (saveHistogram(histogramData, newFilename)) {
          ImageJUtils.log("Saved noise-subtracted histogram to " + newFilename);
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
    final List<Molecule> molecules = analysis.cropToRoi(WindowManager.getCurrentImage());

    if (molecules.size() < 2) {
      error("No results within the crop region");
      return null;
    }

    ImageJUtils.log("Using %d molecules (Density = %s um^-2) @ %s nm", molecules.size(),
        MathUtils.rounded(molecules.size() / analysis.croppedArea),
        MathUtils.rounded(settings.distance));

    final long s1 = System.nanoTime();
    final ClusteringEngine engine =
        new ClusteringEngine(1, clusteringAlgorithm, SimpleImageJTrackProgress.getInstance());
    if (settings.multiThread) {
      engine.setThreadCount(Prefs.getThreads());
    }
    engine.setTracker(SimpleImageJTrackProgress.getInstance());
    IJ.showStatus("Clustering ...");
    final List<Cluster> clusters =
        engine.findClusters(convertToPoint(molecules), settings.distance);
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
    results.setCalibration(CalibrationHelper.create(100, 1, moleculesResults.seconds * 1000));
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

    final HistogramData histogramData = (settings.calibrateHistogram)
        ? new HistogramData(hist, settings.frames, settings.area, Settings.UNITS[settings.units])
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
  private List<ClusterPoint> convertToPoint(List<Molecule> molecules) {
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
  private boolean saveHistogram(HistogramData histogramData) {
    if (!settings.saveHistogram) {
      return false;
    }
    settings.histogramFile = ImageJUtils.getFilename("Histogram_file", settings.histogramFile);
    return saveHistogram(histogramData, settings.histogramFile);
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

    filename = FileUtils.replaceExtension(filename, "tsv");
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
    } catch (final IOException ex) {
      IJ.log("Failed to save histogram to file: " + filename + ". " + ex.getMessage());
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
      gd.addCheckbox("Auto_save noise-subtracted histogram", settings.autoSave);
    }

    // If this is a macro then the dialog will not have Yes or No pressed.
    // Add a checkbox that can be read from the macro arguments by ImageJ.
    final String macroOption = "subtract";
    if (IJ.isMacro()) {
      gd.addCheckbox(macroOption, true);
    }

    gd.addHelp(HelpUrls.getUrl("pc-palm-clusters"));
    gd.showDialog();
    if (!gd.wasOKed()) {
      return null;
    }
    if (allowSave) {
      autoSave = settings.autoSave = gd.getNextBoolean();
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

    settings.noiseFile = ImageJUtils.getFilename("Noise_file", settings.noiseFile);
    if (settings.noiseFile != null) {
      final HistogramData data = loadHistogram(settings.noiseFile);
      // Check the data is calibrated with the same units
      if (data.isCalibrated() && data.units.equalsIgnoreCase(histogramData.units)) {
        return data;
      }
    }
    return null;
  }

  private boolean showDialog() {
    moleculesResults = PcPalmMolecules.getMoleculesResults();
    settings = Settings.load();
    if (moleculesResults.molecules == null || moleculesResults.molecules.size() < 2) {
      ImageJUtils.log(TITLE + " defaulting to File mode");
      fileInput = true;
      // Ensure this gets recorded
      Recorder.recordOption("Method", "File");
    } else {
      final GenericDialog gd = new GenericDialog(TITLE);
      final String[] items = {"Clustering", "File"};

      gd.addMessage("Fit a Binomial distribution to a histogram of cluster sizes.\n \n"
          + "Select the method to generate the histogram:");
      gd.addChoice("Method", items, items[settings.runMode]);
      gd.addHelp(HelpUrls.getUrl("pc-palm-clusters"));
      gd.showDialog();
      if (gd.wasCanceled()) {
        return false;
      }
      settings.runMode = gd.getNextChoiceIndex();
      fileInput = (settings.runMode == 1);
    }

    settings.save();

    if (fileInput) {
      final String newFile = ImageJUtils.getFilename("Histogram_file", settings.histogramFile);
      if (newFile == null) {
        return false;
      }
      settings.histogramFile = newFile;
    }

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);

    // Check if the molecules have weights
    boolean haveWeights = false;
    if (!fileInput) {
      haveWeights = checkForWeights();

      gd.addMessage("Find clusters using centroid-linkage clustering.");

      gd.addNumericField("Distance", settings.distance, 0, 6, "nm");
      final String[] clusteringAlgorithmNames =
          SettingsManager.getNames((Object[]) ClusteringAlgorithm.values());
      gd.addChoice("Algorithm", clusteringAlgorithmNames, settings.clusteringAlgorithm.ordinal());
      gd.addCheckbox("Multi_thread", settings.multiThread);
      if (haveWeights) {
        gd.addCheckbox("Weighted_clustering", settings.weightedClustering);
      }
    }

    gd.addSlider("Min_N", 1, 10, settings.minN);
    gd.addSlider("Max_N", 0, 10, settings.maxN);
    gd.addCheckbox("Show_cumulative_histogram", settings.showCumulativeHistogram);
    gd.addCheckbox("Maximum_likelihood", settings.maximumLikelihood);

    if (!fileInput) {
      gd.addCheckbox("Save_histogram", settings.saveHistogram);
      gd.addMessage("Histogram calibration (optional)");
      gd.addCheckbox("Calibrate_histogram", settings.calibrateHistogram);
      gd.addNumericField("Frames", settings.frames, 0);
      gd.addNumericField("Area", settings.area, 2);
      gd.addChoice("Units", Settings.UNITS, settings.units);
    }

    gd.addHelp(HelpUrls.getUrl("pc-palm-clusters"));
    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    if (!fileInput) {
      settings.distance = gd.getNextNumber();
      clusteringAlgorithm =
          settings.clusteringAlgorithm = ClusteringAlgorithm.values()[gd.getNextChoiceIndex()];
      settings.multiThread = gd.getNextBoolean();
      if (haveWeights) {
        weightedClustering = settings.weightedClustering = gd.getNextBoolean();
      }
    }
    settings.minN = (int) Math.abs(gd.getNextNumber());
    settings.maxN = (int) Math.abs(gd.getNextNumber());
    settings.showCumulativeHistogram = gd.getNextBoolean();
    settings.maximumLikelihood = gd.getNextBoolean();
    if (!fileInput) {
      settings.saveHistogram = gd.getNextBoolean();
      settings.calibrateHistogram = gd.getNextBoolean();
      settings.frames = (int) Math.abs(gd.getNextNumber());
      settings.area = Math.abs(gd.getNextNumber());
      settings.units = gd.getNextChoiceIndex();
    }

    // Check arguments
    try {
      ParameterUtils.isAboveZero("Min N", settings.minN);
      if (!fileInput) {
        ParameterUtils.isAboveZero("Distance", settings.distance);
        ParameterUtils.isAboveZero("Frames", settings.frames);
        ParameterUtils.isAboveZero("Area", settings.area);
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
  private boolean checkForWeights() {
    for (final Molecule m : moleculesResults.molecules) {
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
    final int length = v1.length; // FastMath.max(v1.length, v2.length)
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
    if (settings.showCumulativeHistogram) {
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
    int min = settings.minN;
    final boolean customRange = (settings.minN > 1) || (settings.maxN > 0);
    if (min > countN) {
      min = countN;
    }
    if (settings.maxN > 0 && countN > settings.maxN) {
      countN = settings.maxN;
    }

    ImageJUtils.log("Fitting N from %d to %d%s", min, countN,
        (customRange) ? " (custom-range)" : "");

    // Since varying the N should be done in integer steps do this
    // for n=1,2,3,... until the SS peaks then falls off (is worse then the best
    // score several times in succession)
    final BinomialFitter bf = new BinomialFitter(ImageJPluginLoggerHelper.getLogger(getClass()));
    bf.setMaximumLikelihood(settings.maximumLikelihood);
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

      if (settings.showCumulativeHistogram) {
        addToPlot(n, p, title, plot, new Color((float) n / countN, 0, 1f - (float) n / countN));
      }
    }

    // Add best it in magenta
    if (settings.showCumulativeHistogram && parameters != null) {
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
