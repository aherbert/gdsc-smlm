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
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.text.TextWindow;
import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.concurrent.atomic.AtomicReference;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.MultivariateMatrixFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.exception.TooManyIterationsException;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer.Optimum;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.linear.DiagonalMatrix;
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
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunctionGradient;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.FastMath;
import uk.ac.sussex.gdsc.core.annotation.Nullable;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.utils.FileUtils;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.rng.RandomGeneratorAdapter;
import uk.ac.sussex.gdsc.core.utils.rng.UniformRandomProviders;
import uk.ac.sussex.gdsc.smlm.ij.plugins.HelpUrls;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ParameterUtils;
import uk.ac.sussex.gdsc.smlm.ij.plugins.SmlmUsageTracker;
import uk.ac.sussex.gdsc.smlm.ij.plugins.pcpalm.PcPalmMolecules.MoleculesResults;
import uk.ac.sussex.gdsc.smlm.ij.utils.LoggingOptimiserFunction;
import uk.ac.sussex.gdsc.smlm.math3.optim.nonlinear.scalar.gradient.BfgsOptimizer;

/**
 * Use the PC-PALM protocol to fit correlation curve(s) using the random or clustered model.
 *
 * <p>See Sengupta, et al (2013). Quantifying spatial resolution in point-localisation
 * superresolution images using pair correlation analysis. Nature Protocols 8, pp345-354.
 */
public class PcPalmFitting implements PlugIn {
  /** The title. */
  static final String TITLE = "PC-PALM Fitting";

  private static final String INPUT_FROM_FILE = "Load from file";
  private static final String INPUT_PREVIOUS = "Re-use previous curve";
  private static final String INPUT_ANALYSIS = "Select PC-PALM Analysis results";
  private static final String HEADER_PEAK_DENSITY = "Peak density (um^-2)";
  private static final String HEADER_SPATIAL_DOMAIN = "Spatial domain";

  private RandomModelFunction randomModel;
  private ClusteredModelFunctionGradient clusteredModel;
  private EmulsionModelFunctionGradient emulsionModel;

  private int boundedEvaluations;

  /** The adjusted coefficient of determination (r^2) of the random model. */
  private double randomModelAdjustedR2;
  private boolean valid1;
  private boolean valid2;

  // Used for the results table
  private static AtomicReference<TextWindow> resultsTableRef = new AtomicReference<>();
  private TextWindow resultsTable;

  private boolean doneHeader;

  private int offset;
  private double[][] gr;
  private double peakDensity;
  private boolean spatialDomain;

  // Save the input for the analysis

  /** The latest correlation curve (g(r)). */
  private static AtomicReference<CorrelationCurveResult> latestResult = new AtomicReference<>();

  /** The plugin settings. */
  private Settings settings;

  /**
   * Hold the correlation curve results.
   */
  static class CorrelationCurveResult {
    /** The correlation curve (g(r)). */
    final double[][] gr;
    /** The peak density. */
    final double peakDensity;
    /** The spatial domain. */
    final boolean spatialDomain;

    /**
     * Instantiates a new correlation curve result.
     *
     * @param gr the correlation curve
     * @param peakDensity the peak density
     * @param spatialDomain the spatial domain
     */
    CorrelationCurveResult(double[][] gr, double peakDensity, boolean spatialDomain) {
      this.gr = gr;
      this.peakDensity = peakDensity;
      this.spatialDomain = spatialDomain;
    }
  }

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    String inputOption;
    double correlationDistance;
    double estimatedPrecision;
    double copiedEstimatedPrecision;
    double blinkingRate;
    double copiedBlinkingRate;
    boolean showErrorBars;
    boolean fitClusteredModels;
    boolean saveCorrelationCurve;
    String inputFilename;
    String outputFilename;
    int fitRestarts;
    boolean refitWithGradients;
    double fitAboveEstimatedPrecision;
    double fittingTolerance; // Zero to ignore
    double grProteinThreshold;

    Settings() {
      // Set defaults
      inputOption = "";
      correlationDistance = 800; // nm
      estimatedPrecision = -1;
      copiedEstimatedPrecision = -1;
      blinkingRate = -1;
      copiedBlinkingRate = -1;
      inputFilename = "";
      outputFilename = "";
      fitRestarts = 3;
      grProteinThreshold = 1.5;
    }

    Settings(Settings source) {
      inputOption = source.inputOption;
      correlationDistance = source.correlationDistance;
      estimatedPrecision = source.estimatedPrecision;
      copiedEstimatedPrecision = source.copiedEstimatedPrecision;
      blinkingRate = source.blinkingRate;
      copiedBlinkingRate = source.copiedBlinkingRate;
      showErrorBars = source.showErrorBars;
      fitClusteredModels = source.fitClusteredModels;
      saveCorrelationCurve = source.saveCorrelationCurve;
      inputFilename = source.inputFilename;
      outputFilename = source.outputFilename;
      fitRestarts = source.fitRestarts;
      refitWithGradients = source.refitWithGradients;
      fitAboveEstimatedPrecision = source.fitAboveEstimatedPrecision;
      fittingTolerance = source.fittingTolerance;
      grProteinThreshold = source.grProteinThreshold;
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

  @Override
  public void run(String arg) {
    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    settings = Settings.load();
    if (!showDialog()) {
      return;
    }

    final long start = System.currentTimeMillis();
    header();

    analyse();

    final double seconds = (System.currentTimeMillis() - start) / 1000.0;
    final String msg = TITLE + " complete : " + seconds + "s";
    IJ.showStatus(msg);
    IJ.log(msg);
  }

  private void header() {
    if (!doneHeader) {
      doneHeader = true;
      PcPalmMolecules.logSpacer();
      IJ.log(TITLE);
      PcPalmMolecules.logSpacer();
    }
  }

  private boolean showDialog() {
    // Build a list of results to use in the analysis
    if (!getCorrelationResults()) {
      return false;
    }

    final GenericDialog gd = new GenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("pc-palm-fitting"));
    if (spatialDomain) {
      // Spatial domain results are just combined to a curve
      // Add option to save the results curve
      gd.addMessage("Options:");
      gd.addCheckbox("Save_correlation_curve", settings.saveCorrelationCurve);
      gd.showDialog();
      if (gd.wasCanceled()) {
        return false;
      }
      settings.saveCorrelationCurve = gd.getNextBoolean();
      return true;
    }

    final MoleculesResults moleculesResults = PcPalmMolecules.getMoleculesResults();
    if (settings.estimatedPrecision < 0
        || settings.copiedEstimatedPrecision != moleculesResults.sigmaS) {
      settings.copiedEstimatedPrecision = settings.estimatedPrecision = moleculesResults.sigmaS;
    }
    final double analysisBlinkingRate = PcPalmAnalysis.getBlinkingRate();
    if (settings.blinkingRate < 0 || settings.copiedBlinkingRate != analysisBlinkingRate) {
      settings.copiedBlinkingRate = settings.blinkingRate = analysisBlinkingRate;
    }

    gd.addMessage("Analyse clusters using Pair Correlation.");

    gd.addNumericField("Estimated_precision", settings.estimatedPrecision, 2, 6, "nm");
    gd.addNumericField("Blinking_rate", settings.blinkingRate, 2);
    gd.addCheckbox("Show_error_bars", settings.showErrorBars);
    gd.addSlider("Fit_restarts", 0, 5, settings.fitRestarts);
    gd.addCheckbox("Refit_using_gradients", settings.refitWithGradients);
    gd.addSlider("Fit_above_estimated_precision", 0, 2.5, settings.fitAboveEstimatedPrecision);
    gd.addSlider("Fitting_tolerance", 0, 200, settings.fittingTolerance);
    gd.addSlider("gr_random_threshold", 1, 2.5, settings.grProteinThreshold);
    gd.addCheckbox("Fit_clustered_models", settings.fitClusteredModels);
    gd.addCheckbox("Save_correlation_curve", settings.saveCorrelationCurve);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.estimatedPrecision = gd.getNextNumber();
    settings.blinkingRate = gd.getNextNumber();
    settings.showErrorBars = gd.getNextBoolean();
    settings.fitRestarts = (int) Math.abs(gd.getNextNumber());
    settings.refitWithGradients = gd.getNextBoolean();
    settings.fitAboveEstimatedPrecision = Math.abs(gd.getNextNumber());
    settings.fittingTolerance = Math.abs(gd.getNextNumber());
    settings.grProteinThreshold = gd.getNextNumber();
    settings.fitClusteredModels = gd.getNextBoolean();
    settings.saveCorrelationCurve = gd.getNextBoolean();

    settings.save();

    // Check arguments
    try {
      ParameterUtils.isAbove("Correlation distance", settings.correlationDistance, 1);
      ParameterUtils.isAbove("Estimated precision", settings.estimatedPrecision, 0);
      ParameterUtils.isAbove("Blinking_rate", settings.blinkingRate, 0);
      ParameterUtils.isAbove("gr random threshold", settings.grProteinThreshold, 1);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return true;
  }

  private boolean getCorrelationResults() {
    // Option to:
    // - load a correlation curve
    // - use previous results (if available)
    // - select a set of analysis results (if available)
    String[] options = new String[] {INPUT_FROM_FILE, "", ""};
    int count = 1;
    final CorrelationCurveResult previous = latestResult.get();
    if (previous != null) {
      options[count++] = INPUT_PREVIOUS;
    }
    final List<CorrelationResult> allResults = PcPalmAnalysis.getResults();
    if (!allResults.isEmpty()) {
      options[count++] = INPUT_ANALYSIS;
    }

    options = Arrays.copyOf(options, count);
    final GenericDialog gd = new GenericDialog(TITLE);
    gd.addMessage("Select the source for the correlation curve");
    gd.addChoice("Input", options, settings.inputOption);
    gd.addHelp(HelpUrls.getUrl("pc-palm-fitting"));
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }
    settings.inputOption = gd.getNextChoice();
    settings.save();

    if (settings.inputOption.equals(INPUT_PREVIOUS)) {
      // In the case of a macro the previous results may be null
      if (previous == null) {
        return false;
      }

      gr = previous.gr;
      peakDensity = previous.peakDensity;
      spatialDomain = previous.spatialDomain;
      return true;
    } else if (settings.inputOption.equals(INPUT_FROM_FILE)) {
      return loadCorrelationCurve();
    }

    // Fill the results list with analysis results from PCPALM Analysis
    final ArrayList<CorrelationResult> results = new ArrayList<>();
    if (!selectAnalysisResults(allResults, results)) {
      return false;
    }

    // We have some results. Convert them to the format used for fitting.

    header();
    ImageJUtils.log("Computing combined pair correlation curve (%d datasets)", results.size());

    spatialDomain = results.get(0).spatialDomain;

    // Get average peak density
    peakDensity = 0;

    int size = 0;
    for (final CorrelationResult r : results) {
      peakDensity += r.peakDensity;
      size = Math.max(size, r.gr[0].length);
    }
    peakDensity /= results.size();

    // Combine all selected g(r) curves
    gr = combineCurves(results, size);

    return true;
  }

  private boolean selectAnalysisResults(List<CorrelationResult> allResults,
      ArrayList<CorrelationResult> results) {
    // If no results then fail
    if (allResults.isEmpty()) {
      return false;
    }

    // If only one result then use that
    if (allResults.size() == 1) {
      results.add(allResults.get(0));
      return true;
    }

    // Otherwise build a set of matched analysis results
    try {
      while (selectNextCorrelation(allResults, results)) {
        // All processing done in selectNextCorrelation
      }
    } catch (final Exception ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    // Remove bad results from the dataset.
    final ArrayList<CorrelationResult> newResults = new ArrayList<>(results.size());
    for (int i = 0; i < results.size(); i++) {
      final CorrelationResult r = results.get(i);
      // If the PC-PALM Analysis has been done on too few molecules then the g(r) curve will be bad
      if (r.uniquePoints < 10) {
        header();
        ImageJUtils.log("Excluding dataset ID %d - Too few unique points (%f)", r.id,
            r.uniquePoints);
        continue;
      }
      // If the PC-PALM Analysis has a g(r) curve all below 1 then it is not valid
      final int offset = r.spatialDomain ? 0 : 1;
      double max = 0;
      for (int j = offset; j < r.gr[1].length; j++) {
        if (max < r.gr[1][j]) {
          max = r.gr[1][j];
        }
      }
      if (max < 1) {
        header();
        ImageJUtils.log("Excluding dataset ID %d - g(r) curve is always below 1 (max = %f)", r.id,
            max);
        continue;
      }
      newResults.add(r);
    }

    results = newResults;
    return !results.isEmpty();
  }

  private static boolean selectNextCorrelation(List<CorrelationResult> allResults,
      ArrayList<CorrelationResult> results) {
    final ArrayList<String> titles = buildTitlesList(allResults, results);

    // Show a dialog allowing the user to select an input image
    if (titles.isEmpty()) {
      return false;
    }

    final GenericDialog gd = new GenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("pc-palm-fitting"));
    gd.addMessage(
        "Select the next correlation curve\nFrequency domain curves are identified with *");
    final int n = (results.size() + 1);

    // If in macro mode then we must just use the String input field to allow the macro
    // IJ to return the field values from the macro arguments. Using a Choice input
    // will always return a field value.

    if (IJ.isMacro()) {
      // Use blank default value so bad macro parameters return nothing
      gd.addStringField("R_" + n, "");
    } else {
      gd.addChoice("R_" + n, titles.toArray(new String[0]), "");
    }

    gd.addMessage("Cancel to finish");
    gd.showDialog();
    if (gd.wasCanceled()) {
      return false;
    }

    String title;
    if (IJ.isMacro()) {
      title = gd.getNextString();
    } else {
      title = gd.getNextChoice();
    }

    // Check the correlation exists. If not then exit. This is mainly relevant for Macro mode since
    // the loop will continue otherwise since the titles list is not empty.
    final String[] fields = title.split("\\*?:");
    try {
      final int id = Integer.parseInt(fields[0]);
      for (final CorrelationResult r : allResults) {
        if (r.id == id) {
          results.add(r);
          return true;
        }
      }
    } catch (final NumberFormatException ex) {
      // Ignore
    }
    return false;
  }

  private static ArrayList<String> buildTitlesList(List<CorrelationResult> allResults,
      ArrayList<CorrelationResult> results) {
    // Make all subsequent results match the same nmPerPixel limit
    double nmPerPixel = 0;
    boolean spatialDomain = false;
    boolean filter = false;
    if (!results.isEmpty()) {
      filter = true;
      nmPerPixel = results.get(0).nmPerPixel;
      spatialDomain = results.get(0).spatialDomain;
    }

    final ArrayList<String> titles = new ArrayList<>();
    for (final CorrelationResult r : allResults) {
      if (alreadySelected(results, r)
          || (filter && (r.nmPerPixel != nmPerPixel || r.spatialDomain != spatialDomain))) {
        continue;
      }
      titles.add(String.format("%d%s: %s (%s nm/px)", r.id, (r.spatialDomain) ? "" : "*",
          r.source.getName(), MathUtils.rounded(r.nmPerPixel, 3)));
    }
    return titles;
  }

  private static boolean alreadySelected(ArrayList<CorrelationResult> results,
      CorrelationResult result) {
    for (final CorrelationResult r2 : results) {
      if (result.id == r2.id) {
        return true;
      }
    }
    return false;
  }

  /**
   * Perform the PC Analysis.
   *
   * <p>Spatial domain results can just be combined to an average curve.
   *
   * <p>Frequency domain results can be fit using the g(r) model.
   */
  private void analyse() {
    latestResult.set(new CorrelationCurveResult(gr, peakDensity, spatialDomain));

    String axisTitle;
    if (spatialDomain) {
      offset = 0;
      axisTitle = "molecules/um^2";
    } else {
      // Ignore the r=0 value by starting with an offset if necessary
      offset = (gr[0][0] == 0) ? 1 : 0;
      axisTitle = "g(r)";
    }
    final String title = TITLE + " " + axisTitle;
    final Plot plot = PcPalmAnalysis.plotCorrelation(gr, offset, title, axisTitle, spatialDomain,
        settings.showErrorBars);

    if (spatialDomain) {
      saveCorrelationCurve(gr);
      IJ.log("Created correlation curve from the spatial domain (Plot title = " + title + ")");
      return;
    }

    // -------------
    // Model fitting for g(r) correlation curves
    // -------------
    IJ.log("Fitting g(r) correlation curve from the frequency domain");
    ImageJUtils.log("Average peak density = %s um^-2. Blinking estimate = %s",
        MathUtils.rounded(peakDensity, 4), MathUtils.rounded(settings.blinkingRate, 4));

    resultsTable = createResultsTable();

    // Get the protein density in nm^2.
    final double peakDensityNm2 = peakDensity / 1e6;

    // Use the blinking rate estimate to estimate the density
    // (factors in the over-counting of the same molecules)
    final double proteinDensity = peakDensityNm2 / settings.blinkingRate;

    final ArrayList<double[]> curves = new ArrayList<>();

    // Fit the g(r) curve for r>0 to equation 2
    Color color = Color.red;
    String resultColour = "Red";
    double[] parameters =
        fitRandomModel(gr, settings.estimatedPrecision, proteinDensity, resultColour);
    if (parameters != null) {
      ImageJUtils.log("  Plot %s: Over-counting estimate = %s", randomModel.getName(),
          MathUtils.rounded(peakDensityNm2 / parameters[1], 4));
      ImageJUtils.log("  Plot %s == %s", randomModel.getName(), resultColour);
      plot.setColor(color);
      plot.addPoints(randomModel.getX(), randomModel.value(parameters), Plot.LINE);
      addNonFittedPoints(plot, gr, randomModel, parameters);
      ImageJUtils.display(title, plot);
      if (settings.saveCorrelationCurve) {
        curves.add(extractCurve(gr, randomModel, parameters));
      }
    }

    // Fit the clustered models if the random model fails or if chosen as an option
    if (!valid1 || settings.fitClusteredModels) {
      // Fit the g(r) curve for r>0 to equation 3
      color = Color.blue;
      resultColour = "Blue";
      parameters = fitClusteredModel(gr, settings.estimatedPrecision, proteinDensity, resultColour);

      if (parameters != null) {
        ImageJUtils.log("  Plot %s: Over-counting estimate = %s", clusteredModel.getName(),
            MathUtils.rounded(peakDensityNm2 / parameters[1], 4));
        ImageJUtils.log("  Plot %s == %s", clusteredModel.getName(), resultColour.toString());
        plot.setColor(color);
        plot.addPoints(clusteredModel.getX(), clusteredModel.value(parameters), Plot.LINE);
        addNonFittedPoints(plot, gr, clusteredModel, parameters);
        ImageJUtils.display(title, plot);
        if (settings.saveCorrelationCurve) {
          curves.add(extractCurve(gr, clusteredModel, parameters));
        }
      }

      // Fit to an emulsion model for a distribution confined to circles
      color = Color.magenta;
      resultColour = "Magenta";
      parameters = fitEmulsionModel(gr, settings.estimatedPrecision, proteinDensity, resultColour);

      if (parameters != null) {
        ImageJUtils.log("  Plot %s: Over-counting estimate = %s", emulsionModel.getName(),
            MathUtils.rounded(peakDensityNm2 / parameters[1], 4));
        ImageJUtils.log("  Plot %s == %s", emulsionModel.getName(), resultColour.toString());
        plot.setColor(color);
        plot.addPoints(emulsionModel.getX(), emulsionModel.value(parameters), Plot.LINE);
        addNonFittedPoints(plot, gr, emulsionModel, parameters);
        ImageJUtils.display(title, plot);
        if (settings.saveCorrelationCurve) {
          curves.add(extractCurve(gr, emulsionModel, parameters));
        }
      }
    }

    saveCorrelationCurve(gr, curves.toArray(new double[0][0]));
  }

  private void addNonFittedPoints(Plot plot, double[][] gr, BaseModelFunction model,
      double[] parameters) {
    double[] x = new double[gr[0].length];
    double[] y = new double[x.length];
    int count = 0;
    final double first = randomModel.getX()[0];
    for (int i = offset; i < gr[0].length; i++) {
      // Output points that were not fitted
      if (gr[0][i] < first) {
        x[count] = gr[0][i];
        y[count] = model.evaluate(gr[0][i], parameters);
        count++;
      }
    }
    x = Arrays.copyOf(x, count);
    y = Arrays.copyOf(y, count);
    plot.addPoints(x, y, Plot.CIRCLE);
  }

  private double[] extractCurve(double[][] gr, BaseModelFunction model, double[] parameters) {
    final double[] y = new double[gr[0].length - offset];
    for (int i = offset; i < gr[0].length; i++) {
      y[i - offset] = model.evaluate(gr[0][i], parameters);
    }
    return y;
  }

  private static double[][] combineCurves(ArrayList<CorrelationResult> results, int maxSize) {
    final double[][] gr = new double[3][maxSize];
    final Statistics[] grStats = new Statistics[maxSize];
    for (int i = 0; i < maxSize; i++) {
      grStats[i] = new Statistics();
    }

    for (final CorrelationResult r : results) {
      for (int i = 0; i < r.gr[0].length; i++) {
        gr[0][i] = r.gr[0][i]; // All scales should be the same so over-write is OK

        // Note that sometimes the analysis generates values that are very bad (e.g. if too
        // few points were analysed). Perhaps we should exclude outliers for each distance interval.

        // NaN values can be generated so ignore them
        if (!Double.isNaN(r.gr[1][i])) {
          grStats[i].add(r.gr[1][i]);
        }
      }
    }
    for (int i = 0; i < maxSize; i++) {
      gr[1][i] = grStats[i].getMean();
      gr[2][i] = grStats[i].getStandardError();
    }
    return gr;
  }

  private void saveCorrelationCurve(double[][] gr, double[]... curves) {
    if (!settings.saveCorrelationCurve) {
      return;
    }
    settings.outputFilename =
        ImageJUtils.getFilename("Output_Correlation_File", settings.outputFilename);
    if (settings.outputFilename != null) {
      settings.outputFilename = FileUtils.replaceExtension(settings.outputFilename, "xls");

      try (BufferedWriter output = Files.newBufferedWriter(Paths.get(settings.outputFilename))) {
        writeHeader(output, HEADER_PEAK_DENSITY, Double.toString(peakDensity));
        writeHeader(output, HEADER_SPATIAL_DOMAIN, Boolean.toString(spatialDomain));
        output.write("#r\tg(r)\tS.E.");
        for (int j = 0; j < curves.length; j++) {
          output.write(String.format("\tModel %d", j + 1));
        }
        output.newLine();
        // Ignore the r=0 value by starting with an offset if necessary
        for (int i = offset; i < gr[0].length; i++) {
          output.write(String.format("%f\t%f\t%f", gr[0][i], gr[1][i], gr[2][i]));
          for (int j = 0; j < curves.length; j++) {
            output.write(String.format("\t%f", curves[j][i - offset]));
          }
          output.newLine();
        }
      } catch (final IOException ex) {
        IJ.log("Failed to save correlation curve to file: " + settings.outputFilename + ". "
            + ex.getMessage());
      }
    }
  }

  private static void writeHeader(BufferedWriter output, String header, String value)
      throws IOException {
    output.write("#");
    output.write(header);
    output.write(" = ");
    output.write(value);
    output.newLine();
  }

  /**
   * Load a correlation curve from file. Will set the global gr, peakDensity and spatialDomain
   * variables. If the data fails to be loaded then the method will return false.
   *
   * @return True if loaded
   */
  private boolean loadCorrelationCurve() {
    settings.inputFilename =
        ImageJUtils.getFilename("Input_Correlation_File", settings.inputFilename);
    if (settings.inputFilename == null) {
      return false;
    }

    // Set the analysis variables
    boolean spatialDomainSet = false;
    boolean peakDensitySet = false;

    try (BufferedReader input = Files.newBufferedReader(Paths.get(settings.inputFilename))) {
      String line;
      int count = 0;

      final Pattern pattern = Pattern.compile("#([^=]+) = ([^ ]+)");

      // Read the header
      while ((line = input.readLine()) != null) {
        count++;

        if (line.length() == 0) {
          continue;
        }
        if (line.charAt(0) != '#') {
          // This is the first record
          break;
        }

        // This is a header line. Extract the key-value pair
        final Matcher match = pattern.matcher(line);
        if (match.find()) {
          if (match.group(1).equals(HEADER_SPATIAL_DOMAIN)) {
            // Do not use Boolean.parseBoolean because this will not fail if the field is
            // neither true/false - it only return true for a match to true
            spatialDomainSet = true;
            if (match.group(2).equalsIgnoreCase("true")) {
              spatialDomain = true;
            } else if (match.group(2).equalsIgnoreCase("false")) {
              spatialDomain = false;
            } else {
              // We want to know if the field is not true/false
              spatialDomainSet = false;
            }
          } else if (match.group(1).equals(HEADER_PEAK_DENSITY)) {
            try {
              peakDensity = Double.parseDouble(match.group(2));
              peakDensitySet = true;
            } catch (final NumberFormatException ex) {
              // Ignore this.
            }
          }
        }
      }

      if (!peakDensitySet) {
        IJ.error(TITLE,
            "No valid " + HEADER_PEAK_DENSITY + " record in file " + settings.inputFilename);
        return false;
      }
      if (!spatialDomainSet) {
        IJ.error(TITLE,
            "No valid " + HEADER_SPATIAL_DOMAIN + " record in file " + settings.inputFilename);
        return false;
      }

      // Read the data: gr[0][i], gr[1][i], gr[2][i]
      final ArrayList<double[]> data = new ArrayList<>();
      while (line != null) {
        if (line.length() == 0) {
          continue;
        }
        if (line.charAt(0) == '#') {
          continue;
        }

        // Extract the first 3 fields
        try (Scanner scanner = new Scanner(line)) {
          scanner.useDelimiter("[\t ,]+");

          double radius;
          double gr;
          try {
            radius = scanner.nextDouble();
            gr = scanner.nextDouble();
          } catch (final NoSuchElementException ex) {
            IJ.error(TITLE, "Incorrect fields on line " + count);
            return false;
          }
          // Allow the file to be missing the curve error. This is only used for plotting anyway.
          double error = 0;
          try {
            error = scanner.nextDouble();
          } catch (final NoSuchElementException ex) {
            // Ignore
          }
          data.add(new double[] {radius, gr, error});
        }

        // Read the next line
        line = input.readLine();
        count++;
      }

      if (data.isEmpty()) {
        IJ.error(TITLE, "No data in file " + settings.inputFilename);
        return false;
      }

      gr = new double[3][data.size()];
      for (int i = 0; i < data.size(); i++) {
        final double[] d = data.get(i);
        gr[0][i] = d[0];
        gr[1][i] = d[1];
        gr[2][i] = d[2];
      }
    } catch (final IOException ex) {
      IJ.error(TITLE, "Unable to read from file " + settings.inputFilename);
      return false;
    }

    return true;
  }

  /**
   * Fits the correlation curve with r>0 to the random model using the estimated density and
   * precision. Parameters must be fit within a tolerance of the starting values.
   *
   * @param gr the correlation curve
   * @param sigmaS The estimated precision
   * @param proteinDensity The estimate protein density
   * @param resultColour the result colour
   * @return The fitted parameters [precision, density]
   */
  @Nullable
  private double[] fitRandomModel(double[][] gr, double sigmaS, double proteinDensity,
      String resultColour) {
    final RandomModelFunction function = new RandomModelFunction();
    randomModel = function;

    ImageJUtils.log("Fitting %s: Estimated precision = %f nm, estimated protein density = %g um^-2",
        randomModel.getName(), sigmaS, proteinDensity * 1e6);

    randomModel.setLogging(true);

    for (int i = offset; i < gr[0].length; i++) {
      // Only fit the curve above the estimated resolution (points below it will be subject to
      // error)
      if (gr[0][i] > sigmaS * settings.fitAboveEstimatedPrecision) {
        randomModel.addPoint(gr[0][i], gr[1][i]);
      }
    }
    final LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();

    Optimum optimum;
    try {
      //@formatter:off
      final LeastSquaresProblem problem = new LeastSquaresBuilder()
          .maxEvaluations(Integer.MAX_VALUE)
          .maxIterations(3000)
          .start(new double[] { sigmaS, proteinDensity })
          .target(function.getY())
          .weight(new DiagonalMatrix(function.getWeights()))
          .model(function, function::jacobian)
          .build();
      //@formatter:on

      optimum = optimizer.optimize(problem);
    } catch (final TooManyIterationsException ex) {
      ImageJUtils.log("Failed to fit %s: Too many iterations (%s)", randomModel.getName(),
          ex.getMessage());
      return null;
    } catch (final ConvergenceException ex) {
      ImageJUtils.log("Failed to fit %s: %s", randomModel.getName(), ex.getMessage());
      return null;
    }

    randomModel.setLogging(false);

    final double[] parameters = optimum.getPoint().toArray();
    // Ensure the width is positive
    parameters[0] = Math.abs(parameters[0]);

    final double ss = optimum.getResiduals().dotProduct(optimum.getResiduals());
    final double totalSumSquares = MathUtils.getTotalSumOfSquares(randomModel.getY());
    randomModelAdjustedR2 = MathUtils.getAdjustedCoefficientOfDetermination(ss, totalSumSquares,
        randomModel.size(), parameters.length);

    final double fitSigmaS = parameters[0];
    final double fitProteinDensity = parameters[1];

    // Check the fitted parameters are within tolerance of the initial estimates
    final double e1 = parameterDrift(sigmaS, fitSigmaS);
    final double e2 = parameterDrift(proteinDensity, fitProteinDensity);

    ImageJUtils.log("  %s fit: SS = %f. Adj.R^2 = %f. %d evaluations", randomModel.getName(), ss,
        randomModelAdjustedR2, optimum.getEvaluations());
    ImageJUtils.log("  %s parameters:", randomModel.getName());
    ImageJUtils.log("    Average precision = %s nm (%s%%)", MathUtils.rounded(fitSigmaS, 4),
        MathUtils.rounded(e1, 4));
    ImageJUtils.log("    Average protein density = %s um^-2 (%s%%)",
        MathUtils.rounded(fitProteinDensity * 1e6, 4), MathUtils.rounded(e2, 4));

    valid1 = true;
    if (settings.fittingTolerance > 0
        && (Math.abs(e1) > settings.fittingTolerance || Math.abs(e2) > settings.fittingTolerance)) {
      ImageJUtils.log(
          "  Failed to fit %s within tolerance (%s%%): Average precision = %f nm (%s%%),"
              + " average protein density = %g um^-2 (%s%%)",
          randomModel.getName(), MathUtils.rounded(settings.fittingTolerance, 4), fitSigmaS,
          MathUtils.rounded(e1, 4), fitProteinDensity * 1e6, MathUtils.rounded(e2, 4));
      valid1 = false;
    }

    if (valid1) {
      // ---------
      // TODO - My data does not comply with this criteria.
      // This could be due to the PC-PALM Molecule code limiting the nmPerPixel to fit the images in
      // memory thus removing correlations at small r.
      // It could also be due to the nature of the random simulations being 3D not 2D membranes
      // as per the PC-PALM paper.
      // ---------
      // Evaluate g(r)protein where:
      // g(r)peaks = g(r)protein + g(r)stoch
      // g(r)peaks ~ 1 + g(r)stoch
      // Verify g(r)protein should be <1.5 for all r>0
      final double[] grStoch = randomModel.value(parameters);
      final double[] grPeaks = randomModel.getY();
      final double[] radius = randomModel.getX();

      for (int i = 0; i < grPeaks.length; i++) {
        // Only evaluate above the fitted average precision
        if (radius[i] < fitSigmaS) {
          continue;
        }

        // Note the RandomModelFunction evaluates g(r)stoch + 1
        final double gr_protein_i = grPeaks[i] - (grStoch[i] - 1);

        if (gr_protein_i > settings.grProteinThreshold) {
          // Failed fit
          ImageJUtils.log("  Failed to fit %s: g(r)protein %s > %s @ r=%s", randomModel.getName(),
              MathUtils.rounded(gr_protein_i, 4), MathUtils.rounded(settings.grProteinThreshold, 4),
              MathUtils.rounded(radius[i], 4));
          valid1 = false;
        }
      }
    }

    addResult(randomModel.getName(), resultColour, valid1, fitSigmaS, fitProteinDensity, 0, 0, 0, 0,
        randomModelAdjustedR2);

    return parameters;
  }

  private static double parameterDrift(double start, double end) {
    if (end < start) {
      return -100 * (start - end) / end;
    }
    return 100 * (end - start) / start;
  }

  /**
   * Fits the correlation curve with r>0 to the clustered model using the estimated density and
   * precision. Parameters must be fit within a tolerance of the starting values.
   *
   * @param gr the correlation curve
   * @param sigmaS The estimated precision
   * @param proteinDensity The estimated protein density
   * @param resultColour the result colour
   * @return The fitted parameters [precision, density, clusterRadius, clusterDensity]
   */
  @Nullable
  private double[] fitClusteredModel(double[][] gr, double sigmaS, double proteinDensity,
      String resultColour) {
    final ClusteredModelFunctionGradient function = new ClusteredModelFunctionGradient();
    clusteredModel = function;
    ImageJUtils.log("Fitting %s: Estimated precision = %f nm, estimated protein density = %g um^-2",
        clusteredModel.getName(), sigmaS, proteinDensity * 1e6);

    clusteredModel.setLogging(true);
    for (int i = offset; i < gr[0].length; i++) {
      // Only fit the curve above the estimated resolution (points below it will be subject to
      // error)
      if (gr[0][i] > sigmaS * settings.fitAboveEstimatedPrecision) {
        clusteredModel.addPoint(gr[0][i], gr[1][i]);
      }
    }

    double[] parameters;
    // The model is: sigma, density, range, amplitude
    final double[] initialSolution = new double[] {sigmaS, proteinDensity, sigmaS * 5, 1};

    // Constrain the fitting to be close to the estimated precision (sigmaS) and protein density.
    // LVM fitting does not support constrained fitting so use a bounded optimiser.
    final SumOfSquaresModelFunction clusteredModelMulti =
        new SumOfSquaresModelFunction(clusteredModel);
    final double[] x = clusteredModelMulti.x;

    // Put some bounds around the initial guess. Use the fitting tolerance (in %) if provided.
    final double limit = (settings.fittingTolerance > 0) ? 1 + settings.fittingTolerance / 100 : 2;
    final double[] lB = new double[] {initialSolution[0] / limit, initialSolution[1] / limit, 0, 0};
    // The amplitude and range should not extend beyond the limits of the g(r) curve.
    final double[] uB = new double[] {initialSolution[0] * limit, initialSolution[1] * limit,
        MathUtils.max(x), MathUtils.max(gr[1])};
    ImageJUtils.log("Fitting %s using a bounded search: %s < precision < %s & %s < density < %s",
        clusteredModel.getName(), MathUtils.rounded(lB[0], 4), MathUtils.rounded(uB[0], 4),
        MathUtils.rounded(lB[1] * 1e6, 4), MathUtils.rounded(uB[1] * 1e6, 4));

    final PointValuePair constrainedSolution =
        runBoundedOptimiser(initialSolution, lB, uB, clusteredModelMulti);

    if (constrainedSolution == null) {
      return null;
    }

    parameters = constrainedSolution.getPointRef();
    int evaluations = boundedEvaluations;

    // Refit using a LVM
    if (settings.refitWithGradients) {
      ImageJUtils.log("Re-fitting %s using a gradient optimisation", clusteredModel.getName());
      final LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();
      Optimum lvmSolution;
      try {
        //@formatter:off
        final LeastSquaresProblem problem = new LeastSquaresBuilder()
            .maxEvaluations(Integer.MAX_VALUE)
            .maxIterations(3000)
            .start(parameters)
            .target(function.getY())
            .weight(new DiagonalMatrix(function.getWeights()))
            .model(function, new MultivariateMatrixFunction() {
              @Override
              public double[][] value(double[] point) {
                return function.jacobian(point);
              }} )
            .build();
        //@formatter:on

        lvmSolution = optimizer.optimize(problem);
        evaluations += lvmSolution.getEvaluations();

        final double ss = lvmSolution.getResiduals().dotProduct(lvmSolution.getResiduals());
        if (ss < constrainedSolution.getValue()) {
          ImageJUtils.log("Re-fitting %s improved the SS from %s to %s (-%s%%)",
              clusteredModel.getName(), MathUtils.rounded(constrainedSolution.getValue(), 4),
              MathUtils.rounded(ss, 4), MathUtils.rounded(
                  100 * (constrainedSolution.getValue() - ss) / constrainedSolution.getValue(), 4));
          parameters = lvmSolution.getPoint().toArray();
        }
      } catch (final TooManyIterationsException ex) {
        ImageJUtils.log("Failed to re-fit %s: Too many iterations (%s)", clusteredModel.getName(),
            ex.getMessage());
      } catch (final ConvergenceException ex) {
        ImageJUtils.log("Failed to re-fit %s: %s", clusteredModel.getName(), ex.getMessage());
      }
    }

    clusteredModel.setLogging(false);

    // Ensure the width is positive
    parameters[0] = Math.abs(parameters[0]);

    double ss = 0;
    final double[] obs = clusteredModel.getY();
    final double[] exp = clusteredModel.value(parameters);
    for (int i = 0; i < obs.length; i++) {
      ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);
    }
    final double totalSumSquares = MathUtils.getTotalSumOfSquares(clusteredModel.getY());
    final double adjustedR2 = MathUtils.getAdjustedCoefficientOfDetermination(ss, totalSumSquares,
        clusteredModel.size(), parameters.length);

    final double fitSigmaS = parameters[0];
    final double fitProteinDensity = parameters[1];
    final double domainRadius = parameters[2]; // The radius of the cluster domain
    final double domainDensity = parameters[3]; // The density of the cluster domain

    // This is from the PC-PALM paper. However that paper fits the g(r)protein exponential convolved
    // in 2D with the g(r)PSF. In this method we have just fit the exponential
    final double nCluster =
        2 * domainDensity * Math.PI * domainRadius * domainRadius * fitProteinDensity;

    final double e1 = parameterDrift(sigmaS, fitSigmaS);
    final double e2 = parameterDrift(proteinDensity, fitProteinDensity);

    ImageJUtils.log("  %s fit: SS = %f. Adj.R^2 = %f. %d evaluations", clusteredModel.getName(), ss,
        adjustedR2, evaluations);
    ImageJUtils.log("  %s parameters:", clusteredModel.getName());
    ImageJUtils.log("    Average precision = %s nm (%s%%)", MathUtils.rounded(fitSigmaS, 4),
        MathUtils.rounded(e1, 4));
    ImageJUtils.log("    Average protein density = %s um^-2 (%s%%)",
        MathUtils.rounded(fitProteinDensity * 1e6, 4), MathUtils.rounded(e2, 4));
    ImageJUtils.log("    Domain radius = %s nm", MathUtils.rounded(domainRadius, 4));
    ImageJUtils.log("    Domain density = %s", MathUtils.rounded(domainDensity, 4));
    ImageJUtils.log("    nCluster = %s", MathUtils.rounded(nCluster, 4));

    // Check the fitted parameters are within tolerance of the initial estimates
    valid2 = true;
    if (settings.fittingTolerance > 0
        && (Math.abs(e1) > settings.fittingTolerance || Math.abs(e2) > settings.fittingTolerance)) {
      ImageJUtils.log(
          "  Failed to fit %s within tolerance (%s%%): Average precision = %f nm (%s%%),"
              + " average protein density = %g um^-2 (%s%%)",
          clusteredModel.getName(), MathUtils.rounded(settings.fittingTolerance, 4), fitSigmaS,
          MathUtils.rounded(e1, 4), fitProteinDensity * 1e6, MathUtils.rounded(e2, 4));
      valid2 = false;
    }

    // Check extra parameters. Domain radius should be higher than the precision. Density should be
    // positive
    if (domainRadius < fitSigmaS) {
      ImageJUtils.log(
          "  Failed to fit %s: Domain radius is smaller than the average precision (%s < %s)",
          clusteredModel.getName(), MathUtils.rounded(domainRadius, 4),
          MathUtils.rounded(fitSigmaS, 4));
      valid2 = false;
    }
    if (domainDensity < 0) {
      ImageJUtils.log("  Failed to fit %s: Domain density is negative (%s)",
          clusteredModel.getName(), MathUtils.rounded(domainDensity, 4));
      valid2 = false;
    }

    if (adjustedR2 <= randomModelAdjustedR2) {
      ImageJUtils.log("  Failed to fit %s - Adjusted r^2 has decreased %s%%",
          clusteredModel.getName(), MathUtils
              .rounded((100 * (randomModelAdjustedR2 - adjustedR2) / randomModelAdjustedR2), 4));
      valid2 = false;
    }

    addResult(clusteredModel.getName(), resultColour, valid2, fitSigmaS, fitProteinDensity,
        domainRadius, domainDensity, nCluster, -1, adjustedR2);

    return parameters;
  }

  private PointValuePair runBoundedOptimiser(double[] initialSolution, double[] lowerB,
      double[] upperB, SumOfSquaresModelFunction function) {
    // Create the functions to optimise
    final ObjectiveFunction objective =
        new ObjectiveFunction(new SumOfSquaresMultivariateFunction(function));
    final ObjectiveFunctionGradient gradient =
        new ObjectiveFunctionGradient(new SumOfSquaresMultivariateVectorFunction(function));

    final boolean debug = false;

    // Try a BFGS optimiser since this will produce a deterministic solution and can respect bounds.
    PointValuePair optimum = null;
    boundedEvaluations = 0;
    final MaxEval maxEvaluations = new MaxEval(2000);
    MultivariateOptimizer opt = null;
    for (int iteration = 0; iteration <= settings.fitRestarts; iteration++) {
      try {
        opt = new BfgsOptimizer();
        final double relativeThreshold = 1e-6;

        // Configure maximum step length for each dimension using the bounds
        final double[] stepLength = new double[lowerB.length];
        for (int i = 0; i < stepLength.length; i++) {
          stepLength[i] = (upperB[i] - lowerB[i]) * 0.3333333;
        }

        // The GoalType is always minimise so no need to pass this in
        optimum = opt.optimize(maxEvaluations, gradient, objective,
            new InitialGuess((optimum == null) ? initialSolution : optimum.getPointRef()),
            new SimpleBounds(lowerB, upperB),
            new BfgsOptimizer.GradientTolerance(relativeThreshold),
            new BfgsOptimizer.StepLength(stepLength));
        if (debug) {
          System.out.printf("BFGS Iter %d = %g (%d)\n", iteration, optimum.getValue(),
              opt.getEvaluations());
        }
      } catch (final RuntimeException ex) {
        break; // No need to restart
      } finally {
        if (opt != null) {
          boundedEvaluations += opt.getEvaluations();
        }
      }
    }

    // Try a CMAES optimiser which is non-deterministic. To overcome this we perform restarts.

    // CMAESOptimiser based on Matlab code:
    // https://www.lri.fr/~hansen/cmaes.m
    // Take the defaults from the Matlab documentation
    final double stopFitness = 0;
    final boolean isActiveCma = true;
    final int diagonalOnly = 0;
    final int checkFeasableCount = 1;
    final RandomGenerator random = new RandomGeneratorAdapter(UniformRandomProviders.create());
    final boolean generateStatistics = false;
    final ConvergenceChecker<PointValuePair> checker = new SimpleValueChecker(1e-6, 1e-10);
    // The sigma determines the search range for the variables. It should be 1/3 of the initial
    // search region.
    final double[] range = new double[lowerB.length];
    for (int i = 0; i < lowerB.length; i++) {
      range[i] = (upperB[i] - lowerB[i]) / 3;
    }
    final OptimizationData sigma = new CMAESOptimizer.Sigma(range);
    final OptimizationData popSize = new CMAESOptimizer.PopulationSize(
        (int) (4 + Math.floor(3 * Math.log(initialSolution.length))));
    final SimpleBounds bounds = new SimpleBounds(lowerB, upperB);

    opt = new CMAESOptimizer(maxEvaluations.getMaxEval(), stopFitness, isActiveCma, diagonalOnly,
        checkFeasableCount, random, generateStatistics, checker);
    // Restart the optimiser several times and take the best answer.
    for (int iteration = 0; iteration <= settings.fitRestarts; iteration++) {
      try {
        // Start from the initial solution
        final PointValuePair constrainedSolution = opt.optimize(new InitialGuess(initialSolution),
            objective, GoalType.MINIMIZE, bounds, sigma, popSize, maxEvaluations);
        if (debug) {
          System.out.printf("CMAES Iter %d initial = %g (%d)\n", iteration,
              constrainedSolution.getValue(), opt.getEvaluations());
        }
        boundedEvaluations += opt.getEvaluations();
        if (optimum == null || constrainedSolution.getValue() < optimum.getValue()) {
          optimum = constrainedSolution;
        }
      } catch (final TooManyEvaluationsException | TooManyIterationsException ex) {
        // Ignore
      } finally {
        boundedEvaluations += maxEvaluations.getMaxEval();
      }
      if (optimum == null) {
        continue;
      }
      try {
        // Also restart from the current optimum
        final PointValuePair constrainedSolution =
            opt.optimize(new InitialGuess(optimum.getPointRef()), objective, GoalType.MINIMIZE,
                bounds, sigma, popSize, maxEvaluations);
        if (debug) {
          System.out.printf("CMAES Iter %d restart = %g (%d)\n", iteration,
              constrainedSolution.getValue(), opt.getEvaluations());
        }
        if (constrainedSolution.getValue() < optimum.getValue()) {
          optimum = constrainedSolution;
        }
      } catch (final TooManyEvaluationsException | TooManyIterationsException ex) {
        // Ignore
      } finally {
        boundedEvaluations += maxEvaluations.getMaxEval();
      }
    }
    return optimum;
  }

  /**
   * Fits the correlation curve with r>0 to the clustered model using the estimated density and
   * precision. Parameters must be fit within a tolerance of the starting values.
   *
   * @param gr the correlation curve
   * @param sigmaS The estimated precision
   * @param proteinDensity The estimated protein density
   * @param resultColour the result colour
   * @return The fitted parameters [precision, density, clusterRadius, clusterDensity]
   */
  private double[] fitEmulsionModel(double[][] gr, double sigmaS, double proteinDensity,
      String resultColour) {
    final EmulsionModelFunctionGradient function = new EmulsionModelFunctionGradient();
    emulsionModel = function;
    ImageJUtils.log("Fitting %s: Estimated precision = %f nm, estimated protein density = %g um^-2",
        emulsionModel.getName(), sigmaS, proteinDensity * 1e6);

    emulsionModel.setLogging(true);
    for (int i = offset; i < gr[0].length; i++) {
      // Only fit the curve above the estimated resolution (points below it will be subject to
      // error)
      if (gr[0][i] > sigmaS * settings.fitAboveEstimatedPrecision) {
        emulsionModel.addPoint(gr[0][i], gr[1][i]);
      }
    }

    // The model is: sigma, density, range, amplitude, alpha
    final double[] initialSolution =
        new double[] {sigmaS, proteinDensity, sigmaS * 5, 1, sigmaS * 5};

    // Constrain the fitting to be close to the estimated precision (sigmaS) and protein density.
    // LVM fitting does not support constrained fitting so use a bounded optimiser.
    final SumOfSquaresModelFunction emulsionModelMulti =
        new SumOfSquaresModelFunction(emulsionModel);
    final double[] x = emulsionModelMulti.x;
    final double[] y = emulsionModelMulti.y;

    // Range should be equal to the first time the g(r) curve crosses 1
    for (int i = 0; i < x.length; i++) {
      if (y[i] < 1) {
        initialSolution[4] = initialSolution[2] = (i > 0) ? (x[i - 1] + x[i]) * 0.5 : x[i];
        break;
      }
    }

    // Put some bounds around the initial guess. Use the fitting tolerance (in %) if provided.
    final double limit = (settings.fittingTolerance > 0) ? 1 + settings.fittingTolerance / 100 : 2;
    final double[] lB =
        new double[] {initialSolution[0] / limit, initialSolution[1] / limit, 0, 0, 0};
    // The amplitude and range should not extend beyond the limits of the g(r) curve.
    // TODO - Find out the expected range for the alpha parameter.
    final double[] uB = new double[] {initialSolution[0] * limit, initialSolution[1] * limit,
        MathUtils.max(x), MathUtils.max(gr[1]), MathUtils.max(x) * 2};
    ImageJUtils.log("Fitting %s using a bounded search: %s < precision < %s & %s < density < %s",
        emulsionModel.getName(), MathUtils.rounded(lB[0], 4), MathUtils.rounded(uB[0], 4),
        MathUtils.rounded(lB[1] * 1e6, 4), MathUtils.rounded(uB[1] * 1e6, 4));

    final PointValuePair constrainedSolution =
        runBoundedOptimiser(initialSolution, lB, uB, emulsionModelMulti);

    if (constrainedSolution == null) {
      return null;
    }

    double[] parameters = constrainedSolution.getPointRef();
    int evaluations = boundedEvaluations;

    // Refit using a LVM
    if (settings.refitWithGradients) {
      ImageJUtils.log("Re-fitting %s using a gradient optimisation", emulsionModel.getName());
      final LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();
      Optimum lvmSolution;
      try {
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
        evaluations += lvmSolution.getEvaluations();

        final double ss = lvmSolution.getResiduals().dotProduct(lvmSolution.getResiduals());
        if (ss < constrainedSolution.getValue()) {
          ImageJUtils.log("Re-fitting %s improved the SS from %s to %s (-%s%%)",
              emulsionModel.getName(), MathUtils.rounded(constrainedSolution.getValue(), 4),
              MathUtils.rounded(ss, 4), MathUtils.rounded(
                  100 * (constrainedSolution.getValue() - ss) / constrainedSolution.getValue(), 4));
          parameters = lvmSolution.getPoint().toArray();
        }
      } catch (final TooManyIterationsException ex) {
        ImageJUtils.log("Failed to re-fit %s: Too many iterations (%s)", emulsionModel.getName(),
            ex.getMessage());
      } catch (final ConvergenceException ex) {
        ImageJUtils.log("Failed to re-fit %s: %s", emulsionModel.getName(), ex.getMessage());
      }
    }

    emulsionModel.setLogging(false);

    // Ensure the width is positive
    parameters[0] = Math.abs(parameters[0]);

    double ss = 0;
    final double[] obs = emulsionModel.getY();
    final double[] exp = emulsionModel.value(parameters);
    for (int i = 0; i < obs.length; i++) {
      ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);
    }
    final double totalSumSquares = MathUtils.getTotalSumOfSquares(clusteredModel.getY());
    final double adjustedR2 = MathUtils.getAdjustedCoefficientOfDetermination(ss, totalSumSquares,
        emulsionModel.size(), parameters.length);

    final double fitSigmaS = parameters[0];
    final double fitProteinDensity = parameters[1];
    final double domainRadius = parameters[2]; // The radius of the cluster domain
    final double amplitutde = parameters[3];
    final double coherence = parameters[4]; // The coherence length between circles

    final double e1 = parameterDrift(sigmaS, fitSigmaS);
    final double e2 = parameterDrift(proteinDensity, fitProteinDensity);

    ImageJUtils.log("  %s fit: SS = %f. Adj.R^2 = %f. %d evaluations", emulsionModel.getName(), ss,
        adjustedR2, evaluations);
    ImageJUtils.log("  %s parameters:", emulsionModel.getName());
    ImageJUtils.log("    Average precision = %s nm (%s%%)", MathUtils.rounded(fitSigmaS, 4),
        MathUtils.rounded(e1, 4));
    ImageJUtils.log("    Average protein density = %s um^-2 (%s%%)",
        MathUtils.rounded(fitProteinDensity * 1e6, 4), MathUtils.rounded(e2, 4));
    ImageJUtils.log("    Domain radius = %s nm", MathUtils.rounded(domainRadius, 4));
    ImageJUtils.log("    Domain density = %s", MathUtils.rounded(amplitutde, 4));
    ImageJUtils.log("    Domain coherence = %s", MathUtils.rounded(coherence, 4));

    // Check the fitted parameters are within tolerance of the initial estimates
    valid2 = true;
    if (settings.fittingTolerance > 0
        && (Math.abs(e1) > settings.fittingTolerance || Math.abs(e2) > settings.fittingTolerance)) {
      ImageJUtils.log(
          "  Failed to fit %s within tolerance (%s%%): Average precision = %f nm (%s%%),"
              + " average protein density = %g um^-2 (%s%%)",
          emulsionModel.getName(), MathUtils.rounded(settings.fittingTolerance, 4), fitSigmaS,
          MathUtils.rounded(e1, 4), fitProteinDensity * 1e6, MathUtils.rounded(e2, 4));
      valid2 = false;
    }

    // Check extra parameters. Domain radius should be higher than the precision. Density should be
    // positive
    if (domainRadius < fitSigmaS) {
      ImageJUtils.log(
          "  Failed to fit %s: Domain radius is smaller than the average precision (%s < %s)",
          emulsionModel.getName(), MathUtils.rounded(domainRadius, 4),
          MathUtils.rounded(fitSigmaS, 4));
      valid2 = false;
    }
    if (amplitutde < 0) {
      ImageJUtils.log("  Failed to fit %s: Domain density is negative (%s)",
          emulsionModel.getName(), MathUtils.rounded(amplitutde, 4));
      valid2 = false;
    }

    if (adjustedR2 <= randomModelAdjustedR2) {
      ImageJUtils.log("  Failed to fit %s - Adjusted r^2 has decreased %s%%",
          clusteredModel.getName(), MathUtils
              .rounded((100 * (randomModelAdjustedR2 - adjustedR2) / randomModelAdjustedR2), 4));
      valid2 = false;
    }

    addResult(emulsionModel.getName(), resultColour, valid2, fitSigmaS, fitProteinDensity,
        domainRadius, amplitutde, -1, coherence, adjustedR2);

    return parameters;
  }

  /**
   * Abstract base model function class for common functionality.
   */
  private abstract static class BaseModelFunction extends LoggingOptimiserFunction {
    BaseModelFunction(String name) {
      super(name);
    }

    /**
     * Evaluate the correlation function.
     *
     * @param radius The correlation radius
     * @param parameters The parameters
     * @return the value
     */
    abstract double evaluate(double radius, final double[] parameters);

    /**
     * Evaluate the jacobian of the correlation function for all data points (see
     * {@link #addData(double[], double[])}).
     *
     * @param parameters the parameters
     * @return The jacobian
     */
    abstract double[][] jacobian(double[] parameters);

    /**
     * Get the value of the function for all data points corresponding to the last call to
     * {@link #jacobian(double[])}.
     *
     * @return The corresponding value
     */
    abstract double[] getValue();
  }

  /**
   * Allow optimisation using Apache Commons Math 3 Gradient Optimiser.
   *
   * <p>g(r)peaks = g(r)stoch + 1
   *
   * <p>where
   *
   * <p>g(r)stoch = (1/4*pi*s^2*p) * exp(-r^2/4s^2)
   *
   * <p>s = average single molecule positional uncertainty (precision)
   *
   * <p>p = average protein density
   */
  private static class RandomModelFunction extends BaseModelFunction
      implements MultivariateVectorFunction {
    double[] lastValue;

    RandomModelFunction() {
      super("Random Model");
    }

    @Override
    double[] getValue() {
      return lastValue;
    }

    // Adapted from http://commons.apache.org/proper/commons-math/userguide/optimization.html
    // Use the deprecated API since the new one is not yet documented.

    @Override
    double[][] jacobian(double[] variables) {
      // Compute the gradients using calculus differentiation
      final double sigma = variables[0];
      final double density = variables[1];
      final double[][] jacobian = new double[x.size()][2];
      lastValue = new double[x.size()];

      for (int i = 0; i < jacobian.length; ++i) {
        final double r = this.x.get(i);

        final double a = 1.0 / (4 * Math.PI * density * sigma * sigma);
        final double b = -r * r / (4 * sigma * sigma);
        final double c = FastMath.exp(b);

        // value = a * c
        lastValue[i] = a * c;

        // Differentiate with respect to sigma:
        // value' = a' * c + a * c' [ Product rule ]
        // c = FastMath.exp(b)
        // c' = b' * FastMath.exp(b) [ Chain rule ]
        // value' = a' * c + a * b' * c
        jacobian[i][0] = (-2 * a / sigma) * c + a * (-2 * b / sigma) * c;

        // Differentiate with respect to density:
        // c' = 0 since density does not feature in c
        // => value' = a' * c
        jacobian[i][1] = (-a / density) * c;
      }

      return jacobian;
    }

    @SuppressWarnings("unused")
    private double[][] jacobian2(double[] variables) {
      // Compute the gradients using numerical differentiation
      final double sigma = variables[0];
      final double density = variables[1];
      final double[][] jacobian = new double[x.size()][2];
      lastValue = new double[x.size()];

      final double delta = 0.001;
      final double[][] d = new double[variables.length][variables.length];
      for (int i = 0; i < variables.length; i++) {
        // Should the delta be changed for each parameter?
        d[i][i] = delta * Math.abs(variables[i]);
      }
      for (int i = 0; i < jacobian.length; ++i) {
        final double r = this.x.get(i);
        final double value = lastValue[i] = evaluate(r, sigma, density);
        for (int j = 0; j < variables.length; j++) {
          final double value2 = evaluate(r, sigma + d[0][j], density + d[1][j]);
          jacobian[i][j] = (value2 - value) / d[j][j];
        }
      }
      return jacobian;
    }

    /**
     * Evaluate the correlation function.
     *
     * @param radius The correlation radius
     * @param sigma Average precision
     * @param density Average protein density
     * @return the value
     */
    double evaluate(double radius, final double sigma, final double density) {
      return (1.0 / (4 * Math.PI * density * sigma * sigma))
          * FastMath.exp(-radius * radius / (4 * sigma * sigma)) + 1;
    }

    /**
     * Evaluate the correlation function.
     *
     * @param radius The correlation radius
     * @param parameters The parameters
     * @return the value
     */
    @Override
    double evaluate(double radius, final double[] parameters) {
      return evaluate(radius, parameters[0], parameters[1]);
    }

    @Override
    public double[] value(double[] variables) {
      increment();
      final double[] values = new double[x.size()];
      for (int i = 0; i < values.length; i++) {
        values[i] = evaluate(x.get(i), variables[0], variables[1]);
      }
      return values;
    }
  }

  /**
   * Base implementation of the PC-PALM clustered model. This is used to fit g(r) curves of membrane
   * proteins which appear to be distributed as per a fluctuations model.
   *
   * <p>g(r)peaks = g(r)stoch + g(r)protein
   *
   * <p>where
   *
   * <p>g(r)stoch = (1/4*pi*s^2*p) * exp(-r^2/4s^2)
   *
   * <p>s = average single molecule positional uncertainty (precision)<br> p = average protein
   * density
   *
   * <p>g(r)protein = (A*exp(-r/l)+1) conv g(r)PSF
   *
   * <p>A = proportional to density of proteins in the cluster<br> l = proportional to length of the
   * cluster<br> conv = a convolution operation
   *
   * <p>g(r)PSF = (1/4*pi*s^2) * exp(-r^2/4s^2)
   *
   * <p>Note: The clustered model described in the Veatch PLoS One paper models g(r)protein using
   * the exponential directly, i.e. there is no convolution !!!
   */
  private abstract static class ClusteredModelFunction extends BaseModelFunction {
    double[] lastValue;

    ClusteredModelFunction() {
      super("Clustered Model");
    }

    @Override
    public double[] getValue() {
      return lastValue;
    }

    // Adapted from http://commons.apache.org/proper/commons-math/userguide/optimization.html
    // Use the deprecated API since the new one is not yet documented.

    @Override
    public double[][] jacobian(double[] variables) {
      // Compute the gradients using calculus differentiation
      final double sigma = variables[0];
      final double density = variables[1];
      final double range = variables[2];
      final double amplitude = variables[3];
      final double[][] jacobian = new double[x.size()][variables.length];
      lastValue = new double[x.size()];

      for (int i = 0; i < jacobian.length; ++i) {
        final double r = this.x.get(i);

        final double a = 1.0 / (4 * Math.PI * density * sigma * sigma);
        final double b = -r * r / (4 * sigma * sigma);
        final double c = FastMath.exp(b);

        final double d = -r / range;
        final double e = FastMath.exp(d);

        // value = a * c +
        // amplitude * e + 1
        lastValue[i] = a * c + amplitude * e + 1;

        // Differentiate with respect to sigma:
        // value' = a' * c + a * c' [ Product rule ]
        // c = FastMath.exp(b)
        // c' = b' * FastMath.exp(b) [ Chain rule ]
        // value' = a' * c + a * b' * c
        jacobian[i][0] = (-2 * a / sigma) * c + a * (-2 * b / sigma) * c;

        // Differentiate with respect to density:
        // c' = 0 since density does not feature in c
        // => value' = a' * c
        jacobian[i][1] = (-a / density) * c;

        // Differentiate with respect to range:
        // value' = amplitude * e'
        // e = FastMath.exp(d)
        // e' = d' * FastMath.exp(d) [ Chain rule ]
        jacobian[i][2] = amplitude * (-1 * d / range) * e;

        // Differentiate with respect to amplitude:
        jacobian[i][3] = e;
      }

      return jacobian;
    }

    @SuppressWarnings("unused")
    double[][] jacobian2(double[] variables) {
      // Compute the gradients using numerical differentiation
      final double sigma = variables[0];
      final double density = variables[1];
      final double range = variables[2];
      final double amplitude = variables[3];
      final double[][] jacobian = new double[x.size()][variables.length];
      lastValue = new double[x.size()];

      final double delta = 0.001;
      final double[][] d = new double[variables.length][variables.length];
      for (int i = 0; i < variables.length; i++) {
        // Should the delta be changed for each parameter?
        d[i][i] = delta * Math.abs(variables[i]);
      }
      for (int i = 0; i < jacobian.length; ++i) {
        final double r = this.x.get(i);
        final double value = lastValue[i] = evaluate(r, sigma, density, range, amplitude);
        for (int j = 0; j < variables.length; j++) {
          final double value2 =
              evaluate(r, sigma + d[0][j], density + d[1][j], range + d[2][j], amplitude + d[3][j]);
          jacobian[i][j] = (value2 - value) / d[j][j];
        }
      }
      return jacobian;
    }

    /**
     * Evaluate the correlation function.
     *
     * @param radius The correlation radius
     * @param sigma Average precision
     * @param density Average protein density
     * @param range Range of the cluster
     * @param amplitude Amplitude of the cluster
     * @return the value
     */
    public double evaluate(double radius, final double sigma, final double density,
        final double range, final double amplitude) {
      final double gr_stoch = (1.0 / (4 * Math.PI * density * sigma * sigma))
          * FastMath.exp(-radius * radius / (4 * sigma * sigma));
      final double gr_protein = amplitude * FastMath.exp(-radius / range) + 1;
      return gr_stoch + gr_protein;
    }

    /**
     * Evaluate the correlation function.
     *
     * @param radius The correlation radius
     * @param parameters The parameters
     * @return the value
     */
    @Override
    public double evaluate(double radius, final double[] parameters) {
      return evaluate(radius, parameters[0], parameters[1], parameters[2], parameters[3]);
    }
  }

  /**
   * Allow optimisation using Apache Commons Math 3 Gradient Optimiser.
   */
  private static class ClusteredModelFunctionGradient extends ClusteredModelFunction
      implements MultivariateVectorFunction {
    // Adapted from http://commons.apache.org/proper/commons-math/userguide/optimization.html
    // Use the deprecated API since the new one is not yet documented.

    @Override
    public double[] value(double[] variables) {
      increment();
      final double[] values = new double[x.size()];
      for (int i = 0; i < values.length; i++) {
        values[i] = evaluate(x.get(i), variables[0], variables[1], variables[2], variables[3]);
      }
      return values;
    }
  }

  private static class SumOfSquaresModelFunction {
    BaseModelFunction fun;
    double[] x;
    double[] y;

    // Cache the value
    double[] lastParameters;
    double lastSumSq;

    SumOfSquaresModelFunction(BaseModelFunction fun) {
      this.fun = fun;
      x = fun.getX();
      y = fun.getY();
    }

    double evaluate(double[] parameters) {
      if (sameVariables(parameters)) {
        return lastSumSq;
      }

      lastParameters = null;

      double ss = 0;
      for (int i = x.length; i-- > 0;) {
        final double dx = fun.evaluate(x[i], parameters) - y[i];
        ss += dx * dx;
      }
      return ss;
    }

    /**
     * Check if the variable match those last used for computation of the value.
     *
     * @param parameters the parameters
     * @return True if the variables are the same
     */
    private boolean sameVariables(double[] parameters) {
      if (lastParameters != null) {
        for (int i = 0; i < parameters.length; i++) {
          if (parameters[i] != lastParameters[i]) {
            return false;
          }
        }
        return true;
      }
      return false;
    }

    /**
     * Compute the gradient.
     *
     * @param parameters the parameters
     * @return the gradient
     */
    double[] gradient(double[] parameters) {
      // We can compute the jacobian for all the functions.
      // To get the gradient for the SS we need:
      // f(x) = (g(x) - y)^2
      // f'(x) = 2 * (g(x) - y) * g'(x)

      final double[][] jacobian = fun.jacobian(parameters);
      final double[] gx = fun.getValue();
      lastSumSq = 0;
      lastParameters = parameters.clone();

      final double[] gradient = new double[parameters.length];
      for (int i = 0; i < x.length; i++) {
        final double dx = gx[i] - y[i];
        lastSumSq += dx * dx;
        final double twodx = 2 * dx;
        for (int j = 0; j < gradient.length; j++) {
          final double g1 = twodx * jacobian[i][j];
          gradient[j] += g1;
        }
      }
      return gradient;
    }
  }

  private static class SumOfSquaresMultivariateFunction implements MultivariateFunction {
    SumOfSquaresModelFunction function;

    public SumOfSquaresMultivariateFunction(SumOfSquaresModelFunction function) {
      this.function = function;
    }

    @Override
    public double value(double[] point) {
      return function.evaluate(point);
    }
  }

  private static class SumOfSquaresMultivariateVectorFunction
      implements MultivariateVectorFunction {
    SumOfSquaresModelFunction function;

    public SumOfSquaresMultivariateVectorFunction(SumOfSquaresModelFunction function) {
      this.function = function;
    }

    @Override
    public double[] value(double[] point) {
      return function.gradient(point);
    }
  }

  /**
   * Base implementation of the emulsion clustered model. This model assumes a random distribution
   * of non-overlapping circles in 2D. The molecules can be located at any position within the
   * circles.
   *
   * <p>g(r)peaks = g(r)stoch + g(r)protein
   *
   * <p>where
   *
   * <p>g(r)stoch = (1/4*pi*s^2*p) * exp(-r^2/4s^2)
   *
   * <p>s = average single molecule positional uncertainty (precision)<br> p = average protein
   * density
   *
   * <p>g(r)protein = (A*exp(-r/alpha)*cos(pi*r/(2*r0))+1)
   *
   * <p>A = proportional to density of proteins in the cluster<br> alpha = measure of the coherence
   * length between circles<br> r0 = Average circle radius
   *
   * <p>Note: Described in figure 3 of Veatch, et al (2012) Plos One, e31457
   */
  private abstract static class EmulsionModelFunction extends BaseModelFunction {
    double[] lastValue;

    public EmulsionModelFunction() {
      super("Emulsion Clustered Model");
    }

    @Override
    public double[] getValue() {
      return lastValue;
    }

    // Adapted from http://commons.apache.org/proper/commons-math/userguide/optimization.html
    // Use the deprecated API since the new one is not yet documented.

    @Override
    public double[][] jacobian(double[] variables) {
      // Compute the gradients using calculus differentiation
      final double sigma = variables[0];
      final double density = variables[1];
      final double range = variables[2];
      final double amplitude = variables[3];
      final double alpha = variables[4];
      final double[][] jacobian = new double[x.size()][variables.length];
      lastValue = new double[x.size()];

      for (int i = 0; i < jacobian.length; ++i) {
        final double r = this.x.get(i);

        final double a = 1.0 / (4 * Math.PI * density * sigma * sigma);
        final double b = -r * r / (4 * sigma * sigma);
        final double c = FastMath.exp(b);

        final double d = -r / alpha;
        final double e = FastMath.exp(d);
        final double f = 0.5 * Math.PI * r / range;
        final double g = Math.cos(f);

        // value = a * c +
        // amplitude * e * g + 1
        lastValue[i] = a * c + amplitude * e * g + 1;

        // Differentiate with respect to sigma:
        // value' = a' * c + a * c' [ Product rule ]
        // c = FastMath.exp(b)
        // c' = b' * FastMath.exp(b) [ Chain rule ]
        // value' = a' * c + a * b' * c
        jacobian[i][0] = (-2 * a / sigma) * c + a * (-2 * b / sigma) * c;

        // Differentiate with respect to density:
        // c' = 0 since density does not feature in c
        // => value' = a' * c
        jacobian[i][1] = (-a / density) * c;

        // Differentiate with respect to range:
        // value' = amplitude * e * g'
        // g = Math.cos(f)
        // g' = f' * -Math.sin(f) [ Chain rule ]
        jacobian[i][2] = amplitude * e * (f / range) * Math.sin(f);

        // Differentiate with respect to amplitude:
        jacobian[i][3] = e * g;

        // Differentiate with respect to alpha:
        // value' = amplitude * e' * g
        // e = FastMath.exp(d)
        // e' = d' * FastMath.exp(d) [ Chain rule ]
        // e' = d' * e
        jacobian[i][4] = amplitude * (-1 * d / alpha) * e * g;
      }

      return jacobian;
    }

    @SuppressWarnings("unused")
    double[][] jacobian2(double[] variables) {
      // Compute the gradients using numerical differentiation
      final double sigma = variables[0];
      final double density = variables[1];
      final double range = variables[2];
      final double amplitude = variables[3];
      final double alpha = variables[4];
      final double[][] jacobian = new double[x.size()][variables.length];
      lastValue = new double[x.size()];

      final double delta = 0.001;
      final double[][] d = new double[variables.length][variables.length];
      for (int i = 0; i < variables.length; i++) {
        // Should the delta be changed for each parameter?
        d[i][i] = delta * Math.abs(variables[i]);
      }
      for (int i = 0; i < jacobian.length; ++i) {
        final double r = this.x.get(i);
        final double value = lastValue[i] = evaluate(r, sigma, density, range, amplitude, alpha);
        for (int j = 0; j < variables.length; j++) {
          final double value2 = evaluate(r, sigma + d[0][j], density + d[1][j], range + d[2][j],
              amplitude + d[3][j], alpha + d[4][j]);
          jacobian[i][j] = (value2 - value) / d[j][j];
        }
      }
      return jacobian;
    }

    /**
     * Evaluate the correlation function.
     *
     * @param radius The correlation radius
     * @param sigma Average precision
     * @param density Average protein density
     * @param range Average circle radius
     * @param amplitude Amplitude of the cluster
     * @param alpha Measure of the coherence length between circles
     * @return the value
     */
    public double evaluate(double radius, final double sigma, final double density,
        final double range, final double amplitude, final double alpha) {
      final double gr_stoch = (1.0 / (4 * Math.PI * density * sigma * sigma))
          * FastMath.exp(-radius * radius / (4 * sigma * sigma));
      final double gr_protein =
          amplitude * FastMath.exp(-radius / alpha) * Math.cos(0.5 * Math.PI * radius / range) + 1;
      return gr_stoch + gr_protein;
    }

    /**
     * Evaluate the correlation function.
     *
     * @param radius The correlation radius
     * @param parameters The parameters
     * @return the value
     */
    @Override
    public double evaluate(double radius, final double[] parameters) {
      return evaluate(radius, parameters[0], parameters[1], parameters[2], parameters[3],
          parameters[4]);
    }
  }

  /**
   * Allow optimisation using Apache Commons Math 3 Gradient Optimiser.
   */
  private static class EmulsionModelFunctionGradient extends EmulsionModelFunction
      implements MultivariateVectorFunction {
    @Override
    public double[] value(double[] variables) {
      increment();
      final double[] values = new double[x.size()];
      for (int i = 0; i < values.length; i++) {
        values[i] = evaluate(x.get(i), variables[0], variables[1], variables[2], variables[3],
            variables[4]);
      }
      return values;
    }
  }

  private static TextWindow createResultsTable() {
    return ImageJUtils.refresh(resultsTableRef, () -> {
      final StringBuilder sb = new StringBuilder();
      sb.append("Model\t");
      sb.append("Colour\t");
      sb.append("Valid\t");
      sb.append("Precision (nm)\t");
      sb.append("Density (um^-2)\t");
      sb.append("Domain Radius (nm)\t");
      // TODO - Find out the units of the domain density
      sb.append("Domain Density\t");
      sb.append("N-cluster\t");
      sb.append("Coherence\t");
      sb.append("Adjusted R2\t");
      return new TextWindow(TITLE, sb.toString(), (String) null, 800, 300);
    });
  }

  private void addResult(String model, String resultColour, boolean valid, double precision,
      double density, double domainRadius, double domainDensity, double ncluster, double coherence,
      double adjustedR2) {
    final StringBuilder sb = new StringBuilder();
    sb.append(model).append('\t');
    sb.append(resultColour).append('\t');
    sb.append(valid).append('\t');
    sb.append(MathUtils.rounded(precision, 4)).append('\t');
    sb.append(MathUtils.rounded(density * 1e6, 4)).append('\t');
    sb.append(getString(domainRadius)).append('\t');
    sb.append(getString(domainDensity)).append('\t');
    sb.append(getString(ncluster)).append('\t');
    sb.append(getString(coherence)).append('\t');
    sb.append(MathUtils.rounded(adjustedR2, 4)).append('\t');
    resultsTable.append(sb.toString());
  }

  private static String getString(double value) {
    return (value < 0) ? "-" : MathUtils.rounded(value, 4);
  }

  /**
   * Clear the current results.
   */
  static void clearResults() {
    latestResult.set(null);
  }
}
