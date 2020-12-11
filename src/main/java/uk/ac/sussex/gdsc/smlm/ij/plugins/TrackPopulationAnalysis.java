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

package uk.ac.sussex.gdsc.smlm.ij.plugins;

import ij.IJ;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.process.LUT;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import org.apache.commons.math3.distribution.FDistribution;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyIterationsException;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresFactory;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer.Optimum;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem.Evaluation;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.fitting.leastsquares.ParameterValidator;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.NonPositiveDefiniteMatrixException;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.commons.math3.util.Pair;
import org.apache.commons.rng.sampling.UnitSphereSampler;
import uk.ac.sussex.gdsc.core.data.VisibleForTesting;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot;
import uk.ac.sussex.gdsc.core.ij.HistogramPlot.BinMethod;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.MultiDialog;
import uk.ac.sussex.gdsc.core.ij.plugin.WindowOrganiser;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper.LutColour;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.LocalList;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SortUtils;
import uk.ac.sussex.gdsc.core.utils.Statistics;
import uk.ac.sussex.gdsc.core.utils.StoredData;
import uk.ac.sussex.gdsc.core.utils.rng.UniformRandomProviders;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationReader;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.function.ChiSquaredDistributionTable;
import uk.ac.sussex.gdsc.smlm.math3.distribution.fitting.MultivariateGaussianMixtureExpectationMaximization;
import uk.ac.sussex.gdsc.smlm.math3.distribution.fitting.MultivariateGaussianMixtureExpectationMaximization.DoubleDoubleBiPredicate;
import uk.ac.sussex.gdsc.smlm.math3.distribution.fitting.MultivariateGaussianMixtureExpectationMaximization.MixtureMultivariateGaussianDistribution;
import uk.ac.sussex.gdsc.smlm.math3.distribution.fitting.MultivariateGaussianMixtureExpectationMaximization.MixtureMultivariateGaussianDistribution.MultivariateGaussianDistribution;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.Trace;
import uk.ac.sussex.gdsc.smlm.results.sort.IdFramePeakResultComparator;

/**
 * Analyse tracks using a local sliding window to extract parameters that characterise the current
 * diffusion. Fit the parameters using a multi-variate Gaussian mixture model to identify
 * sub-populations such as bound and unbound particles. Output analysis on the sub-populations.
 *
 * <blockquote>Basu, et al (2020) Live-cell 3D single-molecule tracking reveals how NuRD modulates
 * enhancer dynamics. doi: <a href="https://doi.org/10.1101/2020.04.03.003178">bioRxiv
 * 2020.04.03.003178</a> </blockquote>
 */
public class TrackPopulationAnalysis implements PlugIn {
  private static final String TITLE = "Track Population Analysis";
  private static final String[] FEATURE_NAMES = {"Anomalous exponent",
      "Effective diffusion coefficient", "Length of confinement", "Drift Vector"};
  private static final int SORT_DIMENSION = 1;
  // Limits for fitting the MSD
  private static final double MIN_D = Double.MIN_NORMAL;
  private static final double MIN_SIGMA = Double.MIN_NORMAL;
  private static final double MIN_ALPHA = Math.ulp(1.0);

  /** The plugin settings. */
  private Settings settings;

  /**
   * Contains the settings that are the re-usable state of the plugin.
   */
  private static class Settings {
    /** The last settings used by the plugin. This should be updated after plugin execution. */
    private static final AtomicReference<Settings> lastSettings =
        new AtomicReference<>(new Settings());

    List<String> input;
    int window;
    int maxIterations;
    double relativeError;
    int repeats;
    int seed;
    boolean debug;
    int maxComponents;
    double minWeight;
    double minAlpha;
    double significance;
    boolean ignoreAlpha;

    Settings() {
      // Set defaults
      window = 11;
      maxIterations = 1000;
      relativeError = 1e-6;
      repeats = 20;
      seed = 42;
      maxComponents = 2;
      minWeight = 0.1;
      significance = 0.05;
    }

    Settings(Settings source) {
      this.input = source.input;
      this.window = source.window;
      this.maxIterations = source.maxIterations;
      this.relativeError = source.relativeError;
      this.repeats = source.repeats;
      this.seed = source.seed;
      this.debug = source.debug;
      this.maxComponents = source.maxComponents;
      this.minWeight = source.minWeight;
      this.significance = source.significance;
      this.minAlpha = source.minAlpha;
      this.ignoreAlpha = source.ignoreAlpha;
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

    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "No localisations in memory");
      return;
    }

    settings = Settings.load();
    // Saved by reference so just save now
    settings.save();

    // Read in multiple traced datasets
    // All datasets must have the same pixel pitch and exposure time
    // Get parameters
    // Convert datasets to tracks
    // For each track compute the 4 local track parameters using the configured window
    // Fit a multi-variate Gaussian mixture model to the data
    // (using the configured number of components/populations)
    // Assign each point in the track using the model.
    // Smooth the assignments.
    // Plot histograms of each track parameter, coloured by component
    //
    // Output for the bound component and free components track parameters
    // Compute dwell times
    // Other ...

    final List<MemoryPeakResults> combinedResults = new LocalList<>();

    if (!showInputDialog(combinedResults)) {
      return;
    }

    if (!showDialog()) {
      return;
    }

    if (combinedResults.isEmpty()) {
      return;
    }

    ImageJUtils.log(TITLE + "...");

    final List<Trace> tracks = getTracks(combinedResults, settings.window);
    if (tracks.isEmpty()) {
      return;
    }

    final double[][] data = extractTrackData(tracks);

    // Histogram the raw data.
    final Array2DRowRealMatrix raw = new Array2DRowRealMatrix(data, false);
    final WindowOrganiser wo = new WindowOrganiser();
    // Store the histogram data for plotting the components
    final double[][] columns = new double[FEATURE_NAMES.length][];
    final double[][] limits = new double[FEATURE_NAMES.length][];
    final int[] bins = new int[FEATURE_NAMES.length];
    final Plot[] plots = new Plot[FEATURE_NAMES.length];
    for (int i = 0; i < FEATURE_NAMES.length; i++) {
      columns[i] = raw.getColumn(i);
      limits[i] = MathUtils.limits(columns[i]);
      final StoredData colData = StoredData.create(columns[i]);
      bins[i] = HistogramPlot.getBins(colData, BinMethod.SCOTT);
      final double[][] hist =
          HistogramPlot.calcHistogram(columns[i], limits[i][0], limits[i][1], bins[i]);
      plots[i] = new Plot(TITLE + " " + FEATURE_NAMES[i], FEATURE_NAMES[i], "Frequency");
      plots[i].addPoints(hist[0], hist[1], Plot.BAR);
      ImageJUtils.display(plots[i].getTitle(), plots[i], 0, wo);
    }
    wo.tile();

    // Provide option to not use the anomalous exponent in the population mix.
    int sortDimension = SORT_DIMENSION;
    double[][] fitData = data;
    if (settings.ignoreAlpha) {
      // Remove index 0. This shifts the sort dimension.
      sortDimension--;
      fitData =
          Arrays.stream(data).map(d -> Arrays.copyOfRange(d, 1, d.length)).toArray(double[][]::new);
    }

    final MultivariateGaussianMixtureExpectationMaximization mixed =
        fitGaussianMixture(fitData, sortDimension);

    if (mixed == null) {
      IJ.error(TITLE, "Failed to fit a mixture model");
      return;
    }

    MixtureMultivariateGaussianDistribution model = mixed.getFittedModel();
    model = sortComponents(model, sortDimension);

    // For the best model, assign to the most likely population.
    final int[] component = assignData(fitData, model);
    final int numComponents = mixed.getFittedModel().getWeights().length;

    // Output coloured histograms of the populations.
    final LUT lut = LutHelper.createLut(LutColour.INTENSE);
    for (int i = 0; i < FEATURE_NAMES.length; i++) {
      // Extract the data for each component
      final double[] col = columns[i];
      final Plot plot = plots[i];
      for (int n = 0; n < numComponents; n++) {
        final StoredData feature = new StoredData();
        for (int j = 0; j < component.length; j++) {
          if (component[j] == n) {
            feature.add(col[j]);
          }
        }
        if (feature.size() == 0) {
          continue;
        }
        final double[][] hist =
            HistogramPlot.calcHistogram(feature.values(), limits[i][0], limits[i][1], bins[i]);
        // Colour the points
        plot.setColor(LutHelper.getColour(lut, n));
        plot.addPoints(hist[0], hist[1], Plot.BAR);
      }
      plot.updateImage();
    }
  }

  private boolean showInputDialog(List<MemoryPeakResults> combinedResults) {
    // Show a list box containing all the clustered results.
    // This should remember the last set of chosen items.
    final MultiDialog md = ResultsManager.createMultiDialog(TITLE, MemoryPeakResults::hasId);
    md.setSelected(settings.input);
    md.setHelpUrl(HelpUrls.getUrl("track-population-analysis"));

    md.showDialog();

    if (md.wasCancelled()) {
      return false;
    }

    final List<String> selected = md.getSelectedResults();
    if (selected.isEmpty()) {
      IJ.error(TITLE, "No results were selected");
      return false;
    }
    settings.input = selected;

    for (final String name : selected) {
      final MemoryPeakResults r = MemoryPeakResults.getResults(name);
      if (r != null) {
        combinedResults.add(r);
      }
    }

    if (combinedResults.isEmpty()) {
      return false;
    }

    // Check calibration exists for the first set of results
    if (!checkCalibration(combinedResults.get(0))) {
      return false;
    }

    // Check the calibration is the same for the rest
    final CalibrationReader cal = combinedResults.get(0).getCalibrationReader();
    final double nmPerPixel = cal.getNmPerPixel();
    final double exposureTime = cal.getExposureTime();
    final DistanceUnit distanceUnit = cal.getDistanceUnit();
    for (int i = 1; i < combinedResults.size(); i++) {
      final MemoryPeakResults results = combinedResults.get(i);

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
      gd.addNumericField("Exposure_time", cal.getExposureTime(), 2, 6, "ms");
      gd.addNumericField("Pixel_pitch", cal.getNmPerPixel(), 2, 6, "nm");
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

  private boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(HelpUrls.getUrl("track-population-analysis"));

    gd.addSlider("Window", 5, 20, settings.window);
    gd.addMessage("Anomalous diffusion coefficient");
    gd.addNumericField("Fit_significance", settings.significance, -2);
    gd.addNumericField("Min_alpha", settings.minAlpha, -3);
    gd.addMessage("Multi-variate Gaussian mixture Expectation-Maximisation");
    gd.addCheckbox("Ignore_alpha", settings.ignoreAlpha);
    gd.addSlider("Max_components", 2, 10, settings.maxComponents);
    gd.addSlider("Min_weight", 0.01, 1, settings.minWeight);
    gd.addNumericField("Max_iterations", settings.maxIterations, 0);
    gd.addNumericField("Relative_error", settings.relativeError, -1);
    gd.addNumericField("Repeats", settings.repeats, 0);
    gd.addNumericField("Seed", settings.seed, 0);
    gd.addCheckbox("Debug", settings.debug);

    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    settings.window = (int) gd.getNextNumber();
    settings.significance = gd.getNextNumber();
    settings.minAlpha = gd.getNextNumber();
    settings.ignoreAlpha = gd.getNextBoolean();
    settings.maxComponents = (int) gd.getNextNumber();
    settings.minWeight = gd.getNextNumber();
    settings.maxIterations = (int) gd.getNextNumber();
    settings.relativeError = gd.getNextNumber();
    settings.repeats = (int) gd.getNextNumber();
    settings.seed = (int) gd.getNextNumber();
    settings.debug = gd.getNextBoolean();

    if (gd.invalidNumber()) {
      return false;
    }

    // Check arguments
    try {
      // For fitting the alpha coefficient there should be:
      // (number of parameters) < (number of points)
      // where the number of points is the window size minus 1.
      ParameterUtils.isEqualOrAbove("Window", settings.window, 5);
      ParameterUtils.isEqualOrAbove("Max components", settings.maxComponents, 2);
      ParameterUtils.isPositive("Min_weight", settings.minWeight);
      ParameterUtils.isAboveZero("Max iterations", settings.maxIterations);
      ParameterUtils.isAboveZero("Maximum N", settings.relativeError);
      ParameterUtils.isAboveZero("Repeats", settings.repeats);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return true;
  }

  /**
   * Gets the tracks. Each track has contiguous frames and the length is at least the window size.
   *
   * @param combinedResults the combined results
   * @return the tracks
   */
  private static List<Trace> getTracks(List<MemoryPeakResults> combinedResults, int window) {
    final LocalList<Trace> tracks = new LocalList<>();
    final Statistics stats = new Statistics();
    combinedResults.forEach(results -> {
      final int start = tracks.size();
      // Sort by id then frame
      results = results.copy();
      results.sort(IdFramePeakResultComparator.INSTANCE);
      final int size = results.size();
      // Skip IDs not associated with clustering
      int index = 0;
      while (index < size && results.get(index).getId() < 1) {
        index++;
      }
      // Initialise current id and frame
      int id = results.get(index).getId() - 1;
      int frame = results.get(index).getFrame();
      Trace track = new Trace();
      for (; index < size; index++) {
        final PeakResult result = results.get(index);
        // Same ID and contiguous frames
        if (result.getId() != id || result.getFrame() != frame + 1) {
          addTrack(window, tracks, track);
          track = new Trace();
        }
        id = result.getId();
        frame = result.getFrame();
        track.add(result);
      }
      addTrack(window, tracks, track);

      stats.reset();
      for (int i = start; i < tracks.size(); i++) {
        stats.add(tracks.unsafeGet(i).size());
      }
      ImageJUtils.log("%s tracks=%d, length=%s +/- %s", results.getName(), stats.getN(),
          MathUtils.rounded(stats.getMean(), 3),
          MathUtils.rounded(stats.getStandardDeviation(), 3));
    });
    return tracks;
  }

  /**
   * Adds the track to the list of tracks if it is {@code >= window}.
   *
   * @param window the window
   * @param tracks the tracks
   * @param track the track
   */
  private static void addTrack(int window, final List<Trace> tracks, Trace track) {
    if (track.size() >= window) {
      tracks.add(track);
    }
  }

  /**
   * Extract the track data. This extracts different descriptors of the track using a rolling local
   * window.
   *
   * @param tracks the tracks
   * @return the double[][]
   */
  private double[][] extractTrackData(List<Trace> tracks) {
    final List<double[]> data = new LocalList<>(tracks.size());
    double[] x = new double[0];
    double[] y = new double[0];
    final int window = settings.window;
    final int wm1 = window - 1;
    // Use for fitting
    final double[] s = new double[wm1];
    final LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer();
    final RealVector observed = new ArrayRealVector(s, false);
    final RealMatrix weight = null;
    final ConvergenceChecker<Evaluation> checker = (iteration, previous,
        current) -> DoubleEquality.relativeError(previous.getCost(), current.getCost()) < 1e-6;
    final int maxEvaluations = Integer.MAX_VALUE;
    final int maxIterations = 3000;
    final boolean lazyEvaluation = false;

    // Linear model for Brownian motion
    final MultivariateJacobianFunction model1 = new BrownianDiffusionFunction(wm1);
    final RealVector start1 = new ArrayRealVector(2);
    final ParameterValidator paramValidator1 = point -> {
      // Ensure diffusion coefficient and precision are positive
      final double d = point.getEntry(0);
      final double sigma = point.getEntry(1);
      // Do not use MIN_VALUE here to avoid sub-normal numbers
      if (d < MIN_D) {
        point.setEntry(0, MIN_D);
      }
      if (sigma < MIN_SIGMA) {
        point.setEntry(1, MIN_SIGMA);
      }
      return point;
    };

    // Non-linear model for anomalous diffusion coefficient
    final MultivariateJacobianFunction model2 = new FbmDiffusionFunction(wm1);
    final RealVector start2 = new ArrayRealVector(3);
    final double minAlpha = Math.max(settings.minAlpha, MIN_ALPHA);
    final ParameterValidator paramValidator2 = point -> {
      // Ensure diffusion coefficient and precision are positive.
      paramValidator1.validate(point);
      // Ensure alpha in range (0, 2]
      final double alpha = point.getEntry(2);
      // Since the computations require (alpha + 1) limit this to the ULP of 1.0
      // so the parameter always has an effect.
      if (alpha < minAlpha) {
        point.setEntry(2, MIN_ALPHA);
      } else if (alpha > 2) {
        point.setEntry(2, 2);
      }
      return point;
    };

    // For linear fit estimation
    final SimpleRegression reg = new SimpleRegression(true);
    // For significance test of the least squares fit.
    // numeratorDegreesOfFreedom = numberOfParameters2 - numberOfParameters1
    // denominatorDegreesOfFreedom = numberOfPoints - numberOfParameters2
    final int denominatorDegreesOfFreedom = wm1 - 3;
    final FDistribution distribution = new FDistribution(null, 1, denominatorDegreesOfFreedom);
    final double significance = settings.significance;

    int insignificant = 0;

    // Used for the standard deviations
    final Statistics statsX = new Statistics();
    final Statistics statsY = new Statistics();
    final Ticker ticker = ImageJUtils.createTicker(tracks.size(), 1, "Computing track features...");

    // Process each track
    for (final Trace track : tracks) {
      // Get xy coordinates
      final int size = track.size();
      if (x.length < size) {
        x = new double[size];
        y = new double[size];
      }
      for (int i = 0; i < size; i++) {
        final PeakResult peak = track.get(i);
        x[i] = peak.getXPosition();
        y[i] = peak.getYPosition();
      }
      final int smwm1 = size - wm1;
      for (int k = 0; k < smwm1; k++) {
        final double[] values = new double[4];
        data.add(values);

        // First point in window = k
        // Last point in window = k + w - 1 = k + wm1
        final int end = k + wm1;

        // 1. Anomalous exponent.

        // Compute the MSD for all distances from m=0, to m=window-1
        for (int m = 1; m <= wm1; m++) {
          // Use intermediate points to compute an average
          double msd = 0;
          for (int i = end - m; i >= k; i--) {
            msd += MathUtils.distance2(x[i], y[i], x[i + m], y[i + m]);
          }
          // Number of points = window - m
          s[m - 1] = msd / (window - m);
        }

        // Note:
        // If the alpha cannot significantly improve the fit then its value is likely to be
        // bad and the feature is noise. This can be judged using an F-test.
        // The local MSD curve can be highly variable, especially
        // if the molecule has switched from bound to unbound or vice versa. In this case there
        // is no good value for alpha.
        // Setting alpha to 1.0 for insignificant fits can result in no variance in the feature
        // and failure to invert the covariance matrix during Expectation-Maximisation.
        // Thus the user has the option to remove alpha from the features
        // before fitting the Gaussian mixture. The value is still computed for the results.

        // The data is small so compute the gradient and intercept using a linear regression.
        // This could exploit pre-computation of the regression x components.
        reg.clear();
        for (int n = 1; n < window; n++) {
          reg.addData(n, s[n - 1]);
        }

        // Compute linear fit with the (n-1/3) correction factor:
        // MSD = 4D n - (4D) / 3 + 4 s^2
        final double slope = reg.getSlope();
        if (slope > 4 * MIN_D) {
          start1.setEntry(0, slope / 4);
          start1.setEntry(1, Math.sqrt(Math.max(0, reg.getIntercept() + slope / 3) / 4));
        } else {
          // Do not allow negative slope. Use a flat line and
          // precision is derived from the mean of the MSD.
          // Using D=0 has no gradient and the fit does not explore the parameter space.
          // TODO: Check if any useful values come out of this situation. If not then fitting
          // should be skipped and the alpha set to 1.0 (insignificant).
          start1.setEntry(0, MIN_D);
          start1.setEntry(1, Math.sqrt((MathUtils.sum(s) / wm1) / 4));
        }
        final LeastSquaresProblem problem1 = LeastSquaresFactory.create(model1, observed, start1,
            weight, checker, maxEvaluations, maxIterations, lazyEvaluation, paramValidator1);
        try {
          final Optimum lvmSolution1 = optimizer.optimize(problem1);

          // Fit the anomalous exponent alpha. Start from the linear fit result.
          final RealVector fit1 = lvmSolution1.getPoint();
          start2.setEntry(0, fit1.getEntry(0));
          start2.setEntry(1, fit1.getEntry(1));
          start2.setEntry(2, 1.0);
          final LeastSquaresProblem problem2 = LeastSquaresFactory.create(model2, observed, start2,
              weight, checker, maxEvaluations, maxIterations, lazyEvaluation, paramValidator2);
          final Optimum lvmSolution2 = optimizer.optimize(problem2);

          // Check for model improvement
          final double rss1 = getResidualSumOfSquares(lvmSolution1);
          final double rss2 = getResidualSumOfSquares(lvmSolution2);

          // Compare significance using the method from RegressionUtils
          // F = ((rss1 - rss2) / (p2 - p1)) / (rss2 / (n - p2))
          final double num = rss1 - rss2;
          final double denom = rss2 / denominatorDegreesOfFreedom;
          // Note: If ss2 = 0 then f-statistic is +infinity and the cumulative probability is NaN
          // and the p-value will be accepted
          final double fStatistic = MathUtils.div0(num, denom);
          final double pValue = 1.0 - distribution.cumulativeProbability(fStatistic);

          // If better then this is anomalous diffusion
          final double alpha = lvmSolution2.getPoint().getEntry(2);
          values[0] = alpha;
          if (pValue > significance) {
            insignificant++;
          }

          // // Debug
          // if (pValue < 0.01 /*alpha > 0.0 && alpha < 0.2*/) {
          // final RealVector p = lvmSolution2.getPoint();
          // final String title = "anomalous exponent";
          // final Plot plot = new Plot(title, "time", "MSD");
          // final double[] t = SimpleArrayUtils.newArray(s.length, 1.0, 1.0);
          // plot.addPoints(t, s, Plot.CROSS);
          // plot.addPoints(t, model2.value(p).getFirst().toArray(), Plot.LINE);
          // plot.addLabel(0, 0, lvmSolution2.getPoint().toString() + " p=" + pValue);
          // ImageJUtils.display(title, plot, ImageJUtils.NO_TO_FRONT);
          // System.out.println(lvmSolution2.getPoint());
          // }
        } catch (TooManyIterationsException | ConvergenceException ex) {
          if (settings.debug) {
            ImageJUtils.log("Failed to fit anomalous exponent: " + ex.getMessage());
          }
          // Ignore this and leave as Brownian motion
          values[0] = 1.0;
        }

        // 2. Effective diffusion coefficient.
        // Sum the squared jump distances and normalise by 1 / 2dim, i.e. 1 / 4
        double sum = 0;
        for (int i = k; i < end; i++) {
          sum += MathUtils.distance2(x[i], y[i], x[i + 1], y[i + 1]);
        }
        values[1] = 0.25 * sum;

        // 3. Length of confinement.
        // Compute the average of the standard deviation of the position in each dimension.
        statsX.reset();
        statsY.reset();
        for (int i = k; i <= end; i++) {
          statsX.add(x[i]);
          statsY.add(y[i]);
        }
        values[2] = (statsX.getStandardDeviation() + statsY.getStandardDeviation()) / 2;

        // 4. Magnitude of drift vector.
        // Note: The vector is given as a sum of the distance between successive points.
        // This cancels and the result is the distance between the first and last point.
        values[3] = MathUtils.distance(x[k], y[k], x[end], y[end]);
      }
      ticker.tick();
    }
    ImageJUtils.finished();
    ImageJUtils.log("Insignificant anomalous coefficients: %d / %d", insignificant, data.size());
    return data.toArray(new double[0][0]);
  }

  /**
   * Gets the residual sum of squares.
   *
   * @param lvmSolution1 the lvm solution
   * @return the residual sum of squares
   */
  private static double getResidualSumOfSquares(final Optimum lvmSolution1) {
    final RealVector res = lvmSolution1.getResiduals();
    return res.dotProduct(res);
  }

  /**
   * Fit the Gaussian mixture to the data. The fitter with the highest likelihood from a number of
   * repeats is returned.
   *
   * @param data the data
   * @param sortDimension the sort dimension
   * @return the multivariate gaussian mixture
   */
  private MultivariateGaussianMixtureExpectationMaximization
      fitGaussianMixture(final double[][] data, int sortDimension) {
    // Get the unmixed multivariate Guassian.
    final MultivariateGaussianDistribution unmixed =
        MultivariateGaussianMixtureExpectationMaximization.createUnmixed(data);

    final int dimensions = unmixed.getMeans().length;

    // Record the likelihood of the unmixed model
    double logLikelihood = Arrays.stream(data).mapToDouble(unmixed::density).map(Math::log).sum();
    // x means, x*x covariances
    final int parametersPerGaussian = dimensions + dimensions * dimensions;
    double aic = MathUtils.getAkaikeInformationCriterion(logLikelihood, parametersPerGaussian);
    double bic = MathUtils.getBayesianInformationCriterion(logLikelihood, data.length,
        parametersPerGaussian);
    ImageJUtils.log("1 component log-likelihood=%s. AIC=%s. BIC=%s", logLikelihood, aic, bic);

    // Fit a mixed component model.
    // Increment the number of components up to a maximim or when the model does not improve.
    MultivariateGaussianMixtureExpectationMaximization mixed = null;
    for (int numComponents = 2; numComponents <= settings.maxComponents; numComponents++) {
      final MultivariateGaussianMixtureExpectationMaximization mixed2 =
          createMixed(data, dimensions, numComponents);
      if (mixed2 == null) {
        ImageJUtils.log("Failed to fit a %d component mixture model", numComponents);
        break;
      }
      final double logLikelihood2 = mixed2.getLogLikelihood();
      // n * (means, covariances, 1 weight) - 1
      // (Note: subtract 1 as the weights are constrained by summing to 1)
      final int param2 = numComponents * (parametersPerGaussian + 1) - 1;
      final double aic2 = MathUtils.getAkaikeInformationCriterion(logLikelihood2, param2);
      final double bic2 =
          MathUtils.getBayesianInformationCriterion(logLikelihood2, data.length, param2);

      // Log-likelihood ratio test statistic
      final double lambdaLr = -2 * (logLikelihood - logLikelihood2);
      // DF = difference in dimensionality from previous number of components
      // means, covariances, 1 weight
      final int degreesOfFreedom = parametersPerGaussian + 1;
      final double q = ChiSquaredDistributionTable.computeQValue(lambdaLr, degreesOfFreedom);
      ImageJUtils.log("%d component log-likelihood=%s. AIC=%s. BIC=%s. LLR significance=%s.",
          numComponents, logLikelihood2, aic2, bic2, MathUtils.rounded(q));
      final double[] weights = mixed2.getFittedModel().getWeights();

      // For consistency sort the mixture by the mean of the diffusion coefficient
      final double[] values = Arrays.stream(mixed2.getFittedModel().getDistributions())
          .mapToDouble(d -> d.getMeans()[sortDimension]).toArray();
      SortUtils.sortData(weights, values, false, false);
      ImageJUtils.log("Population weights: " + Arrays.toString(weights));

      if (MathUtils.min(weights) < settings.minWeight) {
        ImageJUtils.log("%d component model has population weight %s under minimum level %s",
            numComponents, MathUtils.min(weights), settings.minWeight);
        break;
      }
      if (aic <= aic2 || bic <= bic2 || q > 0.001) {
        ImageJUtils.log("%d component model is not significant", numComponents);
        break;
      }
      aic = aic2;
      bic = bic2;
      logLikelihood = logLikelihood2;
      mixed = mixed2;
    }

    return mixed;
  }

  /**
   * Creates the multivariate gaussian mixture as the best of many repeats of the expectation
   * maximisation algorithm.
   *
   * @param data the data
   * @param dimensions the dimensions
   * @param numComponents the number of components
   * @return the multivariate gaussian mixture expectation maximization
   */
  private MultivariateGaussianMixtureExpectationMaximization createMixed(final double[][] data,
      int dimensions, int numComponents) {
    // Fit a mixed multivariate Gaussian with different repeats.
    final UnitSphereSampler sampler =
        new UnitSphereSampler(dimensions, UniformRandomProviders.create(settings.seed++));
    final LocalList<CompletableFuture<MultivariateGaussianMixtureExpectationMaximization>> results =
        new LocalList<>(settings.repeats);
    final DoubleDoubleBiPredicate test = createConvergenceTest(settings.relativeError);
    if (settings.debug) {
      ImageJUtils.log("  Fitting %d components", numComponents);
    }
    final Ticker ticker = ImageJUtils.createTicker(settings.repeats, 2, "Fitting...");
    final AtomicInteger failures = new AtomicInteger();
    for (int i = 0; i < settings.repeats; i++) {
      final double[] vector = sampler.nextVector();
      results.add(CompletableFuture.supplyAsync(() -> {
        final MultivariateGaussianMixtureExpectationMaximization fitter =
            new MultivariateGaussianMixtureExpectationMaximization(data);
        try {
          // This may also throw the same exceptions due to inversion of the covariance matrix
          final MixtureMultivariateGaussianDistribution initialMixture =
              MultivariateGaussianMixtureExpectationMaximization.estimate(data, numComponents,
                  point -> {
                    double dot = 0;
                    for (int j = 0; j < dimensions; j++) {
                      dot += vector[j] * point[j];
                    }
                    return dot;
                  });
          final boolean result = fitter.fit(initialMixture, settings.maxIterations, test);
          // Log the result. Note: The ImageJ log is synchronized.
          if (settings.debug) {
            ImageJUtils.log("  Fit: log-likelihood=%s, iter=%d, converged=%b",
                fitter.getLogLikelihood(), fitter.getIterations(), result);
          }
          return result ? fitter : null;
        } catch (NonPositiveDefiniteMatrixException | SingularMatrixException ex) {
          failures.getAndIncrement();
          if (settings.debug) {
            ImageJUtils.log("  Fit failed during iteration %d", fitter.getIterations());
          }
        } finally {
          ticker.tick();
        }
        return null;
      }));
    }
    ImageJUtils.finished();
    if (failures.get() != 0 && settings.debug) {
      ImageJUtils.log("  %d component fit failed %d/%d", numComponents, failures.get(),
          settings.repeats);
    }
    // Collect results and return the best model.
    return results.stream().map(f -> f.join()).filter(f -> f != null)
        .sorted((f1, f2) -> Double.compare(f2.getLogLikelihood(), f1.getLogLikelihood()))
        .findFirst().orElse(null);
  }

  /**
   * Creates the convergence test.
   *
   * @param relativeError the relative error
   * @return the predicate
   */
  private static DoubleDoubleBiPredicate createConvergenceTest(double relativeError) {
    return (v1, v2) -> DoubleEquality.relativeError(v1, v2) < relativeError;
  }

  /**
   * Sort the components by the mean of the given dimension.
   *
   * @param model the model
   * @param dimension the dimension
   * @return the mixture multivariate gaussian distribution
   */
  private static MixtureMultivariateGaussianDistribution
      sortComponents(MixtureMultivariateGaussianDistribution model, int dimension) {
    final double[] weights = model.getWeights();
    final MultivariateGaussianDistribution[] distributions = model.getDistributions();
    final LocalList<Pair<Double, MultivariateGaussianDistribution>> list =
        new LocalList<>(weights.length);
    for (int i = 0; i < weights.length; i++) {
      list.add(Pair.create(weights[i], distributions[i]));
    }
    list.sort((o1, o2) -> Double.compare(o1.getSecond().getMeans()[dimension],
        o2.getSecond().getMeans()[dimension]));
    for (int i = 0; i < weights.length; i++) {
      weights[i] = list.unsafeGet(i).getFirst();
      distributions[i] = list.unsafeGet(i).getSecond();
    }
    return MixtureMultivariateGaussianDistribution.create(weights, distributions);
  }

  /**
   * Assign the data using the mixture model.
   *
   * @param data the data
   * @param model the model
   * @return the assignments
   */
  private static int[] assignData(double[][] data, MixtureMultivariateGaussianDistribution model) {
    final double[] weights = model.getWeights();
    final MultivariateGaussianDistribution[] distributions = model.getDistributions();
    // All initialised as component 0
    final int[] comp = new int[data.length];
    for (int i = 0; i < comp.length; i++) {
      // Assign using the highest probability.
      final double[] x = data[i];
      double max = weights[0] * distributions[0].density(x);
      for (int j = 1; j < weights.length; j++) {
        final double p = weights[j] * distributions[j].density(x);
        if (max < p) {
          max = p;
          comp[i] = j;
        }
      }
    }
    // Median smoothing window of 3 eliminates isolated assignments, e.g. C in between U:
    // U C U => U U U
    // Note: For more than two components this does nothing to isolated assignments
    // between different neighbours:
    // U C X => U C X
    int n = 0;
    if (comp.length >= 3) {
      int p0 = comp[0];
      int p1 = comp[1];
      for (int i = 2; i < comp.length; i++) {
        final int p2 = comp[i];
        if (p0 != p1 && p0 == p2) {
          comp[i - 1] = p2;
          n++;
        }
        p0 = p1;
        p1 = p2;
      }
    }
    ImageJUtils.log("Re-labelled isolated assignments: %d / %d", n, comp.length);
    return comp;
  }

  /**
   * Define a Brownian motion diffusion function for use in fitting MSD for the case of pure
   * Brownian motion. The function is corrected for static (localisation precision) and dynamic
   * error (localisation error due to diffusion):
   *
   * <pre>
   * MSD = 4Dt * (n - 1/3) + 4s^2
   *     = 4Dtn - 4Dt/3 + 4s^2
   * </pre>
   *
   * <p>D is the diffusion coefficient; t is the acquisition time step; n is the number of frames;
   * and s is the static localisation precision (see Backlund, et al (2015), Physical Review E 91:
   * 062716, eq. 1). Note s is composed of the localisation precision of an unmoving particle (s0)
   * and a contribution due to the effect of a moving particle (see eq. 5).
   *
   * <p>This is a special case power function assuming n is an integer in {@code [1, size]} and t is
   * 1 (thus can be ignored).
   */
  @VisibleForTesting
  static final class BrownianDiffusionFunction implements MultivariateJacobianFunction {
    private static final double THIRD = 1 / 3.0;
    private final int size;
    /** The pre-compute 4 * (n - 1/3). */
    private final double[] scale;

    /**
     * Create an instance.
     *
     * @param size the maximum size of n (inclusive)
     */
    BrownianDiffusionFunction(int size) {
      this.size = size;
      // Pre-computation of the scale for n
      scale = new double[size];
      for (int n = 1; n <= size; n++) {
        scale[n - 1] = 4 * (n - THIRD);
      }
    }

    @Override
    public Pair<RealVector, RealMatrix> value(RealVector point) {
      final double[] value = new double[size];
      final double[][] jacobian = new double[size][2];
      final double d = point.getEntry(0);
      final double s = point.getEntry(1);
      // MSD = 4D * (n - 1/3) + 4s^2
      // dy_dD = 4 * (n - 1/3)
      // dy_ds = 8s
      final double ss4 = 4 * s * s;
      final double s8 = 8 * s;
      // Here the variable n represents (n+1) but is used as an index to a pre-computed scale
      for (int n = 0; n < size; n++) {
        value[n] = d * scale[n] + ss4;
        jacobian[n][0] = scale[n];
        jacobian[n][1] = s8;
      }
      return new Pair<>(new ArrayRealVector(value, false),
          new Array2DRowRealMatrix(jacobian, false));
    }
  }


  /**
   * Define a fractional Brownian motion (FBM) diffusion function for use in fitting MSD for the
   * case of anomalous diffusion. The function is corrected for static (localisation precision) and
   * dynamic error (localisation error due to diffusion):
   *
   * <pre>
   * MSD = [4Dt^a / (a+2)(a+1)] * [(n+1)^(a+2) + (n-1)^(a+2) - 2n^(a+2)]
   *       - [8Dt^a / (a+2)(a+1)] + 4s^2
   * </pre>
   *
   * <p>D is the diffusion coefficient; t is the acquisition time step; n is the number of frames; a
   * is the anomalous coefficient alpha; and s is the static localisation precision (see Backlund,
   * et al (2015), Physical Review E 91: 062716, eq. 8). Note s is composed of the localisation
   * precision of an unmoving particle (s0) and a contribution due to the effect of a moving
   * particle (see eq. 5). The factor alpha should be in the range {@code (0, 2]}.
   *
   * <p>This is a special case power function assuming n is an integer in {@code [1, size]} and t is
   * 1 (thus can be ignored).
   */
  @VisibleForTesting
  static final class FbmDiffusionFunction implements MultivariateJacobianFunction {
    private final int size;
    private final double[] log;

    /**
     * Create an instance.
     *
     * @param size the maximum size of n (inclusive)
     */
    FbmDiffusionFunction(int size) {
      this.size = size;
      log = new double[size + 2];
      for (int n = 2; n < log.length; n++) {
        log[n] = Math.log(n);
      }
    }

    @Override
    public Pair<RealVector, RealMatrix> value(RealVector point) {
      final double[] value = new double[size];
      final double[][] jacobian = new double[size][3];
      final double d = point.getEntry(0);
      final double s = point.getEntry(1);
      final double a = point.getEntry(2);
      // MSD = [4D / (a+2)(a+1)] * [(n+1)^(a+2) + (n-1)^(a+2) - 2n^(a+2)]
      // - [8D / (a+2)(a+1)] + 4s^2
      // dy_dD = [4 / (a+2)(a+1)] * [(n+1)^(a+2) + (n-1)^(a+2) - 2n^(a+2)] - 8 / (a+2)(a+1)
      // dy_ds = 8s
      // Use product rule with the following parts:
      // f(a) = 4D / (a+2)(a+1)
      // f'(a) = -4D / (a+2)(a+1)(a+2) -4 / (a+2)(a+1)(a+1)
      // g(a) = (n+1)^(a+2) + (n-1)^(a+2) - 2n^(a+2)
      // g'(a) = (n+1)^(a+2)*log(n+1) + (n-1)^(a+2)*log(n-1) - 2n^(a+2)*log(n)
      // dy_da = (f'(a) * g(a) + f * g'(a)) - 2 * f'(a)
      final double ss4 = 4 * s * s;
      final double s8 = 8 * s;
      final double a1 = a + 1;
      final double a2 = a + 2;
      final double a2a1 = a2 * a1;
      final double fa = 4 / a2a1;

      // Pre-compute powers of n^(a+2)
      final double[] na = new double[size + 2];
      for (int n = 1; n < na.length; n++) {
        na[n] = Math.pow(n, a2);
      }

      for (int n = 1; n <= size; n++) {
        final int i = n - 1;
        final double ga = na[n + 1] + na[i] - 2 * na[n];
        jacobian[i][0] = fa * ga - 2 * fa;
        value[i] = d * jacobian[i][0] + ss4;
        jacobian[i][1] = s8;
        final double fpa = -fa / a2 - fa / a1;
        final double gpa = (na[n + 1] * log[n + 1] + na[i] * log[i] - 2 * na[n] * log[n]);
        jacobian[i][2] = d * ((fpa * ga + fa * gpa) - 2 * fpa);
      }
      return new Pair<>(new ArrayRealVector(value, false),
          new Array2DRowRealMatrix(jacobian, false));
    }
  }

  /**
   * Define a power function for use in fitting.
   *
   * <pre>
   * y = beta * x ^ alpha
   * </pre>
   *
   * <p>This is a special case power function assuming x is an integer in {@code [0, size)}.
   */
  @VisibleForTesting
  static final class PowerFunction implements MultivariateJacobianFunction {
    private final int size;
    private final double[] logx;

    /**
     * Create an instance.
     *
     * @param size the maximum size of x (exclusive)
     */
    PowerFunction(int size) {
      this.size = size;
      logx = new double[size];
      for (int x = 2; x < size; x++) {
        logx[x] = Math.log(x);
      }
    }

    @Override
    public Pair<RealVector, RealMatrix> value(RealVector point) {
      final double[] value = new double[size];
      final double[][] jacobian = new double[size][2];
      final double alpha = point.getEntry(0);
      final double beta = point.getEntry(1);
      // y = beta * x^alpha
      // dy_da = beta * x^alpha * log(x)
      // dy_db = x^alpha
      // Note: for x=0 the value is always zero and the gradient with respect
      // to alpha or beta is zero (i.e. they have no effect on the value).
      // At x=1 the value is always beta and the gradient is fixed.
      value[1] = beta;
      // jacobian[1][0] = 0; // dy_da == 0
      jacobian[1][1] = 1; // dy_db == 1
      // Thus we evaluate for x>1.
      for (int x = 2; x < size; x++) {
        // y = beta * x^alpha
        // dy_da = beta * x^alpha * log(x)
        // dy_db = x^alpha
        final double xa = Math.pow(x, alpha);
        final double betaXa = beta * xa;
        value[x] = betaXa;
        jacobian[x][0] = betaXa * logx[x];
        jacobian[x][1] = xa;
      }
      return new Pair<>(new ArrayRealVector(value, false),
          new Array2DRowRealMatrix(jacobian, false));
    }
  }

  /**
   * Define a power function for use in fitting with an offset C.
   *
   * <pre>
   * y = C + beta * x ^ alpha
   * </pre>
   *
   * <p>This is a special case power function assuming x is an integer in {@code [0, size)}.
   *
   * <p>The constant represents the intercept at t=0 which is composed of the localisation precision
   * (s) and a factor for the MSD correction used to represent the loss of distance observed for
   * representing diffusing particles over the integral of a frame with a single location. This
   * factor is (n-1/3)/n with n the number of frames for the jump distance. Excluding the alpha
   * factor the MSD is represented as:
   *
   * <pre>
   * MSD = 4Dt n * (n - 1/3)/n + 4 s^2
   * MSD = 4Dt n - (4Dt) / 3 + 4 s^2
   * </pre>
   *
   * <p>The intercept can be negative and so C need not be constrained during fitting.
   */
  @VisibleForTesting
  static final class OffsetPowerFunction implements MultivariateJacobianFunction {
    private final int size;
    private final double[] logx;

    /**
     * Create an instance.
     *
     * @param size the maximum size of x (exclusive)
     */
    OffsetPowerFunction(int size) {
      this.size = size;
      logx = new double[size];
      for (int x = 2; x < size; x++) {
        logx[x] = Math.log(x);
      }
    }

    @Override
    public Pair<RealVector, RealMatrix> value(RealVector point) {
      final double[] value = new double[size];
      final double[][] jacobian = new double[size][3];
      final double alpha = point.getEntry(0);
      final double beta = point.getEntry(1);
      final double c = point.getEntry(2);
      // y = C + beta * x^alpha
      // dy_da = beta * x^alpha * log(x)
      // dy_db = x^alpha
      // dy_dc = 1
      // Note: for x=0 the value is always C and the gradient with respect
      // to alpha or beta is zero (i.e. they have no effect on the value).
      // At x=1 the value is always beta+C and the gradient is fixed.
      value[0] = c;
      jacobian[0][2] = 1;
      value[1] = c + beta;
      // jacobian[1][0] = 0; // dy_da == 0
      jacobian[1][1] = 1; // dy_db == 1
      jacobian[1][2] = 1; // dy_dc == 1
      // Thus we evaluate for x>1.
      for (int x = 2; x < size; x++) {
        final double xa = Math.pow(x, alpha);
        final double betaXa = beta * xa;
        value[x] = c + betaXa;
        jacobian[x][0] = betaXa * logx[x];
        jacobian[x][1] = xa;
        jacobian[x][2] = 1;
      }
      return new Pair<>(new ArrayRealVector(value, false),
          new Array2DRowRealMatrix(jacobian, false));
    }
  }

  /**
   * Define a power function for use in fitting with an offset C.
   *
   * <pre>
   * y = C + beta * x ^ alpha
   * </pre>
   *
   * <p>This is a special case power function assuming x is an integer in {@code [1, size)}.
   *
   * <p>The constant represents the intercept at t=0 which is composed of the localisation precision
   * (s) and a factor for the MSD correction used to represent the loss of distance observed for
   * representing diffusing particles over the integral of a frame with a single location. This
   * factor is (n-1/3)/n with n the number of frames for the jump distance. Excluding the alpha
   * factor the MSD is represented as:
   *
   * <pre>
   * MSD = 4Dt n * (n - 1/3)/n + 4 s^2
   * MSD = 4Dt n - (4Dt) / 3 + 4 s^2
   * </pre>
   *
   * <p>The intercept can be negative and so C need not be constrained during fitting.
   *
   * <p>This function is used when the MSD at time point zero is unknown. Fitting with an MSD of
   * zero at t=0 can result in extremely small values of the alpha factor ({@code < 1e-100}) as the
   * curve is constrained by matching a value that is not realistic. Due to localisation error the
   * value at t=0 would be above zero and this may be offset may be significant if the molecule is
   * not diffusing as the MSD curve will be flat.
   */
  @VisibleForTesting
  static final class OffsetPowerFunction1 implements MultivariateJacobianFunction {
    private final int size;
    private final double[] logx;

    /**
     * Create an instance.
     *
     * @param size the maximum size of x (exclusive)
     */
    OffsetPowerFunction1(int size) {
      this.size = size - 1;
      logx = new double[size];
      for (int x = 2; x < size; x++) {
        logx[x] = Math.log(x);
      }
    }

    @Override
    public Pair<RealVector, RealMatrix> value(RealVector point) {
      final double[] value = new double[size];
      final double[][] jacobian = new double[size][3];
      final double alpha = point.getEntry(0);
      final double beta = point.getEntry(1);
      final double c = point.getEntry(2);
      // y = C + beta * x^alpha
      // dy_da = beta * x^alpha * log(x)
      // dy_db = x^alpha
      // dy_dc = 1
      // Note: Skip evaluation for x=0. Store x=1 at index 0, etc.
      // At x=1 the value is always beta+C and the gradient is fixed.
      value[0] = c + beta;
      // jacobian[0][0] = 0; // dy_da == 0
      jacobian[0][1] = 1; // dy_db == 1
      jacobian[0][2] = 1; // dy_dc == 1
      // Thus we evaluate for x>1.
      for (int x = 2; x <= size; x++) {
        final double xa = Math.pow(x, alpha);
        final double betaXa = beta * xa;
        final int i = x - 1;
        value[i] = c + betaXa;
        jacobian[i][0] = betaXa * logx[x];
        jacobian[i][1] = xa;
        jacobian[i][2] = 1;
      }
      return new Pair<>(new ArrayRealVector(value, false),
          new Array2DRowRealMatrix(jacobian, false));
    }
  }
}
