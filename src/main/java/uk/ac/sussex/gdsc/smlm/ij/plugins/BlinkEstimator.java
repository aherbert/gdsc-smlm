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

import java.awt.Color;
import java.util.ArrayList;

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
import org.apache.commons.math3.util.Precision;

import ij.IJ;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.Plot2;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.plugins.pcpalm.Molecule;
import uk.ac.sussex.gdsc.smlm.ij.plugins.pcpalm.PCPALMMolecules;
import uk.ac.sussex.gdsc.smlm.ij.utils.LoggingOptimiserFunction;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.TraceManager;

/**
 * Estimates the flourophore blinking rate from a set of localisations. <p> Uses the method of
 * Annibale, et al (2011). Quantitative Photo Activated Localization Microscopy: Unravelling the
 * Effect of Photoblinking. PLoS ONE 6, e22678.
 */
public class BlinkEstimator implements PlugIn {
  private static String TITLE = "Blink Estimator";

  private static String inputOption = "";
  private static int s_maxDarkTime = 80;
  private static boolean s_relativeDistance = true;
  private static int histogramBins = 50;
  private static boolean showHistogram = true;
  private static double s_searchDistance = 2.5;
  private static int s_nFittedPoints = 5;
  private static int rangeFittedPoints = 0;
  private static boolean fitIntercept = true;
  private static boolean s_timeAtLowerBound = true;

  private BlinkingFunction blinkingModel;
  private double r2;
  private double adjustedR2;

  /** The milliseconds./frame */
  public double msPerFrame;
  /** The max dark time. */
  public int maxDarkTime = s_maxDarkTime;
  /** The relative distance. */
  public boolean relativeDistance = s_relativeDistance;
  /** The search distance. */
  public double searchDistance = s_searchDistance;
  /** The number of fitted points. */
  public int nFittedPoints = s_nFittedPoints;
  /** The time at lower bound flag. */
  public boolean timeAtLowerBound = s_timeAtLowerBound;
  /** The show plots flag. */
  public boolean showPlots = false;

  private double[] parameters = null;
  private boolean increaseNFittedPoints = false;

  /** {@inheritDoc} */
  @Override
  public void run(String arg) {
    SMLMUsageTracker.recordPlugin(this.getClass(), arg);

    // Require some fit results and selected regions
    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "There are no fitting results in memory");
      return;
    }

    if (!showDialog()) {
      return;
    }

    final MemoryPeakResults results =
        ResultsManager.loadInputResults(inputOption, true, null, null);
    if (results == null || results.size() == 0) {
      IJ.error(TITLE, "No results could be loaded");
      IJ.showStatus("");
      return;
    }
    msPerFrame = results.getCalibrationReader().getExposureTime();
    ImageJUtils.log("%s: %d localisations", TITLE, results.size());

    showPlots = true;
    if (rangeFittedPoints > 0) {
      computeFitCurves(results, true);
    } else {
      computeBlinkingRate(results, true);
    }
  }

  private boolean showDialog() {
    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addHelp(About.HELP_URL);

    gd.addMessage(
        "Compute the blinking rate by fitting counts to dark-time.\nSee Annibale et al (2011) PLos ONE 6, e22678.");
    ResultsManager.addInput(gd, inputOption, InputSource.MEMORY);

    gd.addNumericField("Max_dark_time (frames)", s_maxDarkTime, 0);
    gd.addNumericField("Histogram_bins", histogramBins, 0);
    gd.addCheckbox("Show_histogram", showHistogram);
    gd.addSlider("Search_distance", 0.5, 5, s_searchDistance);
    gd.addCheckbox("Relative_distance", s_relativeDistance);
    gd.addSlider("Fitted_points", 4, 15, s_nFittedPoints);
    gd.addSlider("Range_of_fitted_points", 0, 15, rangeFittedPoints);
    gd.addCheckbox("Time_at_lower_bound", s_timeAtLowerBound);
    // gd.addCheckbox("Fit_intercept", fitIntercept);
    gd.showDialog();

    if (gd.wasCanceled()) {
      return false;
    }

    inputOption = gd.getNextChoice();
    maxDarkTime = s_maxDarkTime = (int) gd.getNextNumber();
    histogramBins = (int) gd.getNextNumber();
    showHistogram = gd.getNextBoolean();
    searchDistance = s_searchDistance = gd.getNextNumber();
    relativeDistance = s_relativeDistance = gd.getNextBoolean();
    nFittedPoints = s_nFittedPoints = (int) gd.getNextNumber();
    rangeFittedPoints = (int) gd.getNextNumber();
    timeAtLowerBound = s_timeAtLowerBound = gd.getNextBoolean();
    // fitIntercept = gd.getNextBoolean();

    // Check arguments
    try {
      Parameters.isAbove("Max dark time", maxDarkTime, 3);
      Parameters.isAbove("Histogram bins", histogramBins, 1);
      Parameters.isAboveZero("Search distance", searchDistance);
      Parameters.isAbove("n-Fitted points", nFittedPoints, 3);
      Parameters.isPositive("Range of fitted points", rangeFittedPoints);
    } catch (final IllegalArgumentException ex) {
      IJ.error(TITLE, ex.getMessage());
      return false;
    }

    return true;
  }

  private void computeFitCurves(MemoryPeakResults results, boolean verbose) {
    // Calculate the counts verses dark time curve
    double[] Ntd = calculateCounts(results, maxDarkTime, searchDistance, relativeDistance, verbose);
    double[] td = calculateTd(Ntd);

    Ntd = shift(Ntd);
    td = shift(td);

    // Fit curve
    final double[] nPoints = new double[rangeFittedPoints + 1];
    final double[][] parameters = new double[3][nPoints.length];
    final double[] r2 = new double[rangeFittedPoints + 1];
    final double[] adjustedR2 = new double[rangeFittedPoints + 1];
    for (int n = 0; n <= rangeFittedPoints; n++) {
      nPoints[n] = n + nFittedPoints;
      final double[] p = fit(td, Ntd, (int) nPoints[n], false);
      if (p == null) {
        // Leave as empty in the output plots
        continue;
      }
      for (int i = 0; i < p.length; i++) {
        parameters[i][n] = p[i];
      }
      r2[n] = this.r2;
      adjustedR2[n] = this.adjustedR2;
    }

    // Plot
    plot("Fitted points", "N", nPoints, parameters[0]);
    plot("Fitted points", "nBlinks", nPoints, parameters[1]);
    plot("Fitted points", "tOff", nPoints, parameters[2]);
    plot("Fitted points", "Adjusted R^2", nPoints, adjustedR2);
  }

  /**
   * Remove the first element of the array. Return the rest of the array
   *
   * @param d the d
   * @return the shifted array
   */
  private static double[] shift(double[] d) {
    if (fitIntercept) {
      return d;
    }
    final double[] d2 = new double[d.length - 1];
    System.arraycopy(d, 1, d2, 0, d2.length);
    return d2;
  }

  private static void plot(String xAxisTitle, String yAxisTitle, double[] x, double[] y) {
    final String title = TITLE + " " + yAxisTitle;
    final Plot2 plot = new Plot2(title, xAxisTitle, yAxisTitle, x, y);
    ImageJUtils.display(title, plot);
  }

  /**
   * Compute blinking rate.
   *
   * @param results the results
   * @param verbose the verbose
   * @return the blinking rate
   */
  double computeBlinkingRate(MemoryPeakResults results, boolean verbose) {
    parameters = null;
    increaseNFittedPoints = false;

    // Calculate the counts verses dark time curve
    double[] Ntd = calculateCounts(results, maxDarkTime, searchDistance, relativeDistance, verbose);
    double[] td = calculateTd(Ntd);

    if (verbose) {
      ImageJUtils.log("  Estimate %.0f molecules at td = %.0f ms", Ntd[0], td[0]);
    }

    Ntd = shift(Ntd);
    td = shift(td);

    // Fit curve
    parameters = fit(td, Ntd, nFittedPoints, verbose);
    if (parameters == null) {
      return 0;
    }

    // Display
    if (showPlots) {
      final String title = TITLE + " Molecule Counts";
      final Plot2 plot = new Plot2(title, "td (ms)", "Count", td, Ntd);
      ImageJUtils.display(title, plot);

      plot.setColor(Color.red);
      plot.addPoints(blinkingModel.getX(), blinkingModel.value(parameters), Plot.CIRCLE);

      // Add the rest that is not fitted
      final double[] xOther = new double[td.length - blinkingModel.size()];
      final double[] yOther = new double[xOther.length];
      for (int i = 0, t = blinkingModel.size(); i < xOther.length; i++, t++) {
        xOther[i] = td[t];
        yOther[i] = blinkingModel.evaluate(td[t], parameters);
      }

      plot.setColor(Color.blue);
      plot.addPoints(xOther, yOther, Plot.CROSS);
      ImageJUtils.display(title, plot);
    }

    // Check if the fitted curve asymptotes above the real curve
    if (blinkingModel.evaluate(td[Ntd.length - 1], parameters) < Ntd[Ntd.length - 1]) {
      if (verbose) {
        ImageJUtils.log("  *** Warning ***");
        ImageJUtils.log(
            "  Fitted curve does not asymptote above real curve. Increase the number of fitted points to sample more of the overcounting regime");
        ImageJUtils.log("  ***************");
      }
      increaseNFittedPoints = true;
    }

    // Blinking rate is 1 + nBlinks
    final double blinkingRate = 1 + parameters[1];
    if (verbose) {
      ImageJUtils.log("  Blinking rate = %s", MathUtils.rounded(blinkingRate, 4));
    }
    return blinkingRate;
  }

  /**
   * Compute blinking rate.
   *
   * @param results the results
   * @return the double
   */
  public double computeBlinkingRate(MemoryPeakResults results) {
    return computeBlinkingRate(results, false);
  }

  /**
   * Calculate the counts of molecules using different dark times. The distance threshold for
   * molecule tracing will be absolute or relative. If relative it is set using the average
   * precision multiplied by the search distance. <p> Note that index 0 corresponds to a t-threshold
   * of 1 in the tracing algorithm, i.e. adjacent frames in the sequence. This is equivalent to a
   * dark time of (up to) the frame acquisition rate, i.e. the molecule is not allowed to blink.
   *
   * @param results the results
   * @param maxDarkTime the max dark time
   * @param searchDistance the search distance
   * @param relativeDistance the relative distance
   * @param verbose Output log messages
   * @return the counts of molecules
   */
  private static double[] calculateCounts(MemoryPeakResults results, int maxDarkTime,
      double searchDistance, boolean relativeDistance, boolean verbose) {
    double distanceThreshold;
    if (relativeDistance) {
      final double averagePrecision = calculateAveragePrecision(results, verbose);
      distanceThreshold = averagePrecision * searchDistance / results.getNmPerPixel();
      if (verbose) {
        ImageJUtils.log("Average precision = %f, Distance threshold = %f px", averagePrecision,
            distanceThreshold);
      }
    } else {
      distanceThreshold = searchDistance;
      ImageJUtils.log("Distance threshold = %f px", distanceThreshold);
    }

    final double[] Ntd = new double[maxDarkTime + 1];

    final TraceManager tm = new TraceManager(results);
    IJ.showStatus("Computing counts ...");
    for (int td = 0; td <= maxDarkTime; td++) {
      IJ.showProgress(td, maxDarkTime);
      Ntd[td] = tm.traceMolecules(distanceThreshold, td + 1);
    }
    IJ.showProgress(1);
    IJ.showStatus("");

    return Ntd;
  }

  /**
   * Calculate the dark time corresponding to the molecule counts. <p> Note that index 0 corresponds
   * to a t-threshold of 1 in the tracing algorithm, i.e. adjacent frames in the sequence. This is
   * equivalent to a dark time of (up to) the frame acquisition rate, i.e. the molecule is not
   * allowed to blink. <p> The returned Td values are the lower bounds of the dark time, i.e.
   * t-threshold 1 equals 0 dark frames (0ms), t-threshold 2 equals 1 dark frame (n ms per frame),
   * etc. This behaviour can be changed by setting the {@link #timeAtLowerBound} flag to false. Then
   * the time will reflect the upper bounds of the dark time, i.e. t-threshold 1 equals 1 dark
   * frames (n ms per frame), t-threshold 2 equals 2 dark frames (2n ms per frame), etc.
   *
   * @param Ntd the ntd
   * @return the dark time
   */
  public double[] calculateTd(double[] Ntd) {
    final double[] td = new double[Ntd.length];
    for (int t = 0; t < td.length; t++) {
      if (timeAtLowerBound) {
        // Using the lower bounds of the dark time allows the blink estimator to predict the sampled
        // blinks
        // statistic produced by the Create Data plugin.
        td[t] = t * msPerFrame;
      } else {
        // Adjust for the number of frames that the molecule is allowed to be in a dark state
        td[t] = (t + 1) * msPerFrame;
      }
    }
    return td;
  }

  private static double calculateAveragePrecision(MemoryPeakResults results, boolean verbose) {
    double fittedAverage = 0;

    try {
      final PCPALMMolecules fitter = new PCPALMMolecules();
      final ArrayList<Molecule> molecules = fitter.extractLocalisations(results);
      final String title = (verbose) ? TITLE + " Localisation Precision" : null;

      fittedAverage = fitter.calculateAveragePrecision(molecules, title, histogramBins, true, true);
    } catch (final DataException ex) {
      // This is thrown when the data cannot be converted for precision computation
    }

    // Sense check the precision
    if (fittedAverage < 5 || fittedAverage > 60) {
      final GenericDialog gd = new GenericDialog(TITLE);
      gd.addMessage(
          "Estimated precision is not within expected bounds.\nPlease enter an estimate:");
      gd.addSlider("Precision", 5, 60, fittedAverage);
      gd.showDialog();
      if (!gd.wasCanceled()) {
        fittedAverage = gd.getNextNumber();
      }
    }

    // The fitter does checks for a good fit to the histogram so just return the value
    return fittedAverage;
  }

  /**
   * Fit the dark time to counts of molecules curve. Only use the first n fitted points. <p>
   * Calculates:<br> N = The number of photoblinking molecules in the sample<br> nBlink = The
   * average number of blinks per flourophore<br> tOff = The off-time
   *
   * @param td The dark time
   * @param ntd The counts of molecules
   * @param nFittedPoints the number of fitted points
   * @param log Write the fitting results to the ImageJ log window
   * @return The fitted parameters [N, nBlink, tOff], or null if no fit was possible
   */
  public double[] fit(double[] td, double[] ntd, int nFittedPoints, boolean log) {
    blinkingModel = new BlinkingFunction();
    blinkingModel.setLogging(true);
    for (int i = 0; i < nFittedPoints; i++) {
      blinkingModel.addPoint(td[i], ntd[i]);
    }

    // Different convergence thresholds seem to have no effect on the resulting fit, only the number
    // of
    // iterations for convergence
    final double initialStepBoundFactor = 100;
    final double costRelativeTolerance = 1e-6;
    final double parRelativeTolerance = 1e-6;
    final double orthoTolerance = 1e-6;
    final double threshold = Precision.SAFE_MIN;
    final LevenbergMarquardtOptimizer optimiser =
        new LevenbergMarquardtOptimizer(initialStepBoundFactor, costRelativeTolerance,
            parRelativeTolerance, orthoTolerance, threshold);
    try {
      final double[] obs = blinkingModel.getY();

      //@formatter:off
      final LeastSquaresProblem problem = new LeastSquaresBuilder()
          .maxEvaluations(Integer.MAX_VALUE)
          .maxIterations(1000)
          .start(new double[] { ntd[0], 0.1, td[1] })
          .target(obs)
          .weight(new DiagonalMatrix(blinkingModel.getWeights()))
          .model(blinkingModel, new MultivariateMatrixFunction() {
            @Override
            public double[][] value(double[] point) throws IllegalArgumentException
            {
              return blinkingModel.jacobian(point);
            }} )
          //.checkerPair(checker)
          .build();
      //@formatter:on

      blinkingModel.setLogging(false);

      final Optimum optimum = optimiser.optimize(problem);

      final double[] parameters = optimum.getPoint().toArray();

      // double[] exp = blinkingModel.value(parameters);
      double mean = 0;
      for (final double d : obs) {
        mean += d;
      }
      mean /= obs.length;
      double ssResiduals = 0, ssTotal = 0;
      for (int i = 0; i < obs.length; i++) {
        // ssResiduals += (obs[i] - exp[i]) * (obs[i] - exp[i]);
        ssTotal += (obs[i] - mean) * (obs[i] - mean);
      }
      // This is true if the weights are 1
      ssResiduals = optimum.getResiduals().dotProduct(optimum.getResiduals());

      r2 = 1 - ssResiduals / ssTotal;
      adjustedR2 = getAdjustedCoefficientOfDetermination(ssResiduals, ssTotal, obs.length,
          parameters.length);

      if (log) {
        ImageJUtils.log("  Fit %d points. R^2 = %s. Adjusted R^2 = %s", obs.length,
            MathUtils.rounded(r2, 4), MathUtils.rounded(adjustedR2, 4));
        ImageJUtils.log("  N=%s, nBlink=%s, tOff=%s (%s frames)",
            MathUtils.rounded(parameters[0], 4), MathUtils.rounded(parameters[1], 4),
            MathUtils.rounded(parameters[2], 4), MathUtils.rounded(parameters[2] / msPerFrame, 4));
      }

      return parameters;
    } catch (final TooManyIterationsException ex) {
      if (log) {
        ImageJUtils.log("  Failed to fit %d points: Too many iterations: (%s)",
            blinkingModel.size(), ex.getMessage());
      }
      return null;
    } catch (final ConvergenceException ex) {
      if (log) {
        ImageJUtils.log("  Failed to fit %d points", blinkingModel.size());
      }
      return null;
    }
  }

  /**
   * Gets the n molecules.
   *
   * @return the n molecules
   */
  public double getNMolecules() {
    if (parameters != null) {
      return parameters[0];
    }
    return 0;
  }

  /**
   * Gets the n blinks.
   *
   * @return the n blinks
   */
  public double getNBlinks() {
    if (parameters != null) {
      return parameters[1];
    }
    return 0;
  }

  /**
   * Gets the t off.
   *
   * @return the t off
   */
  public double getTOff() {
    if (parameters != null) {
      return parameters[2];
    }
    return 0;
  }

  /**
   * Gets the adjusted coefficient of determination.
   *
   * @param ssResiduals Sum of squared residuals from the model
   * @param ssTotal SStotal is the sum of the squared differences from the mean of the dependent
   *        variable (total sum of squares)
   * @param n Number of observations
   * @param d Number of parameters in the model
   * @return the adjusted coefficient of determination
   */
  private static double getAdjustedCoefficientOfDetermination(double ssResiduals, double ssTotal,
      int n, int d) {
    if (n - d - 1 <= 0) {
      return 1 - (ssResiduals / ssTotal);
    }
    return 1 - (ssResiduals / ssTotal) * ((n - 1) / (n - d - 1));
  }

  /**
   * @return the coefficient of determination of the previous fit.
   */
  public double getR2() {
    return r2;
  }

  /**
   * @return the adjusted coefficient of determination of the previous fit.
   */
  public double getAdjustedR2() {
    return adjustedR2;
  }

  /**
   * @return the increaseNFittedPoints.
   */
  public boolean isIncreaseNFittedPoints() {
    return increaseNFittedPoints;
  }

  /**
   * Allow optimisation using Apache Commons Math 3 Gradient Optimiser <p> N(td) = N . (1 + nBlink .
   * exp((1-td)/tOff) <p> where <p> N(td) = The number of calculated molecules at different dark
   * times (td)<br> N = The number of photoblinking molecules in the sample<br> nBlink = The average
   * number of blinks per flourophore<br> td = The dark time<br> tOff = The off-time<br>
   */
  public class BlinkingFunction extends LoggingOptimiserFunction
      implements MultivariateVectorFunction // , DifferentiableMultivariateVectorFunction
  {
    /**
     * Instantiates a new blinking function.
     */
    public BlinkingFunction() {
      super("Blinking Model");
    }

    /** {@inheritDoc} */
    @Override
    public double[] getWeights() {
      // Bias the early values
      // double[] w = new double[y.size()];
      // for (int i = 0; i < w.length; i++)
      // w[i] = w.length - i;
      // return w;
      return super.getWeights();
    }

    // Adapted from http://commons.apache.org/proper/commons-math/userguide/optimization.html
    // Use the deprecated API since the new one is not yet documented.

    private double[][] jacobian(double[] variables) {
      // Compute the gradients using calculus differentiation
      final double N = variables[0];
      final double nBlink = variables[1];
      final double tOff = variables[2];
      final double[][] jacobian = new double[x.size()][variables.length];

      for (int i = 0; i < jacobian.length; ++i) {
        final double td = this.x.get(i);

        final double a = (1 - td) / tOff;
        final double b = FastMath.exp(a);

        // value = N * (1 + nBlink * b)
        // = N + N * nBlink * exp(a)

        // Differentiate with respect to N:
        jacobian[i][0] = 1 + nBlink * b;

        // Differentiate with respect to nBlink:
        jacobian[i][1] = N * b;

        // Differentiate with respect to tOff:
        jacobian[i][2] = N * nBlink * b * -a / tOff;
      }

      //// Check numerically ...
      // double[][] jacobian2 = jacobian2(variables);
      // for (int i = 0; i < jacobian.length; i++)
      // {
      // System.out.printf("N = %g : %g = %g. nBlink = %g : %g = %g. tOff = %g : %g = %g\n",
      //// jacobian[i][0],
      // jacobian2[i][0], DoubleEquality.relativeError(jacobian[i][0], jacobian2[i][0]),
      //// jacobian[i][1],
      // jacobian2[i][1], DoubleEquality.relativeError(jacobian[i][1], jacobian2[i][1]),
      //// jacobian[i][2],
      // jacobian2[i][2], DoubleEquality.relativeError(jacobian[i][2], jacobian2[i][2]));
      // }

      return jacobian;
    }

    @SuppressWarnings("unused")
    private double[][] jacobian2(double[] variables) {
      // Compute the gradients using numerical differentiation
      final double N = variables[0];
      final double nBlink = variables[1];
      final double tOff = variables[2];
      final double[][] jacobian = new double[x.size()][variables.length];

      final double delta = 0.001;
      final double[][] d = new double[variables.length][variables.length];
      for (int i = 0; i < variables.length; i++) {
        d[i][i] = delta * Math.abs(variables[i]); // Should the delta be changed for each parameter
                                                  // ?
      }
      for (int i = 0; i < jacobian.length; ++i) {
        final double r = this.x.get(i);
        final double value = evaluate(r, N, nBlink, tOff);
        for (int j = 0; j < variables.length; j++) {
          final double value2 = evaluate(r, N + d[0][j], nBlink + d[1][j], tOff + d[2][j]);
          jacobian[i][j] = (value2 - value) / d[j][j];
        }
      }
      return jacobian;
    }

    /**
     * Evaluate the function.
     *
     * @param td The dark time
     * @param N The number of molecules in the sample
     * @param nBlink The blinking rate
     * @param tOff The off-time
     * @return the value
     */
    public double evaluate(double td, double N, double nBlink, double tOff) {
      return N * (1.0 + nBlink * FastMath.exp((1 - td) / tOff));
    }

    /**
     * Evaluate.
     *
     * @param td the dark time
     * @param parameters the parameters
     * @return the value
     */
    public double evaluate(double td, double[] parameters) {
      return evaluate(td, parameters[0], parameters[1], parameters[2]);
    }

    /** {@inheritDoc} */
    @Override
    public double[] value(double[] variables) {
      increment();
      final double[] values = new double[x.size()];
      for (int i = 0; i < values.length; i++) {
        values[i] = evaluate(x.get(i), variables[0], variables[1], variables[2]);
      }
      return values;
    }

    /**
     * Compute the Jacobian.
     *
     * @return the multivariate matrix function
     */
    public MultivariateMatrixFunction jacobian() {
      return new MultivariateMatrixFunction() {
        @Override
        public double[][] value(double[] variables) {
          return jacobian(variables);
        }
      };
    }
  }
}
