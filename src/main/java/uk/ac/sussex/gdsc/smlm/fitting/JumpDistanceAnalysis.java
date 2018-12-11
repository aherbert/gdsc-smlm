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

package uk.ac.sussex.gdsc.smlm.fitting;

import uk.ac.sussex.gdsc.core.logging.LoggerUtils;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SortUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.function.ChiSquaredDistributionTable;
import uk.ac.sussex.gdsc.smlm.math3.optim.nonlinear.scalar.noderiv.CustomPowellOptimizer;

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
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.util.FastMath;

import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Perform curve fitting on a cumulative histogram of the mean-squared displacement (MSD) per second
 * to calculate the diffusion coefficient of molecules (in um^2/s). The MSD is also known as the
 * Jump Distance, i.e. how far a molecule moves when being tracked.
 *
 * <p>Based on the paper: Weimann, L., Ganzinger, K.A., McColl, J., Irvine, K.L., Davis, S.J., Gay,
 * N.J., Bryant, C.E., Klenerman, D. (2013) A Quantitative Comparison of Single-Dye Tracking
 * Analysis Tools Using Monte Carlo Simulations. PLoS One 8, Issue 5, e64287
 */
public class JumpDistanceAnalysis {
  private static final boolean DEBUG_OPTIMISER = false;
  private static final double THIRD = 1.0 / 3.0;
  private double s2 = 0;
  private boolean msdCorrection = false;
  private int n = 0;
  private double deltaT = 0;
  private boolean calibrated = false;

  /**
   * Interface to logger to record a jump distance curve.
   */
  public interface CurveLogger {
    /**
     * Get the number of points to use for the curve between the minimum and maximum exclusive. The
     * size of the curve arrays will be this value plus 1 (to include the maximum).
     *
     * @return The number of points to use for the curve between the minimum and maximum exclusive
     */
    public int getNumberOfCurvePoints();

    /**
     * Called with the best fit curve using a single population.
     *
     * @param curve the curve
     */
    public void saveSinglePopulationCurve(double[][] curve);

    /**
     * Called with the best fit curve using a mixed population.
     *
     * @param curve the curve
     */
    public void saveMixedPopulationCurve(double[][] curve);
  }

  private final Logger logger;
  private CurveLogger curveLogger;
  private int fitRestarts = 3;
  private double minFraction = 0.1;
  private double minDifference = 2;
  private double minD = 0;
  private int minN = 1;
  private int maxN = 10;
  private double significanceLevel = 0.05;

  // Set by the last call to the doFit functions
  private double ss;
  private double ll;
  private double fitValue;
  // Set by any public fit call
  private double lastFitValue;

  /**
   * Instantiates a new jump distance analysis.
   */
  public JumpDistanceAnalysis() {
    this(null);
  }

  /**
   * Instantiates a new jump distance analysis.
   *
   * @param logger Used to write status messages on the fitting
   */
  public JumpDistanceAnalysis(Logger logger) {
    this.logger = LoggerUtils.createIfNull(logger);
    resetFitResult();
  }

  private void resetFitResult() {
    lastFitValue = Double.NaN;
  }

  /**
   * Fit the jump distances using a fit to a cumulative histogram.
   *
   * <p>The histogram is fit repeatedly using a mixed population model with increasing number of
   * different molecules. Results are sorted by the diffusion coefficient ascending. This process is
   * stopped when: the Adjusted R^2 does not improve; the fraction of one of the populations is
   * below the min fraction; the difference between two consecutive diffusion coefficients is below
   * the min difference; more than one population is below min D.
   *
   * <p>The number of populations must be obtained from the size of the D/fractions arrays.
   *
   * @param jumpDistances The jump distances (in um^2)
   * @return Array containing: { D (um^2), Fractions }. Can be null if no fit was made.
   */
  public double[][] fitJumpDistances(double... jumpDistances) {
    resetFitResult();
    if (jumpDistances == null || jumpDistances.length == 0) {
      return null;
    }
    final double meanJumpDistance = MathUtils.sum(jumpDistances) / jumpDistances.length;
    if (meanJumpDistance == 0) {
      return null;
    }
    final double[][] jdHistogram = cumulativeHistogram(jumpDistances);
    return fitJumpDistanceHistogram(meanJumpDistance, jdHistogram);
  }

  /**
   * Fit the jump distance histogram using a cumulative sum.
   *
   * <p>The histogram is fit repeatedly using a mixed population model with increasing number of
   * different molecules. Results are sorted by the diffusion coefficient ascending. This process is
   * stopped when: the Adjusted R^2 does not improve; the fraction of one of the populations is
   * below the min fraction; the difference between two consecutive diffusion coefficients is below
   * the min difference; more than one population is below min D.
   *
   * <p>The number of populations must be obtained from the size of the D/fractions arrays.
   *
   * @param meanJumpDistance The mean jump distance (in um^2)
   * @param jdHistogram The cumulative jump distance histogram. X-axis is um^2, Y-axis is cumulative
   *        probability. Must be monototic ascending.
   * @return Array containing: { D (um^2), Fractions }. Can be null if no fit was made.
   */
  public double[][] fitJumpDistanceHistogram(double meanJumpDistance, double[][] jdHistogram) {
    resetFitResult();
    // Guess the D
    final double estimatedD = meanJumpDistance / 4;
    if (meanJumpDistance == 0) {
      return null;
    }
    LoggerUtils.log(logger, Level.INFO, "Estimated D = %s um^2", MathUtils.rounded(estimatedD, 4));

    // We use the adjusted R^2 to pick the best model.

    final double[] fitValue = new double[maxN];
    final double[][] coefficients = new double[maxN][];
    final double[][] fractions = new double[maxN][];
    int best = -1;

    if (minN == 1) {
      final double[][] fit = doFitJumpDistanceHistogram(jdHistogram, estimatedD, 1);
      if (fit != null) {
        coefficients[0] = fit[0];
        fractions[0] = fit[1];
        fitValue[0] = this.fitValue;
        saveFitCurve(fit, jdHistogram);
        best = 0;
      }
    }

    // Fit using a mixed population model.
    // Vary n from 2 to N. Stop when the fit fails or the fit is worse.
    int bestMulti = -1;
    for (int n = Math.max(1, minN - 1); n < maxN; n++) {
      final double[][] fit = doFitJumpDistanceHistogram(jdHistogram, estimatedD, n + 1);
      if (fit == null) {
        break;
      }

      coefficients[n] = fit[0];
      fractions[n] = fit[1];
      fitValue[n] = this.fitValue;

      // Store the best multi-model (if none exists)
      if (bestMulti == -1) {
        bestMulti = n;
      }

      // Store the best model (if none exists)
      if (best == -1) {
        best = n;
        continue;
      }

      // Test this model is better using the adjusted R^2 (should be higher)
      final double diff = fitValue[n] - fitValue[best];

      // Stop if not improving
      if (diff < 0) {
        break;
      }

      best = bestMulti = n;
    }

    // Add the best fit to the plot
    if (bestMulti > -1) {
      saveFitCurve(new double[][] {coefficients[bestMulti], fractions[bestMulti]}, jdHistogram);
    }

    if (best > -1) {
      LoggerUtils.log(logger, Level.INFO, "Best fit achieved using %s: %s, %s = %s",
          TextUtils.pleural(best + 1, "population"), formatD(coefficients[best]),
          TextUtils.pleural(best + 1, "Fraction"), format(fractions[best]));
      lastFitValue = fitValue[best];
      return new double[][] {coefficients[best], fractions[best]};
    }
    return null;
  }

  /**
   * Fit the jump distances using a fit to a cumulative histogram with the given number of species.
   *
   * <p>Results are sorted by the diffusion coefficient ascending.
   *
   * @param jumpDistances The jump distances (in um^2)
   * @param n The number of species in the mixed population
   * @return Array containing: { D (um^2), Fractions }. Can be null if no fit was made.
   */
  public double[][] fitJumpDistances(double[] jumpDistances, int n) {
    resetFitResult();
    if (jumpDistances == null || jumpDistances.length == 0) {
      return null;
    }
    final double meanJumpDistance = MathUtils.sum(jumpDistances) / jumpDistances.length;
    if (meanJumpDistance == 0) {
      return null;
    }
    final double[][] jdHistogram = cumulativeHistogram(jumpDistances);
    return fitJumpDistanceHistogram(meanJumpDistance, jdHistogram, n);
  }

  /**
   * Fit the jump distance histogram using a cumulative sum with the given number of species.
   *
   * <p>Results are sorted by the diffusion coefficient ascending.
   *
   * @param meanJumpDistance The mean jump distance (in um^2)
   * @param jdHistogram The cumulative jump distance histogram. X-axis is um^2, Y-axis is cumulative
   *        probability. Must be Monotonic ascending.
   * @param n The number of species in the mixed population
   * @return Array containing: { D (um^2), Fractions }. Can be null if no fit was made.
   */
  public double[][] fitJumpDistanceHistogram(double meanJumpDistance, double[][] jdHistogram,
      int n) {
    resetFitResult();
    // Guess the D
    final double estimatedD = meanJumpDistance / 4;
    if (meanJumpDistance == 0) {
      return null;
    }
    LoggerUtils.log(logger, Level.INFO, "Estimated D = %s um^2", MathUtils.rounded(estimatedD, 4));

    final double[][] fit = doFitJumpDistanceHistogram(jdHistogram, estimatedD, n);
    if (fit != null) {
      saveFitCurve(fit, jdHistogram);
    }
    return fit;
  }

  /**
   * Fit the jump distance histogram using a cumulative sum with the given number of species.
   *
   * <p>Results are sorted by the diffusion coefficient ascending.
   *
   * @param jdHistogram The cumulative jump distance histogram. X-axis is um^2, Y-axis is cumulative
   *        probability. Must be monototic ascending.
   * @param estimatedD The estimated diffusion coefficient
   * @param n The number of species in the mixed population
   * @return Array containing: { D (um^2), Fractions }. Can be null if no fit was made.
   */
  private double[][] doFitJumpDistanceHistogram(double[][] jdHistogram, double estimatedD, int n) {
    calibrated = isCalibrated();

    if (n == 1) {
      // Fit using a single population model
      final LevenbergMarquardtOptimizer lvmOptimizer = new LevenbergMarquardtOptimizer();
      try {
        final JumpDistanceCumulFunction function =
            new JumpDistanceCumulFunction(jdHistogram[0], jdHistogram[1], estimatedD);

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

        final Optimum lvmSolution = lvmOptimizer.optimize(problem);

        final double[] fitParams = lvmSolution.getPoint().toArray();
        // True for an unweighted fit
        ss = lvmSolution.getResiduals().dotProduct(lvmSolution.getResiduals());
        // ss = calculateSumOfSquares(function.getY(), function.value(fitParams));
        lastFitValue = fitValue = MathUtils.getAdjustedCoefficientOfDetermination(ss,
            MathUtils.getTotalSumOfSquares(function.getY()), function.x.length, 1);
        final double[] coefficients = fitParams;
        final double[] fractions = new double[] {1};

        LoggerUtils.log(logger, Level.INFO,
            "Fit Jump distance (N=1) : %s, SS = %s, Adjusted R^2 = %s (%d evaluations)",
            formatD(fitParams[0]), MathUtils.rounded(ss, 4), MathUtils.rounded(fitValue, 4),
            lvmSolution.getEvaluations());

        return new double[][] {coefficients, fractions};
      } catch (final TooManyIterationsException ex) {
        LoggerUtils.log(logger, Level.INFO,
            "LVM optimiser failed to fit (N=1) : Too many iterations : %s", ex.getMessage());
      } catch (final ConvergenceException ex) {
        LoggerUtils.log(logger, Level.INFO, "LVM optimiser failed to fit (N=1) : %s",
            ex.getMessage());
      }
    }

    // Uses a weighted sum of n exponential functions, each function models a fraction of the
    // particles.
    // An LVM fit cannot restrict the parameters so the fractions do not go below zero.
    // Use the CustomPowell/CMEASOptimizer which supports bounded fitting.

    final MixedJumpDistanceCumulFunctionMultivariate function =
        new MixedJumpDistanceCumulFunctionMultivariate(jdHistogram[0], jdHistogram[1], estimatedD,
            n);

    final double[] lB = function.getLowerBounds();

    int evaluations = 0;
    PointValuePair constrainedSolution = null;

    final MaxEval maxEval = new MaxEval(20000);
    final CustomPowellOptimizer powellOptimizer = createCustomPowellOptimizer();
    try {
      // The Powell algorithm can use more general bounds: 0 - Infinity
      constrainedSolution = powellOptimizer.optimize(maxEval, new ObjectiveFunction(function),
          new InitialGuess(function.guess()),
          new SimpleBounds(lB,
              function.getUpperBounds(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY)),
          new CustomPowellOptimizer.BasisStep(function.step()), GoalType.MINIMIZE);

      evaluations = powellOptimizer.getEvaluations();
      LoggerUtils.log(logger, Level.FINE, "Powell optimiser fit (N=%d) : SS = %f (%d evaluations)",
          n, constrainedSolution.getValue(), evaluations);
    } catch (final TooManyEvaluationsException ex) {
      LoggerUtils.log(logger, Level.INFO,
          "Powell optimiser failed to fit (N=%d) : Too many evaluations (%d)", n,
          powellOptimizer.getEvaluations());
    } catch (final TooManyIterationsException ex) {
      LoggerUtils.log(logger, Level.INFO,
          "Powell optimiser failed to fit (N=%d) : Too many iterations (%d)", n,
          powellOptimizer.getIterations());
    } catch (final ConvergenceException ex) {
      LoggerUtils.log(logger, Level.INFO, "Powell optimiser failed to fit (N=%d) : %s", n,
          ex.getMessage());
    }

    if (constrainedSolution == null) {
      LoggerUtils.log(logger, Level.INFO, "Trying CMAES optimiser with restarts ...");

      final double[] uB = function.getUpperBounds();
      final SimpleBounds bounds = new SimpleBounds(lB, uB);

      // The sigma determines the search range for the variables. It should be 1/3 of the initial
      // search region.
      final double[] s = new double[lB.length];
      for (int i = 0; i < s.length; i++) {
        s[i] = (uB[i] - lB[i]) / 3;
      }
      final OptimizationData sigma = new CMAESOptimizer.Sigma(s);
      final OptimizationData popSize = new CMAESOptimizer.PopulationSize(
          (int) (4 + Math.floor(3 * Math.log(function.x.length))));

      // Iterate this for stability in the initial guess
      final CMAESOptimizer cmaesOptimizer = createCMAESOptimizer();

      for (int i = 0; i <= fitRestarts; i++) {
        // Try from the initial guess
        try {
          final PointValuePair solution = cmaesOptimizer.optimize(
              new InitialGuess(function.guess()), new ObjectiveFunction(function),
              GoalType.MINIMIZE, bounds, sigma, popSize, maxEval);
          if (constrainedSolution == null || solution.getValue() < constrainedSolution.getValue()) {
            evaluations = cmaesOptimizer.getEvaluations();
            constrainedSolution = solution;
            LoggerUtils.log(logger, Level.FINE,
                "CMAES optimiser [%da] fit (N=%d) : SS = %f (%d evaluations)", i, n,
                solution.getValue(), evaluations);
          }
        } catch (final TooManyEvaluationsException ex) {
          // No solution
        }

        if (constrainedSolution == null) {
          continue;
        }

        // Try from the current optimum
        try {
          final PointValuePair solution = cmaesOptimizer.optimize(
              new InitialGuess(constrainedSolution.getPointRef()), new ObjectiveFunction(function),
              GoalType.MINIMIZE, bounds, sigma, popSize, maxEval);
          if (solution.getValue() < constrainedSolution.getValue()) {
            evaluations = cmaesOptimizer.getEvaluations();
            constrainedSolution = solution;
            LoggerUtils.log(logger, Level.FINE,
                "CMAES optimiser [%db] fit (N=%d) : SS = %f (%d evaluations)", i, n,
                solution.getValue(), evaluations);
          }
        } catch (final TooManyEvaluationsException ex) {
          // No solution
        }
      }

      if (constrainedSolution != null) {
        // Re-optimise with Powell?
        try {
          final PointValuePair solution = powellOptimizer.optimize(maxEval,
              new ObjectiveFunction(function), new InitialGuess(constrainedSolution.getPointRef()),
              new SimpleBounds(lB,
                  function.getUpperBounds(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY)),
              new CustomPowellOptimizer.BasisStep(function.step()), GoalType.MINIMIZE);
          if (solution.getValue() < constrainedSolution.getValue()) {
            evaluations = cmaesOptimizer.getEvaluations();
            constrainedSolution = solution;
            LoggerUtils.log(logger, Level.INFO,
                "Powell optimiser re-fit (N=%d) : SS = %f (%d evaluations)", n,
                constrainedSolution.getValue(), evaluations);
          }
        } catch (final TooManyEvaluationsException ex) {
          // No solution
        } catch (final TooManyIterationsException ex) {
          // No solution
        } catch (final ConvergenceException ex) {
          // No solution
        }
      }
    }

    if (constrainedSolution == null) {
      LoggerUtils.log(logger, Level.INFO, "Failed to fit N=%d", n);
      return null;
    }

    double[] fitParams = constrainedSolution.getPointRef();
    ss = constrainedSolution.getValue();

    // TODO - Try a bounded BFGS optimiser

    // Try and improve using a LVM fit
    final MixedJumpDistanceCumulFunctionGradient functionGradient =
        new MixedJumpDistanceCumulFunctionGradient(jdHistogram[0], jdHistogram[1], estimatedD, n);

    Optimum lvmSolution;
    final LevenbergMarquardtOptimizer lvmOptimizer = new LevenbergMarquardtOptimizer();
    try {
      //@formatter:off
      final LeastSquaresProblem problem = new LeastSquaresBuilder()
          .maxEvaluations(Integer.MAX_VALUE)
          .maxIterations(3000)
          .start(fitParams)
          .target(functionGradient.getY())
          .weight(new DiagonalMatrix(functionGradient.getWeights()))
          .model(functionGradient, new MultivariateMatrixFunction() {
            @Override
            public double[][] value(double[] point) throws IllegalArgumentException
            {
              return functionGradient.jacobian(point);
            }} )
          .build();
      //@formatter:on

      lvmSolution = lvmOptimizer.optimize(problem);

      // True for an unweighted fit
      final double ss = lvmSolution.getResiduals().dotProduct(lvmSolution.getResiduals());
      // double ss = calculateSumOfSquares(functionGradient.getY(),
      // functionGradient.value(lvmSolution.getPoint().toArray()));

      // All fitted parameters must be above zero
      if (ss < this.ss && MathUtils.min(lvmSolution.getPoint().toArray()) > 0) {
        LoggerUtils.log(logger, Level.INFO, "  Re-fitting improved the SS from %s to %s (-%s%%)",
            MathUtils.rounded(this.ss, 4), MathUtils.rounded(ss, 4),
            MathUtils.rounded(100 * (this.ss - ss) / this.ss, 4));
        fitParams = lvmSolution.getPoint().toArray();
        this.ss = ss;
        evaluations += lvmSolution.getEvaluations();
      }
    } catch (final TooManyIterationsException ex) {
      LoggerUtils.log(logger, Level.WARNING, "Failed to re-fit : Too many iterations : %s",
          ex.getMessage());
    } catch (final ConvergenceException ex) {
      LoggerUtils.log(logger, Level.WARNING, "Failed to re-fit : %s", ex.getMessage());
    }

    // Since the fractions must sum to one we subtract 1 degree of freedom from the number of
    // parameters
    fitValue = MathUtils.getAdjustedCoefficientOfDetermination(ss,
        MathUtils.getTotalSumOfSquares(function.getY()), function.x.length, fitParams.length - 1);

    final double[] d = new double[n];
    final double[] f = new double[n];
    double sum = 0;
    for (int i = 0; i < d.length; i++) {
      f[i] = fitParams[i * 2];
      sum += f[i];
      d[i] = fitParams[i * 2 + 1];
    }
    for (int i = 0; i < f.length; i++) {
      f[i] /= sum;
    }
    // Sort by coefficient size
    sort(d, f);
    final double[] coefficients = d;
    final double[] fractions = f;

    LoggerUtils.log(logger, Level.INFO,
        "Fit Jump distance (N=%d) : %s (%s), SS = %s, Adjusted R^2 = %s (%d evaluations)", n,
        formatD(d), format(f), MathUtils.rounded(ss, 4), MathUtils.rounded(fitValue, 4),
        evaluations);

    if (isValid(d, f)) {
      lastFitValue = fitValue;
      return new double[][] {coefficients, fractions};
    }

    return null;
  }

  private boolean isValid(double[] d, double[] f) {
    int belowMinD = 0;
    for (int i = 0; i < f.length; i++) {
      // Check only one population has a diffusion coefficient below the
      // precision of the experiment
      if (d[i] < minD) {
        if (++belowMinD > 1) {
          LoggerUtils.log(logger, Level.INFO,
              "  Invalid: Multiple populations below minimum D (%s um^2)", MathUtils.rounded(minD));
          return false;
        }
      }
      // Check the fractions and coefficients exist
      if (f[i] <= 0) {
        LoggerUtils.log(logger, Level.INFO, "  Invalid: Fraction is zero");
        return false;
      }
      if (d[i] <= 0) {
        LoggerUtils.log(logger, Level.INFO, "  Invalid: Coefficient is zero");
        return false;
      }
      // Check the fit has fractions above the minimum fraction
      if (f[i] < minFraction) {
        LoggerUtils.log(logger, Level.INFO,
            "  Invalid: Fraction is less than the minimum fraction: %s < %s",
            MathUtils.rounded(f[i]), MathUtils.rounded(minFraction));
        return false;
      }
      // Check the coefficients are different
      if (i > 0 && d[i - 1] / d[i] < minDifference) {
        LoggerUtils.log(logger, Level.INFO,
            "  Invalid: Coefficients are not different: %s / %s = %s < %s",
            MathUtils.rounded(d[i - 1]), MathUtils.rounded(d[i]),
            MathUtils.rounded(d[i - 1] / d[i]), MathUtils.rounded(minDifference));
        return false;
      }
    }
    return true;
  }

  /**
   * Fit the jump distances using a maximum likelihood estimation.
   *
   * <p>The data is fit repeatedly using a mixed population model with increasing number of
   * different molecules. Results are sorted by the diffusion coefficient ascending. This process is
   * stopped when: the log-likelihood ratio does not improve; the fraction of one of the populations
   * is below the min fraction; the difference between two consecutive diffusion coefficients is
   * below the min difference; more than one population is below min D.
   *
   * <p>The number of populations must be obtained from the size of the D/fractions arrays.
   *
   * @param jumpDistances The jump distances (in um^2)
   * @return Array containing: { D (um^2), Fractions }. Can be null if no fit was made.
   */
  public double[][] fitJumpDistancesMLE(double[] jumpDistances) {
    return fitJumpDistancesMLE(jumpDistances, null);
  }

  /**
   * Fit the jump distances using a maximum likelihood estimation.
   *
   * <p>The data is fit repeatedly using a mixed population model with increasing number of
   * different molecules. Results are sorted by the diffusion coefficient ascending. This process is
   * stopped when: the log-likelihood ratio does not improve; the fraction of one of the populations
   * is below the min fraction; the difference between two consecutive diffusion coefficients is
   * below the min difference; more than one population is below min D.
   *
   * <p>The number of populations must be obtained from the size of the D/fractions arrays.
   *
   * @param jumpDistances The jump distances (in um^2)
   * @param jdHistogram The jump distance histogram for the given distances. If null will be
   *        computed using {@link #cumulativeHistogram(double[])}. Only used if the CurveLogger is
   *        not null.
   * @return Array containing: { D (um^2), Fractions }. Can be null if no fit was made.
   */
  public double[][] fitJumpDistancesMLE(double[] jumpDistances, double[][] jdHistogram) {
    resetFitResult();
    if (jumpDistances == null || jumpDistances.length == 0) {
      return null;
    }
    final double meanJumpDistance = MathUtils.sum(jumpDistances) / jumpDistances.length;
    if (meanJumpDistance == 0) {
      return null;
    }

    // Guess the D
    final double estimatedD = meanJumpDistance / 4;
    LoggerUtils.log(logger, Level.INFO, "Estimated D = %s um^2", MathUtils.rounded(estimatedD, 4));

    // Used for saving fitted the curve
    if (curveLogger != null && jdHistogram == null) {
      jdHistogram = cumulativeHistogram(jumpDistances);
    }

    // When performing MLE we can use the Log-Likelihood Ratio (LLR) to do a significance
    // test that the model has improved.

    final double[] fitValue = new double[maxN];
    final double[] ll = new double[maxN];
    Arrays.fill(ll, Double.NaN);
    final double[][] coefficients = new double[maxN][];
    final double[][] fractions = new double[maxN][];
    int best = -1;

    if (minN == 1) {
      final double[][] fit = doFitJumpDistancesMLE(jumpDistances, estimatedD, 1);
      if (fit != null) {
        coefficients[0] = fit[0];
        fractions[0] = fit[1];
        fitValue[0] = this.fitValue;
        ll[0] = this.ll;
        saveFitCurve(fit, jdHistogram);
        best = 0;
      }
    }

    // Fit using a mixed population model.
    // Vary n from 2 to N. Stop when the fit fails or the fit is worse.
    int bestMulti = -1;
    for (int n = Math.max(1, minN - 1); n < maxN; n++) {
      final double[][] fit = doFitJumpDistancesMLE(jumpDistances, estimatedD, n + 1);
      if (fit == null) {
        break;
      }

      coefficients[n] = fit[0];
      fractions[n] = fit[1];
      fitValue[n] = this.fitValue;
      ll[n] = this.ll;

      // Store the best multi-model (if none exists)
      if (bestMulti == -1) {
        bestMulti = n;
      }

      // Store the best model (if none exists)
      if (best == -1) {
        best = n;
        continue;
      }

      // Test this model is better using a log-likelihood ratio test
      final double llr = 2 * (ll[n] - ll[best]);

      // Stop if not improving
      if (llr < 0) {
        break;
      }

      // The difference in the number of fitted parameters will be 2:
      // i.e. 1 extra diffusion coefficient and population fraction
      final double q = ChiSquaredDistributionTable.computeQValue(llr, 2);
      final boolean reject = (q > significanceLevel);
      LoggerUtils.log(logger, Level.INFO,
          "Fit Jump distance (N=%d -> %d) : MLE = %s -> %s, LLR = %s, q-value = %s (Reject=%b)",
          best + 1, n + 1, MathUtils.rounded(ll[best], 4), MathUtils.rounded(ll[n], 4),
          MathUtils.rounded(llr, 4), MathUtils.rounded(q, 4), reject);
      if (reject) {
        break;
      }

      best = bestMulti = n;
    }

    // Add the best fit to the plot
    if (bestMulti > -1) {
      saveFitCurve(new double[][] {coefficients[bestMulti], fractions[bestMulti]}, jdHistogram);
    }

    if (best > -1) {
      LoggerUtils.log(logger, Level.INFO, "Best fit achieved using %s: %s, %s = %s",
          TextUtils.pleural(best + 1, "population"), formatD(coefficients[best]),
          TextUtils.pleural(best + 1, "Fraction"), format(fractions[best]));
      lastFitValue = fitValue[best];
      return new double[][] {coefficients[best], fractions[best]};
    }

    return null;
  }

  /**
   * Fit the jump distances using a maximum likelihood estimation with the given number of species.
   *
   * <p>Results are sorted by the diffusion coefficient ascending.
   *
   * @param jumpDistances The jump distances (in um^2)
   * @param n The number of species in the mixed population
   * @return Array containing: { D (um^2), Fractions }. Can be null if no fit was made.
   */
  public double[][] fitJumpDistancesMLE(double[] jumpDistances, int n) {
    return fitJumpDistancesMLE(jumpDistances, null, n);
  }

  /**
   * Fit the jump distances using a maximum likelihood estimation with the given number of species.
   *
   * <p>Results are sorted by the diffusion coefficient ascending.
   *
   * @param jumpDistances The jump distances (in um^2)
   * @param jdHistogram The jump distance histogram for the given distances. If null will be
   *        computed using {@link #cumulativeHistogram(double[])}. Only used if the CurveLogger is
   *        not null.
   * @param n The number of species in the mixed population
   * @return Array containing: { D (um^2), Fractions }. Can be null if no fit was made.
   */
  public double[][] fitJumpDistancesMLE(double[] jumpDistances, double[][] jdHistogram, int n) {
    resetFitResult();
    if (jumpDistances == null || jumpDistances.length == 0) {
      return null;
    }
    final double meanJumpDistance = MathUtils.sum(jumpDistances) / jumpDistances.length;
    if (meanJumpDistance == 0) {
      return null;
    }

    // Guess the D
    final double estimatedD = meanJumpDistance / 4;
    LoggerUtils.log(logger, Level.INFO, "Estimated D = %s um^2", MathUtils.rounded(estimatedD, 4));

    // Used for saving fitted the curve
    if (curveLogger != null && jdHistogram == null) {
      jdHistogram = cumulativeHistogram(jumpDistances);
    }

    final double[][] fit = doFitJumpDistancesMLE(jumpDistances, estimatedD, n);
    if (fit != null) {
      saveFitCurve(fit, jdHistogram);
    }
    return fit;
  }

  /**
   * Fit the jump distances using a maximum likelihood estimation with the given number of species.
   * | *
   *
   * <p>Results are sorted by the diffusion coefficient ascending.
   *
   * @param jumpDistances The jump distances (in um^2)
   * @param estimatedD The estimated diffusion coefficient
   * @param n The number of species in the mixed population
   * @return Array containing: { D (um^2), Fractions }. Can be null if no fit was made.
   */
  private double[][] doFitJumpDistancesMLE(double[] jumpDistances, double estimatedD, int n) {
    final MaxEval maxEval = new MaxEval(20000);
    final CustomPowellOptimizer powellOptimizer = createCustomPowellOptimizer();
    calibrated = isCalibrated();

    if (n == 1) {
      try {
        final JumpDistanceFunction function = new JumpDistanceFunction(jumpDistances, estimatedD);
        // The Powell algorithm can use more general bounds: 0 - Infinity
        final SimpleBounds bounds = new SimpleBounds(function.getLowerBounds(),
            function.getUpperBounds(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY));
        final PointValuePair solution = powellOptimizer.optimize(maxEval,
            new ObjectiveFunction(function), new InitialGuess(function.guess()), bounds,
            new CustomPowellOptimizer.BasisStep(function.step()), GoalType.MAXIMIZE);

        final double[] fitParams = solution.getPointRef();
        ll = solution.getValue();
        lastFitValue =
            fitValue = MathUtils.getAkaikeInformationCriterion(ll, jumpDistances.length, 1);
        final double[] coefficients = fitParams;
        final double[] fractions = new double[] {1};

        LoggerUtils.log(logger, Level.INFO,
            "Fit Jump distance (N=1) : %s, MLE = %s, Akaike IC = %s (%d evaluations)",
            formatD(fitParams[0]), MathUtils.rounded(ll, 4), MathUtils.rounded(fitValue, 4),
            powellOptimizer.getEvaluations());

        return new double[][] {coefficients, fractions};
      } catch (final TooManyEvaluationsException ex) {
        LoggerUtils.log(logger, Level.INFO,
            "Powell optimiser failed to fit (N=1) : Too many evaluation (%d)",
            powellOptimizer.getEvaluations());
      } catch (final TooManyIterationsException ex) {
        LoggerUtils.log(logger, Level.INFO,
            "Powell optimiser failed to fit (N=1) : Too many iterations (%d)",
            powellOptimizer.getIterations());
      } catch (final ConvergenceException ex) {
        LoggerUtils.log(logger, Level.INFO, "Powell optimiser failed to fit (N=1) : %s",
            ex.getMessage());
      }

      return null;
    }

    final MixedJumpDistanceFunction function =
        new MixedJumpDistanceFunction(jumpDistances, estimatedD, n);

    final double[] lB = function.getLowerBounds();

    int evaluations = 0;
    PointValuePair constrainedSolution = null;

    try {
      // The Powell algorithm can use more general bounds: 0 - Infinity
      constrainedSolution = powellOptimizer.optimize(maxEval, new ObjectiveFunction(function),
          new InitialGuess(function.guess()),
          new SimpleBounds(lB,
              function.getUpperBounds(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY)),
          new CustomPowellOptimizer.BasisStep(function.step()), GoalType.MAXIMIZE);

      evaluations = powellOptimizer.getEvaluations();
      LoggerUtils.log(logger, Level.FINE, "Powell optimiser fit (N=%d) : MLE = %f (%d evaluations)",
          n, constrainedSolution.getValue(), powellOptimizer.getEvaluations());
    } catch (final TooManyEvaluationsException ex) {
      LoggerUtils.log(logger, Level.INFO,
          "Powell optimiser failed to fit (N=%d) : Too many evaluation (%d)", n,
          powellOptimizer.getEvaluations());
    } catch (final TooManyIterationsException ex) {
      LoggerUtils.log(logger, Level.INFO,
          "Powell optimiser failed to fit (N=%d) : Too many iterations (%d)", n,
          powellOptimizer.getIterations());
    } catch (final ConvergenceException ex) {
      LoggerUtils.log(logger, Level.INFO, "Powell optimiser failed to fit (N=%d) : %s", n,
          ex.getMessage());
    }

    if (constrainedSolution == null) {
      LoggerUtils.log(logger, Level.INFO, "Trying CMAES optimiser with restarts ...");

      final double[] uB = function.getUpperBounds();
      final SimpleBounds bounds = new SimpleBounds(lB, uB);

      // Try a bounded CMAES optimiser since the Powell optimiser appears to be
      // sensitive to the order of the parameters. It is not good when the fast particle
      // is the minority fraction. Could this be due to too low an upper bound?

      // The sigma determines the search range for the variables. It should be 1/3 of the initial
      // search region.
      final double[] s = new double[lB.length];
      for (int i = 0; i < s.length; i++) {
        s[i] = (uB[i] - lB[i]) / 3;
      }
      final OptimizationData sigma = new CMAESOptimizer.Sigma(s);
      final OptimizationData popSize = new CMAESOptimizer.PopulationSize(
          (int) (4 + Math.floor(3 * Math.log(function.x.length))));

      // Iterate this for stability in the initial guess
      final CMAESOptimizer cmaesOptimizer = createCMAESOptimizer();

      for (int i = 0; i <= fitRestarts; i++) {
        // Try from the initial guess
        try {
          final PointValuePair solution = cmaesOptimizer.optimize(
              new InitialGuess(function.guess()), new ObjectiveFunction(function),
              GoalType.MAXIMIZE, bounds, sigma, popSize, maxEval);
          if (constrainedSolution == null || solution.getValue() > constrainedSolution.getValue()) {
            evaluations = cmaesOptimizer.getEvaluations();
            constrainedSolution = solution;
            LoggerUtils.log(logger, Level.FINE,
                "CMAES optimiser [%da] fit (N=%d) : MLE = %f (%d evaluations)", i, n,
                solution.getValue(), evaluations);
          }
        } catch (final TooManyEvaluationsException ex) {
          // No solution
        }

        if (constrainedSolution == null) {
          continue;
        }

        // Try from the current optimum
        try {
          final PointValuePair solution = cmaesOptimizer.optimize(
              new InitialGuess(constrainedSolution.getPointRef()), new ObjectiveFunction(function),
              GoalType.MAXIMIZE, bounds, sigma, popSize, maxEval);
          if (solution.getValue() > constrainedSolution.getValue()) {
            evaluations = cmaesOptimizer.getEvaluations();
            constrainedSolution = solution;
            LoggerUtils.log(logger, Level.FINE,
                "CMAES optimiser [%db] fit (N=%d) : MLE = %f (%d evaluations)", i, n,
                solution.getValue(), evaluations);
          }
        } catch (final TooManyEvaluationsException ex) {
          // No solution
        }
      }

      if (constrainedSolution != null) {
        try {
          // Re-optimise with Powell?
          final PointValuePair solution = powellOptimizer.optimize(maxEval,
              new ObjectiveFunction(function), new InitialGuess(constrainedSolution.getPointRef()),
              new SimpleBounds(lB,
                  function.getUpperBounds(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY)),
              new CustomPowellOptimizer.BasisStep(function.step()), GoalType.MAXIMIZE);
          if (solution.getValue() > constrainedSolution.getValue()) {
            evaluations = cmaesOptimizer.getEvaluations();
            constrainedSolution = solution;
            LoggerUtils.log(logger, Level.INFO,
                "Powell optimiser re-fit (N=%d) : MLE = %f (%d evaluations)", n,
                constrainedSolution.getValue(), powellOptimizer.getEvaluations());
          }
        } catch (final TooManyEvaluationsException ex) {
          // No solution
        } catch (final TooManyIterationsException ex) {
          // No solution
        } catch (final ConvergenceException ex) {
          // No solution
        }
      }
    }

    if (constrainedSolution == null) {
      LoggerUtils.log(logger, Level.INFO, "Failed to fit N=%d", n);
      return null;
    }

    final double[] fitParams = constrainedSolution.getPointRef();
    ll = constrainedSolution.getValue();

    // Since the fractions must sum to one we subtract 1 degree of freedom from the number of
    // parameters
    fitValue =
        MathUtils.getAkaikeInformationCriterion(ll, jumpDistances.length, fitParams.length - 1);

    final double[] d = new double[n];
    final double[] f = new double[n];
    double sum = 0;
    for (int i = 0; i < d.length; i++) {
      f[i] = fitParams[i * 2];
      sum += f[i];
      d[i] = fitParams[i * 2 + 1];
    }
    for (int i = 0; i < f.length; i++) {
      f[i] /= sum;
    }
    // Sort by coefficient size
    sort(d, f);
    final double[] coefficients = d;
    final double[] fractions = f;

    LoggerUtils.log(logger, Level.INFO,
        "Fit Jump distance (N=%d) : %s (%s), MLE = %s, Akaike IC = %s (%d evaluations)", n,
        formatD(d), format(f), MathUtils.rounded(ll, 4), MathUtils.rounded(fitValue, 4),
        evaluations);

    if (isValid(d, f)) {
      lastFitValue = fitValue;
      return new double[][] {coefficients, fractions};
    }

    return null;
  }

  private static CustomPowellOptimizer createCustomPowellOptimizer() {
    final double rel = 1e-8;
    final double abs = 1e-10;
    // double lineRel = rel;
    // double lineAbs = abs;
    final ConvergenceChecker<PointValuePair> positionChecker = null;
    // new org.apache.commons.math3.optim.PositionChecker(1e-3, 1e-10);
    final boolean basisConvergence = false;

    // return new CustomPowellOptimizer(rel, abs, lineRel, lineAbs, checker, basisConvergence);
    return new CustomPowellOptimizer(rel, abs, positionChecker, basisConvergence);
  }

  private static CMAESOptimizer createCMAESOptimizer() {
    final double rel = 1e-8;
    final double abs = 1e-10;
    final int maxIterations = 2000;
    final double stopFitness = 0; // Double.NEGATIVE_INFINITY;
    final boolean isActiveCMA = true;
    final int diagonalOnly = 20;
    final int checkFeasableCount = 1;
    final RandomGenerator random = new Well19937c();
    final boolean generateStatistics = false;
    final ConvergenceChecker<PointValuePair> checker = new SimpleValueChecker(rel, abs);

    // Iterate this for stability in the initial guess
    return new CMAESOptimizer(maxIterations, stopFitness, isActiveCMA, diagonalOnly,
        checkFeasableCount, random, generateStatistics, checker);
  }

  /**
   * Format the diffusion coefficients for reporting using the calibration if present.
   *
   * @param jumpD the jump Distances
   * @return The formatted D
   */
  private String formatD(double... jumpD) {
    if (jumpD == null || jumpD.length == 0) {
      return "";
    }
    final StringBuilder sb = new StringBuilder("D = ");
    for (int i = 0; i < jumpD.length; i++) {
      if (i != 0) {
        sb.append(", ");
      }
      sb.append(MathUtils.rounded(jumpD[i], 4));
    }
    sb.append(" um^2");
    if (calibrated) {
      jumpD = calculateApparentDiffusionCoefficient(jumpD);
      sb.append(", D* = ");
      for (int i = 0; i < jumpD.length; i++) {
        if (i != 0) {
          sb.append(", ");
        }
        sb.append(MathUtils.rounded(jumpD[i], 4));
      }
      sb.append(" um^2/s");
    }
    return sb.toString();
  }

  private static String format(double[] data) {
    if (data == null || data.length == 0) {
      return "";
    }
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < data.length; i++) {
      if (i != 0) {
        sb.append(", ");
      }
      sb.append(MathUtils.rounded(data[i], 4));
    }
    return sb.toString();
  }

  /**
   * Sort the arrays by the size of the diffusion coefficient.
   *
   * @param d The diffusion coefficient array
   * @param f The fraction of the population array
   */
  public static void sort(double[] d, double[] f) {
    // Sort by coefficient size
    SortUtils.sortData(f, d, true, false);
  }

  private void saveFitCurve(double[][] fit, double[][] jdHistogram) {
    if (fit[0].length == 1) {
      saveFitCurve(fit[0], jdHistogram);
    } else {
      final double[] params = new double[fit[0].length * 2];
      for (int i = 0; i < fit[0].length; i++) {
        params[i * 2] = fit[1][i];
        params[i * 2 + 1] = fit[0][i];
      }
      saveFitCurve(params, jdHistogram);
    }
  }

  private void saveFitCurve(double[] params, double[][] jdHistogram) {
    if (curveLogger == null) {
      return;
    }
    final int npoints = curveLogger.getNumberOfCurvePoints();
    if (npoints <= 1) {
      return;
    }
    Function function;
    if (params.length == 1) {
      function = new JumpDistanceCumulFunction(null, null, 0);
    } else {
      function = new MixedJumpDistanceCumulFunction(null, null, 0, params.length / 2);
    }

    final double max = jdHistogram[0][jdHistogram[0].length - 1];
    final double interval = max / npoints;
    final double[] x = new double[npoints + 1];
    final double[] y = new double[npoints + 1];

    for (int i = 0; i < npoints; i++) {
      x[i] = i * interval;
      y[i] = function.evaluate(x[i], params);
    }
    x[npoints] = max;
    y[npoints] = function.evaluate(max, params);

    if (params.length == 1) {
      curveLogger.saveSinglePopulationCurve(new double[][] {x, y});
    } else {
      curveLogger.saveMixedPopulationCurve(new double[][] {x, y});
    }
  }

  @SuppressWarnings("unused")
  private static double calculateSumOfSquares(double[] obs, double[] exp) {
    double ss = 0;
    for (int i = 0; i < obs.length; i++) {
      ss += (obs[i] - exp[i]) * (obs[i] - exp[i]);
    }
    return ss;
  }

  /**
   * Function used for least-squares fitting of cumulative histogram of jump distances.
   */
  private abstract static class Function {
    double[] x;
    double[] y;
    double estimatedD;
    int n;

    public Function(double[] x, double[] y, double estimatedD, int n) {
      this.x = x;
      this.y = y;
      this.estimatedD = estimatedD;
      this.n = n;
    }

    /**
     * @return An estimate for the parameters.
     */
    public double[] guess() {
      if (n == 1) {
        return new double[] {estimatedD};
      }

      final double[] guess = new double[n * 2];
      double d = estimatedD;
      for (int i = 0; i < n; i++) {
        // Fraction are all equal
        guess[i * 2] = 1;
        // Diffusion coefficient gets smaller for each fraction
        guess[i * 2 + 1] = d;
        d *= 0.1;
      }
      return guess;
    }

    /**
     * @return An estimate for the initial search step for the parameters.
     */
    public double[] step() {
      if (n == 1) {
        return new double[] {estimatedD * 0.5};
      }

      final double[] step = new double[n * 2];
      double d = estimatedD;
      for (int i = 0; i < n; i++) {
        // Fraction are all equal
        step[i * 2] = 0.1;
        // Diffusion coefficient gets smaller for each fraction
        step[i * 2 + 1] = d * 0.5;
        d *= 0.1;
      }
      return step;
    }

    public double[] getUpperBounds() {
      // Fraction guess is 1 so set the upper limit as 10
      // Diffusion coefficient could be 10x the estimated
      return getUpperBounds(10, estimatedD * 10);
    }

    public double[] getUpperBounds(double fractionLimit, double dLimit) {
      if (n == 1) {
        return new double[] {dLimit};
      }

      final double[] bounds = new double[n * 2];
      for (int i = 0; i < n; i++) {
        bounds[i * 2] = fractionLimit;
        bounds[i * 2 + 1] = dLimit;
      }
      return bounds;
    }

    public double[] getLowerBounds() {
      return getLowerBounds(0, 0);
    }

    public double[] getLowerBounds(double fractionLimit, double dLimit) {
      // Diffusion coefficient could be 0 but this is not practical for
      // testing a mixed population so set to a small value where the optimiser
      // will not be successfully anyway
      fractionLimit = Math.max(fractionLimit, 1e-5);
      dLimit = Math.max(dLimit, 1e-16);

      if (n == 1) {
        return new double[] {dLimit};
      }

      final double[] bounds = new double[n * 2];
      for (int i = 0; i < n; i++) {
        bounds[i * 2] = fractionLimit;
        bounds[i * 2 + 1] = dLimit;
      }
      return bounds;
    }

    public double[] getWeights() {
      final double[] w = new double[x.length];
      Arrays.fill(w, 1);
      return w;
    }

    public double[] getX() {
      return x;
    }

    public double[] getY() {
      return y;
    }

    public abstract double evaluate(double x, double[] parameters);

    public double[][] jacobian(double[] variables) {
      final double[][] jacobian = new double[x.length][variables.length];

      final double delta = 0.001;
      final double[] d = new double[variables.length];
      final double[][] variables2 = new double[variables.length][];
      for (int i = 0; i < variables.length; i++) {
        d[i] = delta * Math.abs(variables[i]); // Should the delta be changed for each parameter ?
        variables2[i] = Arrays.copyOf(variables, variables.length);
        variables2[i][i] += d[i];
      }
      for (int i = 0; i < jacobian.length; ++i) {
        final double value = evaluate(x[i], variables);
        for (int j = 0; j < variables.length; j++) {
          final double value2 = evaluate(x[i], variables2[j]);
          jacobian[i][j] = (value2 - value) / d[j];
        }
      }
      return jacobian;
    }
  }

  /**
   * Compute the probability of mean-squared distance x given a diffusion coefficient.
   *
   * <p>Function used for maximum likelihood fitting.
   */
  static class JumpDistanceFunction extends Function implements MultivariateFunction {
    /**
     * Instantiates a new jump distance function.
     *
     * @param x the x
     * @param estimatedD the estimated D
     */
    public JumpDistanceFunction(double[] x, double estimatedD) {
      super(x, null, estimatedD, 1);
    }

    // Adapted from http://commons.apache.org/proper/commons-math/userguide/optimization.html

    @Override
    public double evaluate(double x, double[] params) {
      // Compute the probability:
      // p = 1/4D * exp(-x/4D)
      // final double fourD = 4 * getD(params[0]);
      final double fourD = 4 * params[0];
      return FastMath.exp(-x / fourD) / fourD;
    }

    /**
     * Evaluate all.
     *
     * @param params the params
     * @return the values
     */
    public double[] evaluateAll(double[] params) {
      // Compute the probability:
      // p = 1/4D * exp(-x/4D)
      final double[] values = new double[x.length];
      // final double one_fourD = 1 / (4 * getD(params[0]));
      final double one_fourD = 1 / (4 * params[0]);
      for (int i = 0; i < values.length; i++) {
        values[i] = one_fourD * FastMath.exp(-x[i] * one_fourD);
      }
      return values;
    }

    /** {@inheritDoc} */
    @Override
    public double value(double[] variables) {
      // Compute the log-likelihood:
      // log(p) = log(1/4D * exp(-x/4D))
      // = log(1/4D) + log(exp(-x/4D))
      // = log(1/4D) + -x/4D
      double l = 0;
      // final double one_fourD = 1 / (4 * getD(variables[0]));
      final double one_fourD = 1 / (4 * variables[0]);
      for (int i = 0; i < x.length; i++) {
        l += -x[i] * one_fourD;
      }
      l += Math.log(one_fourD) * x.length;
      // Debug the call from the optimiser
      if (DEBUG_OPTIMISER) {
        System.out.printf("[1] : [%f] = %f\n", variables[0], l);
      }
      return l;
    }

    // This has not been tested. It could be used for LVM fitting of the p-values. However MLE
    // is less sensitive to outliers of p-values.
    // public double[][] jacobian(double[] variables)
    // {
    // // Compute the gradients using calculus differentiation:
    // // y = 1/4D * exp(-x/4D)
    // // y = aa * a
    // // aa = 1/4D
    // // a = exp(b)
    // // b = -x / 4D
    // //
    // // y' = aa' * a + aa * a'
    // // aa' = -1/4D^2
    // // a' = exp(b) * b'
    // // b' = -1 * -x / 4D^2 = x / 4D^2
    // // y' = -1/4D^2 * exp(-x/4D) + 1/4D * exp(-x/4D) * x / 4D^2
    // // = 1/4D * exp(-x/4D) * (-1/D + x / 4D^2)
    //
    // final double d = variables[0];
    // final double fourD = 4 * d;
    // final double aa = 1 / fourD;
    // final double cc = aa / d;
    // final double c = -1 / d;
    // double[][] jacobian = new double[x.length][variables.length];
    //
    // for (int i = 0; i < jacobian.length; ++i)
    // {
    // jacobian[i][0] = aa * FastMath.exp(-x[i] * aa) * (c + x[i] * cc);
    // }
    //
    // //// Check numerically ...
    // //double[][] jacobian2 = super.jacobian(variables);
    // //for (int i = 0; i < jacobian.length; i++)
    // //{
    // // System.out.printf("dD = %g : %g = %g\n", jacobian[i][0], jacobian2[i][0],
    // // DoubleEquality.relativeError(jacobian[i][0], jacobian2[i][0]));
    // //}
    //
    // return jacobian;
    // }
  }

  /**
   * Compute the probability of mean-squared distance being within x given a diffusion coefficient.
   *
   * <p>Function used for least-squares fitting of cumulative histogram of jump distances.
   */
  static class JumpDistanceCumulFunction extends Function implements MultivariateVectorFunction {
    /**
     * Instantiates a new jump distance cumul function.
     *
     * @param x the x
     * @param y the y
     * @param estimatedD the estimated D
     */
    public JumpDistanceCumulFunction(double[] x, double[] y, double estimatedD) {
      super(x, y, estimatedD, 1);
    }

    // Adapted from http://commons.apache.org/proper/commons-math/userguide/optimization.html

    @Override
    public double evaluate(double x, double[] params) {
      return 1 - FastMath.exp(-x / (4 * params[0]));
    }

    /** {@inheritDoc} */
    @Override
    public double[] value(double[] variables) {
      final double[] values = new double[x.length];
      final double fourD = 4 * variables[0];
      for (int i = 0; i < values.length; i++) {
        values[i] = 1 - FastMath.exp(-x[i] / fourD);
      }
      return values;
    }

    /** {@inheritDoc} */
    @Override
    public double[][] jacobian(double[] variables) {
      // Compute the gradients using calculus differentiation:
      // y = 1 - a
      // a = exp(b)
      // b = -x / 4D
      //
      // y' = -a'
      // a' = exp(b) * b'
      // b' = -1 * -x / 4D^2 = x / 4D^2
      // y' = -exp(b) * x / 4D^2
      // = -a * -b / D
      // = a * b / D
      // = exp(b) * b / D

      final double d = variables[0];
      final double fourD = 4 * d;
      final double[][] jacobian = new double[x.length][variables.length];

      for (int i = 0; i < jacobian.length; ++i) {
        final double b = -x[i] / fourD;
        jacobian[i][0] = FastMath.exp(b) * b / d;
      }

      //// Check numerically ...
      // double[][] jacobian2 = super.jacobian(variables);
      // for (int i = 0; i < jacobian.length; i++)
      // {
      // System.out.printf("dD = %g : %g = %g\n", jacobian[i][0], jacobian2[i][0],
      // DoubleEquality.relativeError(jacobian[i][0], jacobian2[i][0]));
      // }

      return jacobian;
    }
  }

  /**
   * Compute the probability of mean-squared distance x given a mixed population with set fractions
   * and diffusion coefficients.
   *
   * <p>Function used for maximum likelihood fitting.
   */
  static class MixedJumpDistanceFunction extends Function implements MultivariateFunction {
    /**
     * Instantiates a new mixed jump distance function.
     *
     * @param x the x
     * @param estimatedD the estimated D
     * @param n the n
     */
    public MixedJumpDistanceFunction(double[] x, double estimatedD, int n) {
      super(x, null, estimatedD, n);
    }

    @Override
    public double evaluate(double x, double[] params) {
      // Compute the probability:
      // p = sum [ Fj/4Dj * exp(-x/4Dj) ]
      double sum = 0;
      double total = 0;
      for (int i = 0; i < n; i++) {
        // final double f = getF(params[i * 2]);
        // final double fourD = 4 * getD(params[i * 2 + 1]);
        final double f = params[i * 2];
        final double fourD = 4 * params[i * 2 + 1];
        sum += (f / fourD) * FastMath.exp(-x / fourD);
        total += f;
      }
      return sum / total;
    }

    /**
     * Evaluate all.
     *
     * @param params the params
     * @return the values
     */
    public double[] evaluateAll(double[] params) {
      double total = 0;
      final double[] f_d = new double[n];
      for (int i = 0; i < n; i++) {
        // f_d[i] = getF(params[i * 2]);
        f_d[i] = params[i * 2];
        total += f_d[i];
      }

      final double[] fourD = new double[n];
      for (int i = 0; i < n; i++) {
        // fourD[i] = 4 * getD(params[i * 2 + 1]);
        fourD[i] = 4 * params[i * 2 + 1];
        f_d[i] = (f_d[i] / total) / fourD[i];
      }

      // Compute the probability:
      // p = sum [ Fj/4Dj * exp(-x/4Dj) ]
      final double[] values = new double[x.length];
      for (int i = 0; i < x.length; i++) {
        double sum = 0;
        for (int j = 0; j < n; j++) {
          sum += f_d[j] * FastMath.exp(-x[i] / fourD[j]);
        }
        values[i] = sum;
      }
      return values;
    }

    /** {@inheritDoc} */
    @Override
    public double value(double[] params) {
      // Compute the log-likelihood
      double l = 0;
      for (final double p : evaluateAll(params)) {
        l += Math.log(p);
      }
      // Debug the call from the optimiser
      if (DEBUG_OPTIMISER) {
        final double[] F = new double[n];
        final double[] D = new double[n];
        for (int i = 0; i < n; i++) {
          F[i] = params[i * 2];
          D[i] = params[i * 2 + 1];
        }
        System.out.printf("%s : %s = %f\n", Arrays.toString(F), Arrays.toString(D), l);
      }
      return l;
    }
  }

  /**
   * Compute the probability of mean-squared distance being within x given a mixed population with
   * set fractions and diffusion coefficients.
   *
   * <p>Function used for least-squares fitting of cumulative histogram of jump distances.
   */
  static class MixedJumpDistanceCumulFunction extends Function {
    /**
     * Instantiates a new mixed jump distance cumul function.
     *
     * @param x the x
     * @param y the y
     * @param estimatedD the estimated D
     * @param n the n
     */
    public MixedJumpDistanceCumulFunction(double[] x, double[] y, double estimatedD, int n) {
      super(x, y, estimatedD, n);
    }

    @Override
    public double evaluate(double x, double[] params) {
      double sum = 0;
      double total = 0;
      for (int i = 0; i < n; i++) {
        final double f = params[i * 2];
        sum += f * FastMath.exp(-x / (4 * params[i * 2 + 1]));
        total += f;
      }
      return 1 - sum / total;
    }

    /**
     * Gets the value.
     *
     * @param variables the variables
     * @return the values
     */
    public double[] getValue(double[] variables) {
      double total = 0;
      for (int i = 0; i < n; i++) {
        total += variables[i * 2];
      }

      final double[] fourD = new double[n];
      final double[] f = new double[n];
      for (int i = 0; i < n; i++) {
        f[i] = variables[i * 2] / total;
        fourD[i] = 4 * variables[i * 2 + 1];
      }

      final double[] values = new double[x.length];
      for (int i = 0; i < values.length; i++) {
        double sum = 0;
        for (int j = 0; j < n; j++) {
          sum += f[j] * FastMath.exp(-x[i] / fourD[j]);
        }
        values[i] = 1 - sum;
      }
      return values;
    }
  }

  /**
   * Compute the probability of mean-squared distance being within x given a mixed population with
   * set fractions and diffusion coefficients.
   *
   * <p>Function used for least-squares fitting of cumulative histogram of jump distances.
   */
  static class MixedJumpDistanceCumulFunctionGradient extends MixedJumpDistanceCumulFunction
      implements MultivariateVectorFunction {
    /**
     * Instantiates a new mixed jump distance cumul function gradient.
     *
     * @param x the x
     * @param y the y
     * @param estimatedD the estimated D
     * @param n the n
     */
    public MixedJumpDistanceCumulFunctionGradient(double[] x, double[] y, double estimatedD,
        int n) {
      super(x, y, estimatedD, n);
    }

    /** {@inheritDoc} */
    @Override
    public double[] value(double[] point) throws IllegalArgumentException {
      return getValue(point);
    }

    /** {@inheritDoc} */
    @Override
    public double[][] jacobian(double[] variables) {
      // Compute the gradients using calculus differentiation:
      // y = 1 - sum(a)
      // The sum is over n components of the following function
      // a = f * exp(b)
      // b = -x / 4D
      // Each function contributes a fraction f:
      // f = fj / sum_j(f)

      // The gradient is the sum of the individual gradients. The diffusion coefficient is only
      // used per component. The fraction is used in all, either with the fraction as the
      // numerator (A) or part of the denominator (B)
      // E.G.
      // f(A) = A / (A+B+C)
      // Quotient rule: f = g / h => f' = (g'h - gh') / h^2
      // f'(A) = ((A+B+C) - A) / (A+B+C)^2
      // = (B+C) / (A+B+C)^2
      // = (sum(f) - f) / sum(f)^2
      // f'(B) = -A / (A+B+C)^2
      // = -f / sum(f)^2

      // Differentiate with respect to D is easier:
      // y' = -a'
      // a' = f * exp(b) * b'
      // b' = -1 * -x / 4D^2 = x / 4D^2
      // y' = f * -exp(b) * x / 4D^2
      // = f * -a * -b / D
      // = f * a * b / D
      // = f * exp(b) * b / D

      final double[] fourD = new double[n];
      final double[] f = new double[n];
      double total = 0;
      for (int i = 0; i < n; i++) {
        f[i] = variables[i * 2];
        fourD[i] = 4 * variables[i * 2 + 1];
        total += f[i];
      }

      final double[] fraction = new double[n];
      final double[] total_f = new double[n];
      final double[] f_total = new double[n];
      for (int i = 0; i < n; i++) {
        fraction[i] = f[i] / total;
        // Because we use y = 1 - sum(a) all coefficients are inverted
        total_f[i] = -1 * (total - f[i]) / (total * total);
        f_total[i] = -1 * -f[i] / (total * total);
      }

      final double[][] jacobian = new double[x.length][variables.length];

      final double[] b = new double[n];
      for (int i = 0; i < x.length; ++i) {
        for (int j = 0; j < n; j++) {
          b[j] = -x[i] / fourD[j];
        }

        for (int j = 0; j < n; j++) {
          // Gradient for the diffusion coefficient
          jacobian[i][j * 2 + 1] = fraction[j] * FastMath.exp(b[j]) * b[j] / variables[j * 2 + 1];

          // Gradient for the fraction f
          jacobian[i][j * 2] = total_f[j] * FastMath.exp(b[j]);
          for (int k = 0; k < n; k++) {
            if (j == k) {
              continue;
            }
            jacobian[i][j * 2] += f_total[k] * FastMath.exp(b[k]);
          }
        }
      }

      //// Check numerically ...
      // double[][] jacobian2 = super.jacobian(variables);
      // for (int i = 0; i < jacobian.length; i++)
      // {
      // StringBuilder sb = new StringBuilder();
      // for (int j = 0; j < jacobian[i].length; j++)
      // {
      // sb.append(" d").append(j).append(" = ").append(jacobian[i][j]).append(" : ")
      // .append(jacobian2[i][j]).append(" = ")
      // .append(DoubleEquality.relativeError(jacobian[i][j], jacobian2[i][j]));
      // }
      // System.out.println(sb.toString());
      // }

      return jacobian;
    }
  }

  /**
   * Compute the probability of mean-squared distance being within x given a mixed population with
   * set fractions and diffusion coefficients.
   *
   * <p>Function used for least-squares fitting of cumulative histogram of jump distances.
   */
  static class MixedJumpDistanceCumulFunctionMultivariate extends MixedJumpDistanceCumulFunction
      implements MultivariateFunction {
    /**
     * Instantiates a new mixed jump distance cumul function multivariate.
     *
     * @param x the x
     * @param y the y
     * @param estimatedD the estimated D
     * @param n the n
     */
    public MixedJumpDistanceCumulFunctionMultivariate(double[] x, double[] y, double estimatedD,
        int n) {
      super(x, y, estimatedD, n);
    }

    /** {@inheritDoc} */
    @Override
    public double value(double[] parameters) {
      final double[] obs = getValue(parameters);

      // Optimise the sum of squares
      double ss = 0;
      for (int i = x.length; i-- > 0;) {
        final double dx = y[i] - obs[i];
        ss += dx * dx;
      }
      // Debug the call from the optimiser
      if (DEBUG_OPTIMISER) {
        final double[] F = new double[n];
        final double[] D = new double[n];
        for (int i = 0; i < n; i++) {
          F[i] = parameters[i * 2];
          D[i] = parameters[i * 2 + 1];
        }
        System.out.printf("%s : %s = %f\n", Arrays.toString(F), Arrays.toString(D), ss);
      }
      return ss;
    }
  }

  /**
   * @return The number restarts for fitting.
   */
  public int getFitRestarts() {
    return fitRestarts;
  }

  /**
   * @param fitRestarts The number restarts for fitting
   */
  public void setFitRestarts(int fitRestarts) {
    this.fitRestarts = fitRestarts;
  }

  /**
   * @return the min fraction for each population in a mixed population.
   */
  public double getMinFraction() {
    return minFraction;
  }

  /**
   * @param minFraction the min fraction for each population in a mixed population
   */
  public void setMinFraction(double minFraction) {
    this.minFraction = minFraction;
  }

  /**
   * @return the min difference between diffusion coefficients in a mixed population.
   */
  public double getMinDifference() {
    return minDifference;
  }

  /**
   * @param minDifference the min difference between diffusion coefficients in a mixed population
   */
  public void setMinDifference(double minDifference) {
    this.minDifference = minDifference;
  }

  /**
   * @return the minimum number of different molecules to fit in a mixed population model.
   */
  public int getMinN() {
    return minN;
  }

  /**
   * @param n the minimum number of different molecules to fit in a mixed population model
   */
  public void setMinN(int n) {
    if (n < 1) {
      n = 1;
    }
    minN = n;
  }

  /**
   * @return the maximum number of different molecules to fit in a mixed population model.
   */
  public int getMaxN() {
    return maxN;
  }

  /**
   * @param n the maximum number of different molecules to fit in a mixed population model
   */
  public void setMaxN(int n) {
    if (n < 1) {
      n = 1;
    }
    maxN = n;
  }

  /**
   * Sets the curve logger.
   *
   * @param curveLogger the new curve logger
   */
  public void setCurveLogger(CurveLogger curveLogger) {
    this.curveLogger = curveLogger;
  }

  /**
   * Get the cumulative jump distance histogram given a set of jump distance values.
   *
   * @param values The jump distances
   * @return The JD cumulative histogram as two arrays: { MSD, CumulativeProbability }
   */
  public static double[][] cumulativeHistogram(double[] values) {
    return MathUtils.cumulativeHistogram(values, true);
  }

  /**
   * Gets the d.
   *
   * @param d the d
   * @return Return d or minD whichever is larger
   */
  public double getD(double d) {
    return (d < minD) ? minD : d;
  }

  /**
   * Gets the f.
   *
   * @param f the f
   * @return Return f or minFraction whichever is larger
   */
  public double getF(double f) {
    return (f < minFraction) ? minFraction : f;
  }

  /**
   * Gets the minimum diffusion coefficient.
   *
   * @return the minimum diffusion coefficient
   */
  public double getMinD() {
    return minD;
  }

  /**
   * Sets the minimum diffusion coefficient. Only one population is allowed to have a diffusion
   * coefficient below this.
   *
   * <p>This value should be set using an understanding of the localisation precision of the
   * results. Any population with a diffusion coefficient that creates average jumps below the
   * localisation precision is effectively static since the system could not accurately measure the
   * jumps. It is impossible to distinguish populations that move below the localisation precision
   * and so only one population is allowed below this level.
   *
   * @param minD the minimum diffusion coefficient
   */
  public void setMinD(double minD) {
    this.minD = minD;
  }

  /**
   * Get the conversion factor to convert an observed mean-squared distance (MSD) between n frames
   * into the actual MSD.
   *
   * <p>Note that diffusion of a molecule within a frame means that the position of the molecule is
   * an average within the frame. This leads to condensation of the observed distance traveled by
   * the particle between two frames. The start and end frame locations have condensed diffusion
   * within the frame to a single point. This condensation has the effect of reducing the effective
   * time that diffusion occurred in the start and end frame. The observed MSD can be converted to
   * the corrected MSD by applying a factor:
   *
   * <pre>
   * {@code
   * observed = actual * (n - 1/3) / n
   * actual = observed * n / (n - 1/3)
   * }
   * </pre>
   *
   * Note this is only valid for {@code n>=1}.
   *
   * @param n the n
   * @return the conversion factor
   */
  public static double getConversionfactor(int n) {
    return n / (n - THIRD);
  }

  /**
   * Get the conversion factor to convert an observed mean-squared distance (MSD) between n frames
   * into the actual MSD.
   *
   * <p>Note that diffusion of a molecule within a frame means that the position of the molecule is
   * an average within the frame. This leads to condensation of the observed distance traveled by
   * the particle between two frames. The start and end frame locations have condensed diffusion
   * within the frame to a single point. This condensation has the effect of reducing the effective
   * time that diffusion occurred in the start and end frame and the total number of frames should
   * be reduced by 1/3.
   *
   * <p>In the event that n is less than 1 the two frames overlap. Consequently there is
   * interference where some of the same molecule positions have been used to compute the average
   * start and end location within the frames. This must be modelled using a different formula.
   *
   * <p>Simulations using multiple simulation steps within each frame were used to compute the MSD
   * at different frame separation intervals. These curves were compared to the expected MSD for the
   * simulated diffusion coefficient to produce a correction factor curve. This was fitted for
   * {@code n>=1} and {@code n<1}. The observed MSD can be converted to the corrected MSD by
   * applying a factor:
   *
   * <pre>
   * {@code
   * n>=1:
   * observed = actual * (n - 1/3) / n
   * actual = observed * n / (n - 1/3)
   *
   * n<1:
   * observed = actual * (n - n*n / 3)
   * actual = observed / (n - n*n / 3)
   * }
   * </pre>
   *
   * Note this is valid for {@code n>=0}.
   *
   * @param n the n
   * @return the conversion factor
   */
  public static double getConversionfactor(double n) {
    if (n > 1) {
      return n / (n - THIRD);
    }
    if (n > 0) {
      return 1 / (n - n * n / 3.0);
    }
    return 0;
  }

  /**
   * Get the corrected time between n frames for an observed mean-squared distance (MSD).
   *
   * <p>Note that diffusion of a molecule within a frame means that the position of the molecule is
   * an average within the frame. This leads to condensation of the observed distance traveled by
   * the particle between two frames. The start and end frame locations have condensed diffusion
   * within the frame to a single point. This condensation has the effect of reducing the effective
   * time that diffusion occurred in the start and end frame and the total number of frames should
   * be reduced by 1/3.
   *
   * <pre>
   * {@code correctedFrames = n - 1/3}
   * </pre>
   *
   * Note this is only valid for {@code n>=1}.
   *
   * @param n the n
   * @return the corrected time
   */
  public static double getCorrectedTime(int n) {
    return n - THIRD;
  }

  /**
   * Convert an observed mean-squared distance (MSD) between n frames into the actual MSD.
   *
   * <p>Note that diffusion of a molecule within a frame means that the position of the molecule is
   * an average within the frame. This leads to condensation of the observed distance traveled by
   * the particle between two frames. The start and end frame locations have condensed diffusion
   * within the frame to a single point. This condensation has the effect of reducing the effective
   * time that diffusion occurred in the start and end frame and the total number of frames should
   * be reduced by 1/3. The observed MSD can be converted to the corrected MSD by applying a factor:
   *
   * <pre>
   * {@code
   * observed = actual * (n - 1/3) / n
   * actual = observed * n / (n - 1/3)
   * }
   * </pre>
   *
   * Note this is only valid for {@code n>=1}.
   *
   * @param msd The observed MSD
   * @param n The number of frames separating the start and end points for the MSD
   * @return The actual MSD
   */
  public static double convertObservedToActual(double msd, int n) {
    return msd * n / (n - THIRD);
  }

  /**
   * Convert an actual mean-squared distance (MSD) between n frames into the observed MSD.
   *
   * <p>Note that diffusion of a molecule within a frame means that the position of the molecule is
   * an average within the frame. This leads to condensation of the observed distance traveled by
   * the particle between two frames. The start and end frame locations have condensed diffusion
   * within the frame to a single point. This condensation has the effect of reducing the effective
   * time that diffusion occurred in the start and end frame and the total number of frames should
   * be reduced by 1/3. The observed MSD can be converted to the corrected MSD by applying a factor:
   *
   * <pre>
   * {@code
   * observed = actual * (n - 1/3) / n
   * actual = observed * n / (n - 1/3)
   * }
   * </pre>
   *
   * Note this is only valid for {@code n>=1}.
   *
   * @param msd The actual MSD
   * @param n The number of frames separating the start and end points for the MSD
   * @return The actual MSD
   */
  public static double convertActualToObserved(double msd, int n) {
    return msd * (n - THIRD) / n;
  }

  /**
   * @return The localisation error (s) of the start and end coordinates of the jump (in um).
   */
  public double getError() {
    return Math.sqrt(s2);
  }

  /**
   * Set the localisation error (s) of the start and end coordinates of the jump. The error is used
   * to compute the apparent diffusion coefficient: D* = D - s^2
   *
   * @param error The error
   * @param nm True if the error is in nm (default is um)
   */
  public void setError(double error, boolean nm) {
    if (nm) {
      error /= 1000;
    }
    this.s2 = error * error;
  }

  /**
   * @return True if correcting MSD between frames.
   */
  public boolean isMsdCorrection() {
    return msdCorrection;
  }

  /**
   * Set to true to correct the diffusion coefficient to the apparent diffusion coefficient by
   * compensating for the averaging of diffusion with a time frame into a single location.
   *
   * @param msdCorrection True if correcting MSD between frames
   */
  public void setMsdCorrection(boolean msdCorrection) {
    this.msdCorrection = msdCorrection;
  }

  /**
   * @return The number of frames between the start and end coordinates of the jump.
   */
  public int getN() {
    return n;
  }

  /**
   * Set the number of frames between the start and end coordinates of the jump.
   *
   * @param n The number of frames (must be strictly positive)
   */
  public void setN(int n) {
    this.n = n;
  }

  /**
   * Gets the significance level for testing the log-likelihood ratio.
   *
   * @return the significance level
   */
  public double getSignificanceLevel() {
    return significanceLevel;
  }

  /**
   * Sets the significance level. This is used when testing the log-likelihood ratio during maximum
   * likelihood fitting that an increase in model parameters improves the model.
   *
   * @param significanceLevel the new significance level
   */
  public void setSignificanceLevel(double significanceLevel) {
    this.significanceLevel = significanceLevel;
  }

  /**
   * @return The time difference between each frame.
   */
  public double getDeltaT() {
    return deltaT;
  }

  /**
   * Set the time difference between each frame. The total time is {@link #getDeltaT()} *
   * {@link #getN()}
   *
   * @param deltaT The time difference
   */
  public void setDeltaT(double deltaT) {
    this.deltaT = deltaT;
  }

  /**
   * @return True if the number of frames and time delta have been set.
   */
  public boolean isCalibrated() {
    return deltaT > 0 && n > 0;
  }

  /**
   * Convert the diffusion coefficients (um^2/jump) into apparent diffusion coefficients (um^2/s) if
   * the time and frames are calibrated:
   *
   * <pre>
   * D* = factor * (D - s2) / (n * deltaT)
   * </pre>
   *
   * where factor is the conversion factor to increase the MSD to correct for diffusion within the
   * frame
   *
   * @param d The fitted diffusion coefficients D (in um^2)
   * @return The apparent diffusion coefficients D* (in um^2/s)
   */
  public double[] calculateApparentDiffusionCoefficient(double... d) {
    if (isCalibrated()) {
      d = d.clone();
      final double f = ((msdCorrection) ? getConversionfactor(n) : 1) / (n * deltaT);
      for (int i = 0; i < d.length; i++) {
        d[i] = Math.max(0, d[i] - s2) * f;
      }
    }
    return d;
  }

  /**
   * Get the fit value. This will be the Adjusted R^2 for least squares estimation or the
   * information criterion from maximum likelihood estimation.
   *
   * @return The fit value from the last successful fit
   */
  public double getFitValue() {
    return lastFitValue;
  }
}
