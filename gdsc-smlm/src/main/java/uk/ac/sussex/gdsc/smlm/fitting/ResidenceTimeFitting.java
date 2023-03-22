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

package uk.ac.sussex.gdsc.smlm.fitting;

import java.util.Arrays;
import java.util.function.Predicate;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
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
import org.apache.commons.math3.optim.univariate.BracketFinder;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.random.RandomGenerator;
import uk.ac.sussex.gdsc.core.logging.LoggerUtils;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SortUtils;
import uk.ac.sussex.gdsc.core.utils.ValidationUtils;
import uk.ac.sussex.gdsc.core.utils.rng.RandomGeneratorAdapter;
import uk.ac.sussex.gdsc.core.utils.rng.UniformRandomProviders;
import uk.ac.sussex.gdsc.smlm.math3.optim.nonlinear.scalar.noderiv.CustomPowellOptimizer;
import uk.ac.sussex.gdsc.smlm.utils.StdMath;

/**
 * Perform fitting of the residence time. Residence time is modelled using a dissociation rate.
 *
 * <p>Residence times are assumed to be discrete counts at a specified time resolution, e.g.
 * correspond to frame counts taken at a fixed exposure time.
 *
 * <p>Continuous residence time data uses an exponential model. The discrete residence time data is
 * fit by mapping the dissociation rate to the p-value of the geometric distribution (the discrete
 * analogue of the exponential distribution).
 *
 * <p>For analysis of continuous data using an exponential model see the
 * {@link JumpDistanceAnalysis} class (which fits the rate as 1/4D).
 */
public final class ResidenceTimeFitting {
  /** Count at each time point. */
  private final int[] count;
  /** Time resolution (smallest difference between adjacent time points. */
  private final double resolution;

  /**
   * Contains the residence time model for a mixture of populations.
   *
   * <pre>
   * SF(t) = sum_i { f_i exp(-k_i t) }
   * </pre>
   *
   * <p>where {@code k} is the dissociation rate; {@code 1/k} is the mean residence time. The sum of
   * population fractions {@code f_i} equal 1.
   */
  public interface Model {
    /**
     * Gets the number of populations.
     *
     * @return the size
     */
    int getSize();

    /**
     * Gets the dissociation rate for the specified {@code index}.
     *
     * @param index the index
     * @return the k
     */
    double getRate(int index);

    /**
     * Gets the population fraction for the specified {@code index}. The sum of all fractions for
     * the population size should equal 1.
     *
     * @param index the index
     * @return the fraction
     */
    double getFraction(int index);

    /**
     * Gets the mean residence time for the specified {@code index}.
     *
     * @param index the index
     * @return the residence time
     */
    default double getResidenceTime(int index) {
      return 1 / getRate(index);
    }

    /**
     * Return the survival function for the specified time point. This is the complementary
     * cumulative distribution function: {@code P(X > t)}.
     *
     * @param t the time
     * @return sf(t)
     */
    default double sf(double t) {
      double p = 0;
      for (int i = 0; i < getSize(); i++) {
        p += getFraction(i) * StdMath.exp(-getRate(i) * t);
      }
      return p;
    }
  }

  /**
   * Single population model.
   */
  private static class Model1 implements Model {
    /** The rate constant. */
    private final double rate;

    /**
     * Create an instance.
     *
     * @param rate the rate constant
     */
    Model1(double rate) {
      this.rate = rate;
    }

    @Override
    public int getSize() {
      return 1;
    }

    @Override
    public double getRate(int index) {
      assert index == 0;
      return rate;
    }

    @Override
    public double getFraction(int index) {
      assert index == 0;
      return 1;
    }
  }

  /**
   * Dual population model.
   */
  private static class Model2 implements Model {
    /** The rate constant for population 0. */
    private final double k0;
    /** The rate constant for population 1. */
    private final double k1;
    /** Fraction of population 0. */
    private final double f0;

    /**
     * Create and instance.
     *
     * @param k0 the rate constant for 0
     * @param k1 the rate constant for 1
     * @param f0 the fraction of population 0
     */
    Model2(double k0, double k1, double f0) {
      this.k0 = k0;
      this.k1 = k1;
      this.f0 = f0;
    }

    @Override
    public int getSize() {
      return 2;
    }

    @Override
    public double getRate(int index) {
      if (index == 0) {
        return k0;
      }
      assert index == 1;
      return k1;
    }

    @Override
    public double getFraction(int index) {
      if (index == 0) {
        return f0;
      }
      assert index == 1;
      return 1 - f0;
    }
  }

  /**
   * Single population fitting model.
   */
  static class Function1 {
    /** Count at each time point. */
    private final int[] count;
    /** Time resolution. */
    private final double resolution;
    /** Total count. */
    private final int total;

    /**
     * Create an instance.
     *
     * @param n the count at each time point
     * @param resolution the time resolution
     */
    Function1(int[] n, double resolution) {
      this.count = n;
      this.resolution = resolution;
      total = Arrays.stream(n).sum();
    }

    /**
     * Compute the log-likelihood of the rate constant.
     *
     * @param k the rate constant
     * @return the log-likelihood
     */
    double ll(double k) {
      // Use the probability of a geometric distribution:
      // pmf(x) = (1 - p)^x * p
      // == exp(log(1-p) * x) * p
      // log(pmf) = log(1-p) * x + log(p)

      // Exponential X ==> Geometric Y = floor(X)
      // p = 1 - exp(-rate)
      // log(p) = log1p(-exp(-rate))
      // log(1-p) = -rate
      // Note: the exponential has been binned at a set resolution.
      // The effect is to scale the exponential mean by 1/resolution to discrete data.
      // To scale back multiply the rate by the resolution.
      final double logp = Math.log1p(-Math.exp(-k * resolution));
      final double log1mp = -k * resolution;

      // Special case for x=0: n[0] * logp is handled by the final addition of sum(n) * log(p)
      double sum = 0;
      for (int x = 1; x < count.length; x++) {
        sum += count[x] * (log1mp * x);
      }
      // Since the exponential has been scaled by 'resolution' to map to the geometric,
      // all p-values are too small for the equivalent exponential.
      // Scale the p-values back: p / resolution => c * -log(resolution)
      sum += total * (logp - Math.log(resolution));
      return sum;
    }
  }

  /**
   * Dual population fitting model.
   */
  static class Function2 {
    /** Count at each time point. */
    private final int[] count;
    /** Time resolution. */
    private final double resolution;
    /** Total count. */
    private final int total;

    /**
     * Create an instance.
     *
     * @param n the count at each time point
     * @param resolution the time resolution
     */
    Function2(int[] n, double resolution) {
      this.count = n;
      this.resolution = resolution;
      total = Arrays.stream(n).sum();
    }

    /**
     * Compute the log-likelihood of the two rate constants.
     *
     * @param k0 the rate constant for the first population
     * @param k1 the rate constant for the second population
     * @param f0 the fraction of the first population
     * @return the log-likelihood
     */
    double ll(double k0, double k1, double f0) {
      // Compute a weighted sum of geometric distributions:
      // P(x) = f0 (1 - p0)^x * p0 + (1-f0) (1 - p1)^x * p1
      // = f0 exp(log(1-p0)*x) * p0 + (1-f0) exp(log(1-p1)*x) * p1
      // This is logged and summed over all observations.
      final double p0 = -Math.expm1(-k0 * resolution);
      final double log1mp0 = -k0 * resolution;
      final double p1 = -Math.expm1(-k1 * resolution);
      final double log1mp1 = -k1 * resolution;
      final double f0p0 = f0 * p0;
      final double f1p1 = (1 - f0) * p1;
      // Special case for x=0
      double sum = count[0] * Math.log(f0p0 + f1p1);
      for (int x = 1; x < count.length; x++) {
        sum +=
            count[x] * Math.log(f0p0 * StdMath.exp(log1mp0 * x) + f1p1 * StdMath.exp(log1mp1 * x));
      }
      // Since the exponential has been scaled by 'resolution' to map to the geometric,
      // all p-values are too small for the equivalent exponential.
      // Scale the p-values back: p / resolution => c * -log(resolution)
      sum -= total * Math.log(resolution);
      return sum;
    }
  }

  /**
   * Contains results from fitting a model to the residence time data.
   */
  public static class FitResult {
    private final int numberOfPoints;
    private final int numberOfParameters;
    private final double ll;

    /**
     * Create an instance.
     *
     * @param n data points
     * @param p parameters
     * @param ll log-likelihood
     */
    FitResult(int n, int p, double ll) {
      this.numberOfPoints = n;
      this.numberOfParameters = p;
      this.ll = ll;
    }

    /**
     * Gets the number of data points.
     *
     * @return the number of data points
     */
    public int getN() {
      return numberOfPoints;
    }

    /**
     * Gets the number of parameters.
     *
     * @return the number of parameters
     */
    public int getP() {
      return numberOfParameters;
    }

    /**
     * Gets the log likelihood.
     *
     * @return the log likelihood
     */
    public double getLogLikelihood() {
      return ll;
    }
  }

  /**
   * Create an instance.
   *
   * @param n the count at each time point
   * @param resolution the time resolution
   */
  private ResidenceTimeFitting(int[] n, double resolution) {
    this.count = n;
    this.resolution = resolution;
  }

  /**
   * Create an instance.
   *
   * <p>It is assumed the counts are in approximately descending order, otherwise an exponential
   * model is not valid.
   *
   * <p>There must be information in the histogram to allow fitting. This requires at least two
   * non-zero values, otherwise an exception is raised. Support is provided for a total observed
   * count less than a 32-bit signed integer.
   *
   * @param resolution the time resolution (in seconds)
   * @param n the count at each time point
   * @return an instance
   * @throws IllegalArgumentException if the time resolution is not strictly positive; if there are
   *         less than two time points; if any count is negative; if the sum of counts is zero or
   *         {@code >= 2^31}; or if there are less than two non-zero counts.
   */
  public static ResidenceTimeFitting of(double resolution, int[] n) {
    ValidationUtils.checkStrictlyPositive(resolution, "time resolution");
    ValidationUtils.checkArgument(n.length > 1, "not enough time points: %d", n.length);
    int nonZero = 0;
    long sum = 0;
    for (final int v : n) {
      if (v > 0) {
        nonZero++;
        sum += v;
      } else {
        ValidationUtils.checkPositive(v, "count");
      }
    }
    ValidationUtils.checkArgument(nonZero > 1, "too few non-zero counts: %d", nonZero);
    ValidationUtils.checkArgument(sum <= Integer.MAX_VALUE, "sum of counts is too large: %d", sum);
    return new ResidenceTimeFitting(n, resolution);
  }

  /**
   * Gets the observed dissociation rate.
   *
   * <p>This is computed assuming the histogram of residence times follows a geometric distribution.
   *
   * @return the mean
   */
  double getRate() {
    long s = 0;
    int c = 0;
    for (int i = 0; i < count.length; i++) {
      s += i * count[i];
      c += count[i];
    }
    // Note: s and c will be above zero due to validation in the constructor.
    // Mean of observed geometric distribution
    final double mean = (double) s / c;
    // Geometric distribution:
    // mean = (1-p) / p => p = 1 / (1+mean)
    final double p = 1 / (1 + mean);
    // Convert to related exponential
    // p = 1-exp(-rate) => rate = -ln(1 - p)

    // Note: Scale the rate by the resolution
    return -Math.log1p(-p) / resolution;
  }

  /**
   * Fit a model with the specified number of populations to the data.
   *
   * <p>Supports a 1 population or 2 population model.
   *
   * <p>The following options are recognised: <ul>
   *
   * <li>{@link java.util.logging.Logger}: Logger used to record fitting details.
   *
   * <li>{@link org.apache.commons.math3.optim.MaxEval}: maximum evaluations.
   *
   * </ul>
   *
   * @param size the size
   * @param options the options
   * @return the result (or null if fitting failed)
   * @throws IllegalArgumentException if {@code size} is not 1 or 2
   */
  public Pair<FitResult, Model> fit(int size, Object... options) {
    ValidationUtils.checkStrictlyPositive(size);
    ValidationUtils.checkArgument(size <= 2, "Unsupported size: %d", size);
    if (size == 1) {
      return fit1(options);
    }
    return fit2(options);
  }

  /**
   * Fit a 1 population model to the data.
   *
   * @param options the options
   * @return the result (or null if fitting failed)
   */
  private Pair<FitResult, Model> fit1(Object[] options) {
    final Logger logger = orElse(Logger.class, LoggerUtils.createIfNull(null), options);
    final MaxEval maxEval = orElse(new MaxEval(1000), options);

    // Estimate rate constant from the mean for single exponential.
    final double k = getRate();
    final UnivariateFunction f1 = new Function1(count, resolution)::ll;

    try {
      final BracketFinder bf = new BracketFinder();
      bf.search(f1, GoalType.MAXIMIZE, k, k * 1.0001);
      final BrentOptimizer optimizer = new BrentOptimizer(1e-8, Double.MIN_VALUE);
      final UnivariatePointValuePair result =
          optimizer.optimize(new UnivariateObjectiveFunction(f1), GoalType.MAXIMIZE,
              new SearchInterval(bf.getLo(), bf.getHi(), bf.getMid()), maxEval);
      logger.info(() -> String.format("Fit (N=1) : %s, MLE = %f (%d evaluations)",
          formatK(new double[] {k}), result.getValue(), optimizer.getEvaluations()));
      return Pair.of(new FitResult(Arrays.stream(count).sum(), 1, result.getValue()),
          new Model1(result.getPoint()));
    } catch (final TooManyEvaluationsException ex) {
      LoggerUtils.log(logger, Level.INFO, "Failed to fit (N=1) : %s", ex.getMessage());
    }
    return null;
  }

  /**
   * Fit a 2 population model to the data.
   *
   * @param options the options
   * @return the result (or null if fitting failed)
   */
  private Pair<FitResult, Model> fit2(Object[] options) {
    final Logger logger = orElse(Logger.class, LoggerUtils.createIfNull(null), options);
    final MaxEval maxEval = orElse(new MaxEval(20000), options);

    // Create estimates
    final double m = 1 / getRate();
    final double m1 = m * 1.75;
    final double m2 = m * 0.25;
    final InitialGuess guess = new InitialGuess(new double[] {1 / m1, 1 / m2, 0.5});

    // Note:
    // Choice of fitting optimiser has been adapted from the JumpDistanceAnalysis plugin.
    // This is a work-in-progress as it may not be suitable for real data. It works
    // on simulated data where the time resolution allows a histogram with at least
    // 2 bins to be created for the faster dissociation population and there are
    // reasonable counts from each sub-population.

    final CustomPowellOptimizer powellOptimizer = createCustomPowellOptimizer();
    // Initial step size used to bracket the minimum in the line search
    // Use a smaller step size for the fraction component
    final CustomPowellOptimizer.BasisStep step =
        new CustomPowellOptimizer.BasisStep(new double[] {0.1 / m1, 0.1 / m2, 0.05});

    final Function2 f2 = new Function2(count, resolution);
    final ObjectiveFunction fun =
        new ObjectiveFunction(point -> f2.ll(point[0], point[1], point[2]));

    // Bound on rate set using the range of the residence times
    // Limit fraction to the the (0, 1) interval.
    final double hi = 1 / (1e-3 * resolution);
    final double lo = 1 / (1e3 * (resolution * count.length));
    final double[] lB = {lo, lo, 1e-9};
    final double[] uB = {hi, hi, 1 - 1e-9};
    final SimpleBounds bounds = new SimpleBounds(lB, uB);

    int evaluations = 0;
    PointValuePair bestSolution = null;

    try {
      bestSolution = powellOptimizer.optimize(fun, guess, bounds, step, GoalType.MAXIMIZE, maxEval);

      evaluations = powellOptimizer.getEvaluations();
      LoggerUtils.log(logger, Level.FINE, "Powell optimiser fit (N=2) : MLE = %f (%d evaluations)",
          bestSolution.getValue(), powellOptimizer.getEvaluations());
    } catch (final TooManyEvaluationsException | ConvergenceException ex) {
      LoggerUtils.log(logger, Level.INFO, "Powell optimiser failed to fit (N=2) : %s",
          ex.getMessage());
    }

    if (bestSolution == null) {
      LoggerUtils.log(logger, Level.INFO, "Trying CMAES optimiser with restarts ...");

      // Try a bounded CMAES optimiser since the Powell optimiser appears to be
      // sensitive to the order of the parameters. It is not good when the fast particle
      // is the minority fraction. Could this be due to too low an upper bound?

      // The sigma determines the search range for the variables.
      // It should be 1/3 of the initial search region.
      final double[] s = new double[lB.length];
      for (int i = 0; i < s.length; i++) {
        s[i] = (uB[i] - lB[i]) / 3;
      }
      final OptimizationData sigma = new CMAESOptimizer.Sigma(s);
      final OptimizationData popSize =
          new CMAESOptimizer.PopulationSize((int) (4 + Math.floor(3 * Math.log(count.length))));

      // Iterate this for stability in the initial guess.
      // Note: The optimiser does not throw exceptions for too many evaluations; it returns
      // the current optimum. So we manually check the evaluations and ignore the non-converged
      // result.
      final CMAESOptimizer cmaesOptimizer = createCmaesOptimizer();
      final Predicate<CMAESOptimizer> converged =
          opt -> opt.getEvaluations() < maxEval.getMaxEval();

      for (int i = 0; i <= 3; i++) {
        // Try from the initial guess
        PointValuePair solution =
            cmaesOptimizer.optimize(guess, fun, GoalType.MAXIMIZE, bounds, sigma, popSize, maxEval);
        if (converged.test(cmaesOptimizer)
            && (bestSolution == null || solution.getValue() > bestSolution.getValue())) {
          evaluations = cmaesOptimizer.getEvaluations();
          bestSolution = solution;
          LoggerUtils.log(logger, Level.FINE,
              "CMAES optimiser [%da] fit (N=2) : MLE = %f (%d evaluations)", i, solution.getValue(),
              evaluations);
        }

        if (bestSolution == null) {
          continue;
        }

        // Try from the current optimum
        solution = cmaesOptimizer.optimize(new InitialGuess(bestSolution.getPointRef()), fun,
            GoalType.MAXIMIZE, bounds, sigma, popSize, maxEval);
        if (converged.test(cmaesOptimizer) && (solution.getValue() > bestSolution.getValue())) {
          evaluations = cmaesOptimizer.getEvaluations();
          bestSolution = solution;
          LoggerUtils.log(logger, Level.FINE,
              "CMAES optimiser [%db] fit (N=2) : MLE = %f (%d evaluations)", i, solution.getValue(),
              evaluations);
        }
      }

      if (bestSolution != null) {
        try {
          // Re-optimise with Powell?
          final PointValuePair solution =
              powellOptimizer.optimize(fun, new InitialGuess(bestSolution.getPointRef()), bounds,
                  step, GoalType.MAXIMIZE, maxEval);
          if (solution.getValue() > bestSolution.getValue()) {
            evaluations = cmaesOptimizer.getEvaluations();
            bestSolution = solution;
            LoggerUtils.log(logger, Level.INFO,
                "Powell optimiser re-fit (N=2) : MLE = %f (%d evaluations)",
                bestSolution.getValue(), powellOptimizer.getEvaluations());
          }
        } catch (final TooManyEvaluationsException | ConvergenceException ignored) {
          // No solution
        }
      }
    }

    if (bestSolution == null) {
      LoggerUtils.log(logger, Level.INFO, "Failed to fit N=2");
      return null;
    }

    final double[] fitParams = bestSolution.getPointRef();
    final FitResult fr = new FitResult(Arrays.stream(count).sum(), 3, bestSolution.getValue());

    final double[] k = Arrays.copyOf(fitParams, 2);
    final double[] f = new double[] {fitParams[2], 1 - fitParams[2]};

    // Sort by size (ascending)
    SortUtils.sortData(f, k, true, false);

    final int eval = evaluations;
    logger.info(() -> String.format("Fit (N=2) : %s (%s), MLE = %s (%d evaluations)", formatK(k),
        format(f), MathUtils.rounded(fr.getLogLikelihood(), 4), eval));

    return Pair.of(fr, new Model2(k[0], k[1], f[0]));
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

  private static CMAESOptimizer createCmaesOptimizer() {
    final double rel = 1e-8;
    final double abs = 1e-10;
    final int maxIterations = 2000;
    final double stopFitness = 0;
    final boolean isActiveCma = true;
    final int diagonalOnly = 20;
    final int checkFeasableCount = 1;
    final RandomGenerator random = new RandomGeneratorAdapter(UniformRandomProviders.create());
    final boolean generateStatistics = false;
    final ConvergenceChecker<PointValuePair> checker = new SimpleValueChecker(rel, abs);

    // Iterate this for stability in the initial guess
    return new CMAESOptimizer(maxIterations, stopFitness, isActiveCma, diagonalOnly,
        checkFeasableCount, random, generateStatistics, checker);
  }

  private static String formatK(double[] data) {
    return Arrays.stream(data).mapToObj(x -> MathUtils.rounded(x, 4))
        .collect(Collectors.joining(", ", "Rate = ", ""));
  }

  private static String format(double[] data) {
    return Arrays.stream(data).mapToObj(x -> MathUtils.rounded(x, 4))
        .collect(Collectors.joining(", "));
  }

  /**
   * Get the first instance that is assignment compatible with the argument, or else return the
   * argument.
   *
   * @param <T> the generic type
   * @param value the argument
   * @param options the options
   * @return the instance
   */
  @SuppressWarnings("unchecked")
  private static <T> T orElse(T value, Object... options) {
    final Class<?> cls = value.getClass();
    for (final Object o : options) {
      if (cls.isInstance(o)) {
        return (T) o;
      }
    }
    return value;
  }

  /**
   * Get the first instance that is assignment compatible with the class argument, or else return
   * the argument value.
   *
   * @param <T> the generic type
   * @param cls the class
   * @param value the argument value
   * @param options the options
   * @return the instance
   */
  @SuppressWarnings("unchecked")
  private static <T> T orElse(Class<T> cls, T value, Object... options) {
    for (final Object o : options) {
      if (cls.isInstance(o)) {
        return (T) o;
      }
    }
    return value;
  }
}
