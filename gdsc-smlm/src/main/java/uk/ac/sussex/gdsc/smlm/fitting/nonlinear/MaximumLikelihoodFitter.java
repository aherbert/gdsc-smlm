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

package uk.ac.sussex.gdsc.smlm.fitting.nonlinear;

import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.MultivariateVectorFunction;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.exception.TooManyIterationsException;
import org.apache.commons.math3.optim.BaseOptimizer;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.SimpleValueChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateFunctionMappingAdapter;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunctionGradient;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer;
import org.apache.commons.math3.random.RandomGenerator;
import uk.ac.sussex.gdsc.core.utils.rng.RandomGeneratorAdapter;
import uk.ac.sussex.gdsc.core.utils.rng.UniformRandomProviders;
import uk.ac.sussex.gdsc.smlm.fitting.FisherInformationMatrix;
import uk.ac.sussex.gdsc.smlm.fitting.FitStatus;
import uk.ac.sussex.gdsc.smlm.function.FixedNonLinearFunction;
import uk.ac.sussex.gdsc.smlm.function.LikelihoodWrapper;
import uk.ac.sussex.gdsc.smlm.function.NonLinearFunction;
import uk.ac.sussex.gdsc.smlm.function.PoissonGammaGaussianLikelihoodWrapper;
import uk.ac.sussex.gdsc.smlm.function.PoissonGaussianLikelihoodWrapper;
import uk.ac.sussex.gdsc.smlm.function.PoissonLikelihoodWrapper;
import uk.ac.sussex.gdsc.smlm.math3.optim.nonlinear.scalar.gradient.BoundedNonLinearConjugateGradientOptimizer;
import uk.ac.sussex.gdsc.smlm.math3.optim.nonlinear.scalar.gradient.BoundedNonLinearConjugateGradientOptimizer.Formula;
import uk.ac.sussex.gdsc.smlm.math3.optim.nonlinear.scalar.noderiv.CustomPowellOptimizer;

/**
 * Uses Maximum Likelihood Estimation (MLE) to fit a nonlinear model with coefficients (a) for a set
 * of data points (x, y).
 *
 * <p>By default the probability mass function for observed value k is modelled as a Poisson
 * process:<br> pmf = e^-k.(l^k / k!) <br> where: <br> k = Observed number of occurrences <br> l =
 * Expected number of occurrences (the mean)
 *
 * <p>MLE = Max [ sum (ln(e^-k.(l^k / k!)) ] <br> = Max [ sum (k.ln(l) - l) ]
 *
 * <p>The expected number of occurrences can be modelled using any parameterised function, for
 * example the Gaussian 2D function.
 *
 * <p>The probability mass function can be changed to a Poisson-Gaussian or Poisson-Gamma-Gaussian
 * distribution in order to model the counts from a CCD/EMCCD camera.
 */
public class MaximumLikelihoodFitter extends MleBaseFunctionSolver {

  /** The search method. */
  SearchMethod searchMethod = SearchMethod.POWELL;
  private LikelihoodFunction likelihoodFunction = LikelihoodFunction.POISSON;
  private double alpha;
  private double sigma;

  private boolean gradientLineMinimisation = true;
  private double relativeThreshold = 1e-4;
  private double absoluteThreshold = 1e-10;
  private double[] lower;
  private double[] upper;
  private double[] lowerConstraint;
  private double[] upperConstraint;

  /**
   * The function to use for the Powell optimiser (which may have parameters mapped using the sqrt
   * function.
   */
  private MultivariateLikelihood powellFunction;

  /** Maximum iterations for the Powell optimiser. */
  private int maxIterations;

  /**
   * Wrap the LikelihoodFunction with classes that implement the required interfaces.
   */
  private static class Likelihood {
    LikelihoodWrapper fun;

    Likelihood(LikelihoodWrapper fun) {
      this.fun = fun;
    }
  }

  /**
   * Wrap the likelihood function as a MultivariateFunction.
   */
  private static class MultivariateLikelihood extends Likelihood implements MultivariateFunction {
    MultivariateLikelihood(LikelihoodWrapper fun) {
      super(fun);
    }

    @Override
    public double value(double[] point) {
      return fun.likelihood(point);
    }

    boolean isMapped() {
      return false;
    }
  }

  /**
   * Map the specified indices using the sqrt function for use with the Powell optimiser.
   */
  private static class MappedMultivariateLikelihood extends MultivariateLikelihood {
    final int[] map;

    MappedMultivariateLikelihood(LikelihoodWrapper fun, int[] map) {
      super(fun);
      this.map = map;
    }

    @Override
    public double value(double[] point) {
      return fun.likelihood(unmap(point));
    }

    /**
     * Convert the unmapped point to the mapped equivalent. The mapped point is used by the Powell
     * optimiser.
     *
     * <p>This is done by square rooting the value of the mapped indices.
     *
     * @param point the point
     * @return The mapped point
     */
    double[] map(double[] point) {
      point = point.clone();
      for (final int i : map) {
        point[i] = Math.sqrt(Math.abs(point[i])) * Math.signum(point[i]);
      }
      return point;
    }

    /**
     * Convert the mapped point to the unmapped equivalent. The unmapped point is used to evaluate
     * the function.
     *
     * <p>This is done by squaring the value of the mapped indices.
     *
     * @param point the point
     * @return The unmapped point
     */
    double[] unmap(double[] point) {
      point = point.clone();
      for (final int i : map) {
        if (point[i] > 0) {
          point[i] = point[i] * point[i];
        } else {
          point[i] = -(point[i] * point[i]);
        }
      }
      return point;
    }

    @Override
    public boolean isMapped() {
      return true;
    }
  }

  /**
   * Wrap the likelihood function as a MultivariateVectorFunction.
   */
  private static class MultivariateVectorLikelihood extends Likelihood
      implements MultivariateVectorFunction {
    MultivariateVectorLikelihood(LikelihoodWrapper fun) {
      super(fun);
    }

    @Override
    public double[] value(double[] point) {
      final double[] gradient = new double[point.length];
      fun.likelihood(point, gradient);
      return gradient;
    }
  }

  /**
   * The search method.
   */
  public enum SearchMethod {
    /**
     * Search using Powell's conjugate direction method using hard limits to ensure a bounded
     * search.
     */
    POWELL_BOUNDED("Powell (bounded)", false),
    /**
     * Search using Powell's conjugate direction method.
     */
    POWELL("Powell", false),
    /**
     * Search using Powell's conjugate direction method using a mapping adapter to ensure a bounded
     * search.
     */
    POWELL_ADAPTER("Powell (adapter)", false),
    /**
     * Search using Powell's Bound Optimization BY Quadratic Approximation (BOBYQA) algorithm.
     *
     * <p>BOBYQA could also be considered as a replacement of any derivative-based optimizer when
     * the derivatives are approximated by finite differences. This is a bounded search.
     */
    BOBYQA("BOBYQA", false),
    /**
     * Search using active Covariance Matrix Adaptation Evolution Strategy (CMA-ES).
     *
     * <p>The CMA-ES is a reliable stochastic optimization method which should be applied if
     * derivative-based methods, e.g. conjugate gradient, fail due to a rugged search landscape.
     * This is a bounded search.
     */
    CMAES("CMAES", false),
    /**
     * Search using a non-linear conjugate gradient optimiser. Use the Fletcher-Reeves update
     * formulas for the conjugate search directions.
     *
     * <p>This is a bounded search using simple truncation of coordinates at the bounds of the
     * search space.
     */
    CONJUGATE_GRADIENT_FR("Conjugate Gradient Fletcher-Reeves", true),
    /**
     * Search using a non-linear conjugate gradient optimiser. Use the Polak-Ribière update formulas
     * for the conjugate search directions.
     *
     * <p>This is a bounded search using simple truncation of coordinates at the bounds of the
     * search space.
     */
    CONJUGATE_GRADIENT_PR("Conjugate Gradient Polak-Ribière", true);

    private final String name;
    private final boolean usesGradient;

    SearchMethod(String name, boolean usesGradient) {
      this.name = name;
      this.usesGradient = usesGradient;
    }

    @Override
    public String toString() {
      return name;
    }

    /**
     * Check if the search method uses the gradient of the likelihood function.
     *
     * @return True if the search method uses the gradient of the likelihood function.
     */
    public boolean usesGradients() {
      return usesGradient;
    }
  }

  /**
   * The likelihood function.
   */
  public enum LikelihoodFunction {
    /**
     * Use a Poisson likelihood model.
     */
    POISSON("Poisson"),
    /**
     * Use a Poisson likelihood model.
     */
    POISSON_GAUSSIAN("Poisson+Gaussian"),
    /**
     * Use a Poisson likelihood model.
     */
    POISSON_GAMMA_GAUSSIAN("Poisson+Gamma+Gaussian");

    private final String name;

    LikelihoodFunction(String name) {
      this.name = name;
    }

    @Override
    public String toString() {
      return name;
    }
  }

  /**
   * Default constructor.
   *
   * @param function The function
   */
  public MaximumLikelihoodFitter(NonLinearFunction function) {
    super(function);
  }

  @Override
  public FitStatus computeFit(double[] y, double[] fx, double[] a, double[] parametersVariance) {
    final int n = y.length;
    final LikelihoodWrapper maximumLikelihoodFunction =
        createLikelihoodWrapper((NonLinearFunction) function, n, y, a);

    @SuppressWarnings("rawtypes")
    BaseOptimizer baseOptimiser = null;

    try {
      final double[] startPoint = getInitialSolution(a);

      PointValuePair optimum = null;
      if (searchMethod == SearchMethod.POWELL || searchMethod == SearchMethod.POWELL_BOUNDED
          || searchMethod == SearchMethod.POWELL_ADAPTER) {
        // Non-differentiable version using Powell Optimiser

        // Background: see Numerical Recipes 10.5 (Direction Set (Powell's) method).
        // The optimiser could be extended to implement bounds on the directions moved.
        // However the mapping adapter seems to work OK.

        final boolean basisConvergence = false;

        // Perhaps these thresholds should be tighter?
        // The default is to use the sqrt() of the overall tolerance
        // final double lineRel = Math.sqrt(relativeThreshold);
        // final double lineAbs = Math.sqrt(absoluteThreshold);
        // final double lineRel = relativeThreshold * 1e2;
        // final double lineAbs = absoluteThreshold * 1e2;

        // Since we are fitting only a small number of parameters then just use the same tolerance
        // for each search direction
        final double lineRel = relativeThreshold;
        final double lineAbs = absoluteThreshold;

        final CustomPowellOptimizer o = new CustomPowellOptimizer(relativeThreshold,
            absoluteThreshold, lineRel, lineAbs, null, basisConvergence);
        baseOptimiser = o;

        OptimizationData maxIterationData = null;
        if (getMaxIterations() > 0) {
          maxIterationData = new MaxIter(getMaxIterations());
        }

        if (searchMethod == SearchMethod.POWELL_ADAPTER) {
          // Try using the mapping adapter for a bounded Powell search
          final MultivariateFunctionMappingAdapter adapter = new MultivariateFunctionMappingAdapter(
              new MultivariateLikelihood(maximumLikelihoodFunction), lower, upper);
          optimum = o.optimize(maxIterationData, new MaxEval(getMaxEvaluations()),
              new ObjectiveFunction(adapter), GoalType.MINIMIZE,
              new InitialGuess(adapter.boundedToUnbounded(startPoint)));
          final double[] solution = adapter.unboundedToBounded(optimum.getPointRef());
          optimum = new PointValuePair(solution, optimum.getValue());
        } else {
          if (powellFunction == null) {
            powellFunction = new MultivariateLikelihood(maximumLikelihoodFunction);
          }

          // Update the maximum likelihood function in the Powell function wrapper
          powellFunction.fun = maximumLikelihoodFunction;

          final OptimizationData positionChecker = null;
          // new org.apache.commons.math3.optim.PositionChecker(relativeThreshold,
          // absoluteThreshold);
          SimpleBounds simpleBounds = null;
          if (powellFunction.isMapped()) {
            final MappedMultivariateLikelihood adapter =
                (MappedMultivariateLikelihood) powellFunction;
            if (searchMethod == SearchMethod.POWELL_BOUNDED) {
              simpleBounds = new SimpleBounds(adapter.map(lower), adapter.map(upper));
            }
            optimum = o.optimize(maxIterationData, new MaxEval(getMaxEvaluations()),
                new ObjectiveFunction(powellFunction), GoalType.MINIMIZE,
                new InitialGuess(adapter.map(startPoint)), positionChecker, simpleBounds);
            final double[] solution = adapter.unmap(optimum.getPointRef());
            optimum = new PointValuePair(solution, optimum.getValue());
          } else {
            if (searchMethod == SearchMethod.POWELL_BOUNDED) {
              simpleBounds = new SimpleBounds(lower, upper);
            }
            optimum = o.optimize(maxIterationData, new MaxEval(getMaxEvaluations()),
                new ObjectiveFunction(powellFunction), GoalType.MINIMIZE,
                new InitialGuess(startPoint), positionChecker, simpleBounds);
          }
        }
      } else if (searchMethod == SearchMethod.BOBYQA) {
        // Differentiable approximation using Powell's BOBYQA algorithm.
        // This is slower than the Powell optimiser and requires a high number of evaluations.
        final int numberOfInterpolationpoints = this.getNumberOfFittedParameters() + 2;

        final BOBYQAOptimizer o = new BOBYQAOptimizer(numberOfInterpolationpoints);
        baseOptimiser = o;
        optimum = o.optimize(new MaxEval(getMaxEvaluations()),
            new ObjectiveFunction(new MultivariateLikelihood(maximumLikelihoodFunction)),
            GoalType.MINIMIZE, new InitialGuess(startPoint), new SimpleBounds(lower, upper));
      } else if (searchMethod == SearchMethod.CMAES) {
        // TODO - Understand why the CMAES optimiser does not fit very well on test data. It appears
        // to converge too early and the likelihood scores are not as low as the other optimisers.

        // The CMAESOptimiser is based on Matlab code:
        // https://www.lri.fr/~hansen/cmaes.m
        // Take the defaults from the Matlab documentation
        final double stopFitness = 0;
        final boolean isActiveCma = true;
        final int diagonalOnly = 0;
        final int checkFeasableCount = 1;
        final RandomGenerator random = new RandomGeneratorAdapter(UniformRandomProviders.create());
        final boolean generateStatistics = false;
        // The sigma determines the search range for the variables. It should be 1/3 of the initial
        // search region.
        final double[] sigma = new double[lower.length];
        for (int i = 0; i < sigma.length; i++) {
          sigma[i] = (upper[i] - lower[i]) / 3;
        }
        int popSize = (int) (4 + Math.floor(3 * Math.log(sigma.length)));

        // The CMAES optimiser is random and restarting can overcome problems with quick
        // convergence.
        // The Apache commons documentations states that convergence should occur between 30N and
        // 300N^2
        // function evaluations
        final int n30 = Math.min(sigma.length * sigma.length * 30, getMaxEvaluations() / 2);
        evaluations = 0;
        final OptimizationData[] data =
            {new InitialGuess(startPoint), new CMAESOptimizer.PopulationSize(popSize),
                new MaxEval(getMaxEvaluations()), new CMAESOptimizer.Sigma(sigma),
                new ObjectiveFunction(new MultivariateLikelihood(maximumLikelihoodFunction)),
                GoalType.MINIMIZE, new SimpleBounds(lower, upper)};
        // Iterate to prevent early convergence
        int repeat = 0;
        while (evaluations < n30) {
          if (repeat++ > 1) {
            // Update the start point and population size
            if (optimum != null) {
              data[0] = new InitialGuess(optimum.getPointRef());
            }
            popSize *= 2;
            data[1] = new CMAESOptimizer.PopulationSize(popSize);
          }
          final CMAESOptimizer o = new CMAESOptimizer(getMaxIterations(), stopFitness, isActiveCma,
              diagonalOnly, checkFeasableCount, random, generateStatistics,
              new SimpleValueChecker(relativeThreshold, absoluteThreshold));
          baseOptimiser = o;
          final PointValuePair result = o.optimize(data);
          iterations += o.getIterations();
          evaluations += o.getEvaluations();
          if (optimum == null || result.getValue() < optimum.getValue()) {
            optimum = result;
          }
        }
        // Prevent incrementing the iterations again
        baseOptimiser = null;
      } else {
        // The line search algorithm often fails. This is due to searching into a region where the
        // function evaluates to a negative so has been clipped. This means the upper bound of the
        // line cannot be found.
        // Note that running it on an easy problem (200 photons with fixed fitting (no background))
        // the algorithm does sometimes produces results better than the Powell algorithm but it is
        // slower.

        final BoundedNonLinearConjugateGradientOptimizer o =
            new BoundedNonLinearConjugateGradientOptimizer(
                (searchMethod == SearchMethod.CONJUGATE_GRADIENT_FR) ? Formula.FLETCHER_REEVES
                    : Formula.POLAK_RIBIERE,
                new SimpleValueChecker(relativeThreshold, absoluteThreshold));
        baseOptimiser = o;

        // Note: The gradients may become unstable at the edge of the bounds. Or they will not
        // change direction if the true solution is on the bounds since the gradient will always
        // continue towards the bounds. This is key to the conjugate gradient method. It searches
        // along a vector until the direction of the gradient is in the opposite direction (using
        // dot products, i.e. cosine of angle between them)

        // NR 10.7 states there is no advantage of the variable metric DFP or BFGS methods over
        // conjugate gradient methods. So I will try these first.

        // Try this:
        // Adapt the conjugate gradient optimiser to use the gradient to pick the search direction
        // and then for the line minimisation. However if the function is out of bounds then clip
        // the variables at the bounds and continue.
        // If the current point is at the bounds and the gradient is to continue out of bounds then
        // clip the gradient too.

        // Or: just use the gradient for the search direction then use the line minimisation/rest
        // as per the Powell optimiser. The bounds should limit the search.

        // I tried a Bounded conjugate gradient optimiser with clipped variables:
        // This sometimes works. However when the variables go a long way out of the expected range
        // the gradients can have vastly different magnitudes. This results in the algorithm
        // stalling since the gradients can be close to zero and the some of the parameters are no
        // longer adjusted. Perhaps this can be looked for and the algorithm then gives up and
        // resorts to a Powell optimiser from the current point.

        // Changed the bracketing step to very small (default is 1, changed to 0.001). This improves
        // the performance. The gradient direction is very sensitive to small changes in the
        // coordinates so a tighter bracketing of the line search helps.

        // Tried using a non-gradient method for the line search copied from the Powell optimiser:
        // This also works when the bracketing step is small but the number of iterations is higher.

        // 24.10.2014: I have tried to get conjugate gradient to work but the gradient function
        // must not behave suitably for the optimiser. In the current state both methods of using a
        // Bounded Conjugate Gradient Optimiser perform poorly relative to other optimisers:
        // Simulated : n=1000, signal=200, x=0.53, y=0.47
        // LVM : n=1000, signal=171, x=0.537, y=0.471 (1.003s)
        // Powell : n=1000, signal=187, x=0.537, y=0.48 (1.238s)
        // Gradient based PR (constrained): n=858, signal=161, x=0.533, y=0.474 (2.54s)
        // Gradient based PR (bounded): n=948, signal=161, x=0.533, y=0.473 (2.67s)
        // Non-gradient based : n=1000, signal=151.47, x=0.535, y=0.474 (1.626s)
        // The conjugate optimisers are slower, under predict the signal by the most and in the case
        // of the gradient based optimiser, fail to converge on some problems. This is worse when
        // constrained fitting is used and not tightly bounded fitting.
        // I will leave the code in as an option but would not recommend using it. I may remove it
        // in the future.

        // Note: It is strange that the non-gradient based line minimisation is more successful.
        // It may be that the gradient function is not accurate (due to round off error) or that it
        // is simply wrong when far from the optimum. My JUnit tests only evaluate the function
        // within the expected range of the answer.

        // Note the default step size on the Powell optimiser is 1 but the initial directions are
        // unit vectors.
        // So our bracketing step should be a minimum of 1 / average length of the first gradient
        // vector to prevent the first step being too large when bracketing.
        final double[] gradient = new double[startPoint.length];
        maximumLikelihoodFunction.likelihood(startPoint, gradient);
        double length = 0;
        for (final double d : gradient) {
          length += d * d;
        }
        final double bracketingStep = Math.min(0.001, ((length > 1) ? 1.0 / length : 1));

        o.setUseGradientLineSearch(gradientLineMinimisation);

        optimum = o.optimize(new MaxEval(getMaxEvaluations()),
            new ObjectiveFunctionGradient(
                new MultivariateVectorLikelihood(maximumLikelihoodFunction)),
            new ObjectiveFunction(new MultivariateLikelihood(maximumLikelihoodFunction)),
            GoalType.MINIMIZE, new InitialGuess(startPoint),
            new SimpleBounds(lowerConstraint, upperConstraint),
            new BoundedNonLinearConjugateGradientOptimizer.BracketingStep(bracketingStep));
      }

      if (optimum == null) {
        return FitStatus.FAILED_TO_CONVERGE;
      }

      final double[] solution = optimum.getPointRef();

      setSolution(a, solution);

      if (parametersVariance != null) {
        // Compute assuming a Poisson process.
        // Note:
        // If using a Poisson-Gamma-Gaussian model then these will be incorrect.
        // However the variance for the position estimates of a 2D PSF can be
        // scaled by a factor of 2 as in Mortensen, et al (2010) Nature Methods 7, 377-383, SI 4.3.
        // Since the type of function is unknown this cannot be done here.
        final FisherInformationMatrix m =
            new FisherInformationMatrix(maximumLikelihoodFunction.fisherInformation(solution));
        setDeviations(parametersVariance, m);
      }

      // Reverse negative log likelihood for maximum likelihood score
      value = -optimum.getValue();
    } catch (final TooManyIterationsException ex) {
      return FitStatus.TOO_MANY_ITERATIONS;
    } catch (final TooManyEvaluationsException ex) {
      return FitStatus.TOO_MANY_EVALUATIONS;
    } catch (final ConvergenceException ex) {
      // Occurs when QR decomposition fails - mark as a singular non-linear model (no solution)
      return FitStatus.SINGULAR_NON_LINEAR_MODEL;
    } catch (final Exception ex) {
      Logger.getLogger(getClass().getName()).log(Level.SEVERE, "Failed to fit", ex);
      return FitStatus.UNKNOWN;
    } finally {
      if (baseOptimiser != null) {
        iterations += baseOptimiser.getIterations();
        evaluations += baseOptimiser.getEvaluations();
      }
    }

    // Check this as likelihood functions can go wrong
    if (Double.isInfinite(value) || Double.isNaN(value)) {
      return FitStatus.INVALID_LIKELIHOOD;
    }

    return FitStatus.OK;
  }

  private LikelihoodWrapper createLikelihoodWrapper(NonLinearFunction function, int n, double[] y,
      double[] a) {
    LikelihoodWrapper maximumLikelihoodFunction = null;

    final double myAlpha = (this.alpha > 0) ? this.alpha : 1;

    // We can use different likelihood wrapper functions:
    switch (likelihoodFunction) {
      case POISSON_GAMMA_GAUSSIAN:
        // Poisson-Gamma-Gaussian - EM-CCD data
        maximumLikelihoodFunction =
            new PoissonGammaGaussianLikelihoodWrapper(function, a, y, n, myAlpha, sigma);
        break;

      case POISSON_GAUSSIAN:
        // Poisson-Gaussian - CCD data
        // Sigma must be positive, otherwise fall back to a Poisson likelihood function
        if (sigma > 0) {
          maximumLikelihoodFunction =
              new PoissonGaussianLikelihoodWrapper(function, a, y, n, myAlpha, sigma);
        }
        break;

      default:
        // Poisson - most counting data
        break;
    }

    // Check if the method requires the gradient but it cannot be computed
    if (maximumLikelihoodFunction == null
        || (searchMethod.usesGradients() && !maximumLikelihoodFunction.canComputeGradient())) {
      // Ensure no negative data for the Poisson likelihood method.
      // Just truncate the counts for now. These are from noise in the count estimates that we do
      // not model.
      final double[] y2 = new double[n];
      for (int i = 0; i < n; i++) {
        if (y[i] < 0) {
          y2[i] = 0;
        } else {
          y2[i] = y[i];
        }
      }
      final PoissonLikelihoodWrapper mlFunction =
          new PoissonLikelihoodWrapper(function, a, y2, n, myAlpha);
      // This will allow Powell searches. The effect on the gradient search algorithms may be weird
      // so leave alone.
      if (!searchMethod.usesGradients()) {
        mlFunction.setAllowNegativeExpectedValues(true);
      }
      maximumLikelihoodFunction = mlFunction;
    }
    return maximumLikelihoodFunction;
  }

  /**
   * Gets the max iterations for the Powell search method.
   *
   * @return the max iterations for the Powell search method.
   */
  public int getMaxIterations() {
    return maxIterations;
  }

  /**
   * Sets the max iterations for the Powell search method.
   *
   * @param maxIterations the max iterations for the Powell search method
   */
  public void setMaxIterations(int maxIterations) {
    this.maxIterations = maxIterations;
  }

  /**
   * Gets the search method.
   *
   * @return the search method.
   */
  public SearchMethod getSearchMethod() {
    return searchMethod;
  }

  /**
   * Sets the search method.
   *
   * @param searchMethod the search method
   */
  public void setSearchMethod(SearchMethod searchMethod) {
    this.searchMethod = searchMethod;
  }

  /**
   * Gets the likelihood function.
   *
   * @return the likelihood function
   */
  public LikelihoodFunction getLikelihoodFunction() {
    return likelihoodFunction;
  }

  /**
   * Sets the likelihood function.
   *
   * @param likelihoodFunction the likelihood function
   */
  public void setLikelihoodFunction(LikelihoodFunction likelihoodFunction) {
    this.likelihoodFunction = likelihoodFunction;
  }

  /**
   * Gets the alpha for the gamma component of the Poisson-Gamma-Gaussian likelihood function.
   *
   * @return the alpha for the gamma component of the Poisson-Gamma-Gaussian likelihood function.
   */
  public double getAlpha() {
    return alpha;
  }

  /**
   * Sets the alpha for the gamma component of the Poisson-Gamma-Gaussian likelihood function.
   *
   * @param alpha the alpha for the gamma component of the Poisson-Gamma-Gaussian likelihood
   *        function
   */
  public void setAlpha(double alpha) {
    this.alpha = alpha;
  }

  /**
   * Gets the sigma for the Gaussian component of the Poisson-Gaussian/Poisson-Gamma-Gaussian
   * likelihood function.
   *
   * @return the sigma for the Gaussian component of the Poisson-Gaussian/Poisson-Gamma-Gaussian
   *         likelihood function
   */
  public double getSigma() {
    return sigma;
  }

  /**
   * Sets the sigma for the Gaussian component of the Poisson-Gaussian/Poisson-Gamma-Gaussian
   * likelihood function.
   *
   * @param sigma the sigma for the Gaussian component of the
   *        Poisson-Gaussian/Poisson-Gamma-Gaussian likelihood function
   */
  public void setSigma(double sigma) {
    this.sigma = sigma;
  }

  /**
   * This setting applies to the conjugate gradient method of the Maximum Likelihood Estimator.
   *
   * @return the gradientLineMinimisation True if using the gradient for line minimisation
   */
  public boolean isGradientLineMinimisation() {
    return gradientLineMinimisation;
  }

  /**
   * This setting applies to the conjugate gradient method of the Maximum Likelihood Estimator.
   *
   * @param gradientLineMinimisation Set to true to use the gradient for line minimisation
   */
  public void setGradientLineMinimisation(boolean gradientLineMinimisation) {
    this.gradientLineMinimisation = gradientLineMinimisation;
  }

  /**
   * Gets the relative threshold for convergence in the Maximum Likelihood Estimator.
   *
   * @return the relative threshold
   */
  public double getRelativeThreshold() {
    return relativeThreshold;
  }

  /**
   * Sets the relative threshold for convergence in the Maximum Likelihood Estimator.
   *
   * @param relativeThreshold the relative threshold
   */
  public void setRelativeThreshold(double relativeThreshold) {
    this.relativeThreshold = relativeThreshold;
  }

  /**
   * Gets the absolute threshold for convergence in the Maximum Likelihood Estimator.
   *
   * @return the absolute threshold
   */
  public double getAbsoluteThreshold() {
    return absoluteThreshold;
  }

  /**
   * Sets the absolute threshold for convergence in the Maximum Likelihood Estimator.
   *
   * @param absoluteThreshold the absolute threshold
   */
  public void setAbsoluteThreshold(double absoluteThreshold) {
    this.absoluteThreshold = absoluteThreshold;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Note that certain search methods require constraints to function. A null pointer exception
   * will be thrown by the fitter if the constraints are not set for these methods.
   */
  @Override
  public boolean isBounded() {
    switch (searchMethod) {
      case POWELL_ADAPTER:
      case POWELL_BOUNDED:
      case BOBYQA:
      case CMAES:
        return true;
      default:
        return false;
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>Note that certain search methods require constraints to function. A null pointer exception
   * will be thrown by the fitter if the constraints are not set for these methods.
   */
  @Override
  public boolean isConstrained() {
    switch (searchMethod) {
      case CONJUGATE_GRADIENT_FR:
      case CONJUGATE_GRADIENT_PR:
        return true;
      default:
        return false;
    }
  }

  @Override
  public void setBounds(double[] lowerB, double[] upperB) {
    // Extract the bounds for the parameters we are fitting
    final int[] indices = function.gradientIndices();

    lower = new double[indices.length];
    upper = new double[indices.length];
    for (int i = 0; i < indices.length; i++) {
      lower[i] = lowerB[indices[i]];
      upper[i] = upperB[indices[i]];
    }
  }

  @Override
  public void setConstraints(double[] lowerB, double[] upperB) {
    // Extract the bounds for the parameters we are fitting
    final int[] indices = function.gradientIndices();

    lowerConstraint = new double[indices.length];
    upperConstraint = new double[indices.length];
    for (int i = 0; i < indices.length; i++) {
      lowerConstraint[i] = lowerB[indices[i]];
      upperConstraint[i] = upperB[indices[i]];
    }
  }

  @Override
  public boolean computeValue(double[] y, double[] fx, double[] a) {
    final int n = y.length;
    final LikelihoodWrapper maximumLikelihoodFunction =
        createLikelihoodWrapper((NonLinearFunction) function, n, y, a);

    final double l = maximumLikelihoodFunction.likelihood(a);
    if (l == Double.POSITIVE_INFINITY) {
      return false;
    }

    // Reverse negative log likelihood for maximum likelihood score
    value = -l;

    return true;
  }

  @Override
  protected FisherInformationMatrix computeFisherInformationMatrix(double[] y, double[] a) {
    final int n = y.length;
    final LikelihoodWrapper maximumLikelihoodFunction =
        createLikelihoodWrapper((NonLinearFunction) function, n, y, a);
    final double[] solution = getInitialSolution(a);
    // Compute assuming a Poisson process.
    // Note:
    // If using a Poisson-Gamma-Gaussian model then these will be incorrect.
    // However the variance for the position estimates of a 2D PSF can be
    // scaled by a factor of 2 as in Mortensen, et al (2010) Nature Methods 7, 377-383, SI 4.3.
    // Since the type of function is unknown this cannot be done here.
    return new FisherInformationMatrix(maximumLikelihoodFunction.fisherInformation(solution));
  }

  @Override
  protected double computeObservedLogLikelihood(double[] y, double[] a) {
    if (lastY != null) {
      final int n = y.length;
      // The function value must be scaled to the expected value of a Poisson process.
      double[] scaledy;
      if (alpha == 0) {
        scaledy = y;
      } else {
        scaledy = new double[n];
        for (int i = n; i-- > 0;) {
          scaledy[i] = y[i] * alpha;
        }
      }
      final LikelihoodWrapper maximumLikelihoodFunction =
          createLikelihoodWrapper(new FixedNonLinearFunction(scaledy), n, lastY, a);

      final double l = maximumLikelihoodFunction.likelihood(a);
      if (l == Double.POSITIVE_INFINITY) {
        return Double.NEGATIVE_INFINITY;
      }

      // Reverse negative log likelihood for maximum likelihood score
      value = -l;
    }
    throw new IllegalStateException();
  }
}
