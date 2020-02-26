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

package uk.ac.sussex.gdsc.smlm.fitting.nonlinear;

import uk.ac.sussex.gdsc.core.utils.DoubleEquality;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.fitting.FisherInformationMatrix;
import uk.ac.sussex.gdsc.smlm.fitting.FitStatus;
import uk.ac.sussex.gdsc.smlm.fitting.FunctionSolverType;
import uk.ac.sussex.gdsc.smlm.fitting.MleFunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.WLseFunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.linear.EjmlLinearSolver;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.GradientCalculatorUtils;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.stop.ErrorStoppingCriteria;
import uk.ac.sussex.gdsc.smlm.function.ChiSquaredDistributionTable;
import uk.ac.sussex.gdsc.smlm.function.GradientFunction;
import uk.ac.sussex.gdsc.smlm.function.NonLinearFunction;
import uk.ac.sussex.gdsc.smlm.function.PoissonCalculator;

/**
 * Uses Levenberg-Marquardt method to fit a nonlinear model with coefficients (a) for a set of data
 * points (x, y).
 *
 * <p>The MLEFunctionSolver is supported if the flag to use a Poisson MLE model is set. If the
 * function supports weights then the WLSEFunctionSolver is supported. The default implementation
 * supports the LSEFunctionSolver.
 */
public class NonLinearFit extends LseBaseFunctionSolver
    implements MleFunctionSolver, WLseFunctionSolver {
  /**
   * Index for the best sum-of-squares in {@link #sumOfSquaresWorking}.
   */
  protected static final int SUM_OF_SQUARES_BEST = 0;
  /**
   * Index for the new sum-of-squares in {@link #sumOfSquaresWorking}.
   */
  protected static final int SUM_OF_SQUARES_OLD = 1;
  /**
   * Index for the previous sum-of-squares in {@link #sumOfSquaresWorking}.
   */
  protected static final int SUM_OF_SQUARES_NEW = 2;

  /** The solver. */
  protected EjmlLinearSolver solver = new EjmlLinearSolver();
  /** The calculator. */
  protected GradientCalculator calculator;
  /** The stopping criteria. */
  protected StoppingCriteria sc;

  /** The beta. */
  protected double[] beta = new double[0];

  /** Working space for beta. */
  protected double[] da;

  /**
   * The updated parameters a. This is equal to the current parameters a plus the solution x to A x
   * = b with A = {@link #covar} and b = gradient vector (beta).
   */
  protected double[] ap = new double[0];

  /**
   * Working space for the alpha matrix. This is equal to {@link #alpha} with the diagonal scaled by
   * 1 + {@link #lambda}.
   */
  protected double[][] covar;

  /** The alpha matrix. */
  protected double[][] alpha;

  /** The initial lambda. */
  protected double initialLambda = 0.01;

  /** The working lambda. */
  protected double lambda;

  /** The sum of squares results. */
  protected double[] sumOfSquaresWorking;

  /** The initial residual sum of squares. */
  protected double initialResidualSumOfSquares;

  /** The function. */
  protected NonLinearFunction func;

  /** The y data from the last successful fit. */
  protected double[] lastFx;

  /** The log-likelihood. Used for the MLE LVM algorithm. */
  protected double ll = Double.NaN;

  /**
   * Default constructor.
   *
   * @param func The function to fit
   */
  public NonLinearFit(NonLinearFunction func) {
    this(func, null);
  }

  /**
   * Default constructor.
   *
   * @param func The function to fit
   * @param sc The stopping criteria
   */
  public NonLinearFit(NonLinearFunction func, StoppingCriteria sc) {
    this(func, sc, 1e-3, 1e-10);
  }

  /**
   * Default constructor.
   *
   * @param func The function to fit
   * @param sc The stopping criteria
   * @param maxRelativeError Validate the Levenberg-Marquardt fit solution using the specified
   *        maximum relative error
   * @param maxAbsoluteError Validate the Levenberg-Marquardt fit solution using the specified
   *        maximum absolute error
   */
  public NonLinearFit(NonLinearFunction func, StoppingCriteria sc, double maxRelativeError,
      double maxAbsoluteError) {
    super(func);
    this.func = func;
    init(sc, maxRelativeError, maxAbsoluteError);
  }

  @Override
  protected void preProcess() {
    super.preProcess();
    ll = Double.NaN;
  }

  private void init(StoppingCriteria sc, double maxRelativeError, double maxAbsoluteError) {
    setStoppingCriteria(sc);
    solver.setEqual(new DoubleEquality(maxRelativeError, maxAbsoluteError));
  }

  /**
   * Compute a step for the LVM non linear model.
   *
   * @param n the number of data points
   * @param y the data
   * @param a the parameters
   * @param initialStage the initial stage flag
   * @return true, if successful
   */
  protected boolean nonLinearModel(int n, double[] y, double[] a, boolean initialStage) {
    // The NonLinearFunction evaluates a function with parameters a but only computes the gradient
    // for m <= a.length parameters. The parameters can be accessed using the gradientIndices()
    // method.

    final int[] gradientIndices = function.gradientIndices();
    final int m = gradientIndices.length;

    if (initialStage) {
      lambda = initialLambda;
      for (int j = a.length; j-- > 0;) {
        ap[j] = a[j];
      }
      sumOfSquaresWorking[SUM_OF_SQUARES_BEST] =
          calculator.findLinearised(n, y, a, alpha, beta, func);
      initialResidualSumOfSquares = sumOfSquaresWorking[SUM_OF_SQUARES_BEST];
      if (calculator.isNaNGradients()) {
        return false;
      }
    }

    // Set previous using the current best fit result we have
    sumOfSquaresWorking[SUM_OF_SQUARES_OLD] = sumOfSquaresWorking[SUM_OF_SQUARES_BEST];

    // Solve the gradient equation A x = b:
    // A = Hessian matrix (alpha)
    // x = Parameter shift (output da)
    // b = Gradient vector (beta)
    if (!solve(a, m)) {
      return false;
    }

    // Update the parameters. Ensure to use the gradient indices to update the correct parameters
    updateFitParameters(a, gradientIndices, m, da, ap);

    sumOfSquaresWorking[SUM_OF_SQUARES_NEW] = calculator.findLinearised(n, y, ap, covar, da, func);

    if (calculator.isNaNGradients()) {
      return false; // Stop now
    } else if (sumOfSquaresWorking[SUM_OF_SQUARES_NEW] < sumOfSquaresWorking[SUM_OF_SQUARES_OLD]) {
      accepted(a, ap, m);
    } else {
      increaseLambda();
    }

    return true;
  }

  /**
   * Called when there was a successful improvement of the fit. The lambda parameter should be
   * reduced.
   *
   * @param a The old fit parameters
   * @param ap The new fit parameters
   * @param np the number of fitted parameters (matches gradient indicies length)
   */
  protected void accepted(double[] a, double[] ap, int np) {
    decreaseLambda();

    for (int i = 0; i < np; i++) {
      for (int j = np; j-- > 0;) {
        alpha[i][j] = covar[i][j];
      }
    }

    for (int j = np; j-- > 0;) {
      beta[j] = da[j];
    }
    for (int j = a.length; j-- > 0;) {
      a[j] = ap[j];
    }
    sumOfSquaresWorking[SUM_OF_SQUARES_BEST] = sumOfSquaresWorking[SUM_OF_SQUARES_NEW];
  }

  /**
   * Decrease lambda. Call this when the solution to the matrix improved the score.
   */
  protected void decreaseLambda() {
    lambda *= 0.1;
  }

  /**
   * Increase lambda. Call this when the solution to the matrix do not improve the score, or the
   * matrix had no solution.
   */
  protected void increaseLambda() {
    lambda *= 10.0;
  }

  /**
   * Solve the gradient equation A x = b: *
   *
   * <pre>
   * A = Hessian matrix (alpha)
   * x = Parameter shift (output da)
   * b = Gradient vector (beta)
   * </pre>
   *
   * <p>The Hessian and gradient parameter from the current best scoring parameter set are assumed
   * to be in alpha and beta. The lambda parameter is used to weight the diagonal of the Hessian.
   *
   * @param a the current fit parameters
   * @param np the number of fit parameters
   * @return true, if successful
   */
  protected boolean solve(double[] a, final int np) {
    createLinearProblem(np);
    return solve(covar, da);
  }

  /**
   * Solves (one) linear equation, a x = b.
   *
   * <p>On input have a[n][n], b[n]. On output b replaced by x[n].
   *
   * <p>Note: Any zero elements in b are not solved.
   *
   * @param a the a
   * @param b the b
   * @return False if the equation is singular (no solution)
   */
  // CHECKSTYLE.OFF: ParameterName
  protected boolean solve(double[][] a, double[] b) {
    // CHECKSTYLE.ON: ParameterName

    // TODO
    // Q. Do we need a better qr decomposition that uses the largest Eigen column first.
    // There is a version from Apache commons math.
    // We could assess the magnitude of each value in the gradient vector and rearrange.

    return solver.solve(a, b);
  }

  /**
   * Creates the linear problem.
   *
   * <p>The Hessian and gradient parameter from the current best scoring parameter set are assumed
   * to be in alpha and beta. These are copied into the working variables covar and da. The lambda
   * parameter is used to weight the diagonal of the Hessian.
   *
   * @param np the number of fit parameters
   */
  protected void createLinearProblem(final int np) {
    final double weight = (1 + lambda);
    for (int i = np; i-- > 0;) {
      da[i] = beta[i];
      for (int j = np; j-- > 0;) {
        covar[i][j] = alpha[i][j];
      }
      covar[i][i] *= weight;
    }
  }

  /**
   * Update the fit parameters. Note that not all parameters are fit and therefore the gradients
   * indices are used to map the fit parameters to the parameters array.
   *
   * <p>This method can be overridden to provide bounded update to the parameters.
   *
   * @param a the current fit parameters
   * @param gradientIndices the gradient indices (maps the fit parameter index to the parameter
   *        array)
   * @param np the number of fit parameters
   * @param da the parameter shift
   * @param ap the new fit parameters
   */
  protected void updateFitParameters(double[] a, int[] gradientIndices, int np, double[] da,
      double[] ap) {
    for (int j = np; j-- > 0;) {
      ap[gradientIndices[j]] = a[gradientIndices[j]] + da[j];
    }
  }

  private FitStatus doFit(int n, double[] y, double[] fx, double[] a, double[] parametersVariance,
      StoppingCriteria sc) {
    sc.initialise(a);
    if (!nonLinearModel(n, y, a, true)) {
      return (calculator.isNaNGradients()) ? FitStatus.INVALID_GRADIENTS
          : FitStatus.SINGULAR_NON_LINEAR_MODEL;
    }
    sc.evaluate(sumOfSquaresWorking[SUM_OF_SQUARES_OLD], sumOfSquaresWorking[SUM_OF_SQUARES_NEW],
        a);

    while (sc.areNotSatisfied()) {
      if (!nonLinearModel(n, y, a, false)) {
        return (calculator.isNaNGradients()) ? FitStatus.INVALID_GRADIENTS
            : FitStatus.SINGULAR_NON_LINEAR_MODEL;
      }

      sc.evaluate(sumOfSquaresWorking[SUM_OF_SQUARES_OLD], sumOfSquaresWorking[SUM_OF_SQUARES_NEW],
          a);
    }

    if (!sc.areAchieved()) {
      if (sc.getIteration() >= sc.getMaximumIterations()) {
        return FitStatus.TOO_MANY_ITERATIONS;
      }
      return FitStatus.FAILED_TO_CONVERGE;
    }

    if (parametersVariance != null && !computeFitDeviations(n, parametersVariance)) {
      return FitStatus.SINGULAR_NON_LINEAR_SOLUTION;
    }

    value = sumOfSquaresWorking[SUM_OF_SQUARES_BEST];

    computeFitValues(n, fx);

    return FitStatus.OK;
  }

  /**
   * Compute the parameter deviations using the covariance matrix of the solution.
   *
   * @param n the n
   * @param parametersVariance the parameter deviations
   * @return true, if successful
   */
  private boolean computeFitDeviations(int n, double[] parametersVariance) {
    if (isMle()) {
      // The Hessian matrix refers to the log-likelihood ratio.
      // Compute and invert a matrix related to the Poisson log-likelihood.
      // This assumes this does achieve the maximum likelihood estimate for a
      // Poisson process.
      final double[][] fim = calculator.fisherInformationMatrix(n, null, func);
      if (calculator.isNaNGradients()) {
        throw new FunctionSolverException(FitStatus.INVALID_GRADIENTS);
      }

      // Use a dedicated solver optimised for inverting the matrix diagonal.
      final FisherInformationMatrix m = new FisherInformationMatrix(fim);
      setDeviations(parametersVariance, m);
      return true;
    }
    final double[] var = calculator.variance(n, null, func);
    if (var != null) {
      setDeviations(parametersVariance, var);
      return true;
    }
    return false;
  }

  private void computeFitValues(int n, double[] fx) {
    // Compute fitted data points
    if (fx != null) {
      for (int i = 0; i < n; i++) {
        fx[i] = func.eval(i);
      }
    }
  }

  /**
   * Uses Levenberg-Marquardt method to fit a nonlinear model with coefficients (a) for a set of
   * data points (x, y).
   *
   * <p>It is assumed that the data points x[i] corresponding to y[i] are consecutive integers from
   * zero.
   *
   * @param y Set of n data points to fit (input)
   * @param fx Fitted data points (output)
   * @param a Set of m coefficients (input/output)
   * @param parametersVariance Standard deviation of the set of m coefficients (output)
   * @return The fit status
   */
  @Override
  public FitStatus computeFit(double[] y, double[] fx, final double[] a,
      final double[] parametersVariance) {
    final int n = y.length;
    final int nparams = function.gradientIndices().length;

    // Create dynamically for the parameter sizes
    calculator = GradientCalculatorUtils.newCalculator(nparams, isMle());

    // Initialise storage.
    // Note that covar and da are passed to EJMLLinerSolver and so must be the correct size.
    beta = new double[nparams];
    da = new double[nparams];
    covar = new double[nparams][nparams];
    alpha = new double[nparams][nparams];
    ap = new double[a.length];

    // Store the { best, previous, new } sum-of-squares values
    sumOfSquaresWorking = new double[3];

    boolean copyYfit = true;
    if (isMle()) {
      // We must have positive data
      y = ensurePositive(n, y);

      // Store the function values for use in computing the log likelihood
      lastY = y;
      if (fx == null) {
        // Re-use space
        lastFx = SimpleArrayUtils.ensureSize(lastFx, y.length);
        fx = lastFx;
        // We will not need to copy fx later since lastFx is used direct
        copyYfit = false;
      }
    }

    final FitStatus result = doFit(n, y, fx, a, parametersVariance, sc);
    this.evaluations = this.iterations = sc.getIteration();

    // Ensure we have a private copy of fx since the any calling code may modify it
    if (isMle() && copyYfit) {
      lastFx = SimpleArrayUtils.ensureSize(lastFx, y.length);
      System.arraycopy(fx, 0, lastFx, 0, y.length);
    }

    return result;
  }

  /**
   * Sets the the initial lambda for the Levenberg-Marquardt fitting routine.
   *
   * @param initialLambda the initial lambda for the Levenberg-Marquardt fitting routine
   */
  public void setInitialLambda(double initialLambda) {
    this.initialLambda = initialLambda;
  }

  /**
   * Gets the the initial lambda for the Levenberg-Marquardt fitting routine.
   *
   * @return the initial lambda
   */
  public double getInitialLambda() {
    return initialLambda;
  }

  /**
   * Gets the initial residual sum of squares.
   *
   * @return the initialResidualSumOfSquares.
   */
  public double getInitialResidualSumOfSquares() {
    return initialResidualSumOfSquares;
  }

  @Override
  public void setGradientFunction(GradientFunction function) {
    super.setGradientFunction(function);
    if (!(function instanceof NonLinearFunction)) {
      throw new IllegalArgumentException("Function must be a NonLinearFunction");
    }
    func = (NonLinearFunction) function;
  }

  /**
   * Set the stopping criteria for the {@link #fit(double[], double[], double[], double[])} method.
   *
   * @param sc the new stopping criteria
   */
  public void setStoppingCriteria(StoppingCriteria sc) {
    if (sc == null) {
      sc = new ErrorStoppingCriteria();
    }
    this.sc = sc;
  }

  /**
   * Checks if set to perform Maximum Likelihood Estimation assuming Poisson model.
   *
   * @return true if is set to perform MLE
   */
  public boolean isMle() {
    return getType() == FunctionSolverType.MLE;
  }

  /**
   * Sets to true to perform Maximum Likelihood Estimation assuming Poisson model.
   *
   * <p>This modifies the standard LVM as described in Laurence &amp; Chromy (2010) Efficient
   * maximum likelihood estimator. Nature Methods 7, 338-339. The input data must be Poisson
   * distributed for this to be relevant.
   *
   * @param mle true to perform Maximum Likelihood Estimation
   */
  public void setMle(boolean mle) {
    if (mle) {
      setType(FunctionSolverType.MLE);
    } else {
      setType((func.canComputeWeights()) ? FunctionSolverType.WLSE : FunctionSolverType.LSE);
    }
  }

  @Override
  public boolean computeValue(double[] y, double[] fx, double[] a) {
    final int n = y.length;

    // Create dynamically for the parameter sizes
    calculator = GradientCalculatorUtils.newCalculator(function.getNumberOfGradients(), isMle());

    if (isMle()) {
      // We must have positive data
      y = ensurePositive(n, y);

      // Store the function values for use in computing the log likelihood
      lastY = y;
      if (fx == null) {
        // Re-use space
        if (lastFx == null || lastFx.length < y.length) {
          lastFx = new double[y.length];
        }
        fx = lastFx;
      } else {
        lastFx = fx;
      }
    }

    value = calculator.findLinearised(n, y, fx, a, func);

    return true;
  }

  @Override
  public boolean computeDeviations(double[] y, double[] a, double[] parametersVariance) {
    calculator = GradientCalculatorUtils.newCalculator(function.getNumberOfGradients(), isMle());

    if (isMle()) {
      return super.computeDeviations(y, a, parametersVariance);
    }

    // LSE computation
    final double[] var = calculator.variance(y.length, a, func);
    if (var != null) {
      setDeviations(parametersVariance, var);
      return true;
    }
    return false;
  }

  @Override
  protected FisherInformationMatrix computeFisherInformationMatrix(double[] y, double[] a) {
    // Compute and invert a matrix related to the Poisson log-likelihood.
    // This assumes this does achieve the maximum likelihood estimate for a
    // Poisson process.
    final double[][] fim = calculator.fisherInformationMatrix(y.length, a, func);
    if (calculator.isNaNGradients()) {
      throw new FunctionSolverException(FitStatus.INVALID_GRADIENTS);
    }
    return new FisherInformationMatrix(fim);
  }

  @Override
  public double getTotalSumOfSquares() {
    if (getType() == FunctionSolverType.LSE) {
      return super.getTotalSumOfSquares();
    }
    throw new IllegalStateException();
  }

  @Override
  public double getChiSquared() {
    if (getType() == FunctionSolverType.WLSE) {
      // The weighted LSE will produce the chi-squared
      return value;
    }
    throw new IllegalStateException();
  }

  @Override
  public double getLogLikelihood() {
    if (getType() == FunctionSolverType.MLE && lastY != null) {
      // The MLE version directly computes the log-likelihood ratio.
      // We must compute the log likelihood for a Poisson MLE.
      if (Double.isNaN(ll)) {
        ll = PoissonCalculator.fastLogLikelihood(lastFx, lastY);
      }
      return ll;
    }
    throw new IllegalStateException();
  }

  @Override
  public double getLogLikelihoodRatio() {
    if (getType() == FunctionSolverType.MLE) {
      // The MLE version directly computes the log-likelihood ratio
      return value;
    }
    throw new IllegalStateException();
  }

  @Override
  public double getQ() {
    if (getType() == FunctionSolverType.MLE) {
      // Value will be the log-likelihood ratio for the MLE.
      // Wilks theorum states the LLR approaches the chi-squared distribution for large n.
      return ChiSquaredDistributionTable.computeQValue(value,
          getNumberOfFittedPoints() - getNumberOfFittedParameters());
    }
    if (getType() == FunctionSolverType.WLSE) {
      // Value will be the Chi-squared
      return ChiSquaredDistributionTable.computeQValue(value,
          getNumberOfFittedPoints() - getNumberOfFittedParameters());
    }
    throw new IllegalStateException();
  }
}
