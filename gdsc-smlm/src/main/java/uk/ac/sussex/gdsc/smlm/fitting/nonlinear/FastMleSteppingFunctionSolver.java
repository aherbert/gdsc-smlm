/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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

import uk.ac.sussex.gdsc.core.utils.SortUtils;
import uk.ac.sussex.gdsc.smlm.data.NamedObject;
import uk.ac.sussex.gdsc.smlm.fitting.FisherInformationMatrix;
import uk.ac.sussex.gdsc.smlm.fitting.FitStatus;
import uk.ac.sussex.gdsc.smlm.fitting.FunctionSolverType;
import uk.ac.sussex.gdsc.smlm.fitting.MleFunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.FastMleGradient2Procedure;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.FastMleGradient2ProcedureUtils;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.PoissonGradientProcedure;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.PoissonGradientProcedureUtils;
import uk.ac.sussex.gdsc.smlm.function.ChiSquaredDistributionTable;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;
import uk.ac.sussex.gdsc.smlm.function.Gradient2Function;
import uk.ac.sussex.gdsc.smlm.function.Gradient2FunctionValueStore;
import uk.ac.sussex.gdsc.smlm.function.OffsetGradient1Function;
import uk.ac.sussex.gdsc.smlm.function.OffsetGradient2Function;

/**
 * Uses the Fast MLE method to fit a gradient function with coefficients (a).
 *
 * <p>Calculates the Newton-Raphson update vector for a Poisson process using the first and second
 * partial derivatives.
 *
 * <p>Ref: Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically
 * minimum uncertainty. Nature Methods 7, 373-375 (supplementary note), Eq. 12.
 *
 * <p>Ref: Huang et al, (2015). Video-rate nanoscopy using sCMOS camera–specific single-molecule
 * localization algorithms. Nature Methods 10, 653–658.
 */
public class FastMleSteppingFunctionSolver extends SteppingFunctionSolver
    implements MleFunctionSolver {
  /**
   * Define the method to use when the line search direction is not in the same direction as that
   * defined by the first derivative gradient.
   */
  public enum LineSearchMethod implements NamedObject {
    /**
     * Do nothing to handle the incorrect orientation. The default solver action is taken. This may
     * cause the search to take an invalid move or it may error.
     */
    NONE("None"),

    /**
     * Ignore any search direction that is in the opposite direction to the first derivative
     * gradient.
     */
    IGNORE("Ignore"),

    /**
     * Progressively ignore any search direction that is in the opposite direction to the first
     * derivative gradient. Do this in order of the magnitude of the error
     */
    PARTIAL_IGNORE("Partial ignore");

    private final String name;

    LineSearchMethod(String name) {
      this.name = name;
    }

    @Override
    public String toString() {
      return name;
    }

    @Override
    public String getName() {
      return name;
    }
  }

  /** The line search method. */
  protected LineSearchMethod lineSearchMethod = LineSearchMethod.NONE;

  /** The log-likelihood. */
  protected double ll = Double.NaN;
  /** Flag if the log-likelihood is the pseudo log-likelihood. */
  protected boolean isPseudoLogLikelihood;
  /** The log-likelihood ratio. */
  protected double llr = Double.NaN;

  /**
   * The per observation variances. This is not null if fitting using the method of Huang, et al
   * (2015).
   */
  protected double[] obsVariances;
  /**
   * The gradient function used by the procedure. This may be wrapped to add the per observation
   * variances if fitting using the method of Huang, et al (2015).
   */
  protected Gradient2Function f2;
  /** The gradient procedure. */
  protected FastMleGradient2Procedure gradientProcedure;

  /** The old parameters (a). */
  protected double[] oldA;

  /** The search direction. */
  protected double[] searchDirection;

  /** The first evaluation flag. */
  protected boolean firstEvaluation;

  /**
   * Create a new stepping function solver.
   *
   * @param function the function
   * @param maxRelativeError the max relative error
   * @param maxAbsoluteError the max absolute error
   * @throws NullPointerException if the function is null
   */
  public FastMleSteppingFunctionSolver(Gradient2Function function, double maxRelativeError,
      double maxAbsoluteError) {
    this(function, new ToleranceChecker(maxRelativeError, maxAbsoluteError), null);
  }

  /**
   * Create a new stepping function solver.
   *
   * @param function the function
   * @param tc the tolerance checker
   * @param bounds the bounds
   * @throws NullPointerException if the function or tolerance checker is null
   */
  public FastMleSteppingFunctionSolver(Gradient2Function function, ToleranceChecker tc,
      ParameterBounds bounds) {
    super(FunctionSolverType.MLE, function, tc, bounds);
  }

  @Override
  protected void preProcess() {
    ll = llr = Double.NaN;
    isPseudoLogLikelihood = false;
  }

  @Override
  protected double[] prepareFitValue(double[] y, double[] a) {
    firstEvaluation = true;
    // Ensure the gradient procedure is created
    y = prepareY(y);
    gradientProcedure = createGradientProcedure(y);
    // Ensure maximisation
    tc.setMinimiseValue(false);
    return y;
  }

  /**
   * Prepare Y for the gradient procedure by ensuring positive values.
   *
   * @param y the y
   * @return the new y
   */
  protected double[] prepareY(double[] y) {
    // We can handle per-observation variances as detailed in
    // Huang, et al. (2015) by simply adding the variances to the target data.

    final int n = y.length;
    obsVariances = getWeights(n);
    if (obsVariances != null) {
      final double[] x = new double[n];
      for (int i = 0; i < n; i++) {
        // Also ensure the input y is positive
        x[i] = (y[i] > 0) ? y[i] + obsVariances[i] : obsVariances[i];
      }
      return x;
    }
    return ensurePositive(y);
  }

  /**
   * Creates the gradient procedure.
   *
   * @param y the y
   * @return the newton raphson gradient 2 procedure
   */
  protected FastMleGradient2Procedure createGradientProcedure(double[] y) {
    // We can handle per-observation variances as detailed in
    // Huang, et al. (2015) by simply adding the variances to the computed value.
    f2 = (Gradient2Function) function;
    if (obsVariances != null) {
      f2 = OffsetGradient2Function.wrapGradient2Function(f2, obsVariances);
    }
    return FastMleGradient2ProcedureUtils.create(y, f2);
  }

  @Override
  protected double computeFitValue(double[] a) {
    if (lineSearchMethod != LineSearchMethod.NONE) {
      // The code below will adjust the search direction
      if (firstEvaluation) {
        // For the first evaluation we store the old value and initialise
        firstEvaluation = false;
        oldA = a.clone();
        searchDirection = new double[a.length];
      } else {
        // All subsequent calls to computeFitValue() must check the search direction
        for (int i = 0; i < searchDirection.length; i++) {
          // Configure the search direction with the full Newton step
          searchDirection[i] = a[i] - oldA[i];
        }

        final double[] gradient = gradientProcedure.d1;
        final int[] gradientIndices = function.gradientIndices();

        double slope = 0.0;
        for (int i = 0; i < gradient.length; i++) {
          slope += gradient[i] * searchDirection[gradientIndices[i]];
        }

        if (slope <= 0) {
          switch (lineSearchMethod) {
            case IGNORE:
              // Ignore any search direction that is in the opposite direction to the
              // first derivative gradient.
              slope = 0.0;
              for (int i = 0; i < gradient.length; i++) {
                final double slopeComponent = gradient[i] * searchDirection[gradientIndices[i]];
                if (slopeComponent < 0) {
                  // Ignore this component
                  a[gradientIndices[i]] = oldA[gradientIndices[i]];
                } else {
                  slope += slopeComponent;
                }
              }
              if (slope == 0) {
                // No move so just set converged
                tc.setConverged();
                return ll;
                // throw new FunctionSolverException(FitStatus.LINE_SEARCH_ERROR, "No slope");
              }
              break;

            case PARTIAL_IGNORE:
              // Progressively ignore any search direction that is in the opposite direction to
              // the first derivative gradient. Do this in order of the magnitude of the error
              final double[] slopeComponents = new double[gradient.length];
              final int[] indices = new int[slopeComponents.length];
              for (int i = 0; i < slopeComponents.length; i++) {
                slopeComponents[i] = gradient[i] * searchDirection[gradientIndices[i]];
                indices[i] = i;
              }
              SortUtils.sortIndices(indices, slopeComponents, false);
              int count = 0;
              while (slope <= 0 && count < slopeComponents.length
                  && slopeComponents[indices[count]] <= 0) {
                final int i = indices[count];
                // Ignore this component
                slope -= slopeComponents[i];
                a[gradientIndices[i]] = oldA[gradientIndices[i]];
                count++;
              }
              if (count == slopeComponents.length) {
                // No move so just set converged
                tc.setConverged();
                return ll;
                // throw new FunctionSolverException(FitStatus.LINE_SEARCH_ERROR, "No slope");
              }
              break;

            default:
              throw new IllegalStateException("Unknown line search method: " + lineSearchMethod);
          }
        }
      }
    }

    computeGradients(a);

    // Log-likelihood only needs to be computed if the tolerance checker
    // is testing the value. Use the Pseudo log-likelihood for speed.
    if (tc.checkValue) {
      ll = gradientProcedure.computePseudoLogLikelihood();
      isPseudoLogLikelihood = true;
    }

    return ll;
  }

  /**
   * Compute the gradients for the Newton step using the gradient procedure.
   *
   * @param a the funtion parameters
   */
  protected void computeGradients(double[] a) {
    gradientProcedure.computeSecondDerivative(a);
    if (gradientProcedure.isNaNGradients()) {
      throw new FunctionSolverException(FitStatus.INVALID_GRADIENTS);
    }
  }

  @Override
  protected void computeStep(double[] step) {
    final double[] d1 = gradientProcedure.d1;
    final double[] d2 = gradientProcedure.d2;

    // Simple Newton-Raphson update step as per Smith et al, (2010), SI Eq. 13:
    // parameter -> new parameter + delta
    // => new parameter = parameter - delta
    for (int i = 0; i < step.length; i++) {
      step[i] = -d1[i] / d2[i];
    }
  }

  @Override
  protected boolean accept(double currentValue, double[] a, double newValue, double[] newA) {
    // Always accept the step. The Smith, et al (2010) paper used 10 steps until
    // convergence, with no apparent checking of the log-likelihood value or parameters.
    // The Newton-Raphson method converges fast but does require a good initial
    // estimate for the parameters.
    return true;
  }

  @Override
  protected boolean computeValue(double[] y, double[] fx, double[] a) {
    // This is over-ridden since the fx values are computed
    // and stored by the gradient procedure. The super-class SteppingFunctionSolver
    // wraps the function with a Gradient1FunctionStore to store the
    // fx. This is not a gradient 2 function so causes a run-time error
    // in createGradientProcedure(double[])

    gradientIndices = function.gradientIndices();
    lastY = prepareFunctionValue(y, a);
    value = computeFunctionValue(a);
    if (fx != null) {
      copyFunctionValue(fx);
    }
    return true;
  }

  @Override
  protected double[] prepareFunctionValue(double[] y, double[] a) {
    y = prepareY(y);
    gradientProcedure = createGradientProcedure(y);
    return y;
  }

  @Override
  protected double computeFunctionValue(double[] a) {
    ll = gradientProcedure.computeLogLikelihood(a);
    isPseudoLogLikelihood = false;
    return ll;
  }

  /**
   * Copy the function value into the fx array.
   *
   * @param fx the function values
   */
  private void copyFunctionValue(double[] fx) {
    final double[] u = gradientProcedure.u;
    if (obsVariances != null) {
      // The function was wrapped to add the per-observation variances
      // to the computed value, these must be subtracted to get the actual value
      final int n = u.length;
      for (int i = 0; i < n; i++) {
        fx[i] = u[i] - obsVariances[i];
      }
    } else {
      System.arraycopy(u, 0, fx, 0, u.length);
    }
  }

  @Override
  protected void computeValues(double[] fx) {
    copyFunctionValue(fx);
  }

  @Override
  protected FisherInformationMatrix computeLastFisherInformationMatrix(double[] fx) {
    Gradient2Function f2 = (Gradient2Function) function;
    // Capture the y-values if necessary
    if (fx != null && fx.length == f2.size()) {
      f2 = new Gradient2FunctionValueStore(f2, fx);
    }
    // Add the weights if necessary
    if (obsVariances != null) {
      f2 = OffsetGradient2Function.wrapGradient2Function(f2, obsVariances);
    }
    // The fisher information is that for a Poisson process
    final PoissonGradientProcedure p = PoissonGradientProcedureUtils.create(f2);
    initialiseAndRun(p);
    if (p.isNaNGradients()) {
      throw new FunctionSolverException(FitStatus.INVALID_GRADIENTS);
    }
    return new FisherInformationMatrix(p.getLinear(), p.numberOfGradients);
  }

  /**
   * Initialise and run the procedure using the fitted parameters. Provided to allow the
   * a sub-class to initialise the function for gradients.
   *
   * @param procedure the procedure
   */
  protected void initialiseAndRun(PoissonGradientProcedure procedure) {
    procedure.computeFisherInformation(null); // Assume preinitialised function
  }

  @Override
  protected double[] prepareFunctionFisherInformationMatrix(double[] y, double[] a) {
    // Do not create the gradient procedure as we use a different one for the fisher information.
    return prepareY(y);
  }

  @Override
  protected FisherInformationMatrix computeFunctionFisherInformationMatrix(double[] y, double[] a) {
    // The fisher information is that for a Poisson process.
    // We must wrap the gradient function if weights are present.
    Gradient1Function f1 = (Gradient1Function) function;
    if (obsVariances != null) {
      f1 = OffsetGradient1Function.wrapGradient1Function(f1, obsVariances);
    }
    final PoissonGradientProcedure p = PoissonGradientProcedureUtils.create(f1);
    p.computeFisherInformation(a);
    if (p.isNaNGradients()) {
      throw new FunctionSolverException(FitStatus.INVALID_GRADIENTS);
    }
    return new FisherInformationMatrix(p.getLinear(), p.numberOfGradients);
  }

  @Override
  public double getValue() {
    // Override this to return the log likelihood since the value may not
    // actually be computed during computeFitValue(double[])
    return getLogLikelihood();
  }

  @Override
  public double getLogLikelihood() {
    if (Double.isNaN(ll)) {
      ll = gradientProcedure.computeLogLikelihood();
      isPseudoLogLikelihood = false;
    } else if (isPseudoLogLikelihood) {
      isPseudoLogLikelihood = false;
      ll -= gradientProcedure.computeLogXFactorialTerm();
    }
    return ll;
  }

  @Override
  public double getLogLikelihoodRatio() {
    if (Double.isNaN(llr)) {
      llr = gradientProcedure.computeLogLikelihoodRatio(getLogLikelihood());
    }
    return llr;
  }

  @Override
  public double getQ() {
    // Wilks theorum states the LLR approaches the chi-squared distribution for large n.
    return ChiSquaredDistributionTable.computeQValue(getLogLikelihoodRatio(),
        getNumberOfFittedPoints() - getNumberOfFittedParameters());
  }

  /**
   * Gets the line search method.
   *
   * @return the line search method
   */
  public LineSearchMethod getLineSearchMethod() {
    return lineSearchMethod;
  }

  /**
   * Sets the line search method.
   *
   * @param lineSearchMethod the new line search method
   */
  public void setLineSearchMethod(LineSearchMethod lineSearchMethod) {
    this.lineSearchMethod = lineSearchMethod;
  }
}
