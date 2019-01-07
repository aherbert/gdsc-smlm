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

import uk.ac.sussex.gdsc.core.utils.BitFlagUtils;
import uk.ac.sussex.gdsc.smlm.fitting.FisherInformationMatrix;
import uk.ac.sussex.gdsc.smlm.fitting.FitStatus;
import uk.ac.sussex.gdsc.smlm.fitting.FunctionSolverType;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;
import uk.ac.sussex.gdsc.smlm.function.Gradient1FunctionStore;
import uk.ac.sussex.gdsc.smlm.function.GradientFunction;
import uk.ac.sussex.gdsc.smlm.function.ValueFunction;
import uk.ac.sussex.gdsc.smlm.function.ValueProcedure;

/**
 * Abstract class for FunctionSolvers that use update steps to the current parameters.
 */
public abstract class SteppingFunctionSolver extends BaseFunctionSolver {
  /**
   * Simple class to allow the values to be computed.
   */
  private static class SimpleValueProcedure implements ValueProcedure {
    int i;
    double[] fx;

    SimpleValueProcedure(double[] fx) {
      this.fx = fx;
    }

    @Override
    public void execute(double value) {
      fx[i++] = value;
    }
  }

  /** The gradient indices. */
  protected int[] gradientIndices;
  /** The tolerance checker. */
  protected final ToleranceChecker tc;
  /** The bounds. */
  protected final ParameterBounds bounds;
  private double[] weights;

  /**
   * Create a new stepping function solver.
   *
   * @param type the type
   * @param function the function
   * @throws NullPointerException if the function is null
   */
  public SteppingFunctionSolver(FunctionSolverType type, Gradient1Function function) {
    this(type, function, new ToleranceChecker(1e-3, 1e-6), null);
  }

  /**
   * Create a new stepping function solver.
   *
   * @param type the type
   * @param function the function
   * @param tc the tolerance checker
   * @param bounds the bounds
   * @throws NullPointerException if the function or tolerance checker is null
   * @throws IllegalArgumentException if the bounds are not constructed with the same gradient
   *         function
   */
  public SteppingFunctionSolver(FunctionSolverType type, Gradient1Function function,
      ToleranceChecker tc, ParameterBounds bounds) {
    super(type, function);
    if (tc == null) {
      throw new NullPointerException("Null tolerance checker");
    }
    this.tc = tc;
    if (bounds == null) {
      bounds = new ParameterBounds(function);
    } else if (bounds.getGradientFunction() != function) {
      throw new IllegalArgumentException(
          "Bounds must be constructed with the same gradient function");
    }
    this.bounds = bounds;
  }

  /**
   * Compute fit.
   *
   * @param y the y
   * @param fx the fx
   * @param a the a
   * @param parametersVariance the parametersVariance
   * @return the fit status
   */
  @Override
  protected FitStatus computeFit(double[] y, double[] fx, double[] a, double[] parametersVariance) {
    // Lay out a simple iteration loop for a stepping solver.
    // The sub-class must compute the next step.
    // This class handles attenuation of the step.
    // The sub-class determines if the step is accepted or rejected.

    gradientIndices = function.gradientIndices();
    final double[] step = new double[gradientIndices.length];
    final double[] newA = a.clone();

    // Initialise for fitting
    bounds.initialise();
    tc.reset();
    final String name = this.getClass().getSimpleName();

    try {
      lastY = prepareFitValue(y, a);

      // First evaluation
      double currentValue = computeFitValue(a);
      log("%s Value [%s] = %s : %s\n", name, tc.getIterations(), currentValue, a);

      int status = 0;
      while (true) {
        // Compute next step
        computeStep(step);
        log("%s Step [%s] = %s\n", name, tc.getIterations(), step);

        // Apply bounds to the step
        bounds.applyBounds(a, step, newA);

        // Evaluate
        final double newValue = computeFitValue(newA);
        log("%s Value [%s] = %s : %s\n", name, tc.getIterations(), newValue, newA);

        // Check stopping criteria
        status = tc.converged(currentValue, a, newValue, newA);
        log("%s Status [%s] = %s\n", name, tc.getIterations(), status);
        if (status != 0) {
          value = newValue;
          System.arraycopy(newA, 0, a, 0, a.length);
          break;
        }

        // Check if the step was an improvement
        if (accept(currentValue, a, newValue, newA)) {
          log("%s Accepted [%s]\n", name, tc.getIterations());
          currentValue = newValue;
          System.arraycopy(newA, 0, a, 0, a.length);
          bounds.accepted(a, newA);
        }
      }

      log("%s End [%s] = %s\n", name, tc.getIterations(), status);

      if (BitFlagUtils.anySet(status, ToleranceChecker.STATUS_CONVERGED)) {
        log("%s Converged [%s]\n", name, tc.getIterations());
        // A solver may compute both at the same time...
        if (parametersVariance != null) {
          computeDeviationsAndValues(parametersVariance, fx);
        } else if (fx != null) {
          computeValues(fx);
        }
        return FitStatus.OK;
      }

      // Check the iterations
      if (BitFlagUtils.areSet(status, ToleranceChecker.STATUS_MAX_ITERATIONS)) {
        return FitStatus.TOO_MANY_ITERATIONS;
      }

      // We should not reach here unless we missed something
      return FitStatus.FAILED_TO_CONVERGE;
    } catch (final FunctionSolverException ex) {
      // XXX - debugging
      final String msg = ex.getMessage();
      if (msg != null) {
        System.out.printf("%s failed: %s - %s\n", getClass().getSimpleName(),
            ex.fitStatus.getName(), msg);
      } else {
        System.out.printf("%s failed: %s\n", getClass().getSimpleName(), ex.fitStatus.getName());
      }
      return ex.fitStatus;
    } finally {
      iterations = tc.getIterations();
      // Allow subclasses to increment this
      if (evaluations == 0) {
        evaluations = iterations;
      }
    }
  }

  /**
   * Log progress from the solver.
   *
   * @param format the format
   * @param args the args
   */
  private void log(String format, Object... args) {
    // // Convert arrays to a single string
    // for (int i=0; i<args.length; i++)
    // if (args[i] instanceof double[])
    // args[i] = java.util.Arrays.toString((double[])args[i]);
    // System.out.printf(format, args);
  }

  /**
   * Prepare y for fitting, e.g. ensure strictly positive values.
   *
   * @param y the y
   * @param a the parameters
   * @return the new y
   */
  protected abstract double[] prepareFitValue(double[] y, double[] a);

  /**
   * Compute the fit value using the parameters. The y data is the same as that passed to
   * {@link #prepareFitValue(double[], double[])}.
   *
   * <p>This method is followed by a call to {@link #computeStep(double[])} so the step could be
   * pre-computed here.
   *
   * @param a the parameters
   * @return the fit value
   */
  protected abstract double computeFitValue(double[] a);

  /**
   * Compute the update step for the current parameters.
   *
   * @param step the step
   */
  protected abstract void computeStep(double[] step);

  /**
   * Determine if the step should be accepted. If accepted then the current parameters and function
   * value are updated and any bounds on the step size may be updated.
   *
   * <p>Note that although this class handles convergence on the value/parameters it is left to the
   * sub-class to determine if each step should be accepted.
   *
   * @param currentValue the current value
   * @param a the current parameters
   * @param newValue the new value
   * @param newA the new parameters
   * @return true, if successful
   */
  protected abstract boolean accept(double currentValue, double[] a, double newValue,
      double[] newA);

  /**
   * Compute the deviations for the parameters a from the last call to
   * {@link #computeFitValue(double[])}. Optionally store the function values.
   *
   * @param parametersVariance the parameter deviations
   * @param fx the y fit (may be null)
   */
  protected void computeDeviationsAndValues(double[] parametersVariance, double[] fx) {
    // Use a dedicated solver optimised for inverting the matrix diagonal.
    // The last Hessian matrix should be stored in the working alpha.
    final FisherInformationMatrix m = computeFisherInformationMatrix(fx);

    setDeviations(parametersVariance, m);
  }

  /**
   * Compute the Fisher Information matrix for the parameters a from the last call to
   * {@link #computeFitValue(double[])}. This can be used to set the covariances for each of the
   * fitted parameters.
   *
   * <p>Alternatively a sub-class can override
   * {@link #computeDeviationsAndValues(double[], double[])} directly and provide a dummy
   * implementation of this function as it will not be used, e.g. throw an exception.
   *
   * @param fx the y fit (may be null)
   * @return the Fisher Information matrix
   */
  protected abstract FisherInformationMatrix computeFisherInformationMatrix(double[] fx);

  /**
   * Compute the function y-values using the y and parameters a from the last call to
   * {@link #computeFitValue(double[])}.
   *
   * <p>Utility method to compute the function values using the preinitialised function. Sub-classes
   * may override this if they have cached the function values from the last execution of a forEach
   * procedure.
   *
   * <p>The base gradient function is used. If sub-classes wrap the function (e.g. with
   * per-observation weights) then these will be omitted.
   *
   * @param fx the y fit values
   */
  protected void computeValues(double[] fx) {
    final ValueFunction function = (ValueFunction) this.function;
    function.forEach(new SimpleValueProcedure(fx));
  }

  /** {@inheritDoc} */
  @Override
  protected boolean computeValue(double[] y, double[] fx, double[] a) {
    // If the fx array is not null then wrap the gradient function.
    // Compute the value and the wrapper will store the values appropriately.
    // Then reset the gradient function.

    // Note: If a sub class wraps the function with weights
    // then the weights will not be stored in the function value.
    // Only the value produced by the original function is stored:
    // Wrapped (+weights) < FunctionStore < Function

    // However if the base function is already wrapped then this will occur:
    // Wrapped (+weights) < FunctionStore < Wrapped (+precomputed) < Function

    gradientIndices = function.gradientIndices();
    if (fx != null && fx.length == ((Gradient1Function) function).size()) {
      final GradientFunction tmp = function;
      function = new Gradient1FunctionStore((Gradient1Function) function, fx, null);
      lastY = prepareFunctionValue(y, a);
      value = computeFunctionValue(a);
      function = tmp;
    } else {
      lastY = prepareFunctionValue(y, a);
      value = computeFunctionValue(a);
    }
    return true;
  }

  /**
   * Prepare y for computing the function value, e.g. ensure strictly positive values.
   *
   * @param y the y
   * @param a the parameters
   * @return the new y
   */
  protected abstract double[] prepareFunctionValue(double[] y, double[] a);

  /**
   * Compute the function value. The y data is the same as that passed to
   * {@link #prepareFunctionValue(double[], double[])}
   *
   * @param a the parameters
   * @return the function value
   */
  protected abstract double computeFunctionValue(double[] a);

  /** {@inheritDoc} */
  @Override
  protected FisherInformationMatrix computeFisherInformationMatrix(double[] y, double[] a) {
    gradientIndices = function.gradientIndices();
    y = prepareFunctionFisherInformationMatrix(y, a);
    return computeFunctionFisherInformationMatrix(y, a);
  }

  /**
   * Prepare y for computing the Fisher information matrix, e.g. ensure strictly positive values.
   *
   * @param y the y
   * @param a the parameters
   * @return the new y
   */
  protected abstract double[] prepareFunctionFisherInformationMatrix(double[] y, double[] a);

  /**
   * Compute the Fisher information matrix.
   *
   * @param y the y
   * @param a the parameters
   * @return the Fisher Information matrix
   */
  protected abstract FisherInformationMatrix computeFunctionFisherInformationMatrix(double[] y,
      double[] a);

  /** {@inheritDoc} */
  @Override
  public boolean isBounded() {
    // Bounds are tighter than constraints and we support those
    return true;
  }

  /** {@inheritDoc} */
  @Override
  public void setBounds(double[] lower, double[] upper) {
    bounds.setBounds(lower, upper);
  }

  /**
   * Warning: If the function is changed then the clamp values may require updating. However setting
   * a new function does not set the clamp values to null to allow caching when the clamp values are
   * unchanged, e.g. evaluation of a different function in the same parameter space.
   *
   * <p>Setting a new function removes the current bounds.
   *
   * @param function the new gradient function
   */
  @Override
  public void setGradientFunction(GradientFunction function) {
    super.setGradientFunction(function);
    bounds.setGradientFunction(function);
  }

  /** {@inheritDoc} */
  @Override
  public boolean isWeighted() {
    return true;
  }

  /** {@inheritDoc} */
  @Override
  public void setWeights(double[] weights) {
    this.weights = weights;
  }

  /**
   * Gets the weights for observations of size n, e.g. the per observation variance term.
   *
   * @param n the size
   * @return the weights
   */
  public double[] getWeights(int n) {
    return (weights == null || weights.length != n) ? null : weights;
  }
}
