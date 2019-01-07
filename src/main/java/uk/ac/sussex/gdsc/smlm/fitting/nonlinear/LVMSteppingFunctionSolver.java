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
import uk.ac.sussex.gdsc.smlm.fitting.FitStatus;
import uk.ac.sussex.gdsc.smlm.fitting.FunctionSolverType;
import uk.ac.sussex.gdsc.smlm.fitting.linear.EJMLLinearSolver;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.LVMGradientProcedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;
import uk.ac.sussex.gdsc.smlm.function.GradientFunction;

/**
 * Uses the Levenberg-Marquardt method to fit a gradient function with coefficients (a).
 */
public abstract class LVMSteppingFunctionSolver extends SteppingFunctionSolver {
  /**
   * The solver used for solving A x = b to find the update x. <ul> <li>A = Scaled Hessian matrix
   * (Alpha) <li>b = Gradient vector (beta) <li>x = Update vector to modify the parameters </ul>
   */
  protected EJMLLinearSolver solver = new EJMLLinearSolver();

  /** The gradient procedure. */
  protected LVMGradientProcedure gradientProcedure;

  /** The initial lambda value for the LVM algorithm. */
  protected double initialLambda = 0.01;

  /** The current lambda value for the LVM algorithm. */
  protected double lambda;

  // Alpha = Scaled Hessian matrix
  // beta = Gradient vector
  // We want to solve: A x = b to find the update x

  /** Current best alpha (Scaled Hessian matrix). */
  protected double[] alpha;
  /** Current best beta (Gradient vector). */
  protected double[] beta;

  /** Working alpha. */
  protected double[] walpha;
  /** Working beta. */
  protected double[] wbeta;

  // TODO - Determine what a good solution tolerance would be.
  // We may not need to be that strict to accept the solution.

  /** The default max relative error. */
  public static final double DEFAULT_MAX_RELATIVE_ERROR = 1e-3;
  /** The default max absolute error. */
  public static final double DEFAULT_MAX_ABSOLUTE_ERROR = 1e-4;

  /**
   * Create a new stepping function solver.
   *
   * @param type the type
   * @param function the function
   * @throws NullPointerException if the function is null
   */
  public LVMSteppingFunctionSolver(FunctionSolverType type, Gradient1Function function) {
    this(type, function, DEFAULT_MAX_RELATIVE_ERROR, DEFAULT_MAX_ABSOLUTE_ERROR);
  }

  /**
   * Create a new stepping function solver.
   *
   * @param type the type
   * @param function the function
   * @param maxRelativeError Validate the Levenberg-Marquardt fit solution using the specified
   *        maximum relative error
   * @param maxAbsoluteError Validate the Levenberg-Marquardt fit solution using the specified
   *        maximum absolute error
   * @throws NullPointerException if the function is null
   */
  public LVMSteppingFunctionSolver(FunctionSolverType type, Gradient1Function function,
      double maxRelativeError, double maxAbsoluteError) {
    super(type, function);
    solver.setEqual(new DoubleEquality(maxRelativeError, maxAbsoluteError));
  }

  /**
   * Create a new stepping function solver.
   *
   * @param type the type
   * @param function the function
   * @param tc the tolerance checker
   * @param bounds the bounds
   * @throws NullPointerException if the function or tolerance checker is null
   */
  public LVMSteppingFunctionSolver(FunctionSolverType type, Gradient1Function function,
      ToleranceChecker tc, ParameterBounds bounds) {
    this(type, function, tc, bounds, DEFAULT_MAX_RELATIVE_ERROR, DEFAULT_MAX_ABSOLUTE_ERROR);
  }

  /**
   * Create a new stepping function solver.
   *
   * @param type the type
   * @param function the function
   * @param tc the tolerance checker
   * @param bounds the bounds
   * @param maxRelativeError Validate the Levenberg-Marquardt fit solution using the specified
   *        maximum relative error
   * @param maxAbsoluteError Validate the Levenberg-Marquardt fit solution using the specified
   *        maximum absolute error
   * @throws NullPointerException if the function or tolerance checker is null
   */
  public LVMSteppingFunctionSolver(FunctionSolverType type, Gradient1Function function,
      ToleranceChecker tc, ParameterBounds bounds, double maxRelativeError,
      double maxAbsoluteError) {
    super(type, function, tc, bounds);
    solver.setEqual(new DoubleEquality(maxRelativeError, maxAbsoluteError));
  }

  @Override
  protected double[] prepareFitValue(double[] y, double[] a) {
    // Ensure the gradient procedure is created
    y = prepareY(y);
    gradientProcedure = createGradientProcedure(y);

    // Ensure minimisation
    tc.setMinimiseValue(true);

    // Set up the current best Hessian matrix and gradient parameter
    lambda = initialLambda;
    final int n = gradientProcedure.n;
    alpha = null;
    beta = null;
    walpha = new double[n * n];
    wbeta = new double[n];

    return y;
  }

  /**
   * Prepare Y for the gradient procedure, e.g. ensure positive values.
   *
   * @param y the y
   * @return the new y
   */
  protected double[] prepareY(double[] y) {
    return y;
  }

  /**
   * Creates the gradient procedure.
   *
   * @param y the y
   * @return the LVM gradient procedure
   */
  protected abstract LVMGradientProcedure createGradientProcedure(double[] y);

  @Override
  protected double computeFitValue(double[] a) {
    gradientProcedure.gradient(a);

    if (gradientProcedure.isNaNGradients()) {
      throw new FunctionSolverException(FitStatus.INVALID_GRADIENTS);
    }

    if (alpha == null) {
      // This is the first computation:
      // Set the current alpha and beta
      alpha = gradientProcedure.getAlphaLinear();
      beta = gradientProcedure.beta.clone();
    } else {
      // This is a subsequent computation:
      // Store the working alpha and beta which may be accepted
      gradientProcedure.getAlphaLinear(walpha);
      gradientProcedure.getBeta(wbeta);
    }

    return gradientProcedure.value;
  }

  @Override
  protected void computeStep(double[] step) {
    // Alpha = Scaled Hessian matrix
    // beta = Gradient vector
    // We want to solve: A x = b to find the update x

    final int n = gradientProcedure.n;
    System.arraycopy(beta, 0, step, 0, n);
    System.arraycopy(alpha, 0, walpha, 0, alpha.length);
    final double scale = (1.0 + lambda);
    for (int i = 0, j = 0; i < n; i++, j += (n + 1)) {
      // Scale the diagonal of the Hessian to favour direct descent
      walpha[j] *= scale;
    }
    if (!solver.solve(walpha, step)) {
      throw new FunctionSolverException(FitStatus.SINGULAR_NON_LINEAR_MODEL);
    }
  }

  @Override
  protected boolean accept(double currentValue, double[] a, double newValue, double[] newA) {
    if (newValue <= currentValue) {
      // Update the current alpha and beta:
      // We can do this by swapping storage.
      double[] tmp = alpha;
      alpha = walpha;
      walpha = tmp;
      tmp = beta;
      beta = wbeta;
      wbeta = tmp;

      // Decrease Lambda
      lambda *= 0.1;
      return true;
    }

    // Increase Lambda
    lambda *= 10.0;
    return false;
  }

  @Override
  protected double[] prepareFunctionValue(double[] y, double[] a) {
    // Ensure the gradient procedure is created
    y = prepareY(y);
    gradientProcedure = createGradientProcedure(y);
    return y;
  }

  @Override
  protected double computeFunctionValue(double[] a) {
    gradientProcedure.value(a);
    return gradientProcedure.value;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Note: In contrast to {@link #prepareFunctionValue(double[], double[])} this does not create
   * the gradient procedure.
   */
  @Override
  protected double[] prepareFunctionFisherInformationMatrix(double[] y, double[] a) {
    // Do not create the gradient procedure. Sub-classes will create the correct one required.
    return prepareY(y);
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

  @Override
  public void setGradientFunction(GradientFunction function) {
    super.setGradientFunction(function);
    gradientProcedure = null;
  }
}
