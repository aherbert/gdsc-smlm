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

import uk.ac.sussex.gdsc.smlm.fitting.FitStatus;
import uk.ac.sussex.gdsc.smlm.fitting.linear.EjmlLinearSolver;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.FastMleGradient2Procedure;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.FastMleGradient2ProcedureUtils;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.FastMleJacobianGradient2Procedure;
import uk.ac.sussex.gdsc.smlm.function.ExtendedGradient2Function;
import uk.ac.sussex.gdsc.smlm.function.Gradient2Function;
import uk.ac.sussex.gdsc.smlm.function.OffsetExtendedGradient2Function;
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
public class ExtendedFastMleSteppingFunctionSolver extends FastMleSteppingFunctionSolver {
  /** The Jacobian gradient procedure. */
  protected FastMleJacobianGradient2Procedure jacobianGradientProcedure;
  /** The jacobian. */
  protected double[] jacobian;

  /** The solver. */
  protected EjmlLinearSolver solver;

  /** The default max relative error. */
  public static final double DEFAULT_MAX_RELATIVE_ERROR =
      LvmSteppingFunctionSolver.DEFAULT_MAX_RELATIVE_ERROR;
  /** The default max absolute error. */
  public static final double DEFAULT_MAX_ABSOLUTE_ERROR =
      LvmSteppingFunctionSolver.DEFAULT_MAX_ABSOLUTE_ERROR;

  /**
   * Create a new stepping function solver.
   *
   * @param function the function
   * @param maxRelativeError the max relative error
   * @param maxAbsoluteError the max absolute error
   * @throws NullPointerException if the function is null
   */
  public ExtendedFastMleSteppingFunctionSolver(Gradient2Function function, double maxRelativeError,
      double maxAbsoluteError) {
    super(function, maxRelativeError, maxAbsoluteError);
  }

  /**
   * Create a new stepping function solver.
   *
   * @param function the function
   * @param tc the tolerance checker
   * @param bounds the bounds
   * @throws NullPointerException if the function or tolerance checker is null
   */
  public ExtendedFastMleSteppingFunctionSolver(Gradient2Function function, ToleranceChecker tc,
      ParameterBounds bounds) {
    super(function, tc, bounds);
  }

  /**
   * Enable computing the Newton-Raphson step using a full Jacobian solution. This requires
   * computation of the Jacobian of second order partial derivatives with respect to parameters
   * [i,j] and inversion using matrix decomposition.
   *
   * @param enable Set to true to enable
   * @deprecated The computation of the step using the full Jacobian is invalid
   */
  @Deprecated
  void enableJacobianSolution(boolean enable) {
    // This method is defined at the package level for JUnit testing. It is not public
    // as the method does not work.
    enableJacobianSolution(enable, DEFAULT_MAX_RELATIVE_ERROR, DEFAULT_MAX_ABSOLUTE_ERROR);
  }

  /**
   * Enable computing the Newton-Raphson step using a full Jacobian solution. This requires
   * computation of the Jacobian of second order partial derivatives with respect to parameters
   * [i,j] and inversion using matrix decomposition.
   *
   * @param enable Set to true to enable
   * @param maxRelativeError Validate the Jacobian solution using the specified maximum relative
   *        error
   * @param maxAbsoluteError Validate the Jacobian solution using the specified maximum absolute
   *        error
   * @deprecated The computation of the step using the full Jacobian is invalid
   */
  @Deprecated
  void enableJacobianSolution(boolean enable, double maxRelativeError, double maxAbsoluteError) {
    // This method is defined at the package level for JUnit testing. It is not public
    // as the method does not work.
    if (enable) {
      if (!(function instanceof ExtendedGradient2Function)) {
        throw new IllegalStateException(
            "Jacobian requires an " + ExtendedGradient2Function.class.getName());
      }
      solver = new EjmlLinearSolver(maxRelativeError, maxAbsoluteError);
    } else {
      solver = null;
    }
  }

  /**
   * Creates the gradient procedure.
   *
   * @param y the y
   * @return the newton raphson gradient 2 procedure
   */
  @Override
  protected FastMleGradient2Procedure createGradientProcedure(double[] y) {
    // We can handle per-observation variances as detailed in
    // Huang, et al. (2015) by simply adding the variances to the computed value.
    f2 = (Gradient2Function) function;
    if (solver != null) {
      if (obsVariances != null) {
        f2 = OffsetExtendedGradient2Function
            .wrapExtendedGradient2Function((ExtendedGradient2Function) function, obsVariances);
      }
      jacobian = new double[f2.size()];
      return jacobianGradientProcedure =
          new FastMleJacobianGradient2Procedure(y, (ExtendedGradient2Function) f2);
    }
    if (obsVariances != null) {
      f2 = OffsetGradient2Function.wrapGradient2Function(f2, obsVariances);
    }
    return FastMleGradient2ProcedureUtils.create(y, f2);
  }

  /**
   * Compute the gradients for the Newton step using the gradient procedure.
   *
   * @param a the funtion parameters
   */
  @Override
  protected void computeGradients(double[] a) {
    if (solver != null) {
      jacobianGradientProcedure.computeJacobian(a);
    } else {
      gradientProcedure.computeSecondDerivative(a);
    }

    if (gradientProcedure.isNaNGradients()) {
      throw new FunctionSolverException(FitStatus.INVALID_GRADIENTS);
    }
  }

  @Override
  protected void computeStep(double[] step) {
    final double[] d1 = gradientProcedure.d1;
    final double[] d2 = gradientProcedure.d2;

    if (solver != null) {
      // Solve the Jacobian. This is an implementation of the Newton-Raphson method
      // for systems of non-linear equations (see Numerical Recipes in C++, 2nd Ed, section 9.6)
      // XXX This does not work.
      // The the first order derivatives "are not n independent, arbitrary functions,
      // rather they obey so-called integrability conditions that are highly restrictive".
      // This code is deprecated and may be removed.
      for (int i = 0; i < step.length; i++) {
        step[i] = -d1[i];
      }
      jacobianGradientProcedure.getJacobianLinear(jacobian);
      if (solver.solve(jacobian, step)) {
        // XXX - debug the difference
        final double[] step2 = new double[d1.length];
        for (int i = 0; i < step.length; i++) {
          step2[i] = -d1[i] / d2[i];
        }
        return;
      }
    }

    // Simple Newton-Raphson update step as per Smith et al, (2010), SI Eq. 13:
    // parameter -> new parameter + delta
    // => new parameter = parameter - delta
    for (int i = 0; i < step.length; i++) {
      step[i] = -d1[i] / d2[i];
    }
  }
}
