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

import uk.ac.sussex.gdsc.smlm.fitting.FisherInformationMatrix;
import uk.ac.sussex.gdsc.smlm.fitting.FitStatus;
import uk.ac.sussex.gdsc.smlm.fitting.FunctionSolverType;
import uk.ac.sussex.gdsc.smlm.fitting.WLseFunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.LvmGradientProcedure;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.WLsqLvmGradientProcedureUtils;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.WPoissonGradientProcedure;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.WPoissonGradientProcedureUtils;
import uk.ac.sussex.gdsc.smlm.function.ChiSquaredDistributionTable;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;

/**
 * Uses the Levenberg-Marquardt (LVM) method to fit a gradient function with coefficients (a) using
 * weighted least squares estimation.
 *
 * <p>This solver computes a modified Chi-squared expression to perform Weighted Least Squares
 * Estimation (WLSE) assuming a Poisson model with a Gaussian noise component. The weight per
 * observation is equal to 1/[variance + max(y, 0) + 1].
 *
 * <p>See Ruisheng, et al (2017) Algorithmic corrections for localization microscopy with sCMOS
 * cameras - characterisation of a computationally efficient localization approach. Optical Express
 * 25, Issue 10, pp 11701-11716.
 */
public class WLseLvmSteppingFunctionSolver extends LvmSteppingFunctionSolver
    implements WLseFunctionSolver {
  /**
   * Create a new stepping function solver.
   *
   * @param function the function
   * @throws NullPointerException if the function is null
   */
  public WLseLvmSteppingFunctionSolver(Gradient1Function function) {
    super(FunctionSolverType.WLSE, function);
  }

  /**
   * Create a new stepping function solver.
   *
   * @param function the function
   * @param maxRelativeError Validate the Levenberg-Marquardt fit solution using the specified
   *        maximum relative error
   * @param maxAbsoluteError Validate the Levenberg-Marquardt fit solution using the specified
   *        maximum absolute error
   * @throws NullPointerException if the function is null
   */
  public WLseLvmSteppingFunctionSolver(Gradient1Function function, double maxRelativeError,
      double maxAbsoluteError) {
    super(FunctionSolverType.WLSE, function, maxRelativeError, maxAbsoluteError);
  }

  /**
   * Create a new stepping function solver.
   *
   * @param function the function
   * @param tc the tolerance checker
   * @param bounds the bounds
   * @throws NullPointerException if the function or tolerance checker is null
   */
  public WLseLvmSteppingFunctionSolver(Gradient1Function function, ToleranceChecker tc,
      ParameterBounds bounds) {
    super(FunctionSolverType.WLSE, function, tc, bounds);
  }

  /**
   * Create a new stepping function solver.
   *
   * @param function the function
   * @param tc the tolerance checker
   * @param bounds the bounds
   * @param maxRelativeError Validate the Levenberg-Marquardt fit solution using the specified
   *        maximum relative error
   * @param maxAbsoluteError Validate the Levenberg-Marquardt fit solution using the specified
   *        maximum absolute error
   * @throws NullPointerException if the function or tolerance checker is null
   */
  public WLseLvmSteppingFunctionSolver(Gradient1Function function, ToleranceChecker tc,
      ParameterBounds bounds, double maxRelativeError, double maxAbsoluteError) {
    super(FunctionSolverType.WLSE, function, tc, bounds, maxRelativeError, maxAbsoluteError);
  }

  @Override
  protected LvmGradientProcedure createGradientProcedure(double[] y) {
    return WLsqLvmGradientProcedureUtils.create(y, getWeights(y.length),
        (Gradient1Function) function);
  }

  @Override
  protected FisherInformationMatrix computeLastFisherInformationMatrix(double[] fx) {
    // Get the values if necessary
    if (fx != null && fx.length == ((Gradient1Function) function).size()) {
      computeValues(fx);
    }

    // TODO. Check if these deviations are correct.
    // Note: Huang et al (2015) compute:
    // Iab = 1 / (uk + var/g^2) * duda * dudb
    // with uk the expected photon count.
    // This will compute:
    // Iab = 1 / (xk + 1.0 + var/g^2) * duda * dudb
    // with xk the observed photon count.

    // The last Hessian matrix should be stored in the working alpha.
    return new FisherInformationMatrix(walpha, beta.length);
  }

  @Override
  protected FisherInformationMatrix computeFunctionFisherInformationMatrix(double[] y, double[] a) {
    // Compute using the scaled Hessian as per the above method.
    // Use a dedicated procedure that omits computing beta.
    final WPoissonGradientProcedure p = WPoissonGradientProcedureUtils.create(y,
        getWeights(y.length), (Gradient1Function) function);
    p.computeFisherInformation(a);
    if (p.isNaNGradients()) {
      throw new FunctionSolverException(FitStatus.INVALID_GRADIENTS);
    }
    return new FisherInformationMatrix(p.getLinear(), p.numberOfGradients);
  }

  @Override
  public double getChiSquared() {
    return value;
  }

  @Override
  public double getQ() {
    // Value will be the Chi-squared
    return ChiSquaredDistributionTable.computeQValue(value,
        getNumberOfFittedPoints() - getNumberOfFittedParameters());
  }
}
