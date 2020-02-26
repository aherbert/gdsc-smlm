/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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

import uk.ac.sussex.gdsc.core.data.NotImplementedException;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.smlm.fitting.FisherInformationMatrix;
import uk.ac.sussex.gdsc.smlm.fitting.FunctionSolverType;
import uk.ac.sussex.gdsc.smlm.fitting.LseFunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.linear.EjmlLinearSolver;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.LsqLvmGradientProcedureUtils;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.LsqVarianceGradientProcedure;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.LsqVarianceGradientProcedureUtils;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.LvmGradientProcedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;
import uk.ac.sussex.gdsc.smlm.function.Gradient2FunctionValueStore;

/**
 * Uses the Levenberg-Marquardt (LVM) method to fit a gradient function with coefficients (a) using
 * least squares estimation (LSE).
 */
public class LseLvmSteppingFunctionSolver extends LvmSteppingFunctionSolver
    implements LseFunctionSolver {
  /**
   * The solver used to invert the Fisher information matrix to find the Cramér–Rao lower bound
   * (CRLB).
   */
  protected EjmlLinearSolver inversionSolver;

  /** The total sum of squares. */
  protected double totalSumOfSquares = Double.NaN;

  /**
   * Create a new stepping function solver.
   *
   * @param function the function
   * @throws NullPointerException if the function is null
   */
  public LseLvmSteppingFunctionSolver(Gradient1Function function) {
    super(FunctionSolverType.LSE, function);
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
  public LseLvmSteppingFunctionSolver(Gradient1Function function, double maxRelativeError,
      double maxAbsoluteError) {
    super(FunctionSolverType.LSE, function, maxRelativeError, maxAbsoluteError);
  }

  /**
   * Create a new stepping function solver.
   *
   * @param function the function
   * @param tc the tolerance checker
   * @param bounds the bounds
   * @throws NullPointerException if the function or tolerance checker is null
   */
  public LseLvmSteppingFunctionSolver(Gradient1Function function, ToleranceChecker tc,
      ParameterBounds bounds) {
    super(FunctionSolverType.LSE, function, tc, bounds);
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
  public LseLvmSteppingFunctionSolver(Gradient1Function function, ToleranceChecker tc,
      ParameterBounds bounds, double maxRelativeError, double maxAbsoluteError) {
    super(FunctionSolverType.LSE, function, tc, bounds, maxRelativeError, maxAbsoluteError);
  }

  @Override
  protected void preProcess() {
    totalSumOfSquares = Double.NaN;
  }

  @Override
  protected LvmGradientProcedure createGradientProcedure(double[] y) {
    return LsqLvmGradientProcedureUtils.create(y, (Gradient1Function) function);
  }

  @Override
  protected void computeDeviationsAndValues(double[] parametersVariance, double[] fx) {
    Gradient1Function f1 = (Gradient1Function) this.function;
    // Capture the y-values if necessary
    if (fx != null && fx.length == f1.size()) {
      f1 = new Gradient2FunctionValueStore(f1, fx);
    }
    final LsqVarianceGradientProcedure p = createVarianceProcedure(f1);
    if (p.variance(null) == LsqVarianceGradientProcedure.STATUS_OK) {
      setDeviations(parametersVariance, p.variance);
    }
  }

  private LsqVarianceGradientProcedure createVarianceProcedure(Gradient1Function function) {
    if (inversionSolver == null) {
      inversionSolver = EjmlLinearSolver.createForInversion(1e-2);
    }
    return LsqVarianceGradientProcedureUtils.create(function, inversionSolver);
  }

  @Override
  public boolean computeDeviations(double[] y, double[] a, double[] parametersVariance) {
    final LsqVarianceGradientProcedure p = createVarianceProcedure((Gradient1Function) function);
    if (p.variance(a) == LsqVarianceGradientProcedure.STATUS_OK) {
      setDeviations(parametersVariance, p.variance);
      return true;
    }
    return false;
  }

  @Override
  protected FisherInformationMatrix computeLastFisherInformationMatrix(double[] fx) {
    // This solver directly implements computation of the deviations
    throw new NotImplementedException();
  }

  @Override
  protected FisherInformationMatrix computeFunctionFisherInformationMatrix(double[] y, double[] a) {
    // This solver directly implements computation of the deviations
    throw new NotImplementedException();
  }

  @Override
  public double getTotalSumOfSquares() {
    if (Double.isNaN(totalSumOfSquares) && lastY != null) {
      totalSumOfSquares = LseBaseFunctionSolver.computeTotalSumOfSquares(lastY);
    }
    return totalSumOfSquares;
  }

  @Override
  public double getResidualSumOfSquares() {
    return value;
  }

  @Override
  public double getCoefficientOfDetermination() {
    return 1.0 - (value / getTotalSumOfSquares());
  }

  @Override
  public double getAdjustedCoefficientOfDetermination() {
    return MathUtils.getAdjustedCoefficientOfDetermination(value, getTotalSumOfSquares(),
        getNumberOfFittedPoints(), getNumberOfFittedParameters());
  }

  @Override
  public double getMeanSquaredError() {
    return value / (getNumberOfFittedPoints() - getNumberOfFittedParameters());
  }

  @Override
  public boolean isWeighted() {
    // This is a stepping solver that is not weighted
    return false;
  }
}
