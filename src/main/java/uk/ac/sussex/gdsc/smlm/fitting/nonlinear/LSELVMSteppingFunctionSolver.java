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
package uk.ac.sussex.gdsc.smlm.fitting.nonlinear;

import uk.ac.sussex.gdsc.core.data.NotImplementedException;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.smlm.fitting.FisherInformationMatrix;
import uk.ac.sussex.gdsc.smlm.fitting.FunctionSolverType;
import uk.ac.sussex.gdsc.smlm.fitting.LSEFunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.linear.EJMLLinearSolver;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.LSQLVMGradientProcedureFactory;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.LSQVarianceGradientProcedure;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.LSQVarianceGradientProcedureFactory;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.LVMGradientProcedure;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;
import uk.ac.sussex.gdsc.smlm.function.Gradient2FunctionValueStore;

/**
 * Uses the Levenberg-Marquardt method to fit a gradient function with coefficients (a) using least
 * squares estimation.
 */
public class LSELVMSteppingFunctionSolver extends LVMSteppingFunctionSolver
    implements LSEFunctionSolver {
  /**
   * The solver used to invert the Fisher information matrix to find the Cramér–Rao lower bound
   * (CRLB).
   */
  protected EJMLLinearSolver inversionSolver;

  /** The total sum of squares. */
  protected double totalSumOfSquares = Double.NaN;

  /**
   * Create a new stepping function solver.
   *
   * @param f the function
   * @throws NullPointerException if the function is null
   */
  public LSELVMSteppingFunctionSolver(Gradient1Function f) {
    super(FunctionSolverType.LSE, f);
  }

  /**
   * Create a new stepping function solver.
   *
   * @param f the function
   * @param maxRelativeError Validate the Levenberg-Marquardt fit solution using the specified
   *        maximum relative error
   * @param maxAbsoluteError Validate the Levenberg-Marquardt fit solution using the specified
   *        maximum absolute error
   * @throws NullPointerException if the function is null
   */
  public LSELVMSteppingFunctionSolver(Gradient1Function f, double maxRelativeError,
      double maxAbsoluteError) {
    super(FunctionSolverType.LSE, f, maxRelativeError, maxAbsoluteError);
  }

  /**
   * Create a new stepping function solver.
   *
   * @param f the function
   * @param tc the tolerance checker
   * @param bounds the bounds
   * @throws NullPointerException if the function or tolerance checker is null
   */
  public LSELVMSteppingFunctionSolver(Gradient1Function f, ToleranceChecker tc,
      ParameterBounds bounds) {
    super(FunctionSolverType.LSE, f, tc, bounds);
  }

  /**
   * Create a new stepping function solver.
   *
   * @param f the function
   * @param tc the tolerance checker
   * @param bounds the bounds
   * @param maxRelativeError Validate the Levenberg-Marquardt fit solution using the specified
   *        maximum relative error
   * @param maxAbsoluteError Validate the Levenberg-Marquardt fit solution using the specified
   *        maximum absolute error
   * @throws NullPointerException if the function or tolerance checker is null
   */
  public LSELVMSteppingFunctionSolver(Gradient1Function f, ToleranceChecker tc,
      ParameterBounds bounds, double maxRelativeError, double maxAbsoluteError) {
    super(FunctionSolverType.LSE, f, tc, bounds, maxRelativeError, maxAbsoluteError);
  }

  /** {@inheritDoc} */
  @Override
  protected void preProcess() {
    totalSumOfSquares = Double.NaN;
  }

  /** {@inheritDoc} */
  @Override
  protected LVMGradientProcedure createGradientProcedure(double[] y) {
    return LSQLVMGradientProcedureFactory.create(y, (Gradient1Function) f);
  }

  @Override
  protected void computeDeviationsAndValues(double[] aDev, double[] yFit) {
    Gradient1Function f1 = (Gradient1Function) this.f;
    // Capture the y-values if necessary
    if (yFit != null && yFit.length == f1.size()) {
      f1 = new Gradient2FunctionValueStore(f1, yFit);
    }
    final LSQVarianceGradientProcedure p = createVarianceProcedure(f1);
    if (p.variance(null) == LSQVarianceGradientProcedure.STATUS_OK) {
      setDeviations(aDev, p.variance);
    }
  }

  private LSQVarianceGradientProcedure createVarianceProcedure(Gradient1Function f) {
    if (inversionSolver == null) {
      inversionSolver = EJMLLinearSolver.createForInversion(1e-2);
    }
    return LSQVarianceGradientProcedureFactory.create(f, inversionSolver);
  }

  @Override
  public boolean computeDeviations(double[] y, double[] a, double[] aDev) {
    final LSQVarianceGradientProcedure p = createVarianceProcedure((Gradient1Function) f);
    if (p.variance(a) == LSQVarianceGradientProcedure.STATUS_OK) {
      setDeviations(aDev, p.variance);
      return true;
    }
    return false;
  }

  @Override
  protected FisherInformationMatrix computeFisherInformationMatrix(double[] yFit) {
    // This solver directly implements computation of the deviations
    throw new NotImplementedException();
  }

  @Override
  protected FisherInformationMatrix computeFunctionFisherInformationMatrix(double[] y, double[] a) {
    // This solver directly implements computation of the deviations
    throw new NotImplementedException();
  }

  // @Override
  // protected FisherInformationMatrix computeFisherInformationMatrix()
  // {
  // // TODO. Check if these deviations are correct.
  // // The last Hessian matrix should be stored in the working alpha.
  // return new FisherInformationMatrix(walpha, beta.length);
  // }
  //
  // @Override
  // protected FisherInformationMatrix computeFunctionFisherInformationMatrix(double[] y, double[]
  // a)
  // {
  // // Compute using the scaled Hessian as per the above method .
  // // TODO - Use a dedicated procedure that omits computing beta.
  // if (gradientProcedure == null)
  // gradientProcedure = createGradientProcedure(y);
  // gradientProcedure.gradient(a);
  // if (gradientProcedure.isNaNGradients())
  // throw new FunctionSolverException(FitStatus.INVALID_GRADIENTS);
  // return new FisherInformationMatrix(gradientProcedure.getAlphaLinear(),
  // f.getNumberOfGradients());
  // }

  /** {@inheritDoc} */
  @Override
  public double getTotalSumOfSquares() {
    if (Double.isNaN(totalSumOfSquares) && lastY != null) {
      totalSumOfSquares = LSEBaseFunctionSolver.getTotalSumOfSquares(lastY);
    }
    return totalSumOfSquares;
  }

  /** {@inheritDoc} */
  @Override
  public double getResidualSumOfSquares() {
    return value;
  }

  /** {@inheritDoc} */
  @Override
  public double getCoefficientOfDetermination() {
    return 1.0 - (value / getTotalSumOfSquares());
  }

  /** {@inheritDoc} */
  @Override
  public double getAdjustedCoefficientOfDetermination() {
    return MathUtils.getAdjustedCoefficientOfDetermination(value, getTotalSumOfSquares(),
        getNumberOfFittedPoints(), getNumberOfFittedParameters());
  }

  /** {@inheritDoc} */
  @Override
  public double getMeanSquaredError() {
    return value / (getNumberOfFittedPoints() - getNumberOfFittedParameters());
  }

  /** {@inheritDoc} */
  @Override
  public boolean isWeighted() {
    // This is a stepping solver that is not weighted
    return false;
  }
}
