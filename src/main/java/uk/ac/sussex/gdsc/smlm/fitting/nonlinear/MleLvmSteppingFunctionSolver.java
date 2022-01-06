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

import uk.ac.sussex.gdsc.smlm.fitting.FisherInformationMatrix;
import uk.ac.sussex.gdsc.smlm.fitting.FitStatus;
import uk.ac.sussex.gdsc.smlm.fitting.FunctionSolverType;
import uk.ac.sussex.gdsc.smlm.fitting.MleFunctionSolver;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.LvmGradientProcedure;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.MleLvmGradientProcedureUtils;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.PoissonGradientProcedure;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.PoissonGradientProcedureUtils;
import uk.ac.sussex.gdsc.smlm.function.ChiSquaredDistributionTable;
import uk.ac.sussex.gdsc.smlm.function.FastLog;
import uk.ac.sussex.gdsc.smlm.function.FastLogFactory;
import uk.ac.sussex.gdsc.smlm.function.Gradient1Function;
import uk.ac.sussex.gdsc.smlm.function.Gradient2FunctionValueStore;
import uk.ac.sussex.gdsc.smlm.function.OffsetGradient1Function;
import uk.ac.sussex.gdsc.smlm.function.PoissonCalculator;

/**
 * Uses the Levenberg-Marquardt (LVM) method to fit a gradient function with coefficients (a) using
 * maximum likelihood estimation.
 *
 * <p>This solver computes a modified Chi-squared expression to perform Maximum Likelihood
 * Estimation (MLE) assuming Poisson model. See Laurence &amp; Chromy (2010) Efficient maximum
 * likelihood estimator. Nature Methods 7, 338-339. The input data must be Poisson distributed for
 * this to be relevant.
 *
 * <p>Per observation variances can be added to both the target x data and the function value to
 * optimise the LLR sCMOS function for Poisson data. See Huang et al, (2015). Video-rate nanoscopy
 * using sCMOS camera–specific single-molecule localization algorithms. Nature Methods 10, 653–658.
 */
public class MleLvmSteppingFunctionSolver extends LvmSteppingFunctionSolver
    implements MleFunctionSolver {
  /** The last Y fit. */
  protected double[] lastyFit;

  /** The ll. */
  protected double ll = Double.NaN;

  /**
   * The per observation variances. This is not null if fitting using the method of Huang, et al
   * (2015).
   */
  protected double[] weights;
  /**
   * The gradient function used by the procedure. This may be wrapped to add the per observation
   * variances if fitting using the LLR sCMOS method of Huang, et al (2015).
   */
  protected Gradient1Function f1;

  /** The fast log instance for the fast log version of the procedure. */
  private FastLog fastLog;

  /**
   * Create a new stepping function solver.
   *
   * @param function the function
   * @throws NullPointerException if the function is null
   */
  public MleLvmSteppingFunctionSolver(Gradient1Function function) {
    super(FunctionSolverType.MLE, function);
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
  public MleLvmSteppingFunctionSolver(Gradient1Function function, double maxRelativeError,
      double maxAbsoluteError) {
    super(FunctionSolverType.MLE, function, maxRelativeError, maxAbsoluteError);
  }

  /**
   * Create a new stepping function solver.
   *
   * @param function the function
   * @param tc the tolerance checker
   * @param bounds the bounds
   * @throws NullPointerException if the function or tolerance checker is null
   */
  public MleLvmSteppingFunctionSolver(Gradient1Function function, ToleranceChecker tc,
      ParameterBounds bounds) {
    super(FunctionSolverType.MLE, function, tc, bounds);
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
  public MleLvmSteppingFunctionSolver(Gradient1Function function, ToleranceChecker tc,
      ParameterBounds bounds, double maxRelativeError, double maxAbsoluteError) {
    super(FunctionSolverType.MLE, function, tc, bounds, maxRelativeError, maxAbsoluteError);
  }

  @Override
  protected void preProcess() {
    ll = Double.NaN;
    lastyFit = null;
  }

  @Override
  protected double[] prepareY(double[] y) {
    // We can handle per-observation variances as detailed in
    // Huang, et al. (2015) by simply adding the variances to the target data.

    final int n = y.length;
    weights = getWeights(n);
    if (weights != null) {
      final double[] x = new double[n];
      for (int i = 0; i < n; i++) {
        // Also ensure the input y is positive
        x[i] = (y[i] > 0) ? y[i] + weights[i] : weights[i];
      }
      return x;
    }
    return ensurePositive(y);
  }

  @Override
  protected LvmGradientProcedure createGradientProcedure(double[] y) {
    // We can handle per-observation variances as detailed in
    // Huang, et al. (2015) by simply adding the variances to the computed value.
    f1 = (Gradient1Function) function;
    if (weights != null) {
      f1 = OffsetGradient1Function.wrapGradient1Function(f1, weights);
    }
    if (isFastLog()) {
      return MleLvmGradientProcedureUtils.create(y, f1, fastLog);
    }
    return MleLvmGradientProcedureUtils.create(y, f1);
  }

  @Override
  protected double computeFitValue(double[] a) {
    // Cache this
    lastA = a;
    return super.computeFitValue(a);
  }

  @Override
  protected void computeValues(double[] fx) {
    super.computeValues(fx);

    // Cache the values to compute the log-likelihood
    final int size = f1.size();
    if (lastyFit == null) {
      lastyFit = new double[size];
    }
    if (weights != null) {
      // For the log-likelihood we must add the per observation weights
      for (int i = 0; i < size; i++) {
        lastyFit[i] = fx[i] + weights[i];
      }
    } else {
      System.arraycopy(fx, 0, lastyFit, 0, size);
    }
  }

  @Override
  protected FisherInformationMatrix computeLastFisherInformationMatrix(double[] fx) {
    // The Hessian matrix refers to the log-likelihood ratio.
    // Compute and invert a matrix related to the Poisson log-likelihood.
    // This assumes this does achieve the maximum likelihood estimate for a
    // Poisson process.
    Gradient1Function localF1 = (Gradient1Function) function;
    // Capture the y-values if necessary
    if (fx != null && fx.length == localF1.size()) {
      localF1 = new Gradient2FunctionValueStore(localF1, fx);
    }
    // Add the weights if necessary
    if (weights != null) {
      localF1 = OffsetGradient1Function.wrapGradient1Function(localF1, weights);
    }
    final PoissonGradientProcedure p = PoissonGradientProcedureUtils.create(localF1);
    p.computeFisherInformation(lastA);
    if (p.isNaNGradients()) {
      throw new FunctionSolverException(FitStatus.INVALID_GRADIENTS);
    }
    p.getLinear(walpha); // Re-use space
    return new FisherInformationMatrix(walpha, beta.length);
  }

  @Override
  protected FisherInformationMatrix computeFunctionFisherInformationMatrix(double[] y, double[] a) {
    // Compute and invert a matrix related to the Poisson log-likelihood.
    // This assumes this does achieve the maximum likelihood estimate for a
    // Poisson process.
    // We must wrap the gradient function if weights are present.
    Gradient1Function localF1 = (Gradient1Function) function;
    if (weights != null) {
      localF1 = OffsetGradient1Function.wrapGradient1Function(localF1, weights);
    }
    final PoissonGradientProcedure p = PoissonGradientProcedureUtils.create(localF1);
    p.computeFisherInformation(a);
    if (p.isNaNGradients()) {
      throw new FunctionSolverException(FitStatus.INVALID_GRADIENTS);
    }
    return new FisherInformationMatrix(p.getLinear(), function.getNumberOfGradients());
  }

  @Override
  public double getLogLikelihood() {
    if (Double.isNaN(ll)) {
      if (lastyFit == null) {
        // Evaluate the function values if necessary
        final int size = f1.size();
        lastyFit = new double[size];
        super.computeValues(lastyFit);

        if (weights != null) {
          // For the log-likelihood we must add the per observation weights
          for (int i = 0; i < weights.length; i++) {
            lastyFit[i] += weights[i];
          }
        }
      }

      // ll = PoissonCalculator.fastLogLikelihood(lastyFit, lastY)

      // This has a relative error of <1e-4 and is 50% faster than fastLogLikelihood.
      // The value is only used for reporting and so high accuracy is not essential.
      ll = PoissonCalculator.fastLogLikelihood(lastyFit, lastY, FastLogFactory.getFastLog());
    }
    return ll;
  }

  @Override
  public double getLogLikelihoodRatio() {
    // This method computes the log-likelihood ratio directly
    return value;
  }

  @Override
  public double getQ() {
    // Wilks theorum states the LLR approaches the chi-squared distribution for large n.
    return ChiSquaredDistributionTable.computeQValue(getLogLikelihoodRatio(),
        getNumberOfFittedPoints() - getNumberOfFittedParameters());
  }

  /**
   * Gets the fast log instance.
   *
   * @return the fast log
   */
  public FastLog getFastLog() {
    return fastLog;
  }

  /**
   * Checks if is using a fast log instance.
   *
   * @return true, if using a fast log instance
   */
  public boolean isFastLog() {
    return fastLog != null;
  }

  /**
   * Sets the fast log instance to use for the MLE LVM procedure. This may decrease stability on
   * convergence (since the function value has less precision) and should be used with caution.
   *
   * @param fastLog the new fast log instance
   */
  public void setFastLog(FastLog fastLog) {
    this.fastLog = fastLog;
  }
}
