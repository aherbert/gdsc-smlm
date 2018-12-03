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

import uk.ac.sussex.gdsc.smlm.fitting.FunctionSolverType;
import uk.ac.sussex.gdsc.smlm.fitting.MLEFunctionSolver;
import uk.ac.sussex.gdsc.smlm.function.ChiSquaredDistributionTable;
import uk.ac.sussex.gdsc.smlm.function.GradientFunction;

/**
 * Abstract class with utility methods for the MLEFunctionSolver interface.
 */
public abstract class MLEBaseFunctionSolver extends BaseFunctionSolver
    implements MLEFunctionSolver {
  /** The log-likelihood ratio. */
  protected double llr = Double.NaN;

  /**
   * Default constructor.
   *
   * @param f the function
   * @throws NullPointerException if the function is null
   */
  public MLEBaseFunctionSolver(GradientFunction f) {
    super(FunctionSolverType.MLE, f);
  }

  @Override
  protected void preProcess() {
    llr = Double.NaN;
  }

  /** {@inheritDoc} */
  @Override
  public double getLogLikelihood() {
    return value;
  }

  /** {@inheritDoc} */
  @Override
  public double getLogLikelihoodRatio() {
    if (Double.isNaN(llr) && lastY != null) {
      // From https://en.wikipedia.org/wiki/Likelihood-ratio_test#Use:
      // LLR = 2 * [ ln(likelihood for alternative model) - ln(likelihood for null model)]
      // The model with more parameters (here alternative) will always fit at least as well—
      // i.e., have the same or greater log-likelihood—than the model with fewer parameters
      // (here null)

      final double llAlternative = computeObservedLogLikelihood(lastY, lastA);
      final double llNull = getLogLikelihood();

      // The alternative should always fit better (higher value) than the null model
      if (llAlternative < llNull) {
        llr = 0;
      } else {
        llr = 2 * (llAlternative - llNull);
      }
    }
    return llr;
  }

  /**
   * Compute the observed log likelihood (i.e. the log-likelihood with y as the function value).
   *
   * @param y the y
   * @param a the a
   * @return the observed log likelihood
   */
  protected abstract double computeObservedLogLikelihood(double[] y, double[] a);

  /** {@inheritDoc} */
  @Override
  public double getQ() {
    return ChiSquaredDistributionTable.computeQValue(getLogLikelihoodRatio(),
        getNumberOfFittedPoints() - getNumberOfFittedParameters());
  }
}
