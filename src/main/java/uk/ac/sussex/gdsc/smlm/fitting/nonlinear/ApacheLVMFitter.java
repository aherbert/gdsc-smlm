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
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator;
import uk.ac.sussex.gdsc.smlm.fitting.nonlinear.gradient.GradientCalculatorFactory;
import uk.ac.sussex.gdsc.smlm.function.ExtendedNonLinearFunction;
import uk.ac.sussex.gdsc.smlm.function.MultivariateMatrixFunctionWrapper;
import uk.ac.sussex.gdsc.smlm.function.MultivariateVectorFunctionWrapper;
import uk.ac.sussex.gdsc.smlm.function.NonLinearFunction;
import uk.ac.sussex.gdsc.smlm.function.ValueProcedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;

import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.exception.TooManyIterationsException;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer.Optimum;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.ValueAndJacobianFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;
import org.apache.commons.math3.util.Precision;

/**
 * Uses Apache Commons Math Levenberg-Marquardt method to fit a nonlinear model with coefficients
 * (a) for a set of data points (x, y).
 */
public class ApacheLVMFitter extends LSEBaseFunctionSolver {
  /**
   * Default constructor.
   *
   * @param gf the Gaussian function
   */
  public ApacheLVMFitter(Gaussian2DFunction gf) {
    super(gf);
  }

  @Override
  public FitStatus computeFit(double[] y, final double[] fx, double[] a,
      double[] parametersVariance) {
    final int n = y.length;
    try {
      // Different convergence thresholds seem to have no effect on the resulting fit, only the
      // number of
      // iterations for convergence
      final double initialStepBoundFactor = 100;
      final double costRelativeTolerance = 1e-10;
      final double parRelativeTolerance = 1e-10;
      final double orthoTolerance = 1e-10;
      final double threshold = Precision.SAFE_MIN;

      // Extract the parameters to be fitted
      final double[] initialSolution = getInitialSolution(a);

      // TODO - Pass in more advanced stopping criteria.

      // Create the target and weight arrays
      final double[] yd = new double[n];
      // final double[] w = new double[n];
      for (int i = 0; i < n; i++) {
        yd[i] = y[i];
        // w[i] = 1;
      }

      final LevenbergMarquardtOptimizer optimizer =
          new LevenbergMarquardtOptimizer(initialStepBoundFactor, costRelativeTolerance,
              parRelativeTolerance, orthoTolerance, threshold);

      //@formatter:off
      final LeastSquaresBuilder builder = new LeastSquaresBuilder()
          .maxEvaluations(Integer.MAX_VALUE)
          .maxIterations(getMaxEvaluations())
          .start(initialSolution)
          .target(yd);
          // This is not required
          //.weight(new DiagonalMatrix(w));
      //@formatter:on

      if (function instanceof ExtendedNonLinearFunction
          && ((ExtendedNonLinearFunction) function).canComputeValuesAndJacobian()) {
        // Compute together, or each individually
        builder.model(new ValueAndJacobianFunction() {
          final ExtendedNonLinearFunction fun = (ExtendedNonLinearFunction) function;

          @Override
          public Pair<RealVector, RealMatrix> value(RealVector point) {
            final double[] p = point.toArray();
            final uk.ac.sussex.gdsc.smlm.utils.Pair<double[], double[][]> result =
                fun.computeValuesAndJacobian(p);
            return new Pair<>(new ArrayRealVector(result.a, false),
                new Array2DRowRealMatrix(result.b, false));
          }

          @Override
          public RealVector computeValue(double[] params) {
            return new ArrayRealVector(fun.computeValues(params), false);
          }

          @Override
          public RealMatrix computeJacobian(double[] params) {
            return new Array2DRowRealMatrix(fun.computeJacobian(params), false);
          }
        });
      } else {
        // Compute separately
        builder.model(new MultivariateVectorFunctionWrapper((NonLinearFunction) function, a, n),
            new MultivariateMatrixFunctionWrapper((NonLinearFunction) function, a, n));
      }

      final LeastSquaresProblem problem = builder.build();

      final Optimum optimum = optimizer.optimize(problem);

      final double[] parameters = optimum.getPoint().toArray();
      setSolution(a, parameters);
      iterations = optimum.getIterations();
      evaluations = optimum.getEvaluations();
      if (parametersVariance != null) {
        // Set up the Jacobian.
        final RealMatrix j = optimum.getJacobian();

        // Compute transpose(J)J.
        final RealMatrix jTj = j.transpose().multiply(j);

        final double[][] data =
            (jTj instanceof Array2DRowRealMatrix) ? ((Array2DRowRealMatrix) jTj).getDataRef()
                : jTj.getData();
        final FisherInformationMatrix m = new FisherInformationMatrix(data);
        setDeviations(parametersVariance, m);
      }
      // Compute function value
      if (fx != null) {
        final Gaussian2DFunction function = (Gaussian2DFunction) this.function;
        function.initialise0(a);
        function.forEach(new ValueProcedure() {
          int i;

          @Override
          public void execute(double value) {
            fx[i] = value;
          }
        });
      }

      // As this is unweighted then we can do this to get the sum of squared residuals
      // This is the same as optimum.getCost() * optimum.getCost(); The getCost() function
      // just computes the dot product anyway.
      value = optimum.getResiduals().dotProduct(optimum.getResiduals());
    } catch (final TooManyEvaluationsException ex) {
      return FitStatus.TOO_MANY_EVALUATIONS;
    } catch (final TooManyIterationsException ex) {
      return FitStatus.TOO_MANY_ITERATIONS;
    } catch (final ConvergenceException ex) {
      // Occurs when QR decomposition fails - mark as a singular non-linear model (no solution)
      return FitStatus.SINGULAR_NON_LINEAR_MODEL;
    } catch (final Exception ex) {
      // TODO - Find out the other exceptions from the fitter and add return values to match.
      return FitStatus.UNKNOWN;
    }

    return FitStatus.OK;
  }

  @Override
  public boolean computeValue(double[] y, double[] fx, double[] a) {
    final GradientCalculator calculator =
        GradientCalculatorFactory.newCalculator(function.getNumberOfGradients(), false);

    // Since we know the function is a Gaussian2DFunction from the constructor
    value = calculator.findLinearised(y.length, y, fx, a, (NonLinearFunction) function);

    return true;
  }

  @Override
  protected FisherInformationMatrix computeFisherInformationMatrix(double[] y, double[] a) {
    final GradientCalculator c =
        GradientCalculatorFactory.newCalculator(function.getNumberOfGradients(), false);
    // Since we know the function is a Gaussian2DFunction from the constructor
    final double[][] I = c.fisherInformationMatrix(y.length, a, (NonLinearFunction) function);
    if (c.isNaNGradients()) {
      throw new FunctionSolverException(FitStatus.INVALID_GRADIENTS);
    }
    return new FisherInformationMatrix(I);
  }
}
