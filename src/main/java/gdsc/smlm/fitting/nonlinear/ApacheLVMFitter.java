package gdsc.smlm.fitting.nonlinear;

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
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.util.Pair;
import org.apache.commons.math3.util.Precision;

import gdsc.smlm.fitting.FisherInformationMatrix;
import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator;
import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculatorFactory;
import gdsc.smlm.function.ExtendedNonLinearFunction;
import gdsc.smlm.function.MultivariateMatrixFunctionWrapper;
import gdsc.smlm.function.MultivariateVectorFunctionWrapper;
import gdsc.smlm.function.NonLinearFunction;
import gdsc.smlm.function.ValueProcedure;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Uses Apache Commons Math Levenberg-Marquardt method to fit a nonlinear model with coefficients (a) for a
 * set of data points (x, y).
 */
public class ApacheLVMFitter extends LSEBaseFunctionSolver
{
	/**
	 * Default constructor
	 */
	public ApacheLVMFitter(Gaussian2DFunction gf)
	{
		super(gf);
	}

	public FitStatus computeFit(double[] y, final double[] y_fit, double[] a, double[] a_dev)
	{
		int n = y.length;
		try
		{
			// Different convergence thresholds seem to have no effect on the resulting fit, only the number of
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
			final double[] w = new double[n];
			for (int i = 0; i < n; i++)
			{
				yd[i] = y[i];
				w[i] = 1;
			}

			LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer(initialStepBoundFactor,
					costRelativeTolerance, parRelativeTolerance, orthoTolerance, threshold);

			//@formatter:off
			LeastSquaresBuilder builder = new LeastSquaresBuilder()
					.maxEvaluations(Integer.MAX_VALUE)
					.maxIterations(getMaxEvaluations())
					.start(initialSolution)
					.target(yd)
					.weight(new DiagonalMatrix(w));
			//@formatter:on

			if (f instanceof ExtendedNonLinearFunction && ((ExtendedNonLinearFunction) f).canComputeValuesAndJacobian())
			{
				// Compute together, or each individually
				builder.model(new ValueAndJacobianFunction()
				{
					final ExtendedNonLinearFunction fun = (ExtendedNonLinearFunction) f;

					public Pair<RealVector, RealMatrix> value(RealVector point)
					{
						final double[] p = point.toArray();
						final Pair<double[], double[][]> result = fun.computeValuesAndJacobian(p);
						return new Pair<RealVector, RealMatrix>(new ArrayRealVector(result.getFirst(), false),
								new Array2DRowRealMatrix(result.getSecond(), false));
					}

					public RealVector computeValue(double[] params)
					{
						return new ArrayRealVector(fun.computeValues(params), false);
					}

					public RealMatrix computeJacobian(double[] params)
					{
						return new Array2DRowRealMatrix(fun.computeJacobian(params), false);
					}
				});
			}
			else
			{
				// Compute separately
				builder.model(new MultivariateVectorFunctionWrapper((NonLinearFunction) f, a, n),
						new MultivariateMatrixFunctionWrapper((NonLinearFunction) f, a, n));
			}

			LeastSquaresProblem problem = builder.build();

			Optimum optimum = optimizer.optimize(problem);

			final double[] parameters = optimum.getPoint().toArray();
			setSolution(a, parameters);
			iterations = optimum.getIterations();
			evaluations = optimum.getEvaluations();
			if (a_dev != null)
			{
				try
				{
					double[][] covar = optimum.getCovariances(threshold).getData();
					setDeviationsFromMatrix(a_dev, covar);
				}
				catch (SingularMatrixException e)
				{
					// Matrix inversion failed. In order to return a solution 
					// return the reciprocal of the diagonal of the Fisher information 
					// for a loose bound on the limit 
					final int[] gradientIndices = f.gradientIndices();
					final int nparams = gradientIndices.length;
					GradientCalculator calculator = GradientCalculatorFactory.newCalculator(nparams);
					double[][] alpha = new double[nparams][nparams];
					double[] beta = new double[nparams];
					calculator.findLinearised(nparams, y, a, alpha, beta, (NonLinearFunction) f);

					FisherInformationMatrix m = new FisherInformationMatrix(alpha);
					setDeviations(a_dev, m.crlb(true));
				}
			}
			// Compute function value
			if (y_fit != null)
			{
				Gaussian2DFunction f = (Gaussian2DFunction) this.f;
				f.initialise0(a);
				f.forEach(new ValueProcedure()
				{
					int i = 0;

					public void execute(double value)
					{
						y_fit[i] = value;
					}
				});
			}

			// As this is unweighted then we can do this to get the sum of squared residuals
			// This is the same as optimum.getCost() * optimum.getCost(); The getCost() function
			// just computes the dot product anyway.
			value = optimum.getResiduals().dotProduct(optimum.getResiduals());
		}
		catch (TooManyEvaluationsException e)
		{
			return FitStatus.TOO_MANY_EVALUATIONS;
		}
		catch (TooManyIterationsException e)
		{
			return FitStatus.TOO_MANY_ITERATIONS;
		}
		catch (ConvergenceException e)
		{
			// Occurs when QR decomposition fails - mark as a singular non-linear model (no solution)
			return FitStatus.SINGULAR_NON_LINEAR_MODEL;
		}
		catch (Exception e)
		{
			// TODO - Find out the other exceptions from the fitter and add return values to match. 
			return FitStatus.UNKNOWN;
		}

		return FitStatus.OK;
	}

	@Override
	public boolean computeValue(double[] y, double[] y_fit, double[] a)
	{
		final int nparams = f.gradientIndices().length;

		GradientCalculator calculator = GradientCalculatorFactory.newCalculator(nparams, false);

		// Since we know the function is a Gaussian2DFunction
		value = calculator.findLinearised(y.length, y, y_fit, a, (NonLinearFunction) f);

		return true;
	}
}
