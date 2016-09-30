package gdsc.smlm.fitting.nonlinear;

import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.exception.TooManyIterationsException;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer.Optimum;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.util.Precision;

import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator;
import gdsc.smlm.fitting.nonlinear.gradient.GradientCalculatorFactory;
import gdsc.smlm.function.MultivariateMatrixFunctionWrapper;
import gdsc.smlm.function.MultivariateVectorFunctionWrapper;
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
public class ApacheLVMFitter extends BaseFunctionSolver
{
	/**
	 * Default constructor
	 */
	public ApacheLVMFitter(Gaussian2DFunction gf)
	{
		super(gf);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.BaseFunctionSolver#computeFit(int, double[], double[], double[], double[],
	 * double[], double)
	 */
	public FitStatus computeFit(int n, double[] y, double[] y_fit, double[] a, double[] a_dev, double[] error,
			double noise)
	{
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
			LeastSquaresProblem problem = new LeastSquaresBuilder()
					.maxEvaluations(Integer.MAX_VALUE)
					.maxIterations(getMaxEvaluations())
					.start(initialSolution)
					.target(yd)
					.weight(new DiagonalMatrix(w))
					.model(
						new MultivariateVectorFunctionWrapper(f, a, n), 
						new MultivariateMatrixFunctionWrapper(f, a, n))
					.build();
			//@formatter:on

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
					setDeviations(a_dev, covar);
				}
				catch (SingularMatrixException e)
				{
					// Matrix inversion failed. In order to return a solution assume 
					// the fit achieves the Cramer Roa lower bounds and so the covariance can 
					// be obtained from the Fisher Information Matrix. 
					final int[] gradientIndices = f.gradientIndices();
					final int nparams = gradientIndices.length;
					GradientCalculator calculator = GradientCalculatorFactory.newCalculator(nparams);
					final double[] I = calculator.fisherInformationDiagonal(n, a, f);
					for (int i = nparams; i-- > 0;)
						a_dev[gradientIndices[i]] = 1.0 / Math.sqrt(I[i]);
				}
			}
			// Compute sum-of-squares
			if (y_fit != null)
			{
				residualSumOfSquares = 0;
				f.initialise(a);
				for (int i = 0; i < n; i++)
				{
					y_fit[i] = f.eval(i);
					final double dy = y[i] - y_fit[i];
					residualSumOfSquares += dy * dy;
				}
			}
			else
			{
				// As this is unweighted then we can do this
				residualSumOfSquares = optimum.getResiduals().dotProduct(optimum.getResiduals());
			}

			value = residualSumOfSquares;
			error[0] = getError(residualSumOfSquares, noise, n, initialSolution.length);
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
}
