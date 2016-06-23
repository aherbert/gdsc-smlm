package gdsc.smlm.fitting.nonlinear;

import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.function.MultivariateMatrixFunctionWrapper;
import gdsc.smlm.function.MultivariateVectorFunctionWrapper;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;

import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.exception.TooManyIterationsException;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.MaxIter;
import org.apache.commons.math3.optim.PointVectorValuePair;
import org.apache.commons.math3.optim.nonlinear.vector.ModelFunction;
import org.apache.commons.math3.optim.nonlinear.vector.ModelFunctionJacobian;
import org.apache.commons.math3.optim.nonlinear.vector.Target;
import org.apache.commons.math3.optim.nonlinear.vector.Weight;
import org.apache.commons.math3.optim.nonlinear.vector.jacobian.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.util.Precision;

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
			PointVectorValuePair optimum = optimizer.optimize(new MaxIter(getMaxEvaluations()),
					new MaxEval(Integer.MAX_VALUE),
					new ModelFunctionJacobian(new MultivariateMatrixFunctionWrapper(f, a, n)),
					new ModelFunction(new MultivariateVectorFunctionWrapper(f, a, n)), new Target(yd), new Weight(w),
					new InitialGuess(initialSolution));

			final double[] parameters = optimum.getPointRef();
			setSolution(a, parameters);
			iterations = optimizer.getIterations();
			evaluations = optimizer.getEvaluations();
			if (a_dev != null)
			{
				double[][] covar = optimizer.computeCovariances(parameters, threshold);
				setDeviations(a_dev, covar);
			}
			// Compute sum-of-squares
			if (y_fit != null)
			{
				final double[] optimumValue = optimum.getValue();
				System.arraycopy(optimumValue, 0, y_fit, 0, n);
				//f.initialise(a);
				//for (int i = 0; i < n; i++)
				//{
				//	y_fit[i] = f.eval(i);
				//}
			}

			value = residualSumOfSquares = optimizer.getChiSquare();
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
