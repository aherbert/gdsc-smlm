package gdsc.smlm.fitting.nonlinear;

import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.function.DifferentiableGaussian2DFunction;
import gdsc.smlm.fitting.function.Gaussian2DFunction;

import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optimization.PointVectorValuePair;
import org.apache.commons.math3.optimization.general.LevenbergMarquardtOptimizer;
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
	 * @see gdsc.smlm.fitting.FunctionSolver#fit(int, float[], float[], float[], float[], double[], double)
	 */
	public FitStatus fit(int n, float[] y, float[] y_fit, float[] a, float[] a_dev, double[] error, double noise)
	{
		numberOfFittedPoints = n;

		DifferentiableGaussian2DFunction func = new DifferentiableGaussian2DFunction((Gaussian2DFunction) f, a);
		func.addData(n, y);

		try
		{
			// Different convergence thresholds seem to have no effect on the resulting fit, only the number of
			// iterations for convergence
			double initialStepBoundFactor = 100;
			double costRelativeTolerance = 1e-10;
			double parRelativeTolerance = 1e-10;
			double orthoTolerance = 1e-10;
			double threshold = Precision.SAFE_MIN;

			// Extract the parameters to be fitted
			double[] initialSolution = getInitialSolution(a);

			// TODO - Pass in stopping criteria.
			// Do this if the fitting seems to be successful

			LevenbergMarquardtOptimizer optimizer = new LevenbergMarquardtOptimizer(initialStepBoundFactor,
					costRelativeTolerance, parRelativeTolerance, orthoTolerance, threshold);
			PointVectorValuePair optimum = optimizer.optimize(getMaxEvaluations(), func, func.getY(),
					func.getWeights(), initialSolution);

			double[] parameters = optimum.getPoint();
			setSolution(a, parameters);
			iterations = optimizer.getEvaluations();
			if (a_dev != null)
			{
				double[][] covar = optimizer.getCovariances();
				setDeviations(a_dev, covar);
			}
			// Compute fitted function if desired 
			if (y_fit != null)
			{
				f.initialise(a);
				for (int i = 0; i < n; i++)
					y_fit[i] = f.eval(i);
			}

			residualSumOfSquares = error[0] = optimizer.getChiSquare();
			totalSumOfSquares = getSumOfSquares(n, y);
		}
		catch (TooManyEvaluationsException e)
		{
			return FitStatus.FAILED_TO_CONVERGE;
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
