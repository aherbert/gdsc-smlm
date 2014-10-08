package gdsc.smlm.fitting.nonlinear;

import gdsc.smlm.fitting.FitStatus;
import gdsc.smlm.fitting.function.NonLinearFunction;

import java.util.Arrays;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.PointValuePair;
import org.apache.commons.math3.optimization.direct.PowellOptimizer;

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
 * 
 * This is an adaption of the C-code contained in the CcpNmr Analysis Program:
 *   CCPN website (http://www.ccpn.ac.uk/). 
 * The CCPN code was based on Numerical Recipes. 
 *---------------------------------------------------------------------------*/

/**
 * Uses Maximum Likelihood Estimation (MLE) to fit a nonlinear model with coefficients (a) for a
 * set of data points (x, y).
 * <p>
 * The probability mass function for observed value k is modelled as a Poisson process:<br/>
 * pmf = e^-k.(l^k / k!) <br/>
 * where: <br/>
 * k = Observed number of occurrences <br/>
 * l = Expected number of occurrences (the mean)
 * <p>
 * MLE = Max [ sum (ln(e^-k.(l^k / k!)) ] <br/>
 * = Max [ sum (k.ln(l) - l) ]
 * <p>
 * The expected number of occurrences can be modelled using any parameterised function, for example the Gaussian 2D
 * function.
 */
public class MaximumLikelihoodFitter extends BaseFunctionSolver
{
	/**
	 * Wrapper for any function to allow use of the Apache Commons optimiser for Maximum Likelihood Estimation.
	 * <p>
	 * Uses the deprecated API since the new API for version 4.0 is not a fully documented final release.
	 */
	public class PoissonLikelihoodFunction implements MultivariateFunction
	{
		private NonLinearFunction f;
		private float[] a, y;
		int n;

		/**
		 * @param f
		 *            The function to be used to calculated the expected values
		 * @param a
		 *            The initial parameters for the gaussian function
		 * @param y
		 *            The observed values
		 * @param n
		 *            The number of observed values
		 */
		public PoissonLikelihoodFunction(NonLinearFunction f, float[] a, float[] y, int n)
		{
			this.f = f;
			this.a = Arrays.copyOf(a, a.length);
			this.y = y;
			this.n = n;
		}

		/**
		 * Copy the variables into the appropriate parameter positions for the NonLinearFunction
		 * @param variables
		 */
		private void initialiseFunction(double[] variables)
		{
			int[] gradientIndices = f.gradientIndices();
			for (int i = 0; i < gradientIndices.length; i++)
				a[gradientIndices[i]] = (float) variables[i];
			f.initialise(a);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see org.apache.commons.math3.analysis.MultivariateFunction#value(double[])
		 */
		public double value(double[] variables)
		{
			initialiseFunction(variables);

			// Compute the negative log-likelihood to be minimised
			double ll = 0;
			for (int i = 0; i < n; i++)
			{
				final double l = f.eval(i);
				final double k = y[i];
				ll += l - k * Math.log(l);
			}
			return ll;
		}
	}

	/**
	 * Default constructor
	 */
	public MaximumLikelihoodFitter(NonLinearFunction gf)
	{
		super(gf);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.FunctionSolver#fit(int, float[], float[], float[], float[], double[], double)
	 */
	@SuppressWarnings("deprecation")
	public FitStatus fit(int n, float[] y, float[] y_fit, float[] a, float[] a_dev, double[] error, double noise)
	{
		numberOfFittedPoints = n;

		try
		{
			// Non-differentiable version using Powell Optimiser
			double rel = 1e-6; // Relative threshold.
			double abs = 1e-16; // Absolute threshold.
			PowellOptimizer optimizer = new PowellOptimizer(rel, abs);
			double[] startPoint = getInitialSolution(a);

			PointValuePair optimum = optimizer.optimize(getMaxEvaluations(), new PoissonLikelihoodFunction(f, a, y, n),
					GoalType.MINIMIZE, startPoint);

			// TODO - differentiable version using NonLinearConjugateGradientOptimiser
			
			
			
			

			setSolution(a, optimum.getPoint());
			iterations = optimizer.getEvaluations();

			// Compute residuals for the FunctionSolver interface
			if (y_fit == null || y_fit.length < n)
				y_fit = new float[n];
			f.initialise(a);
			residualSumOfSquares = 0;
			for (int i = 0; i < n; i++)
			{
				y_fit[i] = f.eval(i);
				final double residual = y[i] - y_fit[i];
				residualSumOfSquares += residual * residual;
			}

			if (a_dev != null)
			{
				// Covariance not supported for this fitter
			}

			error[0] = NonLinearFit.getError(residualSumOfSquares, noise, n, f.gradientIndices().length);
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
