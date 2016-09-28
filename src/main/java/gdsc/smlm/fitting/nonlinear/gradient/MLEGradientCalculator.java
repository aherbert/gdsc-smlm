package gdsc.smlm.fitting.nonlinear.gradient;

import gdsc.smlm.function.NonLinearFunction;

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
 * Calculates the Hessian matrix (the square matrix of second-order partial derivatives of a function)
 * and the gradient vector of the function's partial first derivatives with respect to the parameters.
 * This is used within the Levenberg-Marquardt method to fit a nonlinear model with coefficients (a) for a
 * set of data points (x, y).
 * <p>
 * This calculator computes a modified Chi-squared expression to perform Maximum Likelihood Estimation assuming Poisson
 * model. See Laurence & Chromy (2010) Efficient maximum likelihood estimator. Nature Methods 7, 338-339. The input data
 * must be Poisson distributed for this to be relevant.
 */
public class MLEGradientCalculator extends GradientCalculator
{
	/**
	 * @param nparams
	 *            The number of gradient parameters
	 */
	public MLEGradientCalculator(final int nparams)
	{
		super(nparams);
	}

	/**
	 * @param y
	 *            Data to fit (must be strictly positive Poisson data)
	 * @return The MLE chi-squared value
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator#findLinearised(int[], double[], double[],
	 *      double[][], double[], gdsc.smlm.function.NonLinearFunction)
	 */
	public double findLinearised(final int[] x, final double[] y, final double[] a, final double[][] alpha,
			final double[] beta, final NonLinearFunction func)
	{
		double chisq = 0;
		final double[] dy_da = new double[a.length];

		for (int i = 0; i < nparams; i++)
		{
			beta[i] = 0;
			for (int j = 0; j <= i; j++)
				alpha[i][j] = 0;
		}

		func.initialise(a);

		for (int i = 0; i < x.length; i++)
		{
			// Function must produce a positive output
			double ymod = func.eval(x[i], dy_da);
			final double dy;
			if (ymod <= 0)
			{
				ymod = Double.MIN_VALUE;
				dy = y[i];
			}
			else
			{
				dy = y[i] - ymod;
			}
			final double sig2i = 1.0 / ymod;
			final double y_ymod = y[i] / ymod;

			// Compute:
			// - the Hessian matrix (the square matrix of second-order partial derivatives of a function; 
			//   that is, it describes the local curvature of a function of many variables.)
			// - the gradient vector of the function's partial first derivatives with respect to the parameters

			for (int j = 0; j < nparams; j++)
			{
				final double wt = dy_da[j] * sig2i;

				for (int k = 0; k <= j; k++)
					// This is the non-optimised version:
					//alpha[j][k] += dy_da[j] * dy_da[k] * y[i] / (ymod * ymod);
					alpha[j][k] += y_ymod * wt * dy_da[k];

				// This is the non-optimised version:
				//beta[j] -= (1 - y[i] / ymod) * dy_da[j];
				beta[j] -= (1 - y_ymod) * dy_da[j];
			}

			if (y[i] == 0)
				chisq += 2 * dy;
			else
				chisq += 2 * (dy - y[i] * Math.log(ymod / y[i]));
		}

		// Generate symmetric matrix
		for (int i = 0; i < nparams - 1; i++)
			for (int j = i + 1; j < nparams; j++)
				alpha[i][j] = alpha[j][i];

		return checkGradients(alpha, beta, nparams, chisq);
	}

	/**
	 * @param y
	 *            Data to fit (must be strictly positive Poisson data)
	 * @return The MLE chi-squared value
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator#findLinearised(int, double[], double[], double[][],
	 *      double[], gdsc.smlm.function.NonLinearFunction)
	 */
	public double findLinearised(final int n, final double[] y, final double[] a, final double[][] alpha,
			final double[] beta, final NonLinearFunction func)
	{
		double chisq = 0;
		final double[] dy_da = new double[a.length];

		for (int i = 0; i < nparams; i++)
		{
			beta[i] = 0;
			for (int j = 0; j <= i; j++)
				alpha[i][j] = 0;
		}

		func.initialise(a);

		for (int i = 0; i < n; i++)
		{
			// Function must produce a positive output
			double ymod = func.eval(i, dy_da);
			final double dy;
			if (ymod <= 0)
			{
				ymod = Double.MIN_VALUE;
				dy = y[i];
			}
			else
			{
				dy = y[i] - ymod;
			}
			final double sig2i = 1.0 / ymod;
			final double y_ymod = y[i] / ymod;

			// Compute:
			// - the Hessian matrix (the square matrix of second-order partial derivatives of a function; 
			//   that is, it describes the local curvature of a function of many variables.)
			// - the gradient vector of the function's partial first derivatives with respect to the parameters

			for (int j = 0; j < nparams; j++)
			{
				final double wt = dy_da[j] * sig2i;

				for (int k = 0; k <= j; k++)
					// This is the non-optimised version:
					//alpha[j][k] += dy_da[j] * dy_da[k] * y[i] / (ymod * ymod);
					alpha[j][k] += y_ymod * wt * dy_da[k];

				// This is the non-optimised version:
				//beta[j] -= (1 - y[i] / ymod) * dy_da[j];
				beta[j] -= (1 - y_ymod) * dy_da[j];
			}

			if (y[i] == 0)
				chisq += 2 * dy;
			else
				chisq += 2 * (dy - y[i] * Math.log(ymod / y[i]));
		}

		// Generate symmetric matrix
		for (int i = 0; i < nparams - 1; i++)
			for (int j = i + 1; j < nparams; j++)
				alpha[i][j] = alpha[j][i];

		return checkGradients(alpha, beta, nparams, chisq);
	}
}
