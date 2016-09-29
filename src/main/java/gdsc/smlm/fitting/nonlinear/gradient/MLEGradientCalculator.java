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
	 * Note: if the function returns a negative value then it is set to zero
	 * 
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
		final double[] dfi_da = new double[nparams];

		zero(alpha, beta);

		func.initialise(a);

		for (int i = 0; i < x.length; i++)
		{
			// Function must produce a positive output.
			// The code provided in Laurence & Chromy (2010) Nature Methods 7, 338-339, SI
			// effectively ignores the any function value below zero. This could lead to a 
			// situation where the best chisq value can be achieved by setting the output
			// function to produce 0 for all evaluations. To cope with this the input 
			// function should be bounded to always produce a positive number. For now
			// we just set it to Double.MIN_VALUE
			final double fi = Math.max(func.eval(x[i], dfi_da), Double.MIN_VALUE);
			final double xi = y[i];

			// We assume y[i] is positive
			if (xi == 0)
				chisq += 2 * fi;
			else
				chisq += 2 * (fi - xi - xi * Math.log(fi / xi));

			compute(alpha, beta, dfi_da, fi, xi);
		}

		symmetric(alpha);

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
		final double[] dfi_da = new double[nparams];

		zero(alpha, beta);

		func.initialise(a);

		for (int i = 0; i < n; i++)
		{
			// Function must produce a positive output.
			// The code provided in Laurence & Chromy (2010) Nature Methods 7, 338-339, SI
			// effectively ignores the any function value below zero. This could lead to a 
			// situation where the best chisq value can be achieved by setting the output
			// function to produce 0 for all evaluations. To cope with this the input 
			// function should be bounded to always produce a positive number. For now
			// we just set it to Double.MIN_VALUE
			final double fi = Math.max(func.eval(i, dfi_da), Double.MIN_VALUE);
			final double xi = y[i];

			// We assume y[i] is positive
			if (xi == 0)
				chisq += 2 * fi;
			else
				chisq += 2 * (fi - xi - xi * Math.log(fi / xi));

			compute(alpha, beta, dfi_da, fi, xi);
		}

		symmetric(alpha);

		return checkGradients(alpha, beta, nparams, chisq);
	}

	/**
	 * Zero the working region of the input matrix alpha and vector beta
	 *
	 * @param alpha
	 *            the alpha
	 * @param beta
	 *            the beta
	 */
	protected void zero(final double[][] alpha, final double[] beta)
	{
		for (int i = 0; i < nparams; i++)
		{
			beta[i] = 0;
			for (int j = 0; j <= i; j++)
				alpha[i][j] = 0;
		}
	}

	/**
	 * Compute the matrix alpha and vector beta
	 *
	 * @param alpha
	 *            the alpha
	 * @param beta
	 *            the beta
	 * @param dfi_da
	 *            the gradient of the function with respect to each parameter a
	 * @param fi
	 *            the function value at index i
	 * @param xi
	 *            the data value at index i
	 */
	protected void compute(final double[][] alpha, final double[] beta, final double[] dfi_da, final double fi,
			final double xi)
	{
		final double xi_fi = xi / fi;
		final double xi_fi2 = xi_fi / fi;
		final double e = 1 - (xi_fi);

		// Compute:
		// Laurence & Chromy (2010) Nature Methods 7, 338-339, SI
		// alpha - the Hessian matrix (the square matrix of second-order partial derivatives of a function; 
		//         that is, it describes the local curvature of a function of many variables.)
		// beta  - the gradient vector of the function's partial first derivatives with respect to the parameters

		for (int k = 0; k < nparams; k++)
		{
			final double w = dfi_da[k] * xi_fi2;

			for (int l = 0; l <= k; l++)
				// This is the non-optimised version:
				//alpha[j][k] += dfi_da[k] * dfi_da[l] * xi / (fi * fi);
				alpha[k][l] += w * dfi_da[l];

			// This is the non-optimised version:
			//beta[j] -= (1 - xi / fi) * dfi_da[k];
			beta[k] -= e * dfi_da[k];
		}
	}

	/**
	 * Generate a symmetric matrix alpha
	 *
	 * @param alpha
	 *            the alpha
	 */
	protected void symmetric(final double[][] alpha)
	{
		for (int i = 0; i < nparams - 1; i++)
			for (int j = i + 1; j < nparams; j++)
				alpha[i][j] = alpha[j][i];
	}
}
