package gdsc.smlm.fitting.nonlinear.gradient;

import gdsc.smlm.fitting.function.NonLinearFunction;

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
 * Calculates the Hessian matrix (the square matrix of second-order partial derivatives of a function)
 * and the gradient vector of the function's partial first derivatives with respect to the parameters.
 * This is used within the Levenberg-Marquardt method to fit a nonlinear model with coefficients (a) for a
 * set of data points (x, y).
 */
public class GradientCalculator
{
	public final int nparams;
	private boolean badGradients;

	/**
	 * @param nparams
	 *            The number of gradient parameters
	 */
	public GradientCalculator(int nparams)
	{
		this.nparams = nparams;
	}

	/**
	 * Evaluate the function and compute the sum-of-squares and the curvature matrix.
	 * 
	 * @param x
	 *            n observations
	 * @param y
	 *            Data to fit
	 * @param a
	 *            Set of m coefficients
	 * @param alpha
	 *            the Hessian curvature matrix (size m*m)
	 * @param beta
	 *            the gradient vector of the function's partial first derivatives with respect to the parameters (size
	 *            m)
	 * @param func
	 *            Non-linear fitting function
	 * @return The sum-of-squares value for the fit
	 */
	public double findLinearised(int[] x, float[] y, float[] a, double[][] alpha, double[] beta, NonLinearFunction func)
	{
		double ssx = 0;
		float[] dy_da = new float[a.length];

		for (int i = 0; i < nparams; i++)
		{
			beta[i] = 0;
			for (int j = 0; j <= i; j++)
				alpha[i][j] = 0;
		}

		func.initialise(a);

		if (func.canComputeWeights())
		{
			float[] w = new float[1];
			for (int i = 0; i < x.length; i++)
			{
				final double dy = y[i] - func.eval(x[i], dy_da, w);
				final double weight = getWeight(w[0]);

				// Compute:
				// - the Hessian matrix (the square matrix of second-order partial derivatives of a function; 
				//   that is, it describes the local curvature of a function of many variables.)
				// - the gradient vector of the function's partial first derivatives with respect to the parameters

				for (int j = 0; j < nparams; j++)
				{
					final double wgt = dy_da[j] * weight;

					for (int k = 0; k <= j; k++)
						alpha[j][k] += wgt * dy_da[k];

					beta[j] += wgt * dy;
				}

				ssx += dy * dy * weight;
			}
		}
		else
		{
			for (int i = 0; i < x.length; i++)
			{
				double dy = y[i] - func.eval(x[i], dy_da);

				// Compute:
				// - the Hessian matrix (the square matrix of second-order partial derivatives of a function; 
				//   that is, it describes the local curvature of a function of many variables.)
				// - the gradient vector of the function's partial first derivatives with respect to the parameters

				for (int j = 0; j < nparams; j++)
				{
					float wgt = dy_da[j];

					for (int k = 0; k <= j; k++)
						alpha[j][k] += wgt * dy_da[k];

					beta[j] += wgt * dy;
				}

				ssx += dy * dy;
			}
		}

		// Generate symmetric matrix
		for (int i = 0; i < nparams - 1; i++)
			for (int j = i + 1; j < nparams; j++)
				alpha[i][j] = alpha[j][i];

		return checkGradients(alpha, beta, nparams, ssx);
	}

	/**
	 * Evaluate the function and compute the sum-of-squares and the curvature matrix.
	 * Assumes the n observations (x) are sequential integers from 0.
	 * <p>
	 * If the function supports weights then these will be used to compute the SS and curvature matrix.
	 * 
	 * @param n
	 *            The number of data points
	 * @param y
	 *            Data to fit
	 * @param a
	 *            Set of m coefficients
	 * @param alpha
	 *            the Hessian curvature matrix (size m*m)
	 * @param beta
	 *            the gradient vector of the function's partial first derivatives with respect to the parameters (size
	 *            m)
	 * @param func
	 *            Non-linear fitting function
	 * @see {@link gdsc.smlm.fitting.function.NonLinearFunction#eval(int, float[])},
	 * @see {@link gdsc.smlm.fitting.function.NonLinearFunction#eval(int, float[], float[])},
	 * @see {@link gdsc.smlm.fitting.function.NonLinearFunction#canComputeWeights()}
	 * @return The sum-of-squares value for the fit. GradientCalculator.FAILED if the gradients contain NaN.
	 */
	public double findLinearised(int n, float[] y, float[] a, double[][] alpha, double[] beta, NonLinearFunction func)
	{
		double ssx = 0;
		float[] dy_da = new float[a.length];

		for (int i = 0; i < nparams; i++)
		{
			beta[i] = 0;
			for (int j = 0; j <= i; j++)
				alpha[i][j] = 0;
		}

		func.initialise(a);

		if (func.canComputeWeights())
		{
			float[] w = new float[1];
			for (int i = 0; i < n; i++)
			{
				final double dy = y[i] - func.eval(i, dy_da, w);
				final double weight = getWeight(w[0]);

				// Compute:
				// - the Hessian matrix (the square matrix of second-order partial derivatives of a function; 
				//   that is, it describes the local curvature of a function of many variables.)
				// - the gradient vector of the function's partial first derivatives with respect to the parameters

				for (int j = 0; j < nparams; j++)
				{
					final double wgt = dy_da[j] * weight;

					for (int k = 0; k <= j; k++)
						alpha[j][k] += wgt * dy_da[k];

					beta[j] += wgt * dy;
				}

				ssx += dy * dy * weight;
			}
		}
		else
		{
			for (int i = 0; i < n; i++)
			{
				final double dy = y[i] - func.eval(i, dy_da);

				// Compute:
				// - the Hessian matrix (the square matrix of second-order partial derivatives of a function; 
				//   that is, it describes the local curvature of a function of many variables.)
				// - the gradient vector of the function's partial first derivatives with respect to the parameters

				for (int j = 0; j < nparams; j++)
				{
					final float wgt = dy_da[j];

					for (int k = 0; k <= j; k++)
						alpha[j][k] += wgt * dy_da[k];

					beta[j] += wgt * dy;
				}

				ssx += dy * dy;
			}
		}

		// Generate symmetric matrix
		for (int i = 0; i < nparams - 1; i++)
			for (int j = i + 1; j < nparams; j++)
				alpha[i][j] = alpha[j][i];

		return checkGradients(alpha, beta, nparams, ssx);
	}

	/**
	 * Get the weight factor using the computed weight
	 * <p>
	 * Check if the weight is below 1 and set to 1 to avoid excessive weights.
	 * 
	 * @param w
	 *            The computed weight
	 * @return The weight factor
	 */
	protected double getWeight(float w)
	{
		// TODO - Check if there is a better way to smooth the weights rather than just truncating them at 1
		return (w < 1) ? 1 : 1.0 / w;
	}

	protected double checkGradients(double[][] alpha, double[] beta, int nparams, double ssx)
	{
		badGradients = checkGradients(alpha, beta, nparams);
		return ssx;
	}

	private boolean checkGradients(double[][] alpha, double[] beta, int nparams)
	{
		for (int i = 0; i < nparams; i++)
		{
			if (Double.isNaN(beta[i]))
				return true;
			for (int j = 0; j <= i; j++)
				if (Double.isNaN(alpha[i][j]))
					return true;
		}

		return false;
	}

	/**
	 * @return True if the last calculation produced gradients with NaN values
	 */
	public boolean isNaNGradients()
	{
		return badGradients;
	}
}
