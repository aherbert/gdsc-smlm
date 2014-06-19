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
public class GradientCalculator3 extends GradientCalculator
{
	public GradientCalculator3()
	{
		super(3);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.fitting.model.GradientCalculator#findLinearised(int[], float[] float[], double[][], double[],
	 * gdsc.fitting.function.NonLinearFunction)
	 */
	public double findLinearised(int[] x, float[] y, float[] a, double[][] alpha, double[] beta, NonLinearFunction func)
	{
		double ssx = 0;
		float[] dy_da = new float[a.length];

		for (int i = 0; i < beta.length; i++)
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
				final double weight = 1 / w[0];

				alpha[0][0] += dy_da[0] * weight * dy_da[0];
				alpha[1][0] += dy_da[1] * weight * dy_da[0];
				alpha[1][1] += dy_da[1] * weight * dy_da[1];
				alpha[2][0] += dy_da[2] * weight * dy_da[0];
				alpha[2][1] += dy_da[2] * weight * dy_da[1];
				alpha[2][2] += dy_da[2] * weight * dy_da[2];

				beta[0] += dy_da[0] * weight * dy;
				beta[1] += dy_da[1] * weight * dy;
				beta[2] += dy_da[2] * weight * dy;

				ssx += dy * dy * weight;
			}
		}
		else
		{
			for (int i = 0; i < x.length; i++)
			{
				double dy = y[i] - func.eval(x[i], dy_da);

				alpha[0][0] += dy_da[0] * dy_da[0];
				alpha[1][0] += dy_da[1] * dy_da[0];
				alpha[1][1] += dy_da[1] * dy_da[1];
				alpha[2][0] += dy_da[2] * dy_da[0];
				alpha[2][1] += dy_da[2] * dy_da[1];
				alpha[2][2] += dy_da[2] * dy_da[2];

				//    		for (int j = beta.length; j-- > 0; )
				//    			beta[j] += dy_da[j] * dy;

				beta[0] += dy_da[0] * dy;
				beta[1] += dy_da[1] * dy;
				beta[2] += dy_da[2] * dy;

				ssx += dy * dy;
			}
		}

		// Generate symmetric matrix
		for (int i = 0; i < beta.length - 1; i++)
			for (int j = i + 1; j < beta.length; j++)
				alpha[i][j] = alpha[j][i];

		//        alpha[0][1] = alpha[1][0];
		//        alpha[0][2] = alpha[2][0];
		//        alpha[1][2] = alpha[2][1];

		return checkGradients(alpha, beta, nparams, ssx);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.fitting.nonlinear.gradient.GradientCalculator#findLinearised(int, float[] float[], double[][],
	 * double[], gdsc.fitting.function.NonLinearFunction)
	 */
	public double findLinearised(int n, float[] y, float[] a, double[][] alpha, double[] beta, NonLinearFunction func)
	{
		double ssx = 0;
		float[] dy_da = new float[a.length];

		for (int i = 0; i < beta.length; i++)
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
				final double weight = 1 / w[0];

				alpha[0][0] += dy_da[0] * weight * dy_da[0];
				alpha[1][0] += dy_da[1] * weight * dy_da[0];
				alpha[1][1] += dy_da[1] * weight * dy_da[1];
				alpha[2][0] += dy_da[2] * weight * dy_da[0];
				alpha[2][1] += dy_da[2] * weight * dy_da[1];
				alpha[2][2] += dy_da[2] * weight * dy_da[2];

				beta[0] += dy_da[0] * weight * dy;
				beta[1] += dy_da[1] * weight * dy;
				beta[2] += dy_da[2] * weight * dy;

				ssx += dy * dy * weight;
			}
		}
		else
		{
			for (int i = 0; i < n; i++)
			{
				double dy = y[i] - func.eval(i, dy_da);

				alpha[0][0] += dy_da[0] * dy_da[0];
				alpha[1][0] += dy_da[1] * dy_da[0];
				alpha[1][1] += dy_da[1] * dy_da[1];
				alpha[2][0] += dy_da[2] * dy_da[0];
				alpha[2][1] += dy_da[2] * dy_da[1];
				alpha[2][2] += dy_da[2] * dy_da[2];

				//    		for (int j = beta.length; j-- > 0; )
				//    			beta[j] += dy_da[j] * dy;

				beta[0] += dy_da[0] * dy;
				beta[1] += dy_da[1] * dy;
				beta[2] += dy_da[2] * dy;

				ssx += dy * dy;
			}
		}

		// Generate symmetric matrix
		for (int i = 0; i < beta.length - 1; i++)
			for (int j = i + 1; j < beta.length; j++)
				alpha[i][j] = alpha[j][i];

		//        alpha[0][1] = alpha[1][0];
		//        alpha[0][2] = alpha[2][0];
		//        alpha[1][2] = alpha[2][1];

		return checkGradients(alpha, beta, nparams, ssx);
	}
}
