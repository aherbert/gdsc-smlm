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
 */
public class GradientCalculator5 extends GradientCalculator
{
	public GradientCalculator5()
	{
		super(5);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.fitting.model.GradientCalculator#findLinearised(int[], double[] double[], double[][], double[],
	 * gdsc.fitting.function.NonLinearFunction)
	 */
	public double findLinearised(int[] x, double[] y, double[] a, double[][] alpha, double[] beta, NonLinearFunction func)
	{
		double ssx = 0;
		final double[] dy_da = new double[5];

		alpha[0][0] = 0;
		alpha[1][0] = 0;
		alpha[1][1] = 0;
		alpha[2][0] = 0;
		alpha[2][1] = 0;
		alpha[2][2] = 0;
		alpha[3][0] = 0;
		alpha[3][1] = 0;
		alpha[3][2] = 0;
		alpha[3][3] = 0;
		alpha[4][0] = 0;
		alpha[4][1] = 0;
		alpha[4][2] = 0;
		alpha[4][3] = 0;
		alpha[4][4] = 0;
		
		beta[0] = 0;
		beta[1] = 0;
		beta[2] = 0;
		beta[3] = 0;
		beta[4] = 0;

		func.initialise(a);

		if (func.canComputeWeights())
		{
			double[] w = new double[1];
			for (int i = 0; i < x.length; i++)
			{
				final double dy = y[i] - func.eval(x[i], dy_da, w);
				final double weight = getWeight(w[0]);

				alpha[0][0] += dy_da[0] * weight * dy_da[0];
				alpha[1][0] += dy_da[1] * weight * dy_da[0];
				alpha[1][1] += dy_da[1] * weight * dy_da[1];
				alpha[2][0] += dy_da[2] * weight * dy_da[0];
				alpha[2][1] += dy_da[2] * weight * dy_da[1];
				alpha[2][2] += dy_da[2] * weight * dy_da[2];
				alpha[3][0] += dy_da[3] * weight * dy_da[0];
				alpha[3][1] += dy_da[3] * weight * dy_da[1];
				alpha[3][2] += dy_da[3] * weight * dy_da[2];
				alpha[3][3] += dy_da[3] * weight * dy_da[3];
				alpha[4][0] += dy_da[4] * weight * dy_da[0];
				alpha[4][1] += dy_da[4] * weight * dy_da[1];
				alpha[4][2] += dy_da[4] * weight * dy_da[2];
				alpha[4][3] += dy_da[4] * weight * dy_da[3];
				alpha[4][4] += dy_da[4] * weight * dy_da[4];

				beta[0] += dy_da[0] * weight * dy;
				beta[1] += dy_da[1] * weight * dy;
				beta[2] += dy_da[2] * weight * dy;
				beta[3] += dy_da[3] * weight * dy;
				beta[4] += dy_da[4] * weight * dy;

				ssx += dy * dy * weight;
			}
		}
		else
		{
			for (int i = 0; i < x.length; i++)
			{
				final double dy = y[i] - func.eval(x[i], dy_da);

				alpha[0][0] += dy_da[0] * dy_da[0];
				alpha[1][0] += dy_da[1] * dy_da[0];
				alpha[1][1] += dy_da[1] * dy_da[1];
				alpha[2][0] += dy_da[2] * dy_da[0];
				alpha[2][1] += dy_da[2] * dy_da[1];
				alpha[2][2] += dy_da[2] * dy_da[2];
				alpha[3][0] += dy_da[3] * dy_da[0];
				alpha[3][1] += dy_da[3] * dy_da[1];
				alpha[3][2] += dy_da[3] * dy_da[2];
				alpha[3][3] += dy_da[3] * dy_da[3];
				alpha[4][0] += dy_da[4] * dy_da[0];
				alpha[4][1] += dy_da[4] * dy_da[1];
				alpha[4][2] += dy_da[4] * dy_da[2];
				alpha[4][3] += dy_da[4] * dy_da[3];
				alpha[4][4] += dy_da[4] * dy_da[4];

				beta[0] += dy_da[0] * dy;
				beta[1] += dy_da[1] * dy;
				beta[2] += dy_da[2] * dy;
				beta[3] += dy_da[3] * dy;
				beta[4] += dy_da[4] * dy;

				ssx += dy * dy;
			}
		}

		// Generate symmetric matrix
		alpha[0][1] = alpha[1][0];
		alpha[0][2] = alpha[2][0];
		alpha[0][3] = alpha[3][0];
		alpha[0][4] = alpha[4][0];
		alpha[1][2] = alpha[2][1];
		alpha[1][3] = alpha[3][1];
		alpha[1][4] = alpha[4][1];
		alpha[2][3] = alpha[3][2];
		alpha[2][4] = alpha[4][2];
		alpha[3][4] = alpha[4][3];

		return checkGradients(alpha, beta, nparams, ssx);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.fitting.nonlinear.gradient.GradientCalculator#findLinearised(int, double[] double[], double[][],
	 * double[], gdsc.fitting.function.NonLinearFunction)
	 */
	public double findLinearised(int n, double[] y, double[] a, double[][] alpha, double[] beta, NonLinearFunction func)
	{
		double ssx = 0;
		final double[] dy_da = new double[5];

		alpha[0][0] = 0;
		alpha[1][0] = 0;
		alpha[1][1] = 0;
		alpha[2][0] = 0;
		alpha[2][1] = 0;
		alpha[2][2] = 0;
		alpha[3][0] = 0;
		alpha[3][1] = 0;
		alpha[3][2] = 0;
		alpha[3][3] = 0;
		alpha[4][0] = 0;
		alpha[4][1] = 0;
		alpha[4][2] = 0;
		alpha[4][3] = 0;
		alpha[4][4] = 0;
		
		beta[0] = 0;
		beta[1] = 0;
		beta[2] = 0;
		beta[3] = 0;
		beta[4] = 0;

		func.initialise(a);

		if (func.canComputeWeights())
		{
			double[] w = new double[1];
			for (int i = 0; i < n; i++)
			{
				final double dy = y[i] - func.eval(i, dy_da, w);
				final double weight = getWeight(w[0]);

				alpha[0][0] += dy_da[0] * weight * dy_da[0];
				alpha[1][0] += dy_da[1] * weight * dy_da[0];
				alpha[1][1] += dy_da[1] * weight * dy_da[1];
				alpha[2][0] += dy_da[2] * weight * dy_da[0];
				alpha[2][1] += dy_da[2] * weight * dy_da[1];
				alpha[2][2] += dy_da[2] * weight * dy_da[2];
				alpha[3][0] += dy_da[3] * weight * dy_da[0];
				alpha[3][1] += dy_da[3] * weight * dy_da[1];
				alpha[3][2] += dy_da[3] * weight * dy_da[2];
				alpha[3][3] += dy_da[3] * weight * dy_da[3];
				alpha[4][0] += dy_da[4] * weight * dy_da[0];
				alpha[4][1] += dy_da[4] * weight * dy_da[1];
				alpha[4][2] += dy_da[4] * weight * dy_da[2];
				alpha[4][3] += dy_da[4] * weight * dy_da[3];
				alpha[4][4] += dy_da[4] * weight * dy_da[4];

				beta[0] += dy_da[0] * weight * dy;
				beta[1] += dy_da[1] * weight * dy;
				beta[2] += dy_da[2] * weight * dy;
				beta[3] += dy_da[3] * weight * dy;
				beta[4] += dy_da[4] * weight * dy;

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
				alpha[3][0] += dy_da[3] * dy_da[0];
				alpha[3][1] += dy_da[3] * dy_da[1];
				alpha[3][2] += dy_da[3] * dy_da[2];
				alpha[3][3] += dy_da[3] * dy_da[3];
				alpha[4][0] += dy_da[4] * dy_da[0];
				alpha[4][1] += dy_da[4] * dy_da[1];
				alpha[4][2] += dy_da[4] * dy_da[2];
				alpha[4][3] += dy_da[4] * dy_da[3];
				alpha[4][4] += dy_da[4] * dy_da[4];

				beta[0] += dy_da[0] * dy;
				beta[1] += dy_da[1] * dy;
				beta[2] += dy_da[2] * dy;
				beta[3] += dy_da[3] * dy;
				beta[4] += dy_da[4] * dy;

				ssx += dy * dy;
			}
		}

		// Generate symmetric matrix
		alpha[0][1] = alpha[1][0];
		alpha[0][2] = alpha[2][0];
		alpha[0][3] = alpha[3][0];
		alpha[0][4] = alpha[4][0];
		alpha[1][2] = alpha[2][1];
		alpha[1][3] = alpha[3][1];
		alpha[1][4] = alpha[4][1];
		alpha[2][3] = alpha[3][2];
		alpha[2][4] = alpha[4][2];
		alpha[3][4] = alpha[4][3];

		return checkGradients(alpha, beta, nparams, ssx);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator#fisherInformationDiagonal(int, double[],
	 * gdsc.smlm.function.NonLinearFunction)
	 */
	public double[] fisherInformationDiagonal(final int n, final double[] a, final NonLinearFunction func)
	{
		final double[] dy_da = new double[a.length];

		final double[] alpha = new double[nparams];

		func.initialise(a);

		for (int i = 0; i < n; i++)
		{
			final double yi = 1.0 / func.eval(i, dy_da);
			alpha[0] += dy_da[0] * dy_da[0] * yi;
			alpha[1] += dy_da[1] * dy_da[1] * yi;
			alpha[2] += dy_da[2] * dy_da[2] * yi;
			alpha[3] += dy_da[3] * dy_da[3] * yi;
			alpha[4] += dy_da[4] * dy_da[4] * yi;
		}

		checkGradients(alpha, nparams);
		return alpha;
	}
}
