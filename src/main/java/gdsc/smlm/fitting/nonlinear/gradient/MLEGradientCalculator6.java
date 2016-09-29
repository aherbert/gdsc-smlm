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
public class MLEGradientCalculator6 extends MLEGradientCalculator
{
	
	/**
	 * Instantiates a new MLE gradient calculator.
	 */
	public MLEGradientCalculator6()
	{
		super(6);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.gradient.MLEGradientCalculator#zero(double[][], double[])
	 */
	@Override
	protected void zero(final double[][] alpha, final double[] beta)
	{
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
		alpha[5][0] = 0;
		alpha[5][1] = 0;
		alpha[5][2] = 0;
		alpha[5][3] = 0;
		alpha[5][4] = 0;
		alpha[5][5] = 0;

		beta[0] = 0;
		beta[1] = 0;
		beta[2] = 0;
		beta[3] = 0;
		beta[4] = 0;
		beta[5] = 0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.gradient.MLEGradientCalculator#compute(double[][], double[], double[], double,
	 * double)
	 */
	@Override
	protected void compute(final double[][] alpha, final double[] beta, final double[] dfi_da, final double fi,
			final double xi)
	{
		final double xi_fi = xi / fi;
		final double xi_fi2 = xi_fi / fi;
		final double e = 1 - (xi_fi);

		alpha[0][0] += dfi_da[0] * xi_fi2 * dfi_da[0];
		double w;
		w = dfi_da[1] * xi_fi2;
		alpha[1][0] += w * dfi_da[0];
		alpha[1][1] += w * dfi_da[1];
		w = dfi_da[2] * xi_fi2;
		alpha[2][0] += w * dfi_da[0];
		alpha[2][1] += w * dfi_da[1];
		alpha[2][2] += w * dfi_da[2];
		w = dfi_da[3] * xi_fi2;
		alpha[3][0] += w * dfi_da[0];
		alpha[3][1] += w * dfi_da[1];
		alpha[3][2] += w * dfi_da[2];
		alpha[3][3] += w * dfi_da[3];
		w = dfi_da[4] * xi_fi2;
		alpha[4][0] += w * dfi_da[0];
		alpha[4][1] += w * dfi_da[1];
		alpha[4][2] += w * dfi_da[2];
		alpha[4][3] += w * dfi_da[3];
		alpha[4][4] += w * dfi_da[4];
		w = dfi_da[5] * xi_fi2;
		alpha[5][0] += w * dfi_da[0];
		alpha[5][1] += w * dfi_da[1];
		alpha[5][2] += w * dfi_da[2];
		alpha[5][3] += w * dfi_da[3];
		alpha[5][4] += w * dfi_da[4];
		alpha[5][5] += w * dfi_da[5];

		beta[0] -= e * dfi_da[0];
		beta[1] -= e * dfi_da[1];
		beta[2] -= e * dfi_da[2];
		beta[3] -= e * dfi_da[3];
		beta[4] -= e * dfi_da[4];
		beta[5] -= e * dfi_da[5];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.gradient.MLEGradientCalculator#symmetric(double[][])
	 */
	@Override
	protected void symmetric(final double[][] alpha)
	{
		alpha[0][1] = alpha[1][0];
		alpha[0][2] = alpha[2][0];
		alpha[0][3] = alpha[3][0];
		alpha[0][4] = alpha[4][0];
		alpha[0][5] = alpha[5][0];
		alpha[1][2] = alpha[2][1];
		alpha[1][3] = alpha[3][1];
		alpha[1][4] = alpha[4][1];
		alpha[1][5] = alpha[5][1];
		alpha[2][3] = alpha[3][2];
		alpha[2][4] = alpha[4][2];
		alpha[2][5] = alpha[5][2];
		alpha[3][4] = alpha[4][3];
		alpha[3][5] = alpha[5][3];
		alpha[4][5] = alpha[5][4];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.fitting.nonlinear.gradient.GradientCalculator#fisherInformationDiagonal(int, double[],
	 * gdsc.smlm.function.NonLinearFunction)
	 */
	@Override
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
			alpha[5] += dy_da[5] * dy_da[5] * yi;
		}

		checkGradients(alpha, nparams);
		return alpha;
	}
}
