package gdsc.smlm.fitting.nonlinear.gradient;

import gdsc.smlm.function.Gradient1Function;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2017 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Calculates the scaled Hessian matrix (the square matrix of second-order partial derivatives of a function)
 * and the scaled gradient vector of the function's partial first derivatives with respect to the parameters.
 * This is used within the Levenberg-Marquardt method to fit a nonlinear model with coefficients (a) for a
 * set of data points (x, y).
 * <p>
 * This procedure computes a modified Chi-squared expression to perform Weighted Least Squares Estimation assuming a
 * Poisson model with a Gaussian noise component. The weight per observation is equal to 1/[variance + max(y, 0) + 1].
 * <p>
 * See Ruisheng, et al (2017) Algorithmic corrections for localization microscopy with sCMOS cameras - characterisation
 * of a computationally efficient localization approach. Optical Express 25, Issue 10, pp 11701-11716.
 */
public class WLSQLVMGradientProcedure extends LSQLVMGradientProcedure
{
	protected final double[] w;

	/**
	 * Instantiates a new WLSQLVM gradient procedure.
	 *
	 * @param y
	 *            Data to fit
	 * @param var
	 *            the base variance of each observation (must be positive)
	 * @param func
	 *            Gradient function
	 */
	public WLSQLVMGradientProcedure(final double[] y, final double[] var, final Gradient1Function func)
	{
		super(y, func);
		final int n = y.length;
		w = new double[n];
		
		// From Ruisheng, et al (2017):
		// Total noise = variance + max(di, 0) + 1 
		
		if (var != null && var.length == n)
		{
			// Include the variance in the weight. Assume variance is positive.
			for (int i = 0; i < n; i++)
				w[i] = (y[i] > 0) ? 1.0 / (var[i] + y[i] + 1.0) : 1.0 / (var[i] + 1.0);
		}
		else
		{
			for (int i = 0; i < n; i++)
				w[i] = (y[i] > 0) ? 1.0 / (y[i] + 1.0) : 1.0;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient1Procedure#execute(double, double[])
	 */
	public void execute(double value, double[] dy_da)
	{
		final double dy = y[++yi] - value;
		final double w = this.w[yi];
		this.value += dy * dy * w;

		// Compute:
		// - the scaled Hessian matrix (the square matrix of second-order partial derivatives of a function; 
		//   that is, it describes the local curvature of a function of many variables.)
		// - the scaled gradient vector of the function's partial first derivatives with respect to the parameters

		for (int j = 0, i = 0; j < n; j++)
		{
			final double wgt = dy_da[j] * w;

			for (int k = 0; k <= j; k++)
			{
				alpha[i++] += wgt * dy_da[k];
			}
			beta[j] += wgt * dy;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.ValueProcedure#execute(double)
	 */
	public void execute(double value)
	{
		// Produce a sum-of-squares
		final double dy = y[++yi] - value;
		this.value += dy * dy * w[yi];
	}
}
