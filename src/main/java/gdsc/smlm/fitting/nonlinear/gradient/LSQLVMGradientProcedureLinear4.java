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
 * Calculates the Hessian matrix (the square matrix of second-order partial derivatives of a function)
 * and the scaled gradient vector of the function's partial first derivatives with respect to the parameters.
 * This is used within the Levenberg-Marquardt method to fit a nonlinear model with coefficients (a) for a
 * set of data points (x, y).
 * <p>
 * Note that the Hessian matrix is scaled by 1/2 and the gradient vector is scaled by -1/2 for convenience in solving
 * the non-linear model. See Numerical Recipes in C++, 2nd Ed. Equation 15.5.8 for Nonlinear Models.
 */
public class LSQLVMGradientProcedureLinear4 extends LSQLVMGradientProcedureLinear
{
	/**
	 * @param y
	 *            Data to fit
	 * @param func
	 *            Gradient function
	 */
	public LSQLVMGradientProcedureLinear4(final double[] y, final Gradient1Function func)
	{
		super(y, func);
		if (n != 4)
			throw new IllegalArgumentException("Function must compute 4 gradients");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient1Procedure#execute(double, double[])
	 */
	public void execute(double value, double[] dy_da)
	{
		final double dy = y[yi++] - value;

		alpha[0] += dy_da[0] * dy_da[0];
		alpha[1] += dy_da[0] * dy_da[1];
		alpha[2] += dy_da[0] * dy_da[2];
		alpha[3] += dy_da[0] * dy_da[3];
		alpha[5] += dy_da[1] * dy_da[1];
		alpha[6] += dy_da[1] * dy_da[2];
		alpha[7] += dy_da[1] * dy_da[3];
		alpha[10] += dy_da[2] * dy_da[2];
		alpha[11] += dy_da[2] * dy_da[3];
		alpha[15] += dy_da[3] * dy_da[3];

		beta[0] += dy_da[0] * dy;
		beta[1] += dy_da[1] * dy;
		beta[2] += dy_da[2] * dy;
		beta[3] += dy_da[3] * dy;

		this.value += dy * dy;
	}

	protected void initialiseGradient()
	{
		alpha[0] = 0;
		alpha[1] = 0;
		alpha[2] = 0;
		alpha[3] = 0;
		alpha[5] = 0;
		alpha[6] = 0;
		alpha[7] = 0;
		alpha[10] = 0;
		alpha[11] = 0;
		alpha[15] = 0;

		beta[0] = 0;
		beta[1] = 0;
		beta[2] = 0;
		beta[3] = 0;
	}

	protected void finishGradient()
	{
		alpha[4] = alpha[1];
		alpha[8] = alpha[2];
		alpha[12] = alpha[3];
		alpha[9] = alpha[6];
		alpha[13] = alpha[7];
		alpha[14] = alpha[11];
	}
}
