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
public class LSQGradientProcedure6 extends LSQGradientProcedure
{
	/**
	 * @param y
	 *            Data to fit
	 * @param func
	 *            Gradient function
	 */
	public LSQGradientProcedure6(final double[] y, final Gradient1Function func)
	{
		super(y, func);
		if (n != 6)
			throw new IllegalArgumentException("Function must compute 6 gradients");
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
		alpha[1] += dy_da[1] * dy_da[0];
		alpha[2] += dy_da[1] * dy_da[1];
		alpha[3] += dy_da[2] * dy_da[0];
		alpha[4] += dy_da[2] * dy_da[1];
		alpha[5] += dy_da[2] * dy_da[2];
		alpha[6] += dy_da[3] * dy_da[0];
		alpha[7] += dy_da[3] * dy_da[1];
		alpha[8] += dy_da[3] * dy_da[2];
		alpha[9] += dy_da[3] * dy_da[3];
		alpha[10] += dy_da[4] * dy_da[0];
		alpha[11] += dy_da[4] * dy_da[1];
		alpha[12] += dy_da[4] * dy_da[2];
		alpha[13] += dy_da[4] * dy_da[3];
		alpha[14] += dy_da[4] * dy_da[4];
		alpha[15] += dy_da[5] * dy_da[0];
		alpha[16] += dy_da[5] * dy_da[1];
		alpha[17] += dy_da[5] * dy_da[2];
		alpha[18] += dy_da[5] * dy_da[3];
		alpha[19] += dy_da[5] * dy_da[4];
		alpha[20] += dy_da[5] * dy_da[5];
		
		beta[0] += dy_da[0] * dy;
		beta[1] += dy_da[1] * dy;
		beta[2] += dy_da[2] * dy;
		beta[3] += dy_da[3] * dy;
		beta[4] += dy_da[4] * dy;
		beta[5] += dy_da[5] * dy;

		ssx += dy * dy;
	}

	protected void initialise()
	{
		alpha[0] = 0;
		alpha[1] = 0;
		alpha[2] = 0;
		alpha[3] = 0;
		alpha[4] = 0;
		alpha[5] = 0;
		alpha[6] = 0;
		alpha[7] = 0;
		alpha[8] = 0;
		alpha[9] = 0;
		alpha[10] = 0;
		alpha[11] = 0;
		alpha[12] = 0;
		alpha[13] = 0;
		alpha[14] = 0;
		alpha[15] = 0;
		alpha[16] = 0;
		alpha[17] = 0;
		alpha[18] = 0;
		alpha[19] = 0;
		alpha[20] = 0;
		beta[0] = 0;
		beta[1] = 0;
		beta[2] = 0;
		beta[3] = 0;
		beta[4] = 0;
		beta[5] = 0;
	}

	@Override
	public void getAlphaMatrix(double[][] alpha)
	{
		// Generate symmetric matrix
		alpha[0][0] = this.alpha[0];
		alpha[1][0] = this.alpha[1];
		alpha[0][1] = this.alpha[1];
		alpha[1][1] = this.alpha[2];
		alpha[2][0] = this.alpha[3];
		alpha[0][2] = this.alpha[3];
		alpha[2][1] = this.alpha[4];
		alpha[1][2] = this.alpha[4];
		alpha[2][2] = this.alpha[5];
		alpha[3][0] = this.alpha[6];
		alpha[0][3] = this.alpha[6];
		alpha[3][1] = this.alpha[7];
		alpha[1][3] = this.alpha[7];
		alpha[3][2] = this.alpha[8];
		alpha[2][3] = this.alpha[8];
		alpha[3][3] = this.alpha[9];
		alpha[4][0] = this.alpha[10];
		alpha[0][4] = this.alpha[10];
		alpha[4][1] = this.alpha[11];
		alpha[1][4] = this.alpha[11];
		alpha[4][2] = this.alpha[12];
		alpha[2][4] = this.alpha[12];
		alpha[4][3] = this.alpha[13];
		alpha[3][4] = this.alpha[13];
		alpha[4][4] = this.alpha[14];
		alpha[5][0] = this.alpha[15];
		alpha[0][5] = this.alpha[15];
		alpha[5][1] = this.alpha[16];
		alpha[1][5] = this.alpha[16];
		alpha[5][2] = this.alpha[17];
		alpha[2][5] = this.alpha[17];
		alpha[5][3] = this.alpha[18];
		alpha[3][5] = this.alpha[18];
		alpha[5][4] = this.alpha[19];
		alpha[4][5] = this.alpha[19];
		alpha[5][5] = this.alpha[20];
	}

	@Override
	public void getAlphaLinear(double[] alpha)
	{
		// Generate symmetric matrix
		alpha[0] = this.alpha[0];
		alpha[6] = this.alpha[1];
		alpha[1] = this.alpha[1];
		alpha[7] = this.alpha[2];
		alpha[12] = this.alpha[3];
		alpha[2] = this.alpha[3];
		alpha[13] = this.alpha[4];
		alpha[8] = this.alpha[4];
		alpha[14] = this.alpha[5];
		alpha[18] = this.alpha[6];
		alpha[3] = this.alpha[6];
		alpha[19] = this.alpha[7];
		alpha[9] = this.alpha[7];
		alpha[20] = this.alpha[8];
		alpha[15] = this.alpha[8];
		alpha[21] = this.alpha[9];
		alpha[24] = this.alpha[10];
		alpha[4] = this.alpha[10];
		alpha[25] = this.alpha[11];
		alpha[10] = this.alpha[11];
		alpha[26] = this.alpha[12];
		alpha[16] = this.alpha[12];
		alpha[27] = this.alpha[13];
		alpha[22] = this.alpha[13];
		alpha[28] = this.alpha[14];
		alpha[30] = this.alpha[15];
		alpha[5] = this.alpha[15];
		alpha[31] = this.alpha[16];
		alpha[11] = this.alpha[16];
		alpha[32] = this.alpha[17];
		alpha[17] = this.alpha[17];
		alpha[33] = this.alpha[18];
		alpha[23] = this.alpha[18];
		alpha[34] = this.alpha[19];
		alpha[29] = this.alpha[19];
		alpha[35] = this.alpha[20];
	}
}
