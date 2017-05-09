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
 * This calculator computes a modified Chi-squared expression to perform Maximum Likelihood Estimation assuming Poisson
 * model. See Laurence & Chromy (2010) Efficient maximum likelihood estimator. Nature Methods 7, 338-339. The input data
 * must be Poisson distributed for this to be relevant.
 */
public class MLELVMGradientProcedure5 extends MLELVMGradientProcedure
{
	/**
	 * @param y
	 *            Data to fit (must be positive)
	 * @param func
	 *            Gradient function
	 */
	public MLELVMGradientProcedure5(final double[] y, final Gradient1Function func)
	{
		super(y, func);
		if (n != 5)
			throw new IllegalArgumentException("Function must compute 5 gradients");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient1Procedure#execute(double, double[])
	 */
	public void execute(double fi, double[] dfi_da)
	{
		if (fi > 0)
		{
			final double xi = y[yi++];
			
			// We assume y[i] is positive
			if (xi == 0)
			{
				value += fi;
				beta[0] -= dfi_da[0];
				beta[1] -= dfi_da[1];
				beta[2] -= dfi_da[2];
				beta[3] -= dfi_da[3];
				beta[4] -= dfi_da[4];
			}
			else
			{
				value += (fi - xi - xi * Math.log(fi / xi));
				
				final double xi_fi2 = xi / fi / fi;
				final double e = 1 - (xi / fi);

				beta[0] -= e * dfi_da[0];
				beta[1] -= e * dfi_da[1];
				beta[2] -= e * dfi_da[2];
				beta[3] -= e * dfi_da[3];
				beta[4] -= e * dfi_da[4];
				
				alpha[0] += dfi_da[0] * xi_fi2 * dfi_da[0];
				double w;
				w = dfi_da[1] * xi_fi2;
				alpha[1] += w * dfi_da[0];
				alpha[2] += w * dfi_da[1];
				w = dfi_da[2] * xi_fi2;
				alpha[3] += w * dfi_da[0];
				alpha[4] += w * dfi_da[1];
				alpha[5] += w * dfi_da[2];
				w = dfi_da[3] * xi_fi2;
				alpha[6] += w * dfi_da[0];
				alpha[7] += w * dfi_da[1];
				alpha[8] += w * dfi_da[2];
				alpha[9] += w * dfi_da[3];
				w = dfi_da[4] * xi_fi2;
				alpha[10] += w * dfi_da[0];
				alpha[11] += w * dfi_da[1];
				alpha[12] += w * dfi_da[2];
				alpha[13] += w * dfi_da[3];
				alpha[14] += w * dfi_da[4];
			}
		}
	}	

	protected void initialiseGradient()
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
		beta[0] = 0;
		beta[1] = 0;
		beta[2] = 0;
		beta[3] = 0;
		beta[4] = 0;
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
	}

	@Override
	public void getAlphaLinear(double[] alpha)
	{
		// Generate symmetric matrix
		alpha[0] = this.alpha[0];
		alpha[5] = this.alpha[1];
		alpha[1] = this.alpha[1];
		alpha[6] = this.alpha[2];
		alpha[10] = this.alpha[3];
		alpha[2] = this.alpha[3];
		alpha[11] = this.alpha[4];
		alpha[7] = this.alpha[4];
		alpha[12] = this.alpha[5];
		alpha[15] = this.alpha[6];
		alpha[3] = this.alpha[6];
		alpha[16] = this.alpha[7];
		alpha[8] = this.alpha[7];
		alpha[17] = this.alpha[8];
		alpha[13] = this.alpha[8];
		alpha[18] = this.alpha[9];
		alpha[20] = this.alpha[10];
		alpha[4] = this.alpha[10];
		alpha[21] = this.alpha[11];
		alpha[9] = this.alpha[11];
		alpha[22] = this.alpha[12];
		alpha[14] = this.alpha[12];
		alpha[23] = this.alpha[13];
		alpha[19] = this.alpha[13];
		alpha[24] = this.alpha[14];
	}
}
