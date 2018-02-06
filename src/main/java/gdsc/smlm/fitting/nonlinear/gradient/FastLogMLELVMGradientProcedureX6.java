package gdsc.smlm.fitting.nonlinear.gradient;

import gdsc.smlm.function.FastLog;
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
 * This procedure computes a modified Chi-squared expression to perform Maximum Likelihood Estimation assuming Poisson
 * model. See Laurence & Chromy (2010) Efficient maximum likelihood estimator. Nature Methods 7, 338-339. The input data
 * must be Poisson distributed for this to be relevant.
 */
public class FastLogMLELVMGradientProcedureX6 extends FastLogMLELVMGradientProcedureX
{
	/**
	 * Instantiates a new fast log MLELVM gradient procedure X 6.
	 *
	 * @param y
	 *            Data to fit (must be strictly positive)
	 * @param func
	 *            Gradient function
	 * @param fastLog
	 *            the fast log
	 */
	public FastLogMLELVMGradientProcedureX6(final double[] y, final Gradient1Function func, FastLog fastLog)
	{
		super(y, func, fastLog);
		if (n != 6)
			throw new IllegalArgumentException("Function must compute 6 gradients");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient1Procedure#execute(double, double[])
	 */
	public void execute(double fi, double[] dfi_da)
	{
		++yi;
		if (fi > 0.0)
		{
			final double xi = y[yi];

			// We assume y[i] is strictly positive
			value += (fi - xi - xi * fastLog.fastLog(fi / xi));

			final double xi_fi2 = xi / fi / fi;
			final double e = 1 - (xi / fi);

			beta[0] -= e * dfi_da[0];
			beta[1] -= e * dfi_da[1];
			beta[2] -= e * dfi_da[2];
			beta[3] -= e * dfi_da[3];
			beta[4] -= e * dfi_da[4];
			beta[5] -= e * dfi_da[5];

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
			w = dfi_da[5] * xi_fi2;
			alpha[15] += w * dfi_da[0];
			alpha[16] += w * dfi_da[1];
			alpha[17] += w * dfi_da[2];
			alpha[18] += w * dfi_da[3];
			alpha[19] += w * dfi_da[4];
			alpha[20] += w * dfi_da[5];
		}
	}

	@Override
	protected void initialiseGradient()
	{
		GradientProcedureHelper.initialiseWorkingMatrix6(alpha);
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
		GradientProcedureHelper.getMatrix6(this.alpha, alpha);
	}

	@Override
	public void getAlphaLinear(double[] alpha)
	{
		GradientProcedureHelper.getMatrix6(this.alpha, alpha);
	}
}
