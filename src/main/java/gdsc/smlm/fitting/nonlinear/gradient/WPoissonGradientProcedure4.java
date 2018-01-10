package gdsc.smlm.fitting.nonlinear.gradient;

import gdsc.smlm.function.Gradient1Function;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2018 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Calculates the Fisher information matrix for a Poisson process.
 * <p>
 * This procedure is based on computation of a modified Chi-squared expression to perform Weighted Least Squares
 * Estimation assuming a Poisson model with a Gaussian noise component. The weight per observation is equal to
 * 1/[variance + max(y, 0) + 1].
 * <p>
 * See Ruisheng, et al (2017) Algorithmic corrections for localization microscopy with sCMOS cameras - characterisation
 * of a computationally efficient localization approach. Optical Express 25, Issue 10, pp 11701-11716.
 */
public class WPoissonGradientProcedure4 extends WPoissonGradientProcedure
{
	/**
	 * @param y
	 *            Data to fit
	 * @param var
	 *            the base variance of each observation (must be positive)
	 * @param func
	 *            Gradient function
	 */
	public WPoissonGradientProcedure4(final double[] y, final double[] var, final Gradient1Function func)
	{
		super(y, var, func);
		if (n != 4)
			throw new IllegalArgumentException("Function must compute 4 gradients");
	}

	@Override
	public void execute(double value, double[] dy_da)
	{
		// Note: Ignore the value
		final double w = this.w[yi++];
		data[0] += dy_da[0] * w * dy_da[0];
		double wgt;
		wgt = dy_da[1] * w;
		data[1] += wgt * dy_da[0];
		data[2] += wgt * dy_da[1];
		wgt = dy_da[2] * w;
		data[3] += wgt * dy_da[0];
		data[4] += wgt * dy_da[1];
		data[5] += wgt * dy_da[2];
		wgt = dy_da[3] * w;
		data[6] += wgt * dy_da[0];
		data[7] += wgt * dy_da[1];
		data[8] += wgt * dy_da[2];
		data[9] += wgt * dy_da[3];
	}

	@Override
	protected void initialiseWorkingMatrix()
	{
		GradientProcedureHelper.initialiseWorkingMatrix5(data);
	}

	@Override
	public void getMatrix(double[][] matrix)
	{
		GradientProcedureHelper.getMatrix5(data, matrix);
	}

	@Override
	public void getLinear(double[] matrix)
	{
		GradientProcedureHelper.getMatrix5(data, matrix);
	}
}
