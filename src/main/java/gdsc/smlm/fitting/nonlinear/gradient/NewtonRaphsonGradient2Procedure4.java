package gdsc.smlm.fitting.nonlinear.gradient;

import gdsc.smlm.function.Gradient2Function;

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
 * Calculates the Newton-Raphson update vector for a Poisson process using the first and second partial derivatives.
 * <p>
 * Ref: Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically minimum uncertainty.
 * Nature Methods 7, 373-375 (supplementary note), Eq. 12.
 */
public class NewtonRaphsonGradient2Procedure4 extends NewtonRaphsonGradient2Procedure
{
	/**
	 * @param x
	 *            Data to fit (must be positive, i.e. the value of a Poisson process)
	 * @param func
	 *            Gradient function (must produce a strictly positive value, i.e. the mean of a Poisson process)
	 */
	public NewtonRaphsonGradient2Procedure4(final double[] x, final Gradient2Function func)
	{
		super(x, func);
		if (n != 4)
			throw new IllegalArgumentException("Function must compute 4 gradients");
	}

	@Override
	protected void reset2()
	{
		d1[0] = 0;
		d1[1] = 0;
		d1[2] = 0;
		d1[3] = 0;
		d2[0] = 0;
		d2[1] = 0;
		d2[2] = 0;
		d2[3] = 0;
	}

	@Override
	public void execute(double uk, double[] duk_dt, double[] d2uk_dt2)
	{
		u[k] = uk;
		final double xk = x[k++];
		final double xk_uk_minus1 = xk / uk - 1.0;
		final double xk_uk2 = xk / (uk * uk);
		d1[0] += duk_dt[0] * xk_uk_minus1;
		d1[1] += duk_dt[1] * xk_uk_minus1;
		d1[2] += duk_dt[2] * xk_uk_minus1;
		d1[3] += duk_dt[3] * xk_uk_minus1;
		d2[0] += d2uk_dt2[0] * xk_uk_minus1 - duk_dt[0] * duk_dt[0] * xk_uk2;
		d2[1] += d2uk_dt2[1] * xk_uk_minus1 - duk_dt[1] * duk_dt[1] * xk_uk2;
		d2[2] += d2uk_dt2[2] * xk_uk_minus1 - duk_dt[2] * duk_dt[2] * xk_uk2;
		d2[3] += d2uk_dt2[3] * xk_uk_minus1 - duk_dt[3] * duk_dt[3] * xk_uk2;
	}

	@Override
	public double[] getUpdate()
	{
		double[] update = new double[n];
		update[0] = d1[0] / d2[0];
		update[1] = d1[1] / d2[1];
		update[2] = d1[2] / d2[2];
		update[3] = d1[3] / d2[3];
		return update;
	}

	/**
	 * Reset the first derivative vector
	 */
	protected void reset1()
	{
		d1[0] = 0;
		d1[1] = 0;
		d1[2] = 0;
		d1[3] = 0;
	}

	@Override
	public void execute(double uk, double[] duk_dt)
	{
		u[k] = uk;
		final double xk = x[k++];
		final double xk_uk_minus1 = xk / uk - 1.0;
		d1[0] += duk_dt[0] * xk_uk_minus1;
		d1[1] += duk_dt[1] * xk_uk_minus1;
		d1[2] += duk_dt[2] * xk_uk_minus1;
		d1[3] += duk_dt[3] * xk_uk_minus1;
	}
}
