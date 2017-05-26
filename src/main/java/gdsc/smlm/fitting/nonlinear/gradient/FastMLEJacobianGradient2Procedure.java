package gdsc.smlm.fitting.nonlinear.gradient;

import java.util.Arrays;

import gdsc.smlm.function.ExtendedGradient2Function;
import gdsc.smlm.function.ExtendedGradient2Procedure;
import gdsc.smlm.function.Gradient2Function;
import gdsc.smlm.function.Gradient2Procedure;

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
 * Computes the Jacobian matrix of the partial derivatives, dFi/dxj, for all n parameters. dFi is the first partial
 * derivative of the log likelihood function with respect to parameter i. dFi/dxj is the first partial derivative of dFi
 * with respect to parameter j.
 * <p>
 * Ref: Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically minimum uncertainty.
 * Nature Methods 7, 373-375 (supplementary note), Eq. 12.
 */
public class FastMLEJacobianGradient2Procedure extends FastMLEGradient2Procedure implements ExtendedGradient2Procedure
{
	protected final ExtendedGradient2Function func;

	/**
	 * The Jacobian matrix of the partial derivatives, dFi/dxj, for all n parameters. dFi is the first partial
	 * derivative of the log likelihood function with respect to parameter i. dFi/dxj is the first partial derivative of
	 * dFi with respect to parameter j.
	 */
	public final double[] J;

	/**
	 * @param x
	 *            Data to fit (must be positive, i.e. the value of a Poisson process)
	 * @param func
	 *            Gradient function (must produce a strictly positive value, i.e. the mean of a Poisson process)
	 */
	public FastMLEJacobianGradient2Procedure(final double[] x, final ExtendedGradient2Function func)
	{
		super(x, func);
		J = new double[n * n];
		this.func = func;
	}

	/**
	 * Calculates the first and second derivative of the Poisson log likelihood with respect to each parameter.
	 *
	 * @param a
	 *            Set of coefficients for the function
	 */
	public void computeJacobian(final double[] a)
	{
		k = 0;
		resetExtended2();
		func.initialiseExtended2(a);
		func.forEach((ExtendedGradient2Procedure) this);
		for (int i = 0, j = 0; i < n; i++, j += n + 1)
			d2[i] = J[j];
	}

	protected void resetExtended2()
	{
		Arrays.fill(d1, 0);
		Arrays.fill(J, 0);
	}

	public void executeExtended(double uk, double[] duk_dt, double[] d2uk_dtds)
	{
		u[k] = uk;
		final double xk = x[k++];

		final double xk_uk_minus1 = xk / uk - 1.0;
		final double xk_uk2 = xk / (uk * uk);
		for (int i = 0, index = 0; i < n; i++)
		{
			d1[i] += duk_dt[i] * xk_uk_minus1;

			for (int j = 0; j < n; j++, index++)
			{
				// This requires the partial second derivative with respect to i and j
				J[index] += d2uk_dtds[index] * xk_uk_minus1 - duk_dt[i] * duk_dt[j] * xk_uk2;
			}
		}
	}
	
	@Override
	public boolean isNaNGradients()
	{
		for (int i = n; i-- > 0;)
		{
			if (Double.isNaN(d1[i]))
				return true;
		}
		for (int i = J.length; i-- > 0;)
		{
			if (Double.isNaN(J[i]))
				return true;
		}
		return false;
	}
}
