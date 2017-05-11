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
public class NewtonRaphsonGradient2ProcedureB extends NewtonRaphsonGradient2Procedure
{
	protected final double[] b;

	/**
	 * @param x
	 *            Data to fit (must be positive, i.e. the value of a Poisson process)
	 * @param b
	 *            Baseline pre-computed y-values
	 * @param func
	 *            Gradient function (must produce a strictly positive value, i.e. the mean of a Poisson process)
	 */
	public NewtonRaphsonGradient2ProcedureB(final double[] x, final double[] b, final Gradient2Function func)
	{
		super(x, func);
		this.b = b;
	}

	@Override
	public void execute(double uk, double[] duk_dt, double[] d2uk_dt2)
	{
		super.execute(uk + b[k], duk_dt, d2uk_dt2);
	}

	@Override
	public void execute(double uk, double[] duk_dt)
	{
		super.execute(uk + b[k], duk_dt);
	}

	@Override
	public void execute(double uk)
	{
		super.execute(uk + b[k]);
	}
}
