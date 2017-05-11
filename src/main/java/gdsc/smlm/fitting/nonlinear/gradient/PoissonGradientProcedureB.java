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
 * Calculates the Fisher information matrix for a Poisson process.
 * <p>
 * Ref: Smith et al, (2010). Fast, single-molecule localisation that achieves theoretically minimum uncertainty.
 * Nature Methods 7, 373-375 (supplementary note), Eq. 9.
 */
public class PoissonGradientProcedureB extends PoissonGradientProcedure
{
	protected final double[] b;
	protected int bi;

	/**
	 * @param b
	 *            Baseline pre-computed y-values
	 * @param func
	 *            Gradient function
	 */
	public PoissonGradientProcedureB(final double[] b, final Gradient1Function func)
	{
		super(func);
		this.b = b;
	}

	@Override
	public void computeFisherInformation(final double[] a)
	{
		bi = -1;
		super.computeFisherInformation(a);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.function.Gradient1Procedure#execute(double, double[])
	 */
	public void execute(double value, double[] dy_da)
	{
		super.execute(value + b[++bi], dy_da);
	}
}
