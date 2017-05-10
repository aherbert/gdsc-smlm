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
 * Nature Methods 7, 373-375 (supplementary note), Eq. 12.
 */
public class PoissonGradientProcedure5 extends PoissonGradientProcedure
{
	/**
	 * @param func
	 *            Gradient function
	 */
	public PoissonGradientProcedure5(final Gradient1Function func)
	{
		super(func);
		if (n != 5)
			throw new IllegalArgumentException("Function must compute 5 gradients");
	}

	@Override
	public void execute(double value, double[] dy_da)
	{
		if (value > 0)
		{
			final double f = 1.0 / value;

			data[0] += dy_da[0] * f * dy_da[0];
			double w;
			w = dy_da[1] * f;
			data[1] += w * dy_da[0];
			data[2] += w * dy_da[1];
			w = dy_da[2] * f;
			data[3] += w * dy_da[0];
			data[4] += w * dy_da[1];
			data[5] += w * dy_da[2];
			w = dy_da[3] * f;
			data[6] += w * dy_da[0];
			data[7] += w * dy_da[1];
			data[8] += w * dy_da[2];
			data[9] += w * dy_da[3];
			w = dy_da[4] * f;
			data[10] += w * dy_da[0];
			data[11] += w * dy_da[1];
			data[12] += w * dy_da[2];
			data[13] += w * dy_da[3];
			data[14] += w * dy_da[4];
		}
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
