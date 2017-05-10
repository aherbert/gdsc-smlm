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
public class PoissonGradientProcedure4 extends PoissonGradientProcedure
{
	/**
	 * @param func
	 *            Gradient function
	 */
	public PoissonGradientProcedure4(final Gradient1Function func)
	{
		super(func);
		if (n != 4)
			throw new IllegalArgumentException("Function must compute 4 gradients");
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
		}
	}

	@Override
	protected void initialiseWorkingMatrix()
	{
		GradientProcedureHelper.initialiseWorkingMatrix4(data);
	}

	@Override
	public void getMatrix(double[][] matrix)
	{
		GradientProcedureHelper.getMatrix4(data, matrix);
	}

	@Override
	public void getLinear(double[] matrix)
	{
		GradientProcedureHelper.getMatrix4(data, matrix);
	}
}
