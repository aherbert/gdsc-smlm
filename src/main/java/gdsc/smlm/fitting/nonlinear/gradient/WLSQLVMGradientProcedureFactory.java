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
 * Create a gradient procedure.
 */
public class WLSQLVMGradientProcedureFactory
{
	/**
	 * Create a new gradient procedure
	 * 
	 * @param y
	 *            Data to fit
	 * @param var
	 *            the base variance of each observation (must be positive)
	 * @param func
	 *            Gradient function
	 * @return the gradient procedure
	 */
	public static WLSQLVMGradientProcedure create(final double[] y, final double[] var, final Gradient1Function func)
	{
		switch (func.getNumberOfGradients())
		{
			case 5:
				return new WLSQLVMGradientProcedure5(y, var, func);
			case 4:
				return new WLSQLVMGradientProcedure4(y, var, func);
			case 6:
				return new WLSQLVMGradientProcedure6(y, var, func);

			default:
				return new WLSQLVMGradientProcedure(y, var, func);
		}
	}
}
