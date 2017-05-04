package gdsc.smlm.fitting.nonlinear.gradient;

import gdsc.smlm.function.GradientFunction;

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
public class LSQGradientProcedureMatrixFactory
{
	/**
	 * Create a new gradient calculator
	 * 
	 * @param y
	 *            Data to fit
	 * @param func
	 *            Gradient function
	 * @return the gradient procedure
	 */
	public static LSQGradientProcedureMatrix create(final double[] y, final GradientFunction func)
	{
		switch (func.getNumberOfGradients())
		{
			case 6:
				// free circular single Gaussian
				return new LSQGradientProcedureMatrix6(y, func);

			default:
				return new LSQGradientProcedureMatrix(y, func);
		}
	}
}
