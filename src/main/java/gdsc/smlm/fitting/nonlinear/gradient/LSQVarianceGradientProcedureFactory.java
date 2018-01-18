package gdsc.smlm.fitting.nonlinear.gradient;

import gdsc.smlm.fitting.linear.EJMLLinearSolver;
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
 * Create a gradient procedure.
 */
public class LSQVarianceGradientProcedureFactory
{
	/**
	 * Create a new gradient procedure.
	 *
	 * @param func
	 *            Gradient function
	 * @return the gradient procedure
	 */
	public static LSQVarianceGradientProcedure create(final Gradient1Function func)
	{
		switch (func.getNumberOfGradients())
		{
			case 5:
				return new LSQVarianceGradientProcedure5(func);
			case 4:
				return new LSQVarianceGradientProcedure4(func);
			case 6:
				return new LSQVarianceGradientProcedure6(func);
			default:
				return new LSQVarianceGradientProcedure(func);
		}
	}

	/**
	 * Create a new gradient procedure.
	 *
	 * @param func
	 *            Gradient function
	 * @return the gradient procedure
	 */
	public static LSQVarianceGradientProcedure create(final Gradient1Function func, EJMLLinearSolver solver)
	{
		switch (func.getNumberOfGradients())
		{
			case 5:
				return new LSQVarianceGradientProcedure5(func, solver);
			case 4:
				return new LSQVarianceGradientProcedure4(func, solver);
			case 6:
				return new LSQVarianceGradientProcedure6(func, solver);
			default:
				return new LSQVarianceGradientProcedure(func, solver);
		}
	}
}
